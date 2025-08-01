import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from bioservices import KEGG
import base64
import itertools

# --- Page and App Configuration ---
st.set_page_config(page_title="Cancer Pathway Analyzer", layout="wide")
st.title("🧬 Advanced Cancer Pathway Analyzer")

# --- Helper Functions for UI ---
def get_icon_img_tag(icon_path, height=20):
    """Encodes a local icon image for embedding in HTML."""
    try:
        with open(icon_path, "rb") as f:
            data = f.read()
            encoded = base64.b64encode(data).decode()
            ext = icon_path.split(".")[-1]
            return f'<img src="data:image/{ext};base64,{encoded}" height="{height}px"/>'
    except FileNotFoundError:
        return "" # Return empty string if icon is not found

def create_gene_row(gene):
    """Generates an HTML table row for a gene with links to external databases."""
    base_urls = {
        "GeneCards": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}",
        "NCBI": f"https://www.ncbi.nlm.nih.gov/gene/?term={gene}",
        "ENSEMBL": f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene}",
        "GEO": f"https://www.ncbi.nlm.nih.gov/gds/?term={gene}"
    }
    icon_folder = "Icons" # Assumes a local folder named "Icons"
    icons = {
        "GeneCards": get_icon_img_tag(f"{icon_folder}/Gene Cards.png"),
        "NCBI": get_icon_img_tag(f"{icon_folder}/NCBI.png"),
        "ENSEMBL": get_icon_img_tag(f"{icon_folder}/ENSEMBL.png"),
        "GEO": get_icon_img_tag(f"{icon_folder}/GEO.png")
    }
    links = " ".join([f'<a href="{base_urls[name]}" target="_blank">{icons[name]}</a>' for name in icons])
    return f"<tr><td>{gene}</td><td>{links}</td></tr>"

def display_gene_table(genes):
    """Creates and displays a markdown table for a list of genes."""
    gene_table_html = """
    <table style="width:100%; border-collapse: collapse;" border="1" cellpadding="6">
      <thead style="background-color:#f2f2f2;">
        <tr><th>Gene</th><th>External Database Links</th></tr>
      </thead>
      <tbody>
    """
    for gene in sorted(genes):
        gene_table_html += create_gene_row(gene)
    gene_table_html += "</tbody></table>"
    st.markdown(gene_table_html, unsafe_allow_html=True)

# --- Data Caching and Retrieval Functions ---
@st.cache_data(show_spinner="Fetching KEGG cancer pathways...")
def get_kegg_cancer_pathways():
    """Fetches a list of human cancer-related pathways from KEGG."""
    k = KEGG()
    raw = k.list("pathway/hsa")
    cancer_pathways = {}
    keywords = ['cancer', 'leukemia', 'lymphoma', 'melanoma', 'glioma', 'carcinoma', 'sarcoma']
    for line in raw.strip().split("\n"):
        pid, desc = line.split("\t")
        name = desc.split(" - ")[0]
        if any(kw in name.lower() for kw in keywords):
            cancer_pathways[name] = pid.replace("path:", "")
    return cancer_pathways

@st.cache_data(show_spinner="Fetching genes for pathway: {pathway_id}...")
def get_genes_from_kegg(pathway_id):
    """Retrieves all genes associated with a specific KEGG pathway ID."""
    k = KEGG()
    k.organism = "hsa"
    data = k.get(pathway_id)
    genes = set()
    if not data:
        return genes
    
    in_gene_section = False
    for line in data.strip().split('\n'):
        if line.startswith('GENE'):
            in_gene_section = True
            line_content = line[12:].strip()
        elif in_gene_section and line.startswith(' '):
            line_content = line.strip()
        else:
            in_gene_section = False
            continue
        
        if line_content:
            parts = line_content.split()
            # Gene symbol is typically the second part, after the entry number
            if len(parts) > 1 and parts[0].isdigit():
                gene_symbol = parts[1].replace(';', '')
                genes.add(gene_symbol)
    return genes

@st.cache_data
def get_compound_gene_targets():
    """Returns a predefined dictionary of compounds and their target genes."""
    return {
        # Chemotherapy Drugs
        "Cisplatin": {"ERCC1", "XPA", "MLH1", "MSH2", "TP53", "ATM", "ATR"},
        "Doxorubicin": {"TOP2A", "TOP2B", "H2AFX", "CBR1", "RAD51"},
        "Paclitaxel": {"TUBB1", "MAP2", "MAP4", "BCL2", "AURKA"},
        "Methotrexate": {"DHFR", "TYMS", "ATIC", "GART"},
        "5-Fluorouracil": {"TYMS", "DPYD", "UPP1", "TP53"},
        "Carboplatin": {"BRCA1", "BRCA2", "ERCC1"},
        "Etoposide": {"TOP2A", "TOP2B", "TP53"},
        "Vincristine": {"TUBB1", "TUBA1B"},
        "Cyclophosphamide": {"TP53", "ATM", "BAX"},
        "Bleomycin": {"TP53", "ATM", "H2AFX"},
        # Natural Products
        "Curcumin": {"NFKB1", "STAT3", "TNF", "IL6", "VEGFA", "EGFR", "PPARG"},
        "Resveratrol": {"SIRT1", "NFKB1", "TP53", "VEGFA", "CASP3", "NFE2L2"},
        "Quercetin": {"PIK3CA", "AKT1", "MTOR", "TP53", "CASP3", "NFKB1"},
        "Genistein": {"ESR1", "EGFR", "NFKB1", "AKT1", "PTK2B"},
        "EGCG": {"EGFR", "VEGFA", "HDAC1", "DNMT1", "FASN"},
        "Berberine": {"MAPK1", "TP53", "BCL2", "CASP3"},
        "Ginsenoside": {"AKT1", "MTOR", "CASP3"},
        "Luteolin": {"NFKB1", "TNF", "IL6"},
        "Sulforaphane": {"NFE2L2", "TP53", "CASP3"},
        "Apigenin": {"AKT1", "CASP3", "BCL2"},
    }

# --- Load Static Data ---
cancer_pathways_dict = get_kegg_cancer_pathways()
compound_targets_dict = get_compound_gene_targets()
chemo_drugs_list = list(compound_targets_dict.keys())[:10]
natural_products_list = list(compound_targets_dict.keys())[10:]

# --- Sidebar UI ---
st.sidebar.header("🔬 Configure Analysis")
analysis_mode = st.sidebar.radio(
    "Select Analysis Mode:",
    ('Chemotherapy vs. Cancer Pathways', 'Natural Products vs. Cancer Pathways', 'Cancer Pathway vs. Cancer Pathway'),
    key="analysis_mode"
)
st.sidebar.markdown("---")


# --- Main Application Logic ---

# MODE 1: Chemotherapy vs. Cancer Pathways
if analysis_mode == 'Chemotherapy vs. Cancer Pathways':
    st.sidebar.subheader("Analysis Selections")
    selected_pathways = st.sidebar.multiselect("Select Cancer Pathways", options=list(cancer_pathways_dict.keys()), default=[list(cancer_pathways_dict.keys())[0]])
    selected_chemo = st.sidebar.multiselect("Select Chemotherapy Compounds", options=chemo_drugs_list, default=[chemo_drugs_list[0]])

    if not selected_pathways or not selected_chemo:
        st.warning("⚠️ Please select at least one cancer pathway and one chemotherapy compound.")
        st.stop()

    # --- Gene Intersection Logic ---
    pathway_genes_map = {p: get_genes_from_kegg(cancer_pathways_dict[p]) for p in selected_pathways}
    compound_genes_map = {c: compound_targets_dict[c] for c in selected_chemo}

    overlap_records = []
    for pathway_name, pathway_genes in pathway_genes_map.items():
        for compound, comp_genes in compound_genes_map.items():
            shared_genes = pathway_genes.intersection(comp_genes)
            if shared_genes:
                overlap_records.append({
                    "Item 1": compound,
                    "Item 2": pathway_name,
                    "Shared Gene Count": len(shared_genes),
                    "Shared Genes": shared_genes
                })
    
    # --- Display Results ---
    st.header("Analysis: Chemotherapy vs. Cancer Pathways")
    if not overlap_records:
        st.warning("No overlapping genes found for the current selection.")
    else:
        st.subheader("Shared Gene Analysis")
        for record in overlap_records:
            st.markdown(f"#### 🧬 **{record['Item 1']}** &harr; **{record['Item 2']}**")
            st.markdown(f"**Shared Genes:** `{record['Shared Gene Count']}`")
            display_gene_table(record['Shared Genes'])
            st.markdown("---")

        # --- Network Graph ---
        st.subheader("Interaction Network Visualization")
        G = nx.Graph()
        # Add nodes with specific types for the legend
        for p in selected_pathways: G.add_node(p, type="pathway", size=25)
        for c in selected_chemo: G.add_node(c, type="chemo", size=20)
        
        all_shared_genes = set().union(*[r['Shared Genes'] for r in overlap_records])
        for gene in all_shared_genes: G.add_node(gene, type="gene", size=15)

        # Add edges from records
        for record in overlap_records:
            for gene in record['Shared Genes']:
                G.add_edge(record['Item 1'], gene)
                G.add_edge(record['Item 2'], gene)

# MODE 2: Natural Products vs. Cancer Pathways
elif analysis_mode == 'Natural Products vs. Cancer Pathways':
    st.sidebar.subheader("Analysis Selections")
    selected_pathways = st.sidebar.multiselect("Select Cancer Pathways", options=list(cancer_pathways_dict.keys()), default=[list(cancer_pathways_dict.keys())[0]])
    selected_natural = st.sidebar.multiselect("Select Natural Products", options=natural_products_list, default=[natural_products_list[0]])

    if not selected_pathways or not selected_natural:
        st.warning("⚠️ Please select at least one cancer pathway and one natural product.")
        st.stop()
        
    # --- Gene Intersection Logic ---
    pathway_genes_map = {p: get_genes_from_kegg(cancer_pathways_dict[p]) for p in selected_pathways}
    compound_genes_map = {c: compound_targets_dict[c] for c in selected_natural}

    overlap_records = []
    for pathway_name, pathway_genes in pathway_genes_map.items():
        for compound, comp_genes in compound_genes_map.items():
            shared_genes = pathway_genes.intersection(comp_genes)
            if shared_genes:
                overlap_records.append({
                    "Item 1": compound,
                    "Item 2": pathway_name,
                    "Shared Gene Count": len(shared_genes),
                    "Shared Genes": shared_genes
                })
    
    # --- Display Results ---
    st.header("Analysis: Natural Products vs. Cancer Pathways")
    if not overlap_records:
        st.warning("No overlapping genes found for the current selection.")
    else:
        st.subheader("Shared Gene Analysis")
        for record in overlap_records:
            st.markdown(f"#### 🌿 **{record['Item 1']}** &harr; **{record['Item 2']}**")
            st.markdown(f"**Shared Genes:** `{record['Shared Gene Count']}`")
            display_gene_table(record['Shared Genes'])
            st.markdown("---")

        # --- Network Graph ---
        st.subheader("Interaction Network Visualization")
        G = nx.Graph()
        # Add nodes with specific types for the legend
        for p in selected_pathways: G.add_node(p, type="pathway", size=25)
        for c in selected_natural: G.add_node(c, type="natural", size=20)
        
        all_shared_genes = set().union(*[r['Shared Genes'] for r in overlap_records])
        for gene in all_shared_genes: G.add_node(gene, type="gene", size=15)
        
        # Add edges
        for record in overlap_records:
            for gene in record['Shared Genes']:
                G.add_edge(record['Item 1'], gene)
                G.add_edge(record['Item 2'], gene)


# MODE 3: Cancer Pathway vs. Cancer Pathway
elif analysis_mode == 'Cancer Pathway vs. Cancer Pathway':
    st.sidebar.subheader("Analysis Selections")
    selected_pathways = st.sidebar.multiselect("Select Cancer Pathways", options=list(cancer_pathways_dict.keys()), default=list(cancer_pathways_dict.keys())[:2])

    if len(selected_pathways) < 2:
        st.warning("⚠️ Please select at least two cancer pathways to compare.")
        st.stop()

    # --- Gene Intersection Logic ---
    pathway_genes_map = {p: get_genes_from_kegg(cancer_pathways_dict[p]) for p in selected_pathways}

    overlap_records = []
    for (p1_name, p2_name) in itertools.combinations(selected_pathways, 2):
        p1_genes = pathway_genes_map[p1_name]
        p2_genes = pathway_genes_map[p2_name]
        shared_genes = p1_genes.intersection(p2_genes)
        if shared_genes:
            overlap_records.append({
                "Item 1": p1_name,
                "Item 2": p2_name,
                "Shared Gene Count": len(shared_genes),
                "Shared Genes": shared_genes
            })

    # --- Display Results ---
    st.header("Analysis: Cancer Pathway vs. Cancer Pathway")
    if not overlap_records:
        st.warning("No overlapping genes found between the selected pathways.")
    else:
        st.subheader("Shared Gene Analysis")
        for record in overlap_records:
            st.markdown(f"#### ↔️ **{record['Item 1']}** &harr; **{record['Item 2']}**")
            st.markdown(f"**Shared Genes:** `{record['Shared Gene Count']}`")
            display_gene_table(record['Shared Genes'])
            st.markdown("---")

        # --- Network Graph ---
        st.subheader("Pathway Overlap Network")
        G = nx.Graph()
        # Add nodes with specific types for the legend
        all_shared_genes = set().union(*[r['Shared Genes'] for r in overlap_records])
        for p in selected_pathways: G.add_node(p, type="pathway", size=25)
        for gene in all_shared_genes: G.add_node(gene, type="gene", size=15)
        
        # Add edges
        for record in overlap_records:
            for gene in record['Shared Genes']:
                G.add_edge(record['Item 1'], gene)
                G.add_edge(record['Item 2'], gene)

# --- Common Plotly Network Drawing Logic (if a graph was generated) ---
if 'G' in locals() and G.number_of_nodes() > 0:
    pos = nx.spring_layout(G, k=0.35, iterations=50, seed=42)
    
    # Define traces for edges and each node type
    edge_trace = go.Scatter(x=[], y=[], line=dict(width=0.7, color='#888'), hoverinfo='none', mode='lines')
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace['x'] += tuple([x0, x1, None])
        edge_trace['y'] += tuple([y0, y1, None])

    # --- Create a trace for each node type to build the legend ---
    node_traces = []
    legend_map = {
        'pathway': {'name': 'Cancer Pathway', 'color': 'purple'},
        'chemo': {'name': 'Chemo Compound', 'color': 'orange'},
        'natural': {'name': 'Natural Product', 'color': 'green'},
        'gene': {'name': 'Shared Gene', 'color': 'magenta'}
    }
    
    # Group nodes by type
    nodes_by_type = {}
    for node, data in G.nodes(data=True):
        node_type = data['type']
        if node_type not in nodes_by_type:
            nodes_by_type[node_type] = []
        nodes_by_type[node_type].append(node)

    # Create a scatter plot for each node type
    for node_type, nodes in nodes_by_type.items():
        if not nodes:
            continue
        
        x, y, text, size = [], [], [], []
        for node in nodes:
            x.append(pos[node][0])
            y.append(pos[node][1])
            text.append(f"{node_type.replace('_', ' ').title()}:<br>{node}")
            size.append(G.nodes[node]['size'])
            
        trace = go.Scatter(
            x=x, y=y, text=text, mode='markers', hoverinfo='text',
            name=legend_map[node_type]['name'], # This name appears in the legend
            marker=dict(
                color=legend_map[node_type]['color'],
                size=size,
                line_width=2
            )
        )
        node_traces.append(trace)

    # --- Create the Figure with CORRECTED layout ---
    fig = go.Figure(data=[edge_trace] + node_traces,
                    layout=go.Layout(
                        title=dict(text=f'Network Graph for {analysis_mode}', x=0.5),
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=80), # MODIFIED: Increased top margin
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        legend=dict(
                            orientation="h",
                            yanchor="bottom",
                            y=1.02,
                            xanchor="right",
                            x=1
                        ))
                    )
    st.plotly_chart(fig, use_container_width=True)
