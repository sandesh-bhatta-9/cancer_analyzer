import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from bioservices import KEGG
import base64
import itertools
from collections import Counter

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
show_all_genes = st.sidebar.checkbox("Show all individual genes (may be slow/cluttered)")
st.sidebar.markdown("---")

# --- Main Application Logic ---
G = nx.Graph()

# Get selections from sidebar
st.sidebar.subheader("Analysis Selections")
if analysis_mode != 'Cancer Pathway vs. Cancer Pathway':
    is_chemo_mode = analysis_mode == 'Chemotherapy vs. Cancer Pathways'
    selected_pathways = st.sidebar.multiselect("Select Cancer Pathways", options=list(cancer_pathways_dict.keys()), default=[list(cancer_pathways_dict.keys())[0]])
    if is_chemo_mode:
        selected_compounds = st.sidebar.multiselect("Select Chemotherapy Compounds", options=chemo_drugs_list, default=[chemo_drugs_list[0]])
    else:
        selected_compounds = st.sidebar.multiselect("Select Natural Products", options=natural_products_list, default=[natural_products_list[0]])
else:
    selected_pathways = st.sidebar.multiselect("Select Cancer Pathways", options=list(cancer_pathways_dict.keys()), default=list(cancer_pathways_dict.keys())[:2])

# Display Headers and Tables first
st.header(f"Analysis: {analysis_mode}")
if analysis_mode != 'Cancer Pathway vs. Cancer Pathway':
    if not selected_pathways or not selected_compounds:
        st.warning("⚠️ Please select at least one pathway and one compound.")
        st.stop()
    
    st.subheader("Shared Gene Analysis")
    pathway_genes_map = {p: get_genes_from_kegg(cancer_pathways_dict[p]) for p in selected_pathways}
    compound_genes_map = {c: compound_targets_dict[c] for c in selected_compounds}
    has_overlap = False
    for p_name, p_genes in pathway_genes_map.items():
        for c_name, c_genes in compound_genes_map.items():
            shared = p_genes.intersection(c_genes)
            if shared:
                has_overlap = True
                emoji = "🧬" if is_chemo_mode else "🌿"
                st.markdown(f"#### {emoji} **{c_name}** &harr; **{p_name}**")
                st.markdown(f"**Shared Genes:** `{len(shared)}`")
                display_gene_table(shared)
                st.markdown("---")
    if not has_overlap:
        st.warning("No overlapping genes found for the current selection.")

else: # Pathway vs Pathway Tables
    if len(selected_pathways) < 2:
        st.warning("⚠️ Please select at least two pathways to compare.")
        st.stop()

    st.subheader("Shared Gene Analysis")
    pathway_genes_map = {p: get_genes_from_kegg(cancer_pathways_dict[p]) for p in selected_pathways}
    has_overlap = False
    for (p1_name, p2_name) in itertools.combinations(selected_pathways, 2):
        shared_genes_pair = pathway_genes_map.get(p1_name, set()).intersection(pathway_genes_map.get(p2_name, set()))
        if shared_genes_pair:
            has_overlap = True
            st.markdown(f"#### ↔️ **{p1_name}** &harr; **{p2_name}**")
            st.markdown(f"**Shared Genes:** `{len(shared_genes_pair)}`")
            display_gene_table(shared_genes_pair)
            st.markdown("---")
    if not has_overlap:
        st.info("No overlapping genes found between any pair of the selected pathways.")


# --- Graph Building Logic ---
if show_all_genes: # DETAILED VIEW
    if analysis_mode != 'Cancer Pathway vs. Cancer Pathway':
        for p_name in selected_pathways:
            for c_name in selected_compounds:
                pathway_genes = pathway_genes_map.get(p_name, set())
                compound_genes = compound_genes_map.get(c_name, set())
                shared_genes = pathway_genes.intersection(compound_genes)

                p_color, c_color = 'purple', 'orange' if is_chemo_mode else 'green'
                G.add_node(p_name, type='pathway', size=20, color=p_color)
                G.add_node(c_name, type='chemo' if is_chemo_mode else 'natural', size=20, color=c_color)

                for gene in pathway_genes:
                    if gene in shared_genes: G.add_node(gene, type='gene', size=15, color='magenta')
                    else: G.add_node(gene, type='unique_pathway_gene', size=8, color='#a6cee3')
                    G.add_edge(p_name, gene)
                
                for gene in compound_genes: # Add edges for compound genes
                    G.add_edge(c_name, gene)
                    if gene not in shared_genes: # Add unique compound genes
                        G.add_node(gene, type='unique_compound_gene', size=8, color='#b2df8a')
    else: # Pathway vs Pathway Detailed View
        light_colors = ['#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']
        pathway_color_map = {name: light_colors[i % len(light_colors)] for i, name in enumerate(selected_pathways)}
        all_genes_flat = [gene for genes in pathway_genes_map.values() for gene in genes]
        gene_counts = Counter(all_genes_flat)
        shared_genes_set = {gene for gene, count in gene_counts.items() if count > 1}

        for pathway_name, genes in pathway_genes_map.items():
            G.add_node(pathway_name, type='pathway', size=20, color='purple')
            for gene in genes:
                if gene in shared_genes_set: G.add_node(gene, type='gene', size=15, color='magenta')
                else: G.add_node(gene, type=f'unique_gene_{pathway_name}', size=8, color=pathway_color_map[pathway_name])
                G.add_edge(pathway_name, gene)

else: # SUMMARY VIEW (DEFAULT)
    if analysis_mode != 'Cancer Pathway vs. Cancer Pathway':
        for p_name in selected_pathways:
            for c_name in selected_compounds:
                pathway_genes = pathway_genes_map.get(p_name, set())
                compound_genes = compound_genes_map.get(c_name, set())
                shared_genes = pathway_genes.intersection(compound_genes)

                if shared_genes:
                    p_color, c_color = 'purple', 'orange' if is_chemo_mode else 'green'
                    G.add_node(p_name, type='pathway', size=30, color=p_color)
                    G.add_node(c_name, type='chemo' if is_chemo_mode else 'natural', size=30, color=c_color)

                    for gene in shared_genes:
                        G.add_node(gene, type='gene', size=15, color='magenta')
                        G.add_edge(p_name, gene); G.add_edge(c_name, gene)

                    unique_p_count = len(pathway_genes - shared_genes)
                    if unique_p_count > 0:
                        summary_p_label = f"{p_name}\n({unique_p_count} Unique Genes)"
                        G.add_node(summary_p_label, type='unique_summary', size=25, color='#a6cee3')
                        G.add_edge(p_name, summary_p_label)
                    
                    unique_c_count = len(compound_genes - shared_genes)
                    if unique_c_count > 0:
                        summary_c_label = f"{c_name}\n({unique_c_count} Unique Targets)"
                        G.add_node(summary_c_label, type='unique_summary', size=25, color='#b2df8a')
                        G.add_edge(c_name, summary_c_label)
    else: # Pathway vs Pathway Summary View
        light_colors = ['#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']
        pathway_color_map = {name: light_colors[i % len(light_colors)] for i, name in enumerate(selected_pathways)}
        all_genes_flat = [gene for genes in pathway_genes_map.values() for gene in genes]
        gene_counts = Counter(all_genes_flat)
        shared_genes_set = {gene for gene, count in gene_counts.items() if count > 1}

        for pathway_name, genes in pathway_genes_map.items():
            G.add_node(pathway_name, type='pathway', size=30, color='purple')
            for gene in genes.intersection(shared_genes_set):
                G.add_node(gene, type='gene', size=15, color='magenta')
                G.add_edge(pathway_name, gene)
            
            unique_genes_count = len(genes - shared_genes_set)
            if unique_genes_count > 0:
                summary_node_label = f"{pathway_name}\n({unique_genes_count} Unique Genes)"
                G.add_node(summary_node_label, type='unique_summary', size=25, color=pathway_color_map[pathway_name])
                G.add_edge(pathway_name, summary_node_label)


# --- Unified Plotting Logic ---
if G.number_of_nodes() > 0:
    st.subheader("Context Network Visualization")
    if show_all_genes: st.info("Showing all individual genes. The network may be dense.")
    else: st.info("Displaying shared genes individually and grouping unique genes into summary nodes for clarity.")
    
    pos = nx.spring_layout(G, k=0.5, iterations=60, seed=42)
    edge_trace = go.Scatter(x=[], y=[], line=dict(width=0.7, color='#888'), hoverinfo='none', mode='lines')
    for edge in G.edges():
        x0, y0 = pos[edge[0]]; x1, y1 = pos[edge[1]]
        edge_trace['x'] += (x0, x1, None); edge_trace['y'] += (y0, y1, None)

    node_traces = []
    nodes_by_type = {}
    for node, data in G.nodes(data=True):
        ntype = data.get('type', 'unknown')
        if ntype not in nodes_by_type: nodes_by_type[ntype] = []
        nodes_by_type[ntype].append(node)
    
    legend_map = {
        'pathway': {'name': 'Cancer Pathway'}, 'chemo': {'name': 'Chemo Compound'}, 'natural': {'name': 'Natural Product'},
        'gene': {'name': 'Shared Gene'}, 'unique_summary': {'name': 'Unique Gene/Target Set'},
        'unique_pathway_gene': {'name': 'Pathway-Unique Gene'}, 'unique_compound_gene': {'name': 'Compound-Unique Target'}
    }

    for node_type, nodes in nodes_by_type.items():
        if nodes:
            x, y, text, size, colors = [], [], [], [], []
            for node in nodes:
                x.append(pos[node][0]); y.append(pos[node][1])
                hover_text = f"{node.replace(chr(10), '<br>')}"
                text.append(hover_text)
                size.append(G.nodes[node]['size'])
                colors.append(G.nodes[node]['color'])
            
            trace_name = "Other"
            if node_type.startswith('unique_gene_'):
                pathway_name = node_type.replace('unique_gene_', '')
                trace_name = f'{pathway_name} Unique Genes'
            elif node_type in legend_map:
                trace_name = legend_map[node_type]['name']
            
            node_traces.append(go.Scatter(
                x=x, y=y, text=text, mode='markers', hoverinfo='text', name=trace_name,
                marker=dict(color=colors, size=size, line_width=2)
            ))

    fig = go.Figure(data=[edge_trace] + node_traces, layout=go.Layout(
        showlegend=True, hovermode='closest', margin=dict(b=5, l=5, r=5, t=5),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="center", x=0.5)
    ))
    st.plotly_chart(fig, use_container_width=True)