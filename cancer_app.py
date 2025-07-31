import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from bioservices import KEGG

st.set_page_config(page_title="Cancer Pathway Analyzer", layout="wide")
st.title("🧬 Advanced Cancer Pathway Analyzer")

@st.cache_data(show_spinner=False)
def get_kegg_cancer_pathways():
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

@st.cache_data(show_spinner=False)
def get_genes_from_kegg(pathway_id):
    k = KEGG()
    k.organism = "hsa"
    data = k.get(pathway_id)
    genes = set()
    if not data:
        return genes
    lines = data.split("\n")
    gene_section = False
    for line in lines:
        if line.startswith("GENE"):
            gene_section = True
            gene_line = line[12:].strip()
        elif gene_section and line.startswith(" "):
            gene_line = line.strip()
        else:
            gene_section = False
            continue
        if gene_line:
            parts = gene_line.split()
            if parts and parts[0].isdigit():
                genes.add(parts[1].replace(";", ""))
    return genes

@st.cache_data(show_spinner=False)
def get_compound_gene_targets():
    # Expanded chemo & natural compounds and their known target genes (curated)
    return {
        # Chemotherapy drugs
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
        # Natural products
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

# Sidebar UI
st.sidebar.header("Configure Analysis")

cancer_pathways = get_kegg_cancer_pathways()
selected_pathways = st.sidebar.multiselect(
    "Select Cancer Pathways (multiple allowed)",
    options=list(cancer_pathways.keys()),
    default=[list(cancer_pathways.keys())[0]]
)

compound_targets = get_compound_gene_targets()
chemo_drugs = list(compound_targets.keys())[:10]
natural_products = list(compound_targets.keys())[10:]

selected_chemo = st.sidebar.multiselect(
    "Select Chemotherapy Compounds",
    options=chemo_drugs,
    default=[chemo_drugs[0]]
)

selected_natural = st.sidebar.multiselect(
    "Select Natural Products",
    options=natural_products,
    default=[natural_products[0]]
)

if not selected_pathways:
    st.warning("Please select at least one cancer pathway to analyze.")
    st.stop()

if not (selected_chemo or selected_natural):
    st.warning("Please select at least one chemotherapy compound or natural product.")
    st.stop()

# Fetch genes for selected cancer pathways
cancer_genes_map = {}
for p in selected_pathways:
    pid = cancer_pathways[p]
    cancer_genes_map[p] = get_genes_from_kegg(pid)

# Combine selected compounds genes
selected_compounds = selected_chemo + selected_natural
compound_genes_map = {c: compound_targets[c] for c in selected_compounds}

# Display summary
st.header("Analysis Summary")
st.write(f"**Selected cancer pathways:** {', '.join(selected_pathways)}")
st.write(f"**Selected chemotherapy compounds:** {', '.join(selected_chemo)}")
st.write(f"**Selected natural products:** {', '.join(selected_natural)}")

# Build overlap table
overlap_records = []
for pathway_name, pathway_genes in cancer_genes_map.items():
    for compound, comp_genes in compound_genes_map.items():
        shared = pathway_genes.intersection(comp_genes)
        if shared:
            overlap_records.append({
                "Cancer Pathway": pathway_name,
                "Compound": compound,
                "Type": "Chemotherapy" if compound in selected_chemo else "Natural Product",
                "Shared Gene Count": len(shared),
                "Shared Genes": ", ".join(sorted(shared))
            })

if overlap_records:
    df_overlap = pd.DataFrame(overlap_records)
    st.subheader("Gene Overlap Table")
    st.dataframe(df_overlap)

    csv_data = df_overlap.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download Overlap Table as CSV",
        data=csv_data,
        file_name="cancer_compound_gene_overlap.csv",
        mime="text/csv",
    )
else:
    st.warning("No overlapping genes found between selected cancer pathways and compounds.")

# Legend
# Corrected Node color legend
# Corrected Node color legend (all circles)
st.markdown("""
### 🗺️ Node Color Legend

<style>
.color-circle {
    display: inline-block;
    width: 15px;
    height: 15px;
    border-radius: 50%;
    margin-right: 10px;
    vertical-align: middle;
}
.legend-item {
    margin-bottom: 6px;
    font-size: 16px;
}
</style>

<div class="legend-item"><span class="color-circle" style="background-color:#800080;"></span> <b>Cancer Pathway</b></div>
<div class="legend-item"><span class="color-circle" style="background-color:#1f77b4;"></span> Gene (in pathway only)</div>
<div class="legend-item"><span class="color-circle" style="background-color:#ff69b4;"></span> Overlapping Gene (in both pathway & compound)</div>
<div class="legend-item"><span class="color-circle" style="background-color:#ff0000;"></span> Gene (only in compound, not in pathway)</div>
<div class="legend-item"><span class="color-circle" style="background-color:#f39c12;"></span> Chemotherapy Compound</div>
<div class="legend-item"><span class="color-circle" style="background-color:#27ae60;"></span> Natural Product</div>
""", unsafe_allow_html=True)


# Build network graph
st.subheader("Gene Interaction Network Visualization")


G = nx.Graph()

# Add cancer pathways nodes
for p in selected_pathways:
    G.add_node(p, type="pathway", color="purple", size=25)

# Add compound nodes
for c in selected_compounds:
    node_color = "orange" if c in selected_chemo else "green"
    G.add_node(c, type="compound", color=node_color, size=20)

# Add gene nodes and edges
all_genes = set()
for pathway_name, pathway_genes in cancer_genes_map.items():
    for gene in pathway_genes:
        G.add_node(gene, type="gene", color="blue", size=10)
        G.add_edge(pathway_name, gene)
        all_genes.add(gene)

for compound, comp_genes in compound_genes_map.items():
    for gene in comp_genes:
        # If gene exists, update color if overlaps
        if gene in all_genes:
            G.nodes[gene]['color'] = "magenta"  # overlapping gene
            G.nodes[gene]['size'] = 15
        else:
            G.add_node(gene, type="gene", color="red", size=10)
        G.add_edge(compound, gene)

# Network layout and plotting
# Network layout and plotting
pos = nx.spring_layout(G, k=0.25, iterations=50, seed=42)

edge_x = []
edge_y = []
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_x += [x0, x1, None]
    edge_y += [y0, y1, None]

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line=dict(width=0.5, color="#888"),
    hoverinfo='none',
    mode='lines'
)

node_x = []
node_y = []
node_text = []
node_color = []
node_size = []

for node in G.nodes():
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    node_text.append(node)  # Show name on hover only
    node_color.append(G.nodes[node]['color'])
    node_size.append(G.nodes[node]['size'])

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    text=node_text,  # tooltip on hover
    marker=dict(
        showscale=False,
        color=node_color,
        size=node_size,
        line_width=1
    )
)

fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title="Cancer Pathways, Compounds, and Gene Network",
                    title_x=0.5,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                )
)

st.plotly_chart(fig, use_container_width=True)
