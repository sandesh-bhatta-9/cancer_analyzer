import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from collections import Counter
from bioservices import KEGG
import base64
from pathlib import Path

# --- Page Configuration ---
st.set_page_config(page_title="Cancer Pathway Analyzer", layout="wide")
st.title("🔬 Cancer-Specific Pathway Analyzer")
st.write("Compare cancer pathways against chemotherapy and natural product pathways.")

# --- Icon Utilities ---
ICON_NAMES = ["Gene Cards", "NCBI", "ENSEMBL", "KEGG", "GEO"]
ICONS_DIR = Path(__file__).parent / "Icons"

@st.cache_data
def load_icons_b64():
    icon_map = {}
    for name in ICON_NAMES:
        fp = ICONS_DIR / f"{name}.png"
        if fp.is_file():
            with open(fp, "rb") as f:
                b64 = base64.b64encode(f.read()).decode("utf-8")
            icon_map[name] = f"data:image/png;base64,{b64}"
        else:
            icon_map[name] = None
    return icon_map

ICON_B64 = load_icons_b64()

# --- KEGG Helpers ---
@st.cache_data
def get_all_kegg_pathways():
    k = KEGG(); pathways_raw = k.list("pathway/hsa")
    d = {};
    for line in pathways_raw.strip().split('\n'):
        pid, name_desc = line.split('\t'); name = name_desc.split(' - ')[0]
        d[name] = pid.replace('path:', '')
    return d

@st.cache_data
def get_genes(pathway_id: str):
    k = KEGG(); k.organism = "hsa"; data = k.get(pathway_id) or ""
    genes = set(); rec = False
    for line in data.split('\n'):
        if line.startswith('GENE'): rec = True; line_cont = line[12:].strip()
        elif rec and line.startswith(' '): line_cont = line.strip()
        else: rec = False; continue
        if line_cont:
            parts = line_cont.split()
            if parts and parts[0].isdigit(): genes.add(parts[1])
    return genes

def generate_icon_links(gene_name: str) -> dict:
    urls = {
        "Gene Cards": f"https://www.genecards.org/Search/Keyword?queryString={gene_name}",
        "NCBI": f"https://www.ncbi.nlm.nih.gov/gene/?term={gene_name}",
        "ENSEMBL": f"https://useast.ensembl.org/Search/Results?q={gene_name}",
        "GEO": f"https://www.ncbi.nlm.nih.gov/gds/?term={gene_name}",
    }
    html = {}
    for db, url in urls.items():
        b64src = ICON_B64.get(db)
        if b64src:
            html[db] = (f'<a href="{url}" target="_blank"><img src="{b64src}" width="24" height="24" style="margin:2px; border-radius:4px; border:1px solid #ccc;"></a>')
        else:
            html[db] = f'<a href="{url}" target="_blank">{db}</a>'
    return html

# --- Sidebar ---
all_paths = get_all_kegg_pathways()
st.sidebar.header("Select Pathways for Comparison")

CANCER_KEYWORDS = ['cancer', 'glioma', 'leukemia', 'lymphoma', 'melanoma']
CHEMO_KEYWORDS = ['fluoropyrimidine', 'folate', 'platinum']
NATURAL_PRODUCT_KEYWORDS = ['flavonoid', 'terpenoid', 'stilbenoid']

cancer_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in CANCER_KEYWORDS)}
chemo_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in CHEMO_KEYWORDS)}
natural_options = {n:pid for n,pid in all_paths.items() if any(k in n.lower() for k in NATURAL_PRODUCT_KEYWORDS)}

st.sidebar.subheader("Cancer Pathways")
sel_cancer = st.sidebar.multiselect("Select cancer types:", list(cancer_options.keys()), default=list(cancer_options.keys())[:3])

st.sidebar.subheader("Therapeutic Pathways")
sel_chemo = st.sidebar.multiselect("Select chemotherapy pathways:", list(chemo_options.keys()), default=list(chemo_options.keys()))
sel_natural = st.sidebar.multiselect("Select natural product pathways:", list(natural_options.keys()), default=list(natural_options.keys()))

combined_options = {**cancer_options, **chemo_options, **natural_options}
selected_pathways = sel_cancer + sel_chemo + sel_natural

# --- Main App Logic ---
if selected_pathways:
    with st.spinner("Fetching and analyzing gene data..."):
        p2g = {n: get_genes(combined_options[n]) for n in selected_pathways}

    st.header("📊 Shared Gene Analysis")
    flat_genes = [g for genes in p2g.values() for g in genes]
    gene_counts = Counter(flat_genes)
    genes_by_freq = {cnt: [g for g, c in gene_counts.items() if c == cnt] for cnt in range(2, len(selected_pathways) + 1)}

    if any(genes_by_freq.values()):
        for freq in sorted(genes_by_freq.keys(), reverse=True):
            gene_list = genes_by_freq[freq]
            if not gene_list: continue
            
            st.subheader(f"Genes in {freq} pathways")
            page_data = [{'Gene': g, 'Pathways': ', '.join([n for n, gs in p2g.items() if g in gs]), **generate_icon_links(g)} for g in gene_list]
            df_page = pd.DataFrame(page_data).sort_values(by="Gene")
            
            if len(df_page) <= 20:
                st.write(f"Showing all {len(df_page)} shared genes:")
                display_df = df_page
            else:
                st.write(f"Showing a preview of the top 20 of {len(df_page)} total genes:")
                display_df = df_page.head(20)
            
            cols = ['Gene','Pathways','Gene Cards','NCBI','ENSEMBL','GEO']
            st.write(display_df[cols].to_html(escape=False,index=False), unsafe_allow_html=True)

            csv_string = df_page.to_csv(index=False).encode('utf-8')
            st.download_button(label=f"Download Full List ({len(df_page)} genes)", data=csv_string, file_name=f"shared_genes_{freq}_pathways.csv", mime="text/csv", key=f"dl_gene_{freq}")
            st.markdown("---")
    else:
        st.write("⚠️ No shared genes found.")

    st.header("🕸️ Pathway Network")
    G = nx.Graph()
    for pathway in p2g: G.add_node(pathway, type="pathway")
    for gene, cnt in gene_counts.items(): G.add_node(gene, type="gene", count=cnt)
    for pathway, genes in p2g.items():
        for gene in genes: G.add_edge(pathway, gene)
        
    pos = nx.spring_layout(G, k=0.5, seed=42)
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]; x1, y1 = pos[v]
        edge_x += [x0, x1, None]; edge_y += [y0, y1, None]
    edge_trace = go.Scatter(x=edge_x, y=edge_y, mode='lines', line=dict(color='#888', width=0.5), hoverinfo='none')

    node_x, node_y, text, color = [], [], [], []
    for node, attr in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x); node_y.append(y)
        if attr["type"] == "pathway":
            if node in sel_cancer: color.append("crimson")
            elif node in sel_chemo: color.append("mediumblue")
            elif node in sel_natural: color.append("forestgreen")
            else: color.append("grey")
            text.append(node)
        else:
            count = attr.get("count", 1)
            text.append(f"{node} (in {count} pathways)")
            color.append("orange" if count > 1 else "lightgray")
            
    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers', marker=dict(size=10, color=color, line_width=2), text=text, hoverinfo='text')
    fig = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(title="Cancer, Chemotherapy, and Natural Product Pathway Network", xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), margin=dict(b=20, l=5, r=5, t=40), hovermode='closest'))
    st.plotly_chart(fig, use_container_width=True)

else:
    st.info("☝️ Please select at least one pathway from the sidebar to begin.")