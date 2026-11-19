"""Load the human PPI network from the BioGRID tab3 file."""

import zipfile

import networkx as nx
import pandas as pd
import requests

from cmn.paths import DATA

RELEASE = "5.0.262"
BIOGRID_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
    f"BIOGRID-{RELEASE}/BIOGRID-ORGANISM-{RELEASE}.tab3.zip"
)
BIOGRID_FILE = DATA / f"BIOGRID-ORGANISM-Homo_sapiens-{RELEASE}.tab3.txt"
EDGES_FILE = DATA / "human_ppi_edges.tsv"
SUBGRAPH_FILE = DATA / "bglap_subgraph_edges.tsv"
SEED = 632  # BGLAP (Osteocalcin), Entrez Gene ID
SYSTEMS = ("Two-hybrid", "Affinity Capture-MS")
HUMAN = 9606

COLUMNS = {
    "Entrez Gene Interactor A": "a",
    "Entrez Gene Interactor B": "b",
    "Official Symbol Interactor A": "symbol_a",
    "Official Symbol Interactor B": "symbol_b",
    "Experimental System": "system",
    "Experimental System Type": "system_type",
    "Organism ID Interactor A": "org_a",
    "Organism ID Interactor B": "org_b",
}


def download_biogrid(path=BIOGRID_FILE):
    """Download the BioGRID organism archive and extract only the human file (skipped if it already exists)."""
    if path.exists():
        return path
    DATA.mkdir(parents=True, exist_ok=True)
    archive = DATA / BIOGRID_URL.rsplit("/", 1)[1]
    with requests.get(BIOGRID_URL, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(archive, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    with zipfile.ZipFile(archive) as z:
        z.extract(path.name, DATA)
    archive.unlink()
    return path


def read_interactions(path=BIOGRID_FILE, systems=SYSTEMS):
    """Physical human-human interactions from the given experimental systems, one row per BioGRID record."""
    download_biogrid(path)
    df = pd.read_csv(path, sep="\t", usecols=list(COLUMNS), dtype=str).rename(columns=COLUMNS)
    df = df[
        df["system"].isin(systems)
        & (df["system_type"] == "physical")
        & (df["org_a"] == str(HUMAN))
        & (df["org_b"] == str(HUMAN))
    ]
    df = df[df["a"].str.isdigit() & df["b"].str.isdigit()]
    return df.astype({"a": int, "b": int})


def unique_edges(df):
    """One row per undirected edge (self-loops removed) with the set of systems supporting it."""
    df = df[df["a"] != df["b"]].copy()
    df["u"] = df[["a", "b"]].min(axis=1)
    df["v"] = df[["a", "b"]].max(axis=1)
    edges = df.groupby(["u", "v"])["system"].agg(lambda s: ";".join(sorted(set(s)))).reset_index()
    return edges.rename(columns={"system": "systems"})


def symbols(df):
    """Entrez Gene ID -> official symbol."""
    a = df[["a", "symbol_a"]].set_axis(["id", "symbol"], axis=1)
    b = df[["b", "symbol_b"]].set_axis(["id", "symbol"], axis=1)
    return pd.concat([a, b]).drop_duplicates("id").set_index("id")["symbol"].to_dict()


def build_graph(edges, symbol_of):
    """Undirected graph on Entrez Gene IDs, with node attributes entrez/symbol and edge attribute systems."""
    G = nx.Graph()
    for u, v, s in edges[["u", "v", "systems"]].itertuples(index=False):
        G.add_edge(u, v, systems=s)
    nx.set_node_attributes(G, {n: n for n in G}, "entrez")
    nx.set_node_attributes(G, {n: symbol_of.get(n, str(n)) for n in G}, "symbol")
    return G


def save_edges(edges, symbol_of, path=EDGES_FILE):
    """Write the processed edge list (Entrez IDs, symbols, systems) for reuse in Part B."""
    out = edges.assign(symbol_u=edges["u"].map(symbol_of), symbol_v=edges["v"].map(symbol_of))
    out[["u", "v", "symbol_u", "symbol_v", "systems"]].to_csv(path, sep="\t", index=False)


def save_graph(G, path):
    """Write a graph built by build_graph as an edge list in the same format as save_edges."""
    rows = [(u, v, G.nodes[u]["symbol"], G.nodes[v]["symbol"], s) for u, v, s in G.edges(data="systems")]
    pd.DataFrame(rows, columns=["u", "v", "symbol_u", "symbol_v", "systems"]).to_csv(path, sep="\t", index=False)


def load_graph(path=EDGES_FILE):
    """Load the processed human PPI network written by save_edges."""
    edges = pd.read_csv(path, sep="\t")
    symbol_of = dict(zip(edges["u"], edges["symbol_u"])) | dict(zip(edges["v"], edges["symbol_v"]))
    return build_graph(edges, symbol_of)
