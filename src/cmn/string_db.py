"""Download STRING protein.links files and turn them into thresholded, deduplicated edge lists."""

import gzip

import pandas as pd
import requests

from cmn.paths import DATA

VERSION = "12.5"
STRING_DIR = DATA / "string"
SPECIES = {
    "human": (9606, "Homo sapiens"),
    "yeast": (4932, "Saccharomyces cerevisiae"),
    "worm": (6239, "Caenorhabditis elegans"),
    "fly": (7227, "Drosophila melanogaster"),
    "cress": (3702, "Arabidopsis thaliana"),
}
MIN_SCORE = 700  # combined_score is stored as an integer, 0.7 * 1000


def links_url(taxid):
    return (f"https://stringdb-downloads.org/download/protein.links.v{VERSION}/"
            f"{taxid}.protein.links.v{VERSION}.txt.gz")


def download_links(taxid):
    """Download the protein.links file of a species (skipped if it already exists)."""
    STRING_DIR.mkdir(parents=True, exist_ok=True)
    path = STRING_DIR / f"{taxid}.protein.links.v{VERSION}.txt.gz"
    if not path.exists():
        with requests.get(links_url(taxid), stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
    return path


def edge_list(taxid, min_score=MIN_SCORE):
    """Undirected edges with combined_score >= min_score, each pair once, no self-loops."""
    with gzip.open(download_links(taxid), "rt") as f:
        df = pd.read_csv(f, sep=" ")
    df = df[(df["combined_score"] >= min_score) & (df["protein1"] != df["protein2"])]
    u = df[["protein1", "protein2"]].min(axis=1)
    v = df[["protein1", "protein2"]].max(axis=1)
    return pd.DataFrame({"u": u, "v": v}).drop_duplicates().reset_index(drop=True)


def edges_file(name):
    return STRING_DIR / f"{name}.txt"


def write_edges(edges, name):
    """Edge-list file with one "node1 node2" pair per line (the format MAGNA++ reads)."""
    path = edges_file(name)
    edges.to_csv(path, sep=" ", header=False, index=False)
    return path
