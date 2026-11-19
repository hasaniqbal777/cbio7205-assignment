"""Q3 Part B (c) - UPGMA phenetic tree from MAGNA++ S3 scores, and the NCBI Taxonomy tree of the five species."""

from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

from cmn.apis import ncbi_lineages
from cmn.paths import DATA, save_fig
from cmn.plots import INK, MUTED, SERIES
from cmn.string_db import SPECIES

ALIGNMENTS_FILE = DATA / "q3_B_b_alignments.tsv"
LABEL = {name: f"{latin} ({name})" for name, (_, latin) in SPECIES.items()}


def newick(Z, names):
    """Newick string of a SciPy linkage, with branch lengths (UPGMA heights)."""
    def rec(node, parent_h):
        h = node.dist / 2
        if node.is_leaf():
            return f"{names[node.id]}:{parent_h - h:.4f}"
        return f"({rec(node.left, h)},{rec(node.right, h)}):{parent_h - h:.4f}"
    root = to_tree(Z)
    return f"({rec(root.left, root.dist / 2)},{rec(root.right, root.dist / 2)});"


def draw(Z, names, xlabel, internal=None, xlim=None, hide_x=False):
    """Horizontal dendrogram; `internal` optionally labels merges (by linkage row) with a name."""
    fig, ax = plt.subplots(figsize=(6.5, 2.8))
    d = dendrogram(Z, labels=[LABEL[n] for n in names], orientation="left", ax=ax,
                   color_threshold=0, above_threshold_color=SERIES[0], leaf_font_size=9)
    for coll in ax.collections:
        coll.set_linewidth(2)
    if internal:
        # each U-shaped link is one merge; match it to its linkage row by height (distinct here)
        for x, y in zip(d["dcoord"], d["icoord"]):
            row = int(np.argmin(np.abs(Z[:, 2] - x[1])))
            ax.annotate(internal[row], (x[1], (y[1] + y[2]) / 2), xytext=(-4, 3), textcoords="offset points",
                        ha="right", va="bottom", fontsize=7.5, color=MUTED)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="x", colors=MUTED, labelcolor=INK, labelsize=8)
    ax.tick_params(axis="y", length=0, labelcolor=INK)
    ax.set_xlabel(xlabel, color=INK, fontsize=9)
    if xlim:
        ax.set_xlim(*xlim)
    if hide_x:
        ax.xaxis.set_visible(False)
        ax.spines["bottom"].set_visible(False)
    return fig


def main():
    names = list(SPECIES)
    al = pd.read_csv(ALIGNMENTS_FILE, sep="\t")
    S3 = pd.DataFrame(1.0, index=names, columns=names)
    for r in al.itertuples(index=False):
        S3.loc[r[0], r[1]] = S3.loc[r[1], r[0]] = r.S3
    D = 1 - S3
    print("I. dissimilarity 1 - S3")
    print(D.round(4).to_string())

    Z = linkage(squareform(D.values, checks=False), method="average")
    print("\n  UPGMA merges (cluster, cluster, height, size)")
    clusters = {i: [n] for i, n in enumerate(names)}
    for i, (a, b, h, k) in enumerate(Z):
        clusters[len(names) + i] = clusters[int(a)] + clusters[int(b)]
        print(f"    {clusters[int(a)]} + {clusters[int(b)]} at 1 - S3 = {h:.4f} (S3 = {1 - h:.4f})")
    phenetic = newick(Z, names)
    print(f"  Newick: {phenetic}")
    save_fig(draw(Z, names, "dissimilarity $1 - S^3$ (UPGMA, axis starts at 0.85)", xlim=(Z[:, 2].max() + 0.01, 0.85)),
             "q3_B_c_phenetic_tree.png")

    print("\nII. NCBI Taxonomy")
    lineage = ncbi_lineages([SPECIES[n][0] for n in names])
    lin = {n: lineage[SPECIES[n][0]] for n in names}
    for n in names:
        print(f"  {n}: {'; '.join(lin[n][:8])} ...")
    shared = pd.DataFrame(0, index=names, columns=names)
    lca = {}
    for a, b in combinations(names, 2):
        k = next((i for i, (x, y) in enumerate(zip(lin[a], lin[b])) if x != y), min(len(lin[a]), len(lin[b])))
        shared.loc[a, b] = shared.loc[b, a] = k
        lca[frozenset((a, b))] = lin[a][k - 1]
        print(f"  LCA({a}, {b}) = {lin[a][k - 1]} ({k} shared taxa)")
    depth = max(len(l) for l in lin.values())
    T = depth - shared.to_numpy(dtype=float)
    np.fill_diagonal(T, 0)
    Zt = linkage(squareform(T, checks=False), method="average")
    clusters = {i: [n] for i, n in enumerate(names)}
    internal = {}
    for i, (a, b, h, k) in enumerate(Zt):
        clusters[len(names) + i] = clusters[int(a)] + clusters[int(b)]
        internal[i] = lca[frozenset((clusters[int(a)][0], clusters[int(b)][0]))]
        print(f"    {clusters[int(a)]} + {clusters[int(b)]}  ({internal[i]})")
    print(f"  Newick (topology): {newick(Zt, names)}")
    save_fig(draw(Zt, names, "taxonomic levels below the common ancestor", internal, xlim=(Zt[:, 2].max() + 6, 25), hide_x=True),
             "q3_B_c_taxonomy_tree.png")


if __name__ == "__main__":
    main()
