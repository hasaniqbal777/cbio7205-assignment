"""Q2 Part B (b) - network properties of the Osteocalcin induced subgraph."""

from collections import Counter

import networkx as nx
import numpy as np
import pandas as pd

from cmn.biogrid import SEED, SUBGRAPH_FILE, load_graph
from cmn.paths import DATA, save_fig
from cmn.plots import INK, SERIES, axes, bars

NODES_FILE = DATA / "bglap_subgraph_nodes.tsv"


def top(series, symbol, k=5):
    return ", ".join(f"{symbol[v]} ({series[v]:.4g})" for v in series.sort_values(ascending=False).index[:k])


def main():
    H = load_graph(SUBGRAPH_FILE)
    symbol = nx.get_node_attributes(H, "symbol")
    n, m = H.number_of_nodes(), H.number_of_edges()
    print(f"induced subgraph: n = {n}, m = {m}")

    deg = pd.Series(dict(H.degree()))
    print("\nI. Degree distribution")
    print(f"  mean = {deg.mean():.2f}, median = {deg.median():.0f}, min = {deg.min()}, max = {deg.max()}")
    print(f"  degree-1 nodes = {(deg == 1).sum()} ({(deg == 1).mean():.3f})")
    print(f"  BGLAP degree = {deg[SEED]}")
    print(f"  top: {top(deg, symbol, 10)}")
    pk = deg.value_counts().sort_index() / n

    fig, ax = axes("degree $k$", "$P(k)$")
    ax.scatter(pk.index, pk.values, s=12, color=SERIES[0], edgecolors="white", linewidths=0.4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    save_fig(fig, "q2_B_b_degree.png")

    print("\nII. Average clustering coefficient")
    C = pd.Series(nx.clustering(H))
    print(f"  average C = {C.mean():.4f}")
    print(f"  average C over nodes with degree >= 2 = {C[deg >= 2].mean():.4f}")
    print(f"  global transitivity = {nx.transitivity(H):.4f}")
    print(f"  C(BGLAP) = {C[SEED]:.4f}")

    print("\nIII. Clustering spectrum C(k)")
    Ck = C.groupby(deg).mean()
    for k in (1, 2, 5, 10, 20, 50, 100, 200, 500):
        near = Ck[Ck.index <= k]
        if len(near):
            print(f"  C({near.index[-1]}) = {near.iloc[-1]:.4f}")
    slope = np.polyfit(np.log(Ck[Ck > 0].index), np.log(Ck[Ck > 0].values), 1)[0]
    print(f"  log-log slope of C(k) (k >= 2): {slope:.3f}")

    fig, ax = axes("degree $k$", "$C(k)$")
    ax.scatter(Ck.index[Ck.index >= 2], Ck[Ck.index >= 2], s=12, color=SERIES[0],
               edgecolors="white", linewidths=0.4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    save_fig(fig, "q2_B_b_clustering_spectrum.png")

    print("\nIV/V. Shortest paths")
    lengths = Counter()
    for _, d in nx.all_pairs_shortest_path_length(H):
        lengths.update(d.values())
    del lengths[0]
    pairs = sum(lengths.values()) // 2
    sp = pd.Series({d: c // 2 for d, c in sorted(lengths.items())})
    print(f"  pairs = {pairs}")
    for d, c in sp.items():
        print(f"  d = {d}: {c} pairs ({c / pairs:.4f})")
    avg = sum(d * c for d, c in sp.items()) / pairs
    print(f"  diameter = {sp.index.max()}, average shortest path length = {avg:.4f}")

    fig, ax = axes("shortest path length $d$", "fraction of node pairs")
    bars(ax, sp.index, sp.values / pairs)
    ax.set_xticks(sp.index)
    ax.set_ylim(0, (sp.values / pairs).max() * 1.15)
    save_fig(fig, "q2_B_b_shortest_paths.png")

    print("\nVI. Betweenness centrality")
    B = pd.Series(nx.betweenness_centrality(H, normalized=True))
    print(f"  mean = {B.mean():.3e}, median = {B.median():.3e}, max = {B.max():.4f}")
    print(f"  nodes with B = 0: {(B == 0).sum()} ({(B == 0).mean():.3f})")
    print(f"  B(BGLAP) = {B[SEED]:.3e}")
    print(f"  top: {top(B, symbol, 10)}")
    print(f"  share of total betweenness held by top node: {B.max() / B.sum():.3f}")

    fig, ax = axes("normalised betweenness (log scale)", "number of nodes")
    pos = B[B > 0]
    bins = np.logspace(np.log10(pos.min()), np.log10(pos.max()), 40)
    ax.hist(pos, bins=bins, color=SERIES[0], edgecolor="white", linewidth=0.6)
    ax.set_xscale("log")
    hub = B.idxmax()
    ax.annotate(symbol[hub], (B[hub], 1), xytext=(0, 14), textcoords="offset points", ha="center",
                fontsize=8, color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
    save_fig(fig, "q2_B_b_betweenness.png")

    print("\nVII. Eccentricity")
    E = pd.Series(nx.eccentricity(H))
    ec = E.value_counts().sort_index()
    for e, c in ec.items():
        print(f"  ecc = {e}: {c} nodes ({c / n:.4f})")
    print(f"  radius = {E.min()}, diameter = {E.max()}, E(BGLAP) = {E[SEED]}")
    print(f"  centre: {', '.join(symbol[v] for v in E[E == E.min()].index)}")

    fig, ax = axes("eccentricity", "fraction of nodes")
    bars(ax, ec.index, ec.values / n, fmt="{:.4f}")
    ax.set_xticks(ec.index)
    ax.set_ylim(0, (ec.values / n).max() * 1.15)
    save_fig(fig, "q2_B_b_eccentricity.png")

    nodes = pd.DataFrame({"symbol": pd.Series(symbol), "degree": deg, "clustering": C,
                          "betweenness": B, "eccentricity": E})
    nodes.index.name = "entrez"
    nodes.sort_values("degree", ascending=False).to_csv(NODES_FILE, sep="\t")
    print(f"\nper-node metrics saved to {NODES_FILE}")


if __name__ == "__main__":
    main()
