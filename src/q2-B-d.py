"""Q2 Part B (d) - graphlet orbit counts (orbits 0-14) and graphlet counts of the Osteocalcin induced subgraph."""

import networkx as nx
import numpy as np

from cmn.biogrid import SEED, SUBGRAPH_FILE, load_graph
from cmn.graphlets import GRAPHLET_FROM_ORBIT, GRAPHLET_NAMES, graphlet_counts, orbit_counts
from cmn.paths import DATA, save_fig
from cmn.plots import INK, SERIES, axes

ORBITS_FILE = DATA / "bglap_subgraph_orbits.tsv"

# Every other orbit of a graphlet and the number of its nodes in that orbit: a second, independent way to count it.
CROSS_CHECKS = {"G1": (1, 2), "G3": (4, 2), "G4": (6, 3), "G6": (9, 1), "G6'": (10, 2), "G7": (12, 2)}


def main():
    H = load_graph(SUBGRAPH_FILE)
    nodes = sorted(H)
    O = orbit_counts(H, nodes)
    O.insert(0, "symbol", [H.nodes[v]["symbol"] for v in nodes])
    O.index.name = "entrez"
    O.to_csv(ORBITS_FILE, sep="\t")
    O = O.drop(columns="symbol")
    print(f"orbit counts: {O.shape[0]} nodes x {O.shape[1]} orbits, saved to {ORBITS_FILE}")

    print("\norbit totals over all nodes")
    totals = O.sum()
    for i, t in totals.items():
        print(f"  orbit {i:2d}: {t}")

    counts = graphlet_counts(O)
    print("\ngraphlet counts (orbit sum / nodes per orbit)")
    for g, c in counts.items():
        orbit, k = GRAPHLET_FROM_ORBIT[g]
        print(f"  {g} {GRAPHLET_NAMES[g]:20s} orbit {orbit:2d}: {totals[orbit]:>14d} / {k} = {c:>12d}")

    print("\nchecks")
    m = H.number_of_edges()
    print(f"  orbit 0 / 2 = {totals[0] // 2}  vs  edges = {m}  ->  {totals[0] // 2 == m}")
    tri = sum(nx.triangles(H).values()) // 3
    print(f"  orbit 3 / 3 = {counts['G2']}  vs  networkx triangles = {tri}  ->  {counts['G2'] == tri}")
    for g, (orbit, k) in CROSS_CHECKS.items():
        print(f"  {g[:2]} via orbit {orbit:2d}: {totals[orbit]} / {k} = {totals[orbit] // k}  "
              f"->  {totals[orbit] // k == counts[g[:2]]}")
    paths3 = sum(d * (d - 1) // 2 for _, d in H.degree()) - 3 * tri
    print(f"  G1 from degrees: sum C(k,2) - 3 triangles = {paths3}  ->  {paths3 == counts['G1']}")

    print(f"\nBGLAP orbit signature: {O.loc[SEED].tolist()}")

    centre = O[7].sort_values(ascending=False)
    print("\nlargest centres of 4-node stars (orbit 7)")
    for v in centre.index[:5]:
        print(f"  {H.nodes[v]['symbol']:8s} degree {H.degree(v):5d}  stars {centre[v]:>11d}  "
              f"({centre[v] / counts['G4']:.3f} of all G4)")
    print(f"  top 10 centres hold {centre.iloc[:10].sum() / counts['G4']:.3f} of all G4")

    names = list(counts)
    fig, ax = axes("", "number of graphlets (log scale)", figsize=(8, 3.8))
    ax.bar(range(len(names)), [counts[g] for g in names], width=0.6, color=SERIES[0], edgecolor="white", linewidth=1)
    for i, g in enumerate(names):
        ax.annotate(f"{counts[g]:.2e}".replace("e+0", "e"), (i, counts[g]), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=7.5, color=INK)
    ax.set_yscale("log")
    ax.set_ylim(top=max(counts.values()) * 8)
    ax.set_xticks(range(len(names)), [f"{g}\n" + GRAPHLET_NAMES[g].replace(" + ", "\n+ ").replace("node ", "node\n")
                                      for g in names], fontsize=7.5)
    ax.grid(axis="x", visible=False)
    save_fig(fig, "q2_B_d_graphlets.png")


if __name__ == "__main__":
    main()
