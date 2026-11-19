"""Q2 Part B (c) - clustering of the Osteocalcin subgraph compared with ER and degree-preserving (ER-DD) random graphs."""

import networkx as nx
import numpy as np
import pandas as pd

from cmn.biogrid import SUBGRAPH_FILE, load_graph
from cmn.paths import DATA, save_fig
from cmn.plots import INK, MUTED, SERIES, axes

RUNS = 10
SWAPS_PER_EDGE = 3
RESULTS_FILE = DATA / "q2_B_c_random_clustering.tsv"


def er_graphs(n, m, runs):
    """G(n, m) Erdos-Renyi graphs with the same number of nodes and edges."""
    return [nx.gnm_random_graph(n, m, seed=r) for r in range(runs)]


def erdd_graphs(H, runs, swaps_per_edge):
    """Degree-preserving randomisations by repeated double edge swaps."""
    out = []
    m = H.number_of_edges()
    for r in range(runs):
        R = nx.Graph(H.edges())
        nx.double_edge_swap(R, nswap=swaps_per_edge * m, max_tries=100 * m, seed=r)
        out.append(R)
    return out


def summary(name, graphs):
    rows = [(name, r, nx.average_clustering(R), nx.transitivity(R)) for r, R in enumerate(graphs)]
    return pd.DataFrame(rows, columns=["model", "run", "avg_clustering", "transitivity"])


def main():
    H = load_graph(SUBGRAPH_FILE)
    n, m = H.number_of_nodes(), H.number_of_edges()
    C_obs, T_obs = nx.average_clustering(H), nx.transitivity(H)
    p = 2 * m / (n * (n - 1))
    print(f"induced subgraph: n = {n}, m = {m}")
    print(f"I. observed average clustering = {C_obs:.4f}, transitivity = {T_obs:.4f}")

    print("\nII. ER model")
    print(f"  expected C = p = edge density = {p:.4f}")
    er = summary("ER", er_graphs(n, m, RUNS))

    print("  ER-DD model: double_edge_swap, "
          f"{SWAPS_PER_EDGE} swaps per edge = {SWAPS_PER_EDGE * m} swaps, {RUNS} runs")
    erdd = summary("ER-DD", erdd_graphs(H, RUNS, SWAPS_PER_EDGE))
    df = pd.concat([er, erdd])
    df.to_csv(RESULTS_FILE, sep="\t", index=False)
    print(df.to_string(index=False))

    print("\n  model     mean C    sd C      ratio   z(C)      mean T    sd T      ratio   z(T)")
    stats = df.groupby("model").agg(["mean", "std"])
    for model in ("ER", "ER-DD"):
        mc, sc = stats.loc[model, ("avg_clustering", "mean")], stats.loc[model, ("avg_clustering", "std")]
        mt, st = stats.loc[model, ("transitivity", "mean")], stats.loc[model, ("transitivity", "std")]
        print(f"  {model:8s} {mc:.4f}   {sc:.5f}   {C_obs / mc:5.2f}   {(C_obs - mc) / sc:7.1f}   "
              f"{mt:.4f}   {st:.5f}   {T_obs / mt:5.2f}   {(T_obs - mt) / st:7.1f}")

    k = np.array([d for _, d in H.degree()], dtype=float)
    k1, k2 = k.mean(), (k ** 2).mean()
    print("\n  configuration-model expectation (k2 - k1)^2 / (n k1^3) "
          f"= {(k2 - k1) ** 2 / (n * k1 ** 3):.4f}  (<k> = {k1:.1f}, <k^2> = {k2:.0f}, k_max/n = {k.max() / n:.3f})")
    hubs = sorted(H, key=H.degree, reverse=True)[:10]
    covered = set(hubs).union(*(H[h] for h in hubs))
    print(f"  nodes that are, or are adjacent to, one of the 10 largest hubs: {len(covered)} ({len(covered) / n:.3f})")

    print("\n  average clustering by degree class, observed vs ER-DD run 0")
    R = erdd_graphs(H, 1, SWAPS_PER_EDGE)[0]
    deg = pd.Series(dict(H.degree()))
    obs = pd.Series(nx.clustering(H))
    rnd = pd.Series(nx.clustering(R))
    bins = [1, 2, 10, 50, 200, 2000]
    cls = pd.cut(deg, bins, right=False)
    table = pd.DataFrame({"nodes": deg.groupby(cls, observed=True).size(),
                          "C observed": obs.groupby(cls, observed=True).mean(),
                          "C ER-DD": rnd.reindex(deg.index).groupby(cls, observed=True).mean()})
    print(table.to_string(float_format=lambda x: f"{x:.3f}"))

    fig, ax = axes("", "average clustering coefficient", figsize=(5, 3.4))
    labels = ["observed", "ER-DD\n(degree-preserving)", "ER"]
    means = [C_obs, stats.loc["ER-DD", ("avg_clustering", "mean")], stats.loc["ER", ("avg_clustering", "mean")]]
    ax.bar(range(3), means, width=0.48, color=[SERIES[0], SERIES[1], SERIES[2]], edgecolor="white", linewidth=1)
    tops = list(means)
    for i, model in ((1, "ER-DD"), (2, "ER")):
        vals = df.loc[df["model"] == model, "avg_clustering"]
        ax.scatter(np.full(len(vals), i) + np.linspace(-0.12, 0.12, len(vals)), vals, s=10, color=INK, zorder=3)
        tops[i] = max(tops[i], vals.max())
    for i, (v, t) in enumerate(zip(means, tops)):
        ax.annotate(f"{v:.3f}", (i, t), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=9, color=INK)
    ax.axhline(p, color=MUTED, lw=0.8, ls="--")
    ax.annotate(f"$p$ = {p:.4f}", (1.5, p), xytext=(0, 4), textcoords="offset points",
                ha="center", fontsize=8, color=MUTED)
    ax.set_xticks(range(3), labels)
    ax.set_ylim(0, C_obs * 1.18)
    ax.grid(axis="x", visible=False)
    save_fig(fig, "q2_B_c_clustering.png")
    print(f"\nresults saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
