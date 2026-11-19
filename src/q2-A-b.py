"""Q2 Part A (b) - size, edge density and edge sources of the human Y2H + AP-MS PPI network from BioGRID."""

import networkx as nx

from cmn.biogrid import SYSTEMS, build_graph, read_interactions, save_edges, symbols, unique_edges


def main():
    df = read_interactions()
    print("(b) Human PPI network (Two-hybrid + Affinity Capture-MS, both interactors human)")
    print(f"  BioGRID records kept: {len(df)}")
    for system, n in df["system"].value_counts().items():
        print(f"    {system}: {n}")

    n_self = int((df["a"] == df["b"]).sum())
    edges = unique_edges(df)
    symbol_of = symbols(df)
    G = build_graph(edges, symbol_of)
    save_edges(edges, symbol_of)

    n, m = G.number_of_nodes(), G.number_of_edges()
    print(f"  self-loop records removed: {n_self}")
    print(f"  duplicate records collapsed: {len(df) - n_self - m}")

    print("\n  I. Size")
    print(f"    nodes = {n}")
    print(f"    edges = {m}")

    print("\n  II. Edge density")
    density = 2 * m / (n * (n - 1))
    print(f"    2m / (n(n-1)) = {density:.6e}  (networkx: {nx.density(G):.6e})")
    print(f"    average degree = {2 * m / n:.2f}")

    print("\n  III. Edge sources")
    counts = edges["systems"].value_counts()
    for s, c in counts.items():
        print(f"    {s:35s} {c:8d}  {c / m:.4f}")
    for system in SYSTEMS:
        c = edges["systems"].str.contains(system, regex=False).sum()
        print(f"    supported by {system:22s} {c:8d}  {c / m:.4f}")

    print(f"\n  components = {nx.number_connected_components(G)}, "
          f"largest = {len(max(nx.connected_components(G), key=len))} nodes")
    print(f"  BGLAP (632) in network: {632 in G}, degree = {G.degree(632) if 632 in G else 0}")



if __name__ == "__main__":
    main()
