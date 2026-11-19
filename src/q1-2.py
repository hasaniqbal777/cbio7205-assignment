"""Q1 Question 2 - cross-check the by-hand graph measures on the 9-node graph."""

from collections import Counter
from itertools import combinations

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

from cmn.graphlets import enumerate_graphlets, graphlet_counts, orbit_counts

# Edges read off the figure: a 4-cycle 1-2-3-4, bridge 4-5, and the K4 {6,7,8,9} with 5 joined to 6 and 7.
EDGES = [(1, 2), (2, 3), (3, 4), (1, 4), (4, 5), (5, 6), (5, 7),
         (6, 7), (6, 8), (6, 9), (7, 8), (7, 9), (8, 9)]


def distribution(values):
    """Count and fraction of nodes (or pairs) taking each value."""
    c = Counter(values)
    n = sum(c.values())
    return {v: (c[v], c[v] / n) for v in sorted(c)}


def main():
    G = nx.Graph(EDGES)
    n = G.number_of_nodes()
    print(f"n = {n}, m = {G.number_of_edges()}")

    print("\nI. Degree distribution")
    deg = dict(G.degree())
    print("  degrees:", deg)
    for k, (c, p) in distribution(deg.values()).items():
        print(f"  P({k}) = {c}/{n} = {p:.3f}")

    print("\nII. Clustering")
    C = nx.clustering(G)
    for v in G:
        print(f"  C({v}) = {C[v]:.4f}  (triangles = {nx.triangles(G, v)}, deg = {deg[v]})")
    print(f"  average C = {nx.average_clustering(G):.4f}  (= {sum(C.values())}/{n})")
    print("  clustering spectrum C(k):")
    for k in sorted(set(deg.values())):
        nodes = [v for v in G if deg[v] == k]
        print(f"    C({k}) = {sum(C[v] for v in nodes) / len(nodes):.4f}  over nodes {nodes}")

    print("\nIII. Shortest path lengths")
    pairs = {(u, v): nx.shortest_path_length(G, u, v) for u, v in combinations(sorted(G), 2)}
    for d, (c, p) in distribution(pairs.values()).items():
        print(f"  length {d}: {c}/{len(pairs)} = {p:.3f}  {[uv for uv, l in pairs.items() if l == d]}")
    print(f"  average = {sum(pairs.values()) / len(pairs):.4f}, diameter = {nx.diameter(G)}")

    print("\nIV. Betweenness centrality")
    B = nx.betweenness_centrality(G, normalized=False)
    Bn = nx.betweenness_centrality(G, normalized=True)
    for v in G:
        print(f"  B({v}) = {B[v]:.4f}   normalized = {Bn[v]:.4f}")
    for b, (c, p) in distribution(round(x, 4) for x in B.values()).items():
        print(f"  B = {b}: {c}/{n}")

    print("\nV. Graphlets G5 and G6")
    found = enumerate_graphlets(G)
    for g in ("G5", "G6"):
        print(f"  {g} ({len(found[g])}):")
        for nodes in found[g]:
            H = G.subgraph(nodes)
            print(f"    {nodes}  edges {sorted(tuple(sorted(e)) for e in H.edges())}")
    print("  all graphlet counts (enumeration):", {g: len(v) for g, v in found.items()})

    print("\nVI. Automorphisms and orbits")
    autos = list(GraphMatcher(G, G).isomorphisms_iter())
    print(f"  {len(autos)} automorphisms")
    for a in autos:
        moved = {v: a[v] for v in sorted(a) if a[v] != v}
        print("   ", moved or "identity")
    orbits = {frozenset(a[v] for a in autos) for v in G}
    print("  orbits:", sorted(sorted(o) for o in orbits))

    print("\nVII. Graphlet signatures (orbits 0-14)")
    O = orbit_counts(G, sorted(G))
    print(O.to_string())
    print("  graphlet counts (orbits):", graphlet_counts(O))
    print("  nodes 2, 4, 7, 9:")
    print(O.loc[[2, 4, 7, 9]].to_string())


if __name__ == "__main__":
    main()
