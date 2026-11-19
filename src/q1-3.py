"""Q1 Question 3 - brute-force the maximum number of G6 graphlets in a connected graph with 10 nodes and 10 edges."""

from itertools import combinations

import networkx as nx

from cmn.graphlets import graphlet_counts, orbit_counts

N = 10


def main():
    # A connected graph with n nodes and n edges is a spanning tree plus one extra edge,
    # so every such graph arises from some non-isomorphic tree on 10 nodes.
    trees = list(nx.nonisomorphic_trees(N))
    best, best_graphs, checked = -1, [], 0
    for T in trees:
        for u, v in combinations(range(N), 2):
            if T.has_edge(u, v):
                continue
            G = T.copy()
            G.add_edge(u, v)
            checked += 1
            g6 = graphlet_counts(orbit_counts(G))["G6"]
            if g6 > best:
                best, best_graphs = g6, [G]
            elif g6 == best:
                best_graphs.append(G)

    print(f"{len(trees)} trees, {checked} graphs checked")
    print(f"maximum number of G6 graphlets: {best}")

    unique = []
    for G in best_graphs:
        if not any(nx.is_isomorphic(G, H) for H in unique):
            unique.append(G)
    print(f"{len(unique)} non-isomorphic graphs reach it:")
    for G in unique:
        tri = next(c for c in nx.enumerate_all_cliques(G) if len(c) == 3)
        pendants = sorted((G.degree(v) - 2 for v in tri), reverse=True)
        print(f"  triangle with {pendants} pendants on its three nodes")


if __name__ == "__main__":
    main()
