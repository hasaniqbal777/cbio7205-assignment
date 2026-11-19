"""Q1 Question 4 - maximum number of 4-node cliques in a graph with 10 triangles."""

from math import comb

import networkx as nx


def cliques_by_size(G):
    counts = {}
    for c in nx.enumerate_all_cliques(G):
        counts[len(c)] = counts.get(len(c), 0) + 1
    return counts


def kk_shadow(q, k):
    """Minimum number of (k-1)-sets covered by q k-sets (Kruskal-Katona lower bound on the shadow)."""
    shadow, j = 0, k
    while q > 0 and j > 0:
        a = j
        while comb(a + 1, j) <= q:
            a += 1
        shadow += comb(a, j - 1)
        q -= comb(a, j)
        j -= 1
    return shadow


def main():
    K5 = nx.complete_graph(5)
    print("K5:", cliques_by_size(K5))

    print("\nKruskal-Katona: fewest triangles forced by q 4-cliques")
    for q in range(1, 8):
        print(f"  q = {q}: at least {kk_shadow(q, 4)} triangles")

    # Exhaustive check over every graph with at most 7 nodes (the networkx graph atlas).
    counts = [cliques_by_size(G) for G in nx.graph_atlas_g()]
    best = max(c.get(4, 0) for c in counts if c.get(3, 0) == 10)
    print(f"\nAtlas (all graphs up to 7 nodes) with exactly 10 triangles: max 4-cliques = {best}")


if __name__ == "__main__":
    main()
