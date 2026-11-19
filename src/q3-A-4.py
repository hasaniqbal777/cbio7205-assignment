"""Q3 Part A Question 4 - edge correctness of the Figure 2 alignment and annotation transfer for Figure 3."""

import networkx as nx

# Figure 2 (edges read off the figure); the alignment maps x -> x'.
FIG2_N1 = [("a", "b"), ("b", "c"), ("c", "d"), ("d", "e"), ("d", "f")]
FIG2_N2 = [("a'", "b'"), ("b'", "c'"), ("c'", "d'"), ("c'", "e'"), ("d'", "f'"), ("e'", "f'"), ("f'", "g'")]
FIG2_MAP = {x: f"{x}'" for x in "abcdef"}

# Figure 3: networks, alignment and annotation table.
FIG3_N1 = [("a", "b"), ("b", "c"), ("c", "d"), ("c", "e")]
FIG3_N2 = [("a'", "b'"), ("b'", "c'"), ("c'", "d'"), ("c'", "e'"), ("d'", "e'"), ("e'", "f'")]
FIG3_MAP = {x: f"{x}'" for x in "abcde"}
FIG3_ANN = {
    "a": {1, 2}, "b": {3}, "c": set(), "d": {1, 4}, "e": {5},
    "a'": {2, 3}, "b'": {2}, "c'": {1}, "d'": {1, 2, 4}, "e'": {5}, "f'": {3, 4},
}


def alignment_scores(G1, G2, f):
    """Edge correctness, induced conserved substructure and symmetric substructure score of alignment f: V1 -> V2."""
    conserved = [(u, v) for u, v in G1.edges() if G2.has_edge(f[u], f[v])]
    induced = G2.subgraph(f.values()).number_of_edges()
    m1, c = G1.number_of_edges(), len(conserved)
    return conserved, c / m1, c / induced, c / (m1 + induced - c), induced


def main():
    G1, G2 = nx.Graph(FIG2_N1), nx.Graph(FIG2_N2)
    print("Figure 2")
    print(f"  N1: {G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges; "
          f"N2: {G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges")
    print(f"  aligned N1 nodes: {len(FIG2_MAP)}/{G1.number_of_nodes()}, injective: "
          f"{len(set(FIG2_MAP.values())) == len(FIG2_MAP)}, unaligned N2 nodes: {sorted(set(G2) - set(FIG2_MAP.values()))}")
    for u, v in G1.edges():
        print(f"    {u}-{v} -> {FIG2_MAP[u]}-{FIG2_MAP[v]}: {'conserved' if G2.has_edge(FIG2_MAP[u], FIG2_MAP[v]) else 'not an edge'}")
    conserved, ec, ics, s3, induced = alignment_scores(G1, G2, FIG2_MAP)
    print(f"  EC  = {len(conserved)}/{G1.number_of_edges()} = {ec:.3f}")
    print(f"  ICS = {len(conserved)}/{induced} = {ics:.3f}")
    print(f"  S3  = {len(conserved)}/({G1.number_of_edges()} + {induced} - {len(conserved)}) = {s3:.3f}")

    print("\nFigure 3: annotation transfer between aligned nodes")
    new = {n: set(a) for n, a in FIG3_ANN.items()}
    for x, y in FIG3_MAP.items():
        union = FIG3_ANN[x] | FIG3_ANN[y]
        new[x], new[y] = set(union), set(union)
    for n in FIG3_ANN:
        added = new[n] - FIG3_ANN[n]
        print(f"  {n:3s} before {sorted(FIG3_ANN[n])!s:12s} after {sorted(new[n])!s:12s} added {sorted(added)}")
    conserved, ec, _, s3, _ = alignment_scores(nx.Graph(FIG3_N1), nx.Graph(FIG3_N2), FIG3_MAP)
    print(f"  (Figure 3 alignment: EC = {ec:.3f}, S3 = {s3:.3f})")


if __name__ == "__main__":
    main()
