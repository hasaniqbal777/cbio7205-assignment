"""Graphlet helpers for 2-to-4 node graphlets (G0-G8, orbits 0-14), numbered as in Figure 1 of the homework."""

from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd
from orbit_count import node_orbit_counts

GRAPHLET_NAMES = {
    "G0": "edge",
    "G1": "3-node path",
    "G2": "triangle",
    "G3": "4-node path",
    "G4": "4-node star",
    "G5": "4-cycle",
    "G6": "triangle + pendant",
    "G7": "diamond",
    "G8": "4-clique",
}

# For each graphlet: one of its orbits and how many of its nodes sit in that orbit.
# Summing the orbit over all nodes and dividing by this number gives the graphlet count.
GRAPHLET_FROM_ORBIT = {
    "G0": (0, 2),
    "G1": (2, 1),
    "G2": (3, 3),
    "G3": (5, 2),
    "G4": (7, 1),
    "G5": (8, 4),
    "G6": (11, 1),
    "G7": (13, 2),
    "G8": (14, 4),
}

ORBITS = list(range(15))


def graphlet_type(G, nodes):
    """Graphlet type (G0-G8) of the subgraph induced on `nodes`, or None if it is disconnected."""
    H = G.subgraph(nodes)
    if not nx.is_connected(H):
        return None
    n, m = H.number_of_nodes(), H.number_of_edges()
    degrees = sorted(d for _, d in H.degree())
    if n == 2:
        return "G0"
    if n == 3:
        return "G2" if m == 3 else "G1"
    if m == 3:
        return "G4" if degrees[-1] == 3 else "G3"
    if m == 4:
        return "G5" if degrees == [2, 2, 2, 2] else "G6"
    return "G7" if m == 5 else "G8"


def enumerate_graphlets(G, sizes=(2, 3, 4)):
    """All induced connected graphlets of G, as {type: [sorted node tuples]}."""
    found = {g: [] for g in GRAPHLET_NAMES}
    for k in sizes:
        for nodes in combinations(sorted(G.nodes()), k):
            g = graphlet_type(G, nodes)
            if g is not None:
                found[g].append(nodes)
    return found


def orbit_counts(G, nodes=None):
    """Orbit counts (orbits 0-14) per node, as a DataFrame indexed by node."""
    nodes = list(G.nodes()) if nodes is None else list(nodes)
    O = node_orbit_counts(G, graphlet_size=4, node_list=nodes)
    return pd.DataFrame(O, index=nodes, columns=ORBITS)


def graphlet_counts(O):
    """Number of graphlets of each type from a per-node orbit-count DataFrame."""
    return {g: int(O[orbit].sum() // k) for g, (orbit, k) in GRAPHLET_FROM_ORBIT.items()}


def signature_similarity(O, u):
    """Signature similarity of Milenkovic and Przulj (2008), equal orbit weights, between node u and every node."""
    X = np.asarray(O, dtype=float)
    x = np.asarray(O.loc[u], dtype=float)
    D = np.abs(np.log(X + 1) - np.log(x + 1)) / np.log(np.maximum(X, x) + 2)
    return pd.Series(1 - D.mean(axis=1), index=O.index)
