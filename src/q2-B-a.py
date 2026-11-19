"""Q2 Part B (a) - first and second neighbourhoods of Osteocalcin (BGLAP) and their induced subgraph."""

import matplotlib.pyplot as plt
import networkx as nx
from adjustText import adjust_text

from cmn.biogrid import SEED, SUBGRAPH_FILE, load_graph, save_graph
from cmn.paths import save_fig

# Categorical slots 1-3 of the dataviz reference palette (validated for CVD separation).
COLORS = {0: "#2a78d6", 1: "#eb6834", 2: "#1baf7a"}
LAYER_NAMES = {0: "Osteocalcin (BGLAP)", 1: "first neighbours", 2: "second neighbours"}
INK = "#2b2b2b"


def neighbourhoods(G, seed, radius=2):
    """Distance from the seed for every node within `radius` hops."""
    return nx.single_source_shortest_path_length(G, seed, cutoff=radius)


def draw_subgraph(H, layer, symbol):
    """Force-directed drawing of the induced subgraph, nodes coloured by distance from the seed."""
    pos = nx.spring_layout(H, k=0.6 / H.number_of_nodes() ** 0.5, iterations=200, seed=7)
    fig, ax = plt.subplots(figsize=(9, 9))
    nx.draw_networkx_edges(H, pos, ax=ax, width=0.15, alpha=0.08, edge_color="#7f7f7f")
    core = [(SEED, v) for v in H[SEED]]
    nx.draw_networkx_edges(H, pos, edgelist=core, ax=ax, width=1.2, edge_color=COLORS[1])
    handles = {}
    for d in (2, 1, 0):
        nodes = [v for v in H if layer[v] == d]
        handles[d] = nx.draw_networkx_nodes(H, pos, nodelist=nodes, ax=ax, node_size={0: 400, 1: 160, 2: 8}[d],
                                            node_color=COLORS[d], edgecolors="white",
                                            linewidths={0: 1.5, 1: 1, 2: 0}[d],
                                            label=f"{LAYER_NAMES[d]} ({len(nodes)})")
    texts = [ax.text(*pos[v], symbol[v], fontsize=16, color=INK, fontweight="bold", ha="center", va="center",
                     bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
             for v in H if layer[v] <= 1]
    adjust_text(texts, ax=ax, expand=(1.6, 2.0), arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
    ax.legend(handles=[handles[d] for d in (0, 1, 2)], loc="lower left", frameon=False, labelcolor=INK, fontsize=15,
              markerscale=1.5)
    ax.set_axis_off()
    return fig


def draw_core(H, symbol):
    """BGLAP, its first neighbours and the edges among them."""
    core = H.subgraph([SEED, *H[SEED]])
    pos = nx.shell_layout(core, nlist=[[SEED], [v for v in core if v != SEED]])
    fig, ax = plt.subplots(figsize=(5, 5))
    systems = nx.get_edge_attributes(core, "systems")
    y2h = [e for e, s in systems.items() if "Two-hybrid" in s]
    ms = [e for e, s in systems.items() if "Two-hybrid" not in s]
    nx.draw_networkx_edges(core, pos, edgelist=ms, ax=ax, width=2, edge_color="#7f7f7f")
    nx.draw_networkx_edges(core, pos, edgelist=y2h, ax=ax, width=2, edge_color="#7f7f7f", style="dashed")
    nx.draw_networkx_nodes(core, pos, nodelist=[SEED], ax=ax, node_size=1500, node_color=COLORS[0])
    nx.draw_networkx_nodes(core, pos, nodelist=[v for v in core if v != SEED], ax=ax, node_size=1300,
                           node_color=COLORS[1])
    nx.draw_networkx_labels(core, pos, labels={v: symbol[v] for v in core}, ax=ax, font_size=8,
                            font_color="white", font_weight="bold")
    ax.plot([], [], color="#7f7f7f", lw=2, label="Affinity Capture-MS")
    ax.plot([], [], color="#7f7f7f", lw=2, ls="--", label="Two-hybrid")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False, labelcolor=INK)
    ax.set_axis_off()
    ax.margins(0.12)
    return fig


def main():
    G = load_graph()
    symbol = nx.get_node_attributes(G, "symbol")

    layer = neighbourhoods(G, SEED)
    N1 = sorted((v for v, d in layer.items() if d == 1), key=G.degree, reverse=True)
    N2 = [v for v, d in layer.items() if d == 2]

    print("I. Neighbourhoods of BGLAP (Entrez 632)")
    print(f"  first neighbours: {len(N1)}")
    for v in N1:
        own = sum(1 for w in G[v] if layer.get(w) == 2)
        only = sum(1 for w in G[v] if layer.get(w) == 2 and sum(1 for x in G[w] if x in N1) == 1)
        print(f"    {symbol[v]:8s} entrez {v:<7d} degree {G.degree(v):5d}  "
              f"system {G[SEED][v]['systems']:20s} second neighbours reached {own:5d} (only via it: {only})")
    print(f"  second neighbours: {len(N2)}")

    H = G.subgraph(layer).copy()
    save_graph(H, SUBGRAPH_FILE)
    print("\nII/IV. Induced subgraph on BGLAP + first + second neighbours")
    print(f"  nodes = {H.number_of_nodes()}, edges = {H.number_of_edges()}, connected = {nx.is_connected(H)}")
    print(f"  edges among first neighbours: {G.subgraph(N1).number_of_edges()}")
    print(f"  edge density = {nx.density(H):.4f}")
    print(f"  saved to {SUBGRAPH_FILE}")

    print("\nIII. Figures")
    save_fig(draw_subgraph(H, layer, symbol), "q2_B_a_subgraph.png")
    save_fig(draw_core(H, symbol), "q2_B_a_core.png")
    print("  q2_B_a_subgraph.png, q2_B_a_core.png")


if __name__ == "__main__":
    main()
