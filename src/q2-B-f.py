"""Q2 Part B (f) - the three highest-degree proteins of the Osteocalcin subgraph and their KEGG pathways."""

import networkx as nx
import pandas as pd

from cmn.apis import kegg, kegg_id_from_entrez, kegg_pairs, uniprot_by_entrez
from cmn.biogrid import SEED, SUBGRAPH_FILE, load_graph
from cmn.paths import DATA

SIMILAR_FILE = DATA / "q2_B_e_similar.tsv"
RESULTS_FILE = DATA / "q2_B_f_hub_pathways.tsv"
TOP = 3


def brite(entry):
    """ORTHOLOGY line and the BRITE hierarchy of a KEGG gene entry."""
    out, inside = [], False
    for line in entry.splitlines():
        field = line[:12].strip()
        if field == "ORTHOLOGY":
            out.append(line[12:].strip())
        if field:
            inside = field == "BRITE"
        if inside:
            out.append(line[12:].rstrip())
    return out


def main():
    H = load_graph(SUBGRAPH_FILE)
    symbol = nx.get_node_attributes(H, "symbol")
    hubs = sorted(H, key=H.degree, reverse=True)[:TOP]

    print("I/II. highest-degree proteins")
    for v in hubs:
        names = "; ".join(e["Protein names"].split(" (")[0] for e in uniprot_by_entrez(v))
        print(f"  {symbol[v]:8s} Entrez {v:<7d} degree {H.degree(v):5d}  "
              f"neighbour of BGLAP: {H.has_edge(SEED, v)!s:5s}  {names}")
    print(f"  next: " + ", ".join(f"{symbol[v]} ({H.degree(v)})"
                                  for v in sorted(H, key=H.degree, reverse=True)[TOP:TOP + 3]))

    print("\nIII/IV. KEGG pathways")
    pathway_names = dict(kegg_pairs("list", "pathway/hsa"))
    rows = []
    for v in hubs:
        kid = kegg_id_from_entrez(v)
        name_line = next((l for l in kegg("get", kid).splitlines() if l.startswith("NAME")), "")
        paths = [p for _, p in kegg_pairs("link", f"pathway/{kid}")]
        print(f"  {symbol[v]} ({kid}, {name_line[12:].strip()}): {len(paths)} pathways")
        if not paths:
            print("    no KEGG pathway map; KEGG orthology and BRITE classification:")
            for line in brite(kegg("get", kid)):
                print(f"      {line}")
        for p in paths:
            name = pathway_names.get(p.replace("path:", ""), "")
            print(f"    {p.replace('path:', '')}  {name}")
            rows.append({"symbol": symbol[v], "entrez": v, "kegg_gene": kid,
                         "pathway": p.replace("path:", ""), "name": name})
    pd.DataFrame(rows).to_csv(RESULTS_FILE, sep="\t", index=False)

    print("\nV. intersection with the 10 proteins most similar to BGLAP (part e)")
    similar = pd.read_csv(SIMILAR_FILE, sep="\t")
    common = set(similar["entrez"]) & set(hubs)
    print(f"  common proteins: {sorted(symbol[v] for v in common) or 'none'}")
    for v in hubs:
        linked = similar.loc[similar["neighbours"].str.split(", ").apply(lambda ns: symbol[v] in ns), "symbol"]
        print(f"  {symbol[v]:8s} interacts with {len(linked)} of the 10: {', '.join(linked) or '-'}")
    print(f"  degrees of the 10 similar proteins: {similar['degree'].min()}-{similar['degree'].max()}, "
          f"hubs: {', '.join(str(H.degree(v)) for v in hubs)}")
    print(f"\nresults saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
