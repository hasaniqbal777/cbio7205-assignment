"""Q2 Part B (e) - proteins with the most similar graphlet signature to Osteocalcin, and their GO annotations."""

import networkx as nx
import pandas as pd

from cmn.apis import go_terms, uniprot_by_entrez
from cmn.biogrid import SEED, SUBGRAPH_FILE, load_graph
from cmn.graphlets import ORBITS, signature_similarity
from cmn.paths import DATA

ORBITS_FILE = DATA / "bglap_subgraph_orbits.tsv"
RESULTS_FILE = DATA / "q2_B_e_similar.tsv"
TOP = 10
# GO cellular component: extracellular region, extracellular space.
EXTRACELLULAR = {"GO:0005576", "GO:0005615"}


def go_of(entrez):
    """Protein names and GO biological process terms (name [GO:id]) of a gene, merged over its reviewed entries."""
    entries = uniprot_by_entrez(entrez, fields=("accession", "gene_primary", "protein_name", "go_p", "go_c"))
    names = "; ".join(e["Protein names"].split(" (")[0] for e in entries)
    terms, places = {}, set()
    for e in entries:
        for t in (e["Gene Ontology (biological process)"] or "").split("; "):
            if t:
                terms[go_terms(t).pop()] = t.rsplit(" [", 1)[0]
        places |= go_terms(e["Gene Ontology (cellular component)"])
    return names, terms, places


def main():
    H = load_graph(SUBGRAPH_FILE)
    table = pd.read_csv(ORBITS_FILE, sep="\t", index_col="entrez")
    symbol = table["symbol"]
    O = table[[str(i) for i in ORBITS]]
    O.columns = ORBITS

    S = signature_similarity(O, SEED).drop(SEED).sort_values(ascending=False)
    print(f"I. signature similarity to BGLAP over {len(S)} proteins: "
          f"max {S.max():.4f}, median {S.median():.4f}, min {S.min():.4f}")

    print("\nII. identical signatures")
    same = (O == O.loc[SEED]).all(axis=1).drop(SEED)
    print(f"  proteins with the same 15 orbit counts as BGLAP: {int(same.sum())}")
    print(f"  proteins with similarity 1: {int((S >= 1 - 1e-12).sum())}")

    seed_names, seed_go, seed_places = go_of(SEED)
    print(f"\n  BGLAP: {len(seed_go)} GO biological process terms")

    rows = []
    for rank, v in enumerate(S.index[:TOP], 1):
        names, terms, places = go_of(v)
        shared = sorted(set(terms) & set(seed_go))
        rows.append({
            "rank": rank, "entrez": v, "symbol": symbol[v], "protein": names, "similarity": S[v],
            "degree": H.degree(v), "distance": nx.shortest_path_length(H, SEED, v),
            "neighbour_of_bglap": H.has_edge(SEED, v), "go_terms": len(terms),
            "shared_go": "; ".join(f"{seed_go[t]} [{t}]" for t in shared),
            "signature": " ".join(str(x) for x in O.loc[v]),
            "neighbours": ", ".join(sorted(symbol[w] for w in H[v])),
            "extracellular": bool(places & EXTRACELLULAR), "er_lumen": "GO:0005788" in places,
        })
    top = pd.DataFrame(rows)
    top.to_csv(RESULTS_FILE, sep="\t", index=False)

    print(f"\nI/V. top {TOP} proteins by signature similarity")
    print(top[["rank", "symbol", "entrez", "similarity", "degree", "distance", "neighbour_of_bglap", "go_terms"]]
          .to_string(index=False))
    print(f"\n  BGLAP signature: {' '.join(str(x) for x in O.loc[SEED])}")
    for r in rows:
        print(f"  {r['symbol']:10s} {r['signature']}")

    print("\nIII. GO biological process terms shared with BGLAP")
    for r in rows:
        print(f"  {r['symbol']:10s} {r['protein'][:60]:60s} shared: {r['shared_go'] or '-'}")

    print("\n  location (GO cellular component) and interaction partners")
    print(f"  {'BGLAP':10s} extracellular {bool(seed_places & EXTRACELLULAR)!s:5s} ER lumen {'GO:0005788' in seed_places!s:5s}"
          f" partners: {', '.join(sorted(symbol[w] for w in H[SEED]))}")
    for r in rows:
        print(f"  {r['symbol']:10s} extracellular {r['extracellular']!s:5s} ER lumen {r['er_lumen']!s:5s}"
              f" partners: {r['neighbours']}")
    via = pd.Series([w for r in rows for w in r["neighbours"].split(", ")]).value_counts()
    print(f"  most common partners of the top {TOP}: " + ", ".join(f"{p} ({c})" for p, c in via.head(5).items()))

    print("\nV. first neighbours of BGLAP and their rank")
    rank = {v: i for i, v in enumerate(S.index, 1)}
    for v in H[SEED]:
        print(f"  {symbol[v]:8s} rank {rank[v]:5d}  similarity {S[v]:.4f}")
    print(f"\nresults saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
