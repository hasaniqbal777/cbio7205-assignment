"""Q3 Part B (a) - STRING v12.5 PPI networks of five species (combined score >= 0.7) as edge lists."""

import pandas as pd

from cmn.paths import DATA
from cmn.string_db import MIN_SCORE, SPECIES, VERSION, edge_list, links_url, write_edges

SUMMARY_FILE = DATA / "q3_B_a_networks.tsv"


def main():
    rows = []
    for name, (taxid, latin) in SPECIES.items():
        print(f"{name} ({latin}, {taxid}): {links_url(taxid)}")
        edges = edge_list(taxid)
        path = write_edges(edges, name)
        n = pd.unique(edges[["u", "v"]].values.ravel()).size
        m = len(edges)
        rows.append({"species": name, "latin": latin, "taxid": taxid, "nodes": n, "edges": m,
                     "density": 2 * m / (n * (n - 1)), "avg_degree": 2 * m / n})
        print(f"  -> {path.name}: {n} nodes, {m} edges")

    df = pd.DataFrame(rows)
    df.to_csv(SUMMARY_FILE, sep="\t", index=False)
    print(f"\nSTRING v{VERSION}, combined_score >= {MIN_SCORE}")
    print(df.to_string(index=False, formatters={"density": "{:.3e}".format, "avg_degree": "{:.1f}".format}))


if __name__ == "__main__":
    main()
