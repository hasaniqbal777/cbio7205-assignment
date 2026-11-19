"""Q3 Part B (b) - pairwise MAGNA++ alignments of the five STRING networks, maximising S3 only (alpha = 1)."""

import time
from itertools import combinations

import pandas as pd

from cmn.magna import run, write_gw
from cmn.paths import DATA
from cmn.string_db import SPECIES

POPULATION = 1500
GENERATIONS = 200
THREADS = 14
RESULTS_FILE = DATA / "q3_B_b_alignments.tsv"


def main():
    size = {}
    for name in SPECIES:
        _, n, m = write_gw(name)
        size[name] = n
        print(f"{name}: {n} nodes, {m} edges -> {name}.gw", flush=True)

    rows = []
    for a, b in combinations(SPECIES, 2):
        small, large = sorted((a, b), key=size.get)   # MAGNA++ needs the smaller network as G
        start = time.time()
        stats = run(small, large, POPULATION, GENERATIONS, measure="S3", threads=THREADS)
        rows.append({"G (smaller)": small, "H (larger)": large, "S3": float(stats["S3_score"]),
                     "EC": float(stats["EC_score"]), "ICS": float(stats["ICS_score"]),
                     "seconds": round(time.time() - start)})
        print(f"{small:6s} -> {large:6s}  S3 = {rows[-1]['S3']:.4f}  EC = {rows[-1]['EC']:.4f}  "
              f"ICS = {rows[-1]['ICS']:.4f}  ({rows[-1]['seconds']} s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS_FILE, sep="\t", index=False)
    S3 = pd.DataFrame(1.0, index=list(SPECIES), columns=list(SPECIES))
    for r in rows:
        S3.loc[r["G (smaller)"], r["H (larger)"]] = S3.loc[r["H (larger)"], r["G (smaller)"]] = r["S3"]
    print(f"\npopulation {POPULATION}, generations {GENERATIONS}, measure S3, alpha 1 (no node similarity)")
    print("\nS3 matrix")
    print(S3.round(4).to_string())
    print(f"\nresults saved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
