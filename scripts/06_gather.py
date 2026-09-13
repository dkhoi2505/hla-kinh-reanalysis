#!/usr/bin/env python3
"""Gop 303 ket qua -> bang genotype 101 mau kem nhan on dinh."""
import re
from pathlib import Path
import pandas as pd

LOCI = {"A": ["A1","A2"], "B": ["B1","B2"], "C": ["C1","C2"]}

rows = []
for f in sorted(Path("results/main").glob("*_result.tsv")):
    m = re.match(r"(SRR\d+)_(\d+)_s(\d+)_result", f.stem)
    if not m:
        continue
    d = pd.read_csv(f, sep="\t").iloc[0]
    r = {"sample": m.group(1), "seed": int(m.group(3)),
         "read_dung": int(d["Reads"]), "objective": float(d["Objective"])}
    for loc, cols in LOCI.items():
        # thu tu 2 allele trong 1 locus khong co nghia -> sap xep
        r[loc] = " ".join(sorted(str(d[c]) for c in cols))
    rows.append(r)

raw = pd.DataFrame(rows)
raw.to_csv("results/genotypes_by_seed.csv", index=False)
print(f"Doc {len(raw)} ket qua, {raw['sample'].nunique()} mau")

# --- gop theo mau, danh gia on dinh tung locus ---
out = []
for s, g in raw.groupby("sample"):
    r = {"sample": s, "n_seed": len(g),
         "read_min": g.read_dung.min(), "read_max": g.read_dung.max()}
    n_bat_on = 0
    for loc in LOCI:
        vals = g[loc].value_counts()
        r[loc] = vals.index[0]                    # genotype da so
        r[f"{loc}_on_dinh"] = (len(vals) == 1)
        if len(vals) > 1:
            n_bat_on += 1
            r[f"{loc}_khac"] = " | ".join(vals.index[1:])
    r["so_locus_bat_on"] = n_bat_on
    r["tin_cay"] = "cao" if n_bat_on == 0 else ("trung binh" if n_bat_on == 1 else "thap")
    out.append(r)

t = pd.DataFrame(out).sort_values(["so_locus_bat_on","sample"], ascending=[False,True])
t.to_csv("results/genotypes_consensus.csv", index=False)

print("\n=== DO ON DINH THEO LOCUS ===")
for loc in LOCI:
    n = t[f"{loc}_on_dinh"].sum()
    print(f"HLA-{loc}: {n}/{len(t)} mau on dinh ({100*n/len(t):.1f}%)")

print("\n=== PHAN LOAI MAU ===")
print(t["tin_cay"].value_counts().to_string())

print("\n=== MAU BAT ON ===")
bo = t[t.so_locus_bat_on > 0]
if len(bo):
    cols = ["sample","so_locus_bat_on","A","B","C","read_min","read_max"]
    print(bo[cols].to_string(index=False))
else:
    print("khong co")

print("\n-> results/genotypes_consensus.csv")
print("-> results/genotypes_by_seed.csv")
