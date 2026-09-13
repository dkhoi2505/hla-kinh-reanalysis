#!/usr/bin/env python3
"""Tong hop thi nghiem beta: so ket qua o cac muc beta khac nhau."""
import re
from pathlib import Path
import pandas as pd

LOCI = {"A": ["A1","A2"], "B": ["B1","B2"], "C": ["C1","C2"]}
D = Path("results/beta_experiment")

rows = []
for f in sorted(D.glob("*_beta*_result.tsv")):
    m = re.match(r"(SRR\d+)_beta([\d.]+)_result", f.stem)
    if not m:
        continue
    d = pd.read_csv(f, sep="\t").iloc[0]
    r = {"sample": m.group(1), "beta": float(m.group(2))}
    for loc, cols in LOCI.items():
        v = sorted(str(d[c]) for c in cols)
        r[loc] = " ".join(v)
        r[f"{loc}_homo"] = (v[0] == v[1])
    rows.append(r)

t = pd.DataFrame(rows).sort_values(["sample","beta"])
t.to_csv("results/beta_experiment_summary.csv", index=False)

print("=== KET QUA THEO MAU ===")
for s, g in t.groupby("sample"):
    print(f"\n{s}")
    for _, r in g.iterrows():
        print(f"  beta={r.beta:<6} A: {r.A:<22} B: {r.B:<22} C: {r.C}")

# CHI so tren nhung mau chay du MOI muc beta (mau so phai bang nhau)
du = t.groupby("sample")["beta"].nunique()
chung = du[du == t.beta.nunique()].index
print(f"\n=== DONG HOP TU THEO BETA (chi {len(chung)} mau chay du moi muc) ===")
if len(chung):
    sub = t[t["sample"].isin(chung)]
    for b, g in sub.groupby("beta"):
        n = sum(g[f"{loc}_homo"].sum() for loc in LOCI)
        print(f"beta={b:<6}: {n}/{3*len(g)} locus dong hop tu")
else:
    print("khong co mau nao chay du moi muc beta")

# Chieu 2: so voi ket qua goc (beta mac dinh 0.009) trong genotypes_consensus.csv
print("\n=== CHIEU 2: kiem tra nguoc (so voi beta mac dinh) ===")
goc_f = Path("results/genotypes_consensus.csv")
rieng = t[~t["sample"].isin(chung)]
if goc_f.exists() and len(rieng):
    goc = pd.read_csv(goc_f).set_index("sample")
    n_hong = 0
    for _, r in rieng.iterrows():
        s_ = r["sample"]
        if s_ not in goc.index:
            continue
        for loc in LOCI:
            cu_ = " ".join(sorted(str(goc.loc[s_, loc]).split()))
            moi_ = r[loc]
            if cu_ != moi_:
                n_hong += 1
                print(f"{s_} {loc}: {cu_}  ->  {moi_}   <<< DOI")
    print(f"\nTong: {n_hong} locus bi doi khi ha beta xuong {rieng.beta.iloc[0]}")
else:
    print("thieu results/genotypes_consensus.csv de so sanh")

print("\n=== KET LUAN ===")
print("beta cao  -> false HOMOzygote (allele hiem bi goi thanh dong hop tu)")
print("beta thap -> false HETEROzygote (dong hop tu dung bi tach nham)")
print("=> Khong co gia tri beta nao dung cho moi mau.")

print("\n-> results/beta_experiment_summary.csv")
