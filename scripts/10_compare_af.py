#!/usr/bin/env python3
"""So tan suat allele: ket qua OptiType vs Bang 1 (Do et al. 2020)."""
import pandas as pd

def to2(a):
    """Gop ve 2-field: A*11:01:01 -> A*11:01 ; A*11:04 giu nguyen."""
    p = a.split(":")
    return ":".join(p[:2]) if len(p) >= 2 else a

# --- tan suat tu ket qua cua minh ---
t = pd.read_csv("results/genotypes_consensus.csv")
dem = {}
for loc in ["A","B","C"]:
    for g in t[loc]:
        for a in g.split():
            dem[a] = dem.get(a, 0) + 1

mine = pd.Series(dem, name="count_toi").to_frame()
mine["locus"] = mine.index.str[0]
mine["af_toi"] = mine.groupby("locus")["count_toi"].transform(lambda x: x/x.sum())

# --- Bang 1, gop ve 2-field ---
ref = pd.read_csv("data/kinh_af.tsv", sep="\t", names=["allele","af"])
ref["allele2"] = ref["allele"].map(to2)
ref = ref.groupby("allele2")["af"].sum().to_frame("af_bai")
ref["locus"] = ref.index.str[0]

# --- ghep ---
m = mine[["af_toi","count_toi"]].join(ref[["af_bai"]], how="outer")
m["locus"] = m.index.str[0]
m = m.fillna({"af_toi":0, "af_bai":0, "count_toi":0})
m["chenh"] = (m.af_toi - m.af_bai).round(4)
m = m.sort_values("chenh", key=abs, ascending=False)

m.to_csv("results/af_comparison.csv")

print(f"Allele trong ket qua cua toi : {(m.af_toi>0).sum()}")
print(f"Allele trong Bang 1 (2-field): {(m.af_bai>0).sum()}")
print(f"Co o ca hai                  : {((m.af_toi>0)&(m.af_bai>0)).sum()}")
print(f"CHI co o ket qua toi         : {((m.af_toi>0)&(m.af_bai==0)).sum()}")
print(f"CHI co o Bang 1              : {((m.af_toi==0)&(m.af_bai>0)).sum()}")

print("\n=== TUONG QUAN (chi allele co o ca hai) ===")
both = m[(m.af_toi>0)&(m.af_bai>0)]
for loc in ["A","B","C"]:
    s = both[both.locus==loc]
    if len(s) > 2:
        print(f"HLA-{loc}: n={len(s)}, Spearman={s.af_toi.corr(s.af_bai, method='spearman'):.3f}, "
              f"Pearson={s.af_toi.corr(s.af_bai):.3f}")

print("\n=== 20 ALLELE LECH NHIEU NHAT ===")
print(m.head(20)[["locus","count_toi","af_toi","af_bai","chenh"]].to_string())

print("\n-> results/af_comparison.csv")
