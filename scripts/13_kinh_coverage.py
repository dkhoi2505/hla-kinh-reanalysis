#!/usr/bin/env python3
"""
Ghep tan suat allele nguoi Kinh (Bang 1, Do et al. 2020) voi so luong
peptide trong du lieu huan luyen cua MHCflurry.

Cau hoi: cac allele pho bien o quan the nay duoc dac trung toi dau trong
du lieu huan luyen cua cong cu du doan gan ket peptide-MHC?

LUU Y: MHCflurry la pan-allele. "It du lieu" KHONG dong nghia "du doan kem"
- mo hinh khai quat hoa tu allele co pseudosequence tuong tu. Buoc nay do
du lieu TRUC TIEP; khoang cach pseudosequence do o buoc sau.

Dau vao : data/kinh_af.tsv                      (allele 3-field + AF)
          results/mhcflurry_training_counts.csv (dem theo allele 2-field)
Dau ra  : results/kinh_training_coverage.csv
"""
import sys
from pathlib import Path
import pandas as pd

NGUONG = 50   # nguong "duoc dac trung tot", theo MHCflurry 2.0

def to2(a):
    """A*11:01:01 -> A*11:01 ; A*11:04 giu nguyen"""
    p = a.split(":")
    return ":".join(p[:2]) if len(p) >= 2 else a

for f in ["data/kinh_af.tsv", "results/mhcflurry_training_counts.csv"]:
    if not Path(f).exists():
        sys.exit(f"Khong thay {f}")

# --- tan suat Kinh, gop ve 2-field ---
kinh = pd.read_csv("data/kinh_af.tsv", sep="\t", names=["allele", "af"])
kinh["allele2"] = kinh["allele"].map(to2)
kinh = kinh.groupby("allele2")["af"].sum().to_frame("af_kinh")
kinh["locus"] = kinh.index.str[0]
kinh = kinh[kinh.locus.isin(["A", "B", "C"])]

# --- dem du lieu huan luyen ---
cnt = pd.read_csv("results/mhcflurry_training_counts.csv", index_col=0)

t = kinh.join(cnt[["BA", "EL", "tong"]], how="left").fillna(0)
for c in ["BA", "EL", "tong"]:
    t[c] = t[c].astype(int)
t["ti_le_BA"] = (100 * t.BA / t.tong).round(1).where(t.tong > 0, 0)
t = t.sort_values("af_kinh", ascending=False)
t.to_csv("results/kinh_training_coverage.csv")

print(f"Allele class I trong Bang 1 (2-field): {len(t)}")
print("Tong AF: " + ", ".join(f"{l}={g.af_kinh.sum():.3f}"
                              for l, g in t.groupby("locus")))

print("\n=== 20 ALLELE PHO BIEN NHAT ===")
print(t.head(20)[["locus","af_kinh","BA","EL","tong","ti_le_BA"]].to_string())

print("\n=== GANH NANG ALLELE THEO MUC DU LIEU (trong so tan suat) ===")
for loc, g in t.groupby("locus"):
    s = g.af_kinh.sum()
    print(f"HLA-{loc}  (tong AF {s:.3f})")
    print(f"   khong co du lieu : {100*g[g.tong==0].af_kinh.sum()/s:5.1f}%")
    print(f"   < {NGUONG} peptide     : {100*g[(g.tong>0)&(g.tong<NGUONG)].af_kinh.sum()/s:5.1f}%")
    print(f"   khong co BA      : {100*g[g.BA==0].af_kinh.sum()/s:5.1f}%")
    print(f"   BA < 100         : {100*g[g.BA<100].af_kinh.sum()/s:5.1f}%")

print("\n-> results/kinh_training_coverage.csv")
