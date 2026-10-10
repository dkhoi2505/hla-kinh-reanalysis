#!/usr/bin/env python3
"""
Khoang cach pseudosequence tu moi allele Kinh toi allele gan nhat
co du du lieu huan luyen.

Ly do: MHCflurry la pan-allele - no ma hoa allele bang pseudosequence
34 acid amin (cac vi tri tao pocket cua ranh gan peptide) va hoc quy luat
chung pocket-peptide. Mot allele it du lieu rieng van co the duoc du doan
tot NEU pseudosequence cua no gan mot allele giau du lieu.

NGUONG_BA = 50 la NGUONG TUY TIEN, khong co can cu nao.
Mot ban truoc cua docstring nay quy nguong 50 cho O'Donnell et al. 2020
(MHCflurry 2.0). Da tra y van 2026-10-10: KHONG dung duoc.
  - MHCflurry 1.x (O'Donnell et al. 2018) dung nguong 200 phep do moi allele,
    va noi ro "no information is shared between alleles" -- tuc khong co
    suy luan qua allele ho hang.
  - MHCflurry 2.0 la pan-allele, nen nguong du lieu toi thieu moi allele
    khong ap dung theo kieu do.
  - Khong truy cap duoc toan van 2.0 (ScienceDirect, PMC, EuropePMC deu chan),
    nen trang thai dung la "khong xac minh duoc", khong phai "da bac bo".
=> Doc ket qua nhu thang do TUONG DOI (allele nay gan ho hang giau du lieu
   hon allele kia), KHONG nhu nguong co y nghia sinh hoc. Chua kiem do nhay
   cua cac phan tram khoang cach theo nguong nay.

Dau vao : results/kinh_training_coverage.csv
          class1_pseudosequences.csv (tu mhcflurry-downloads)
          results/mhcflurry_training_counts.csv
Dau ra  : results/kinh_pseudoseq_distance.csv
"""
import subprocess, sys
from pathlib import Path
import pandas as pd

NGUONG_BA = 50

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

# --- doc pseudosequence ---
D = Path(subprocess.run(["mhcflurry-downloads","path","allele_sequences"],
                        capture_output=True, text=True).stdout.strip())
ps = pd.read_csv(D / "class1_pseudosequences.csv")
ps = ps[ps.allele.str.match(r"^HLA-[ABC]\d+:\d+$")].copy()
ps["allele2"] = ps.allele.str.replace("HLA-", "", regex=False)
ps["allele2"] = ps.allele2.str[0] + "*" + ps.allele2.str[1:]
ps = ps.drop_duplicates("allele2").set_index("allele2")["pseudosequence"]

# --- allele giau du lieu (nguon tham chieu) ---
cnt = pd.read_csv("results/mhcflurry_training_counts.csv", index_col=0)
giau = cnt[cnt.BA >= NGUONG_BA].index.intersection(ps.index)
print(f"Allele co >= {NGUONG_BA} phep do BA va co pseudosequence: {len(giau)}")
print("  theo locus: " + ", ".join(
    f"{l}={sum(1 for a in giau if a[0]==l)}" for l in "ABC"))

# --- allele Kinh ---
kinh = pd.read_csv("results/kinh_training_coverage.csv", index_col=0)

rows = []
for a in kinh.index:
    seq = ps[a]
    # trong cung locus
    ung_vien = [g for g in giau if g[0] == a[0]]
    d_loc = sorted((hamming(seq, ps[g]), g) for g in ung_vien)
    # xuyen locus
    d_all = sorted((hamming(seq, ps[g]), g) for g in giau)
    rows.append({
        "allele": a,
        "af_kinh": kinh.loc[a, "af_kinh"],
        "BA": kinh.loc[a, "BA"],
        "tong": kinh.loc[a, "tong"],
        "gan_nhat_cung_locus": d_loc[0][1] if d_loc else None,
        "kc_cung_locus":       d_loc[0][0] if d_loc else None,
        "BA_cua_no":           cnt.loc[d_loc[0][1], "BA"] if d_loc else None,
        "kc_xuyen_locus":      d_all[0][0] if d_all else None,
    })

t = pd.DataFrame(rows).set_index("allele").sort_values("af_kinh", ascending=False)
t["locus"] = t.index.str[0]
t.to_csv("results/kinh_pseudoseq_distance.csv")

print("\n=== 20 ALLELE PHO BIEN NHAT ===")
print(t.head(20)[["af_kinh","BA","gan_nhat_cung_locus","kc_cung_locus","BA_cua_no"]]
      .to_string())

print("\n=== ALLELE NGHEO DU LIEU (BA < 50) ===")
ngheo = t[t.BA < NGUONG_BA].sort_values("af_kinh", ascending=False)
print(f"So allele: {len(ngheo)}  |  tong AF: {ngheo.af_kinh.sum():.3f}")
print(ngheo[["af_kinh","BA","tong","gan_nhat_cung_locus","kc_cung_locus","BA_cua_no"]]
      .to_string())

print("\n=== PHAN BO KHOANG CACH (trong so AF, chi allele ngheo) ===")
for loc, g in ngheo.groupby("locus"):
    s = g.af_kinh.sum()
    if s == 0: continue
    print(f"HLA-{loc}  (AF ngheo {s:.3f})")
    for lo, hi, nhan in [(0,2,"kc 0-2  (rat gan)"),
                         (3,5,"kc 3-5  (gan)"),
                         (6,9,"kc 6-9  (xa)"),
                         (10,99,"kc 10+  (rat xa)")]:
        w = g[(g.kc_cung_locus>=lo)&(g.kc_cung_locus<=hi)].af_kinh.sum()
        print(f"   {nhan}: {100*w/s:5.1f}%")

print("\n-> results/kinh_pseudoseq_distance.csv")
