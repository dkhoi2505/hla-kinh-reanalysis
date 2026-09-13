#!/usr/bin/env python3
"""Dem so peptide huan luyen cho tung allele HLA class I."""
import subprocess
from pathlib import Path
import pandas as pd

D = Path(subprocess.run(["mhcflurry-downloads","path","data_curated"],
                        capture_output=True, text=True).stdout.strip())
f = D / "curated_training_data.csv.bz2"
print(f"Doc: {f}")

df = pd.read_csv(f, compression="bz2",
                 usecols=["allele","peptide","measurement_kind"])
print(f"Tong dong: {len(df):,}")

# chi giu HLA class I nguoi
df = df[df["allele"].str.match(r"^HLA-[ABC]\*")]
print(f"HLA class I: {len(df):,}")

# chuan hoa ve 2-field
df["allele2"] = df["allele"].str.replace("HLA-","",regex=False).str.extract(r"^([ABC]\*\d+:\d+)")

# dem peptide DUY NHAT, khong dem so phep do
g = df.groupby(["allele2","measurement_kind"])["peptide"].nunique().unstack(fill_value=0)
g.columns.name = None
g = g.rename(columns={"affinity":"BA","mass_spec":"EL"})
for c in ["BA","EL"]:
    if c not in g: g[c] = 0
g["tong"] = g["BA"] + g["EL"]
g = g.sort_values("tong", ascending=False)

g.to_csv("results/mhcflurry_training_counts.csv")
print(f"\nSo allele 2-field co du lieu: {len(g)}")
print("\n=== 15 ALLELE GIAU DU LIEU NHAT ===")
print(g.head(15).to_string())

print("\n=== PHAN BO ===")
print(f"co >= 50 peptide : {(g['tong']>=50).sum()}")
print(f"co <  50 peptide : {(g['tong']<50).sum()}")
print(f"KHONG co BA nao  : {(g['BA']==0).sum()}")
print(f"KHONG co EL nao  : {(g['EL']==0).sum()}")
print("\n-> results/mhcflurry_training_counts.csv")
