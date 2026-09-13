#!/usr/bin/env python3
"""
Kiem toan tham chieu gen cua OptiType.

Cau hoi: bao nhieu muc trong tham chieu dung trinh tu intron MUON
tu allele ho hang, thay vi trinh tu gen that cua chinh allele do?

Dau ra:
  results/reference_audit.csv       - tung muc trong tham chieu
  results/reference_donors.csv      - allele nao cho muon nhieu nhat
"""
import os
import re
import sys
from pathlib import Path
import pandas as pd

def tim_reference():
    p = os.environ.get("CONDA_PREFIX")
    if not p:
        sys.exit("Chua kich hoat conda env. Chay: conda activate hla")
    f = Path(p) / "share/optitype/data/hla_reference_dna.fasta"
    if not f.exists():
        sys.exit(f"Khong thay tham chieu tai {f}")
    return f

def doc_header(f):
    """Moi header dang:
       >HLA00430 HLA-C*06:02:01:01
       >HLA00431_HLA00430 HLA-C*06:03:01 (introns from HLA-C*06:02:01:01)
    """
    rows = []
    for line in f.open():
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        parts = h.split()
        if len(parts) < 2 or not parts[1].startswith("HLA-"):
            continue                      # bo qua pseudogene E/F/G/H/J/K/L/V
        dich = parts[1].replace("HLA-", "")
        m = re.search(r"introns from (\S+)", h)
        donor = m.group(1).replace("HLA-", "").rstrip(")") if m else None
        rows.append({
            "muc_id": parts[0],
            "allele_dich": dich,
            "dich_2field": ":".join(dich.split(":")[:2]),
            "donor": donor,
            "donor_2field": ":".join(donor.split(":")[:2]) if donor else None,
            "intron_muon": donor is not None,
            "locus": dich[0],
        })
    return pd.DataFrame(rows)

def main():
    ref = tim_reference()
    print(f"Tham chieu: {ref}")
    df = doc_header(ref)
    df = df[df.locus.isin(["A", "B", "C"])]

    n_muc = len(df)
    n_muon = df.intron_muon.sum()
    n_dich = df.allele_dich.nunique()

    print(f"\n=== QUY MO ===")
    print(f"So muc trong tham chieu : {n_muc:,}")
    print(f"So allele dich khac nhau: {n_dich:,}")
    print(f"He so no (muc/allele)   : {n_muc/n_dich:.2f}")

    print(f"\n=== INTRON MUON ===")
    print(f"Muc dung intron muon : {n_muon:,} ({100*n_muon/n_muc:.1f}%)")
    print(f"Muc co gen that      : {n_muc-n_muon:,} ({100*(n_muc-n_muon)/n_muc:.1f}%)")

    print(f"\n=== THEO LOCUS ===")
    t = df.groupby("locus").agg(
        muc=("muc_id", "count"),
        allele_dich=("allele_dich", "nunique"),
        muon=("intron_muon", "sum"),
    )
    t["ti_le_muon_%"] = (100*t.muon/t.muc).round(1)
    print(t.to_string())

    # allele nao cho muon nhieu nhat
    donors = (df[df.intron_muon]
              .groupby("donor_2field")["dich_2field"].nunique()
              .sort_values(ascending=False)
              .to_frame("so_allele_2field_muon"))
    print(f"\n=== 10 DONOR LON NHAT ===")
    print(donors.head(10).to_string())

    Path("results").mkdir(exist_ok=True)
    df.to_csv("results/reference_audit.csv", index=False)
    donors.to_csv("results/reference_donors.csv")
    print("\n-> results/reference_audit.csv")
    print("-> results/reference_donors.csv")

if __name__ == "__main__":
    main()
