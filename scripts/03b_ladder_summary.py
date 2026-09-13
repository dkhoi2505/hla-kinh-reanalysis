#!/usr/bin/env python3
"""Gop ket qua thang do phu thanh mot bang."""
import re
from pathlib import Path
import pandas as pd

rows = []
for f in sorted(Path("results/ladder").glob("*_result.tsv")):
    # ten file dang: SRR11212911_20000_s100_result.tsv
    m = re.match(r"(SRR\d+)_(\d+)_s(\d+)_result", f.stem)
    if not m:
        continue
    sample, level, seed = m.group(1), int(m.group(2)), int(m.group(3))

    df = pd.read_csv(f, sep="\t")
    best = df.iloc[0]           # dap an tot nhat
    obj  = df["Objective"]

    # khoang cach diem giua dap an 1 va 2 (%)
    gap = (obj.iloc[0] - obj.iloc[1]) / obj.iloc[0] * 100 if len(obj) > 1 else float("nan")

    rows.append({
        "sample": sample,
        "cap_read": level,
        "genotype": " ".join(str(best[c]) for c in ["A1","A2","B1","B2","C1","C2"]),
        "read_dung": int(best["Reads"]),
        "objective": round(best["Objective"], 1),
        "gap_pct": round(gap, 2),
    })

t = pd.DataFrame(rows).sort_values(["sample", "cap_read"])

# danh dau genotype co doi so voi muc CAO NHAT cua cung mau
t["khop_muc_cao_nhat"] = t.groupby("sample")["genotype"].transform(
    lambda g: g == g.iloc[-1]
)

pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 60)
print(t.to_string(index=False))

t.to_csv("results/ladder_summary.csv", index=False)
print("\n-> results/ladder_summary.csv")
