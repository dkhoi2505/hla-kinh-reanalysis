#!/usr/bin/env python3
"""
An danh hoa bang genotype truoc khi cong bo.

Thay SRR accession bang ma sample_NNN theo thu tu NGAU NHIEN (co seed).
Bang anh xa luu local, KHONG len repo (xem .gitignore).

Ly do: du lieu SRA la cong khai, nhung gan genotype HLA voi accession
cu the la cong bo mot kieu gen co y nghia lam sang gan voi mot ID.
Phan tich thong ke van tai lap duoc ma khong can lien ket do.

Dau vao : results/genotypes_by_seed.csv
          results/genotypes_consensus.csv
Dau ra  : public/genotypes_by_seed_anon.csv
          public/genotypes_consensus_anon.csv
          data/sample_mapping.csv          <- GIU LOCAL
"""
import sys
from pathlib import Path
import pandas as pd

SEED = 20260913   # doi seed nay thi ma sample doi theo

def main():
    src = Path("results")
    for f in ["genotypes_by_seed.csv", "genotypes_consensus.csv"]:
        if not (src / f).exists():
            sys.exit(f"Khong thay {src/f}")

    raw = pd.read_csv(src / "genotypes_by_seed.csv")
    con = pd.read_csv(src / "genotypes_consensus.csv")

    # xao tron accession roi danh so
    accs = sorted(raw["sample"].unique())
    shuffled = pd.Series(accs).sample(frac=1, random_state=SEED).tolist()
    mapping = {a: f"sample_{i:03d}" for i, a in enumerate(shuffled, start=1)}

    Path("public").mkdir(exist_ok=True)

    for name, df in [("genotypes_by_seed", raw), ("genotypes_consensus", con)]:
        d = df.copy()
        d["sample"] = d["sample"].map(mapping)
        d = d.sort_values("sample")
        out = Path("public") / f"{name}_anon.csv"
        d.to_csv(out, index=False)
        print(f"-> {out}  ({len(d)} dong)")

    m = pd.DataFrame(sorted(mapping.items()), columns=["accession", "ma_an_danh"])
    m.to_csv("data/sample_mapping.csv", index=False)
    print(f"-> data/sample_mapping.csv  (GIU LOCAL, khong commit)")

    # kiem tra khong con accession nao sot lai
    for f in Path("public").glob("*_anon.csv"):
        txt = f.read_text()
        if "SRR" in txt:
            sys.exit(f"LOI: {f} van con chuoi 'SRR'")
    print("\nKiem tra: khong con accession trong file cong bo. OK")

if __name__ == "__main__":
    main()
