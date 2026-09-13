#!/usr/bin/env python3
"""
Phan tich do on dinh cua ket qua goi HLA.

Ba cau hoi:
  1. Moi mau chay 3 seed - locus nao on dinh, locus nao khong?
  2. Allele nao GAY bat on (co o seed nay, vang o seed khac)?
  3. So dong hop tu quan sat co lech so voi ky vong Hardy-Weinberg?

LUU Y: day la thuoc do PRECISION (lap lai duoc), KHONG phai ACCURACY.
Mot cong cu co thien lech he thong se on dinh o moi seed ma van sai.

Dau vao : results/genotypes_by_seed.csv  (tu 05_gather.py)
Dau ra  : results/stability_by_allele.csv
          results/hwe_homozygosity.csv
"""
import sys
from collections import Counter
from pathlib import Path
import pandas as pd

LOCI = ["A", "B", "C"]
NGUONG_HIEM = 3          # allele xuat hien <= 3 lan coi la hiem


def doc_du_lieu():
    f = Path("results/genotypes_by_seed.csv")
    if not f.exists():
        sys.exit(f"Khong thay {f}. Chay 05_gather.py truoc.")
    return pd.read_csv(f)


def on_dinh_theo_locus(raw):
    """Moi (mau, locus): 3 seed co cho cung genotype khong?"""
    out = []
    for s, g in raw.groupby("sample"):
        r = {"sample": s, "n_seed": len(g)}
        n_bat_on = 0
        for loc in LOCI:
            r[f"{loc}_on_dinh"] = (g[loc].nunique() == 1)
            n_bat_on += 0 if g[loc].nunique() == 1 else 1
        r["so_locus_bat_on"] = n_bat_on
        out.append(r)
    return pd.DataFrame(out)


def allele_gay_bat_on(raw):
    """
    Allele GAY bat on = co mat o mot so seed nhung vang o seed khac.

    Phan biet voi allele chi TINH CO nam trong mau bat on:
      seed 100: {C*06:02, C*08:01}
      seed 300: {C*06:03, C*08:01}
      -> C*06:02, C*06:03 gay bat on
      -> C*08:01 co o moi seed, KHONG gay bat on
    """
    xuat_hien, gay_bat_on = Counter(), Counter()
    for _, g in raw.groupby("sample"):
        for loc in LOCI:
            tap = [set(x.split()) for x in g[loc]]
            hop, giao = set().union(*tap), set.intersection(*tap)
            for a in hop:
                xuat_hien[a] += 1
                if a not in giao:
                    gay_bat_on[a] += 1

    df = pd.DataFrame({"xuat_hien": pd.Series(xuat_hien),
                       "gay_bat_on": pd.Series(gay_bat_on)}).fillna(0)
    df["gay_bat_on"] = df.gay_bat_on.astype(int)
    df["ti_le_%"] = (100 * df.gay_bat_on / df.xuat_hien).round(1)
    df["nhom"] = ["hiem" if x <= NGUONG_HIEM else "pho_bien" for x in df.xuat_hien]
    return df.sort_values(["gay_bat_on", "ti_le_%"], ascending=False)


def hwe_dong_hop_tu(raw):
    """
    So nguoi dong hop tu quan sat vs ky vong theo Hardy-Weinberg.
    Ky vong = n * sum(p_i^2), voi p_i la tan suat allele i.
    Ty le > 1 = thua dong hop tu.
    Dung genotype cua seed dau tien (seed nho nhat) cho moi mau.
    """
    seed0 = raw.seed.min()
    d = raw[raw.seed == seed0]
    out = []
    for loc in LOCI:
        geno = [g.split() for g in d[loc]]
        n = len(geno)
        quan_sat = sum(1 for g in geno if g[0] == g[1])
        dem = Counter(a for g in geno for a in g)
        p = {a: c / (2 * n) for a, c in dem.items()}
        ky_vong = n * sum(v ** 2 for v in p.values())
        out.append({"locus": loc, "n_mau": n,
                    "dong_hop_quan_sat": quan_sat,
                    "dong_hop_ky_vong": round(ky_vong, 1),
                    "ty_le": round(quan_sat / ky_vong, 2),
                    "chenh_tuyet_doi": round(quan_sat - ky_vong, 1)})
    return pd.DataFrame(out)


def main():
    raw = doc_du_lieu()
    print(f"Doc {len(raw)} ket qua, {raw['sample'].nunique()} mau, "
          f"{raw['seed'].nunique()} seed\n")

    t = on_dinh_theo_locus(raw)
    print("=== ON DINH THEO LOCUS ===")
    for loc in LOCI:
        n = t[f"{loc}_on_dinh"].sum()
        print(f"HLA-{loc}: {n}/{len(t)} mau on dinh ({100*n/len(t):.1f}%)")
    print("\nSo locus bat on moi mau:")
    print(t.so_locus_bat_on.value_counts().sort_index().to_string())

    df = allele_gay_bat_on(raw)
    print("\n=== ALLELE GAY BAT ON: HIEM vs PHO BIEN ===")
    for nhom in ["hiem", "pho_bien"]:
        s = df[df.nhom == nhom]
        if len(s):
            print(f"{nhom:9s} (xuat hien {'<=' if nhom=='hiem' else '>'}{NGUONG_HIEM} lan): "
                  f"{s.gay_bat_on.sum()}/{s.xuat_hien.sum()} = "
                  f"{100*s.gay_bat_on.sum()/s.xuat_hien.sum():.1f}%")

    print("\n=== 10 ALLELE GAY BAT ON NHIEU NHAT ===")
    print(df[df.gay_bat_on > 0].head(10).to_string())

    h = hwe_dong_hop_tu(raw)
    print("\n=== DONG HOP TU: QUAN SAT vs KY VONG (HWE) ===")
    print(h.to_string(index=False))
    print("\nLuu y: chenh tuyet doi vai ca voi n~100 CHUA du y nghia thong ke.")

    Path("results").mkdir(exist_ok=True)
    df.to_csv("results/stability_by_allele.csv")
    h.to_csv("results/hwe_homozygosity.csv", index=False)
    print("\n-> results/stability_by_allele.csv")
    print("-> results/hwe_homozygosity.csv")


if __name__ == "__main__":
    main()
