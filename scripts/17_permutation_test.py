#!/usr/bin/env python3
"""
Chang 9.5 -- Kiem dinh hoan vi cho cac tuong quan Spearman o Result 2.

VAN DE
  README bao chin gia tri Spearman (+0.573 den +0.815) nhung khong co uoc
  luong bat dinh nao. Muc han che ghi: "No uncertainty estimates on the
  correlations. Group sizes are 13 to 41 per locus." Mot he so +0.8 tren 19
  allele co the hoan toan do ngau nhien -- khong noi duoc neu khong do.

CACH LAM
  Kiem dinh hoan vi (permutation test), khong gia dinh phan phoi:
    1. Tinh rho quan sat duoc.
    2. Xao tron NGAU NHIEN mot trong hai bien N_HOAN_VI lan, moi lan tinh
       lai rho. Day la phan phoi rho khi KHONG co lien he nao -- phan phoi
       rong (null distribution).
    3. p = ti le lan hoan vi cho |rho| >= |rho quan sat|.
  Vi sao hop o day: n nho (19 den 41), co khoi dong hang lon (28/85 allele
  co BA = 0), nen cong thuc tiem can cua Spearman khong dang tin. Hoan vi
  khong gia dinh gi ve phan phoi.

UOC LUONG p THEO KIEU CONG MOT
  p = (1 + so lan vuot) / (1 + N_HOAN_VI).
  Khi khong lan nao vuot, p = 1/10001 = 0.0001 -- KHONG phai 0. Mot kiem
  dinh 10.000 lan khong the chung minh p nho hon 1/10001; viet "p = 0"
  hoac "p < 0.0001" ma khong noi ro la noi qua.

TU KIEM
  Script nay tinh lai rho tu dau. No PHAI tai tao dung bang Spearman cua
  README (script 15). Neu lech => mot trong hai sai, DUNG LAI.

PHU: in co nhom theo locus o ca ba nguong, de kiem cau "Group sizes are
  13 to 41 per locus" trong muc han che.

Dau vao : results/locus_vs_population.csv
Dau ra  : results/permutation_spearman.csv

Chay trong env mhc:  python scripts/17_permutation_test.py
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"

N_CHROM = 202
NGUONG_K = [1, 2, 3]
N_HOAN_VI = 10000
HAT_GIONG = 100          # cung ho seed voi 05_run_all.sh
LOCI = ["A", "B", "C"]

# Bang README (script 15) -- script nay phai tai tao dung
README_RHO = {
    1: {"gop": 0.712, "A": 0.752, "B": 0.762, "C": 0.815},
    2: {"gop": 0.593, "A": 0.666, "B": 0.689, "C": 0.770},
    3: {"gop": 0.573, "A": 0.605, "B": 0.578, "C": 0.775},
}


def xep_hang(v):
    """Hang trung binh cho gia tri dong hang (giong phuong phap 'average')."""
    v = np.asarray(v, dtype=float)
    n = len(v)
    tt = np.argsort(v, kind="mergesort")
    h = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[tt[j + 1]] == v[tt[i]]:
            j += 1
        h[tt[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return h


def pearson(a, b):
    a = a - a.mean()
    b = b - b.mean()
    da, db = np.sqrt((a * a).sum()), np.sqrt((b * b).sum())
    if da == 0 or db == 0:
        return float("nan")
    return float((a * b).sum() / (da * db))


def spearman(x, y):
    return pearson(xep_hang(x), xep_hang(y))


def khoi_dong_hang(v):
    """Kich thuoc khoi dong hang lon nhat, tinh theo % mau."""
    v = pd.Series(np.asarray(v))
    if len(v) == 0:
        return 0, 0.0
    n = int(v.value_counts().iloc[0])
    return n, 100.0 * n / len(v)


def hoan_vi(x, y, n_lan, rng):
    """Tra ve (rho quan sat, mang rho cua phan phoi rong).

    Chi xao tron VECTO HANG, khong xao tron gia tri goc: hang la bat bien
    qua phep hoan vi cua chinh tap gia tri do, nen hai cach tuong duong
    nhau ve mat toan hoc nhung cach nay nhanh hon nhieu.
    """
    rx, ry = xep_hang(x), xep_hang(y)
    rho = pearson(rx, ry)
    rx_c = rx - rx.mean()
    chuan_x = np.sqrt((rx_c * rx_c).sum())
    ry2 = ry.copy()
    out = np.empty(n_lan)
    for i in range(n_lan):
        rng.shuffle(ry2)
        b = ry2 - ry2.mean()
        db = np.sqrt((b * b).sum())
        out[i] = (rx_c * b).sum() / (chuan_x * db) if chuan_x and db else np.nan
    return rho, out


def doc():
    p = RES / "locus_vs_population.csv"
    if not p.exists():
        sys.exit(f"Thieu {p} -- chay scripts/15_locus_vs_population.py truoc")
    df = pd.read_csv(p)
    c = {x.lower().strip(): x for x in df.columns}
    need = {}
    for ten, ung in (("locus", ("locus",)),
                     ("af_duc", ("af_duc", "af_duc_chinh_xac", "af_germany")),
                     ("BA", ("ba", "n_ba")),
                     ("k_kinh", ("k_kinh", "k")),
                     ("af_kinh", ("af_kinh",))):
        for u in ung:
            if u in c:
                need[ten] = c[u]
                break
    for ten in ("locus", "af_duc", "BA"):
        if ten not in need:
            sys.exit(f"Khong tim thay cot {ten}. Cac cot: {list(df.columns)}")

    out = pd.DataFrame({
        "locus": df[need["locus"]].astype(str).str.strip(),
        "af_duc": pd.to_numeric(df[need["af_duc"]], errors="coerce").fillna(0.0),
        "BA": pd.to_numeric(df[need["BA"]], errors="coerce").fillna(0).astype(int),
    })
    if "k_kinh" in need:
        out["k"] = pd.to_numeric(df[need["k_kinh"]], errors="coerce").fillna(0).astype(int)
    elif "af_kinh" in need:
        out["k"] = (pd.to_numeric(df[need["af_kinh"]], errors="coerce").fillna(0)
                    * N_CHROM).round().astype(int)
    else:
        sys.exit("Khong tim thay cot k_kinh hay af_kinh")
    out = out[out.locus.isin(LOCI)]
    print(f"doc {len(out)} dong, {len(out[out.k >= 1])} allele co k >= 1")
    return out


def main():
    print("=== NGUON DU LIEU ===")
    d = doc()
    rng = np.random.default_rng(HAT_GIONG)
    print(f"\n{N_HOAN_VI} hoan vi moi o, hat giong {HAT_GIONG}")
    print("(chay mat khoang 30-60 giay)\n")

    hang_ra, lech = [], []
    print("=" * 78)
    print(f"{'nguong':>7s} {'nhom':>5s} {'n':>4s} {'rho':>8s} {'README':>8s} "
          f"{'p (2 phia)':>11s} {'95% phan phoi rong':>22s} {'dong hang':>10s}")
    print("=" * 78)

    for K in NGUONG_K:
        s = d[d.k >= K]
        for nhan in ["gop"] + LOCI:
            g = s if nhan == "gop" else s[s.locus == nhan]
            if len(g) < 4:
                print(f"{K:>7d} {nhan:>5s} {len(g):>4d}   (qua it, bo qua)")
                continue
            rho, rong = hoan_vi(g.af_duc.values, g.BA.values, N_HOAN_VI, rng)
            vuot = int(np.sum(np.abs(rong) >= abs(rho) - 1e-12))
            p = (1 + vuot) / (1 + N_HOAN_VI)
            lo, hi = np.percentile(rong, [2.5, 97.5])
            nd, npct = khoi_dong_hang(g.BA.values)

            kv = README_RHO[K][nhan]
            if abs(round(rho, 3) - kv) > 0.0015:
                lech.append((K, nhan, round(rho, 3), kv))
            dau = "" if abs(round(rho, 3) - kv) <= 0.0015 else "  <-- LECH"

            print(f"{K:>7d} {nhan:>5s} {len(g):>4d} {rho:>+8.3f} {kv:>+8.3f} "
                  f"{p:>11.4f} [{lo:>+7.3f},{hi:>+7.3f}] {npct:>9.1f}%{dau}")
            hang_ra.append({
                "nguong_k": K, "nhom": nhan, "n": len(g),
                "rho": round(rho, 4), "rho_readme": kv,
                "p_hai_phia": round(p, 5), "so_lan_vuot": vuot,
                "n_hoan_vi": N_HOAN_VI,
                "rong_2.5pct": round(float(lo), 4),
                "rong_97.5pct": round(float(hi), 4),
                "khoi_dong_hang_pct": round(npct, 1),
            })
        print("-" * 78)

    # ------------------------------------------------- co nhom theo locus
    print("\n=== CO NHOM THEO LOCUS (de kiem cau 'Group sizes are 13 to 41') ===")
    print(f"{'nguong':>7s} " + " ".join(f"{'HLA-'+l:>7s}" for l in LOCI)
          + f" {'gop':>6s}")
    tat_ca_n = []
    for K in NGUONG_K:
        s = d[d.k >= K]
        ns = [int((s.locus == l).sum()) for l in LOCI]
        tat_ca_n += ns
        print(f"{K:>7d} " + " ".join(f"{x:>7d}" for x in ns) + f" {len(s):>6d}")
    print(f"\n  co nhom theo locus nho nhat = {min(tat_ca_n)}, "
          f"lon nhat = {max(tat_ca_n)}")
    print(f"  README dang ghi: 'Group sizes are 13 to 41 per locus'")
    if min(tat_ca_n) == 13 and max(tat_ca_n) == 41:
        print("  -> KHOP")
    else:
        print(f"  -> LECH: phai la '{min(tat_ca_n)} to {max(tat_ca_n)}'")

    # ------------------------------------------------- doi chieu
    print("\n=== DOI CHIEU voi bang Spearman cua README (script 15) ===")
    if lech:
        print("  >>> LECH, DUNG LAI -- mot trong hai script sai:")
        for K, nhan, co, can in lech:
            print(f"      nguong {K}, {nhan}: tinh ra {co:+.3f}, README {can:+.3f}")
        print("  Khong ghi file ket qua.")
        sys.exit(1)
    print("  -> KHOP ca 12 o. Logic script nay trung script 15, nen cac gia")
    print("     tri p va khoang phan phoi rong dang tin.")

    dest = RES / "permutation_spearman.csv"
    pd.DataFrame(hang_ra).to_csv(dest, index=False)
    print(f"\n-> {dest.relative_to(ROOT)}")

    print("\nDOC KET QUA:")
    print("  p nho = rho quan sat nam ngoai pham vi ma ngau nhien tao ra duoc.")
    print("  Khoang 95% cua phan phoi rong cho biet ngau nhien 'du suc' tao ra")
    print("  rho lon den dau o co nhom DO. Neu rho quan sat nam xa ngoai khoang")
    print("  do thi lien he khong phai do co nhom nho.")
    print("  LUU Y: hoan vi kiem GIA THUYET RONG, khong cho khoang tin cay cua")
    print("  chinh rho. No tra loi 'co phai ngau nhien khong', KHONG tra loi")
    print("  'rho that nam trong khoang nao'.")


if __name__ == "__main__":
    main()
