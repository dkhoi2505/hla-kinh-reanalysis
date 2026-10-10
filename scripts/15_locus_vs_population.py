#!/usr/bin/env python3
"""
Chang 8.4 -- Tach hai gia thuyet canh tranh.  (ban 3)

Cau hoi: do day du lieu huan luyen A >> B >> C la do
  (a) LOCUS BIAS          -- HLA-C it duoc nghien cuu bat ke quan the, hay
  (b) POPULATION INEQUITY -- allele pho bien ngoai chau Au bi thieu du lieu.

Hai gia thuyet long vao nhau trong du lieu Kinh, vi mot allele Kinh dong thoi
(i) thuoc mot locus nao do va (ii) hiem o quan the chau Au. Can quan the doi
chieu de thao go.

Doi chieu: Germany DKMS (AFND pop 3767, n=3,456,066, NGS, 2021).

THIET KE
  Phan 1  Khoa quan the, doc theo locus. Chi allele PHO BIEN O DUC -- day la
          allele ma cac quan the sinh ra du lieu huan luyen dang mang, nen moi
          chenh lech con lai giua ba locus KHONG the quy cho quan the.
  Phan 2  Khoa locus, doc theo quan the. Trong tung locus: pho bien o Duc
          vs pho bien o Kinh ma hiem o Duc.
  Phan 3  Tuong quan Spearman AF Duc <-> so phep do BA, trong tung locus.
  Phan 4  Do nhay theo nguong Kinh (k >= 1, 2, 3).
  Phan 5  Chan doan: tach "BA=0 vi vang khoi bang" khoi "BA=0 that".

================================================================= NGUONG
Hai quan the duoc cat bang HAI LOAI phep toan khac nhau, co chu dich.

KINH -- cat bang SO LAN QUAN SAT k, khong bang AF.
  n=101 => 202 nhiem sac the. AF = k/202, k nguyen.
  => AF bi LUONG TU HOA: chi 203 gia tri kha di, cac bac cach nhau 0.00495.
     Khong co AF nao nam giua 0.00495 (k=1) va 0.00990 (k=2).
  => Mot nguong AF tro thanh HAM BAC THANG: ca mot dai nguong cho ket qua
     y het nhau, roi nhay bac.
         nguong 0.0001 - 0.00495  ->  k>=1  ->  85 allele
         nguong 0.0050 - 0.00990  ->  k>=2  ->  63 allele
         nguong 0.0099 - 0.01485  ->  k>=3  ->  50 allele

  BAY (da sap HAI lan trong du an nay):
    Ban 1 dat nguong 0.005 va goi do la "gioi han phat hien = 1/202".
    Nhung 1/202 = 0.004950495, con bai bao ghi 0.00495 -- NHO HON 1/202.
    Nen nguong 0.005 loai sach 22 allele singleton. Sai 1% ve so, 26% ve mau.
    Ban 2 dat nguong do nhay = 2/202 va dan nhan "(loai singleton)".
    Cung bay: 2/202 = 0.009900990 > 0.00990 trong file, nen no loai ca k=2.
    Thuc te no do k>=3, khong phai k>=2. Nhan SAI.

  SUA: cat bang k = round(af_kinh * N_CHROM). So nguyen so voi so nguyen,
  khong co epsilon, mien nhiem voi lam tron cua bai bao, va tu giai thich.

DUC -- cat bang AF 0.005, va goi dung ten: QUY UOC 0.5%.
  n=3,456,066 => 6,912,132 nhiem sac the, bac cach nhau 1.4e-7.
  O quy mo do AF thuc te la LIEN TUC, nen nguong AF la phep toan dung.
  0.005 KHONG phai gia tri dan xuat tu thiet ke nghien cuu -- no la quy uoc
  pho bien trong di truyen quan the HLA de dinh nghia "allele thuong gap".
  Script in do phu cua tap nay de dua vao phan han che.
=============================================================================

HAI FILE AF DUC
  data/germany_af.tsv        chi allele AF >= 0.005 (65 dong). Doi chieu.
  data/germany_af_full.tsv   moi allele AF > 0 (1,143 dong). Nguon chinh.
  Nguong KHONG duoc nuong vao file; script tu cat. Gan af_duc = 0 cho allele
  duoi nguong la bia so lieu (A*33:03 that la 0.0037, khong phai 0).

Chay trong env co pandas (vd: conda activate mhc):
    python scripts/15_locus_vs_population.py
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"

N_MAU = 101                  # so ca the trong Do et al. 2020
N_CHROM = 2 * N_MAU          # 202 nhiem sac the moi locus
K_MIN = 1                    # Kinh: allele quan sat duoc it nhat 1 lan
K_DO_NHAY = [2, 3]           # Phan 4
NGUONG_DUC = 0.005           # quy uoc 0.5%, KHONG phai gia tri dan xuat
LOCI = ["A", "B", "C"]


# ----------------------------------------------------------------- doc du lieu
def doc_af(path, ten, bat_buoc=True):
    if not path.exists():
        if bat_buoc:
            sys.exit(f"Thieu file bat buoc: {path}")
        return None
    df = pd.read_csv(path, sep="\t", header=None, names=["allele", "af"])
    df["allele"] = df.allele.str.strip()
    t = df.allele.str.split(":")
    df["allele2"] = t.str[0] + ":" + t.str[1]
    out = df.groupby("allele2", as_index=False)["af"].sum()
    out["locus"] = out.allele2.str[0]
    out = out[out.locus.isin(LOCI)].copy()
    print(f"{ten:10s}: {len(df):5d} dong -> {len(out):5d} allele 2-truong"
          f"  (gop {len(df)-len(out)})  tong AF {out.af.sum():.5f}")
    return out.rename(columns={"af": f"af_{ten}"})


def doc_counts():
    path = RES / "mhcflurry_training_counts.csv"
    if not path.exists():
        sys.exit(f"Thieu file: {path}\nChay scripts/11_mhcflurry_counts.py truoc.")
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}

    def tim(*uv):
        for v in uv:
            if v in cols:
                return cols[v]
        return None

    c_al = tim("allele", "allele2", "allele_2field") or df.columns[0]
    c_ba = tim("ba", "n_ba", "ba_count", "count_ba", "ba_peptides")
    c_el = tim("el", "n_el", "el_count", "count_el", "el_peptides")
    if c_ba is None:
        sys.exit(f"Khong tim thay cot BA. Cac cot: {list(df.columns)}")
    out = pd.DataFrame({
        "allele2": df[c_al].astype(str).str.strip(),
        "BA": pd.to_numeric(df[c_ba], errors="coerce").fillna(0).astype(int),
        "EL": (pd.to_numeric(df[c_el], errors="coerce").fillna(0).astype(int)
               if c_el else 0),
    })
    out = out.groupby("allele2", as_index=False)[["BA", "EL"]].sum()
    print(f"{'counts':10s}: {len(out):5d} allele co ban ghi MHCflurry "
          f"(BA='{c_ba}', EL='{c_el}')")
    return out


# ----------------------------------------------------------------- thong ke
def spearman(x, y):
    """Tuong quan hang Spearman, khong can scipy."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def khoi_dong_hang(v):
    """Kich thuoc khoi dong hang lon nhat, tinh theo % mau.
    Dong hang lam giam phuong sai cot hang -> keo |Spearman| ve 0."""
    v = pd.Series(np.asarray(v))
    if len(v) == 0:
        return 0, 0.0
    n = int(v.value_counts().iloc[0])
    return n, 100.0 * n / len(v)


def mo_ta(v):
    v = np.asarray(v)
    if len(v) == 0:
        return "n=  0"
    return (f"n={len(v):3d}  trung vi={np.median(v):8.0f}  "
            f"p25={np.percentile(v,25):7.0f}  p75={np.percentile(v,75):8.0f}  "
            f"BA=0: {int((v==0).sum())}/{len(v)}")


def ty_le(a, b):
    ma, mb = float(np.median(a)), float(np.median(b))
    if ma > 0 and mb > 0:
        return f"ty le trung vi {ma/mb:.1f}x"
    if ma == 0 and mb == 0:
        return "ca hai trung vi = 0 (khong xac dinh)"
    if mb == 0:
        return f"trung vi nhom sau = 0 (nhom truoc {ma:.0f}); ty le vo cuc"
    return f"trung vi nhom truoc = 0 (nhom sau {mb:.0f}); NGHICH huong"


# ----------------------------------------------------------------- chinh
def main():
    print("=== NGUON DU LIEU ===")
    kinh = doc_af(DATA / "kinh_af.tsv", "kinh")
    duc = doc_af(DATA / "germany_af.tsv", "duc")
    duc_full = doc_af(DATA / "germany_af_full.tsv", "duc_full")
    counts = doc_counts()

    df = kinh.merge(duc_full, on=["allele2", "locus"], how="outer")
    df["locus"] = df.allele2.str[0]
    df["af_kinh"] = df.af_kinh.fillna(0.0)
    df["af_duc"] = df.af_duc_full.fillna(0.0)   # vang trong file day du = AF 0 that

    # Kinh: cat bang SO LAN QUAN SAT, khong bang AF (xem ghi chu NGUONG)
    df["k_kinh"] = (df.af_kinh * N_CHROM).round().astype(int)
    df["pb_kinh"] = df.k_kinh >= K_MIN
    df["pb_duc"] = df.af_duc >= NGUONG_DUC

    df = df.merge(counts, on="allele2", how="left")
    df["co_trong_bang"] = df.BA.notna()
    df[["BA", "EL"]] = df[["BA", "EL"]].fillna(0).astype(int)

    print(f"\n=== NGUONG ===")
    print(f"Kinh: k >= {K_MIN}   (k = round(AF x {N_CHROM}); n={N_MAU} ca the)")
    print(f"      AF bi luong tu hoa, bac cach nhau {1/N_CHROM:.5f}"
          f" -> cat bang k, khong bang AF")
    print(f"Duc : AF >= {NGUONG_DUC}   (QUY UOC 0.5%, khong phai gia tri dan xuat)")
    print(f"      6,912,132 nhiem sac the -> bac cach nhau 1.4e-07,"
          f" thuc te lien tuc")

    # kiem tra nhat quan hai file Duc
    n_thr = int((df.af_duc >= NGUONG_DUC).sum())
    print(f"\nKiem tra: {n_thr} allele Duc >= {NGUONG_DUC}"
          f" (germany_af.tsv co {len(duc)})"
          f"  -> {'KHOP' if n_thr == len(duc) else 'LECH -- DUNG LAI'}")

    # Do phu cua tap Duc da cat -- cho phan han che.
    # MAU SO LA 1.0, khong phai tong AFND bao: moi nhiem sac the deu mang mot
    # allele, nen tong THAT phai la 1.0. Chia cho tong AFND bao (~0.995) se am
    # tham hap thu phan hut cua AFND vao phan "giu lai".
    print("\nDo phu cua tap Duc da cat (mau so = 1.0 = toan bo nhiem sac the):")
    for loc in LOCI:
        s = df[df.locus == loc]
        giu = s[s.pb_duc].af_duc.sum()
        bao = s.af_duc.sum()
        print(f"  HLA-{loc}: giu {giu:.5f} -> {100*giu:.1f}% nhiem sac the"
              f"   | AFND bao tong {bao:.5f}, hut {100*(1-bao):.2f}%"
              f" (nguyen nhan chua ro)"
              f"   | bo qua do cat nguong {100*(bao-giu):.1f}%")

    lq = df[df.pb_duc | df.pb_kinh]
    print(f"\nKhung du lieu (AF>0 o it nhat 1 ben): {len(df)}")
    print(f"Tap LIEN QUAN (pho bien o >= 1 quan the): {len(lq)}")
    print(f"  pho bien o ca hai         : {(lq.pb_duc & lq.pb_kinh).sum()}")
    print(f"  chi pho bien o Duc        : {(lq.pb_duc & ~lq.pb_kinh).sum()}")
    print(f"  chi pho bien o Kinh       : {(~lq.pb_duc & lq.pb_kinh).sum()}")
    print(f"  vang khoi bang MHCflurry  : {(~lq.co_trong_bang).sum()}"
          f"   (tren ca khung: {(~df.co_trong_bang).sum()}, phan lon la"
          f" allele sieu hiem khong vao thong ke nao)")

    print(f"\nCai thang AF Kinh (chung minh viec luong tu hoa):")
    kk = df[df.k_kinh > 0]
    for k in sorted(kk.k_kinh.unique())[:6]:
        print(f"  k={k:2d}  AF = {k}/{N_CHROM} = {k/N_CHROM:.6f}"
              f"  ->  {int((kk.k_kinh==k).sum()):3d} allele")
    print(f"  k>=7                                 ->"
          f"  {int((kk.k_kinh>=7).sum()):3d} allele")
    for k in [1] + K_DO_NHAY:
        print(f"  cat o k>={k}: {int((df.k_kinh>=k).sum()):3d} allele")

    # ---------------------------------------------------------------- Phan 1
    print("\n" + "=" * 74)
    print("PHAN 1 -- KHOA QUAN THE, DOC THEO LOCUS  (chi allele pho bien o Duc)")
    print("=" * 74)
    g1 = df[df.pb_duc]
    for loc in LOCI:
        s = g1[g1.locus == loc]
        ba, el = int(s.BA.sum()), int(s.EL.sum())
        ratio = f"{el/ba:.1f}" if ba else "n/a"
        print(f"HLA-{loc}  {mo_ta(s.BA.to_numpy())}")
        print(f"        tong BA = {ba:>9,}   tong EL = {el:>9,}   EL/BA = {ratio}")
    tba, tel = int(g1.BA.sum()), int(g1.EL.sum())
    if tba and tel:
        print("\nTy trong theo locus (BA vs EL -- cung allele, cung quan the,"
              " chi khac cong nghe do):")
        print(f"  {'locus':6s} {'BA':>7s} {'EL':>7s}")
        for loc in LOCI:
            s = g1[g1.locus == loc]
            print(f"  HLA-{loc}  {100*s.BA.sum()/tba:6.1f}%"
                  f" {100*s.EL.sum()/tel:6.1f}%")

    # ---------------------------------------------------------------- Phan 2+3
    def phan_2_3(k_min):
        nhan = {1: "moi allele quan sat duoc",
                2: "loai k=1 (singleton)",
                3: "loai k=1 va k=2"}.get(k_min, "")
        pk = df.k_kinh >= k_min
        print(f"\n--- Kinh: k >= {k_min}  [{nhan}]  "
              f"(AF >= {k_min}/{N_CHROM} = {k_min/N_CHROM:.6f})  "
              f"-> {int(pk.sum())} allele ---")
        for loc in LOCI:
            sub = df[df.locus == loc]
            a = sub[sub.pb_duc].BA.to_numpy()
            b = sub[(sub.k_kinh >= k_min) & ~sub.pb_duc].BA.to_numpy()
            print(f"HLA-{loc}")
            print(f"  pho bien o Duc           : {mo_ta(a)}")
            print(f"  pho bien Kinh / hiem Duc : {mo_ta(b)}")
            if len(a) and len(b):
                print(f"  -> {ty_le(a, b)}")
        s = df[pk]
        print(f"  Spearman AF_duc <-> BA   chung (n={len(s):3d}): "
              f"{spearman(s.af_duc, s.BA):+.3f}")
        for loc in LOCI:
            t = s[s.locus == loc]
            print(f"                           HLA-{loc} (n={len(t):3d}): "
                  f"{spearman(t.af_duc, t.BA):+.3f}")
        nd, pd_ = khoi_dong_hang(s.af_duc)
        nb, pb_ = khoi_dong_hang(s.BA)
        print(f"  Khoi dong hang lon nhat: af_duc {nd} ({pd_:.0f}% mau),"
              f" BA {nb} ({pb_:.0f}% mau)")
        print(f"    -> dong hang keo |Spearman| ve 0, nen cac so tren co the la"
              f" DANH GIA THAP")

    print("\n" + "=" * 74)
    print("PHAN 2 + 3 -- KHOA LOCUS, DOC THEO QUAN THE  (nguong chinh)")
    print("=" * 74)
    phan_2_3(K_MIN)

    print("\n" + "=" * 74)
    print("PHAN 4 -- DO NHAY THEO NGUONG KINH")
    print("Ket luan co doi khi thay nguong khong? Neu khong -> nguong vo hai.")
    print("=" * 74)
    for k in K_DO_NHAY:
        phan_2_3(k)

    # ---------------------------------------------------------------- Phan 5
    print("\n" + "=" * 74)
    print("PHAN 5 -- CHAN DOAN BA = 0")
    print("'vang khoi bang'   = allele khong co trong MHCflurry chut nao.")
    print("'trong bang, BA=0' = co du lieu EL nhung khong co phep do BA.")
    print("=" * 74)
    for ten, mask in [("pho bien o Duc", df.pb_duc),
                      ("pho bien o Kinh", df.pb_kinh)]:
        s = df[mask & (df.BA == 0)]
        vang = int((~s.co_trong_bang).sum())
        print(f"{ten:18s}: {len(s):3d} allele BA=0  "
              f"| vang khoi bang {vang:3d}  | trong bang nhung BA=0 "
              f"{len(s)-vang:3d}")
        trong = s[s.co_trong_bang]
        if len(trong):
            print(f"  (trong bang, BA=0) -> EL trung vi {np.median(trong.EL):.0f}")
        if len(s):
            nhan = []
            for _, r in s.sort_values("af_kinh", ascending=False).head(12).iterrows():
                tag = "vang" if not r.co_trong_bang else f"EL={r.EL}"
                nhan.append(f"{r.allele2}({tag})")
            print("  " + ", ".join(nhan))

    # ---------------------------------------------------------------- xuat
    cot = ["allele2", "locus", "k_kinh", "af_kinh", "af_duc",
           "pb_kinh", "pb_duc", "co_trong_bang", "BA", "EL"]
    out = df.sort_values(["locus", "k_kinh"], ascending=[True, False])[cot]
    dest = RES / "locus_vs_population.csv"
    out.to_csv(dest, index=False)
    print(f"\n-> {dest.relative_to(ROOT)}")

    # ---------------------------------------------------------------- bang
    print("\n" + "=" * 74)
    print("TOP 20 ALLELE KINH")
    print("=" * 74)
    print(f"{'allele':9s} {'k':>3s} {'af_kinh':>8s} {'af_duc':>11s} "
          f"{'ty le':>9s} {'BA':>8s} {'EL':>8s}")
    for _, r in df[df.pb_kinh].nlargest(20, "k_kinh").iterrows():
        ty = f"{r.af_kinh/r.af_duc:,.0f}x" if r.af_duc > 0 else "inf"
        af = f"{r.af_duc:.6f}" if r.af_duc > 0 else "0"
        print(f"{r.allele2:9s} {r.k_kinh:3d} {r.af_kinh:8.5f} "
              f"{af:>11s} {ty:>9s} {r.BA:8,} {r.EL:8,}")


if __name__ == "__main__":
    main()
