#!/usr/bin/env python3
"""
Chang 8.6 -- Do nhay cua khoang cach pseudosequence theo nguong NGUONG_BA.

VAN DE
  scripts/14_pseudoseq_distance.py coi mot allele la "da duoc dac trung" khi
  no co >= 50 phep do BA. Con so 50 la TUY TIEN -- khong co can cu nao
  (xem docstring cua script 14). README trich ba phan tram khoang cach
  (100% / 87.5% / 73.6% nam trong vong 2 thay the) ma ba so do PHU THUOC
  vao cai 50 do.

CAU HOI
  Neu doi nguong thanh 25 / 100 / 200 thi ket luan "phep ngoai suy la ngan"
  con dung khong?

LUU Y VE THIET KE
  Khi nguong T tang, HAI tap deu doi cung luc:
    - tap "ngheo"      (BA <  T) RONG ra
    - tap "da dac trung"(BA >= T) HEP lai
  Ca hai deu lam khoang cach TANG. Nen T=200 la phep thu KHAC NGHIET nhat.
  Neu ket luan van dung o T=200 thi nguong vo hai.

TU KIEM
  Script nay KHONG dung lai logic cua script 14; no tinh lai tu dau. Nhung
  o T=50 no PHAI tai tao dung ba con so cua script 14. Dong "T=50" duoc danh
  dau [DOI CHIEU]. Neu lech => mot trong hai script sai, DUNG LAI.

Dau vao : data/kinh_af.tsv
          results/mhcflurry_training_counts.csv
          class1_pseudosequences.csv (tu mhcflurry-downloads)
Dau ra  : results/pseudoseq_sensitivity.csv

Chay trong env mhc:  python scripts/16_pseudoseq_sensitivity.py
"""

from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"

N_CHROM = 202            # 101 ca the x 2 nhiem sac the
K_MIN = 1                # allele Kinh quan sat duoc it nhat 1 lan
NGUONG = [25, 50, 100, 200]
DOI_CHIEU = 50           # nguong script 14 dang dung
LOCI = ["A", "B", "C"]
BIN = [(0, 2, "rat gan"), (3, 5, "gan"), (6, 9, "xa"), (10, 99, "rat xa")]


def doc_pseudoseq():
    """Doc class1_pseudosequences.csv tu mhcflurry-downloads.
    Tu phat hien ten cot va dang ten allele (HLA-A02:01 hoac HLA-A*02:01)."""
    try:
        d = subprocess.run(
            ["mhcflurry-downloads", "path", "allele_sequences"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception as e:
        sys.exit(f"Khong chay duoc mhcflurry-downloads: {e}\n"
                 f"Da 'conda activate mhc' chua?")
    p = Path(d) / "class1_pseudosequences.csv"
    if not p.exists():
        sys.exit(f"Khong thay {p}")
    df = pd.read_csv(p)
    cols = {c.lower().strip(): c for c in df.columns}
    c_al = cols.get("allele") or df.columns[0]
    c_sq = (cols.get("sequence") or cols.get("pseudosequence")
            or cols.get("seq") or df.columns[1])

    out = pd.DataFrame({"ten": df[c_al].astype(str).str.strip(),
                        "seq": df[c_sq].astype(str).str.strip()})
    # chi giu HLA-A/B/C class I co dien, 2 truong
    m = out.ten.str.match(r"^HLA-[ABC]\*?\d+:\d+$")
    out = out[m].copy()
    if out.empty:
        sys.exit("Khong khop duoc ten allele nao trong file pseudosequence.\n"
                 f"Vi du 5 ten dau: {df[c_al].head().tolist()}")
    # chuan hoa ve dang A*02:01
    t = out.ten.str.replace("HLA-", "", regex=False).str.replace("*", "", regex=False)
    out["allele2"] = t.str[0] + "*" + t.str[1:]
    out["locus"] = out.allele2.str[0]

    dai = out.seq.str.len().value_counts()
    if len(dai) > 1:
        print(f"  CANH BAO: pseudosequence co nhieu do dai khac nhau: "
              f"{dai.to_dict()} -- chi giu do dai pho bien nhat")
        out = out[out.seq.str.len() == dai.index[0]]
    print(f"{'pseudoseq':12s}: {len(out):5d} allele A/B/C, do dai {dai.index[0]} aa")
    return out.drop_duplicates("allele2").set_index("allele2")


def doc_counts():
    p = RES / "mhcflurry_training_counts.csv"
    if not p.exists():
        sys.exit(f"Thieu {p}")
    df = pd.read_csv(p)
    cols = {c.lower().strip(): c for c in df.columns}
    c_al = cols.get("allele") or cols.get("allele2") or df.columns[0]
    c_ba = next((cols[k] for k in ("ba", "n_ba", "ba_count", "count_ba")
                 if k in cols), None)
    if c_ba is None:
        sys.exit(f"Khong tim thay cot BA. Cac cot: {list(df.columns)}")
    out = pd.DataFrame({
        "allele2": df[c_al].astype(str).str.strip(),
        "BA": pd.to_numeric(df[c_ba], errors="coerce").fillna(0).astype(int)})
    out = out.groupby("allele2", as_index=False)["BA"].sum()
    print(f"{'counts':12s}: {len(out):5d} allele co ban ghi MHCflurry")
    return out.set_index("allele2")["BA"]


def doc_kinh():
    p = DATA / "kinh_af.tsv"
    if not p.exists():
        sys.exit(f"Thieu {p}")
    df = pd.read_csv(p, sep="\t", header=None, names=["allele", "af"])
    t = df.allele.str.strip().str.split(":")
    df["allele2"] = t.str[0] + ":" + t.str[1]
    out = df.groupby("allele2", as_index=False)["af"].sum()
    out["locus"] = out.allele2.str[0]
    out = out[out.locus.isin(LOCI)]
    out["k"] = (out.af * N_CHROM).round().astype(int)
    out = out[out.k >= K_MIN]
    print(f"{'kinh':12s}: {len(out):5d} allele 2-truong (k >= {K_MIN}), "
          f"tong AF {out.af.sum():.5f}")
    return out.set_index("allele2")


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def chay_mot_nguong(T, kinh, ba, ps):
    """Tra ve (bang phan bo theo locus, so allele khong do duoc)."""
    giau = {loc: [] for loc in LOCI}
    for al, seq in ps.seq.items():
        if ba.get(al, 0) >= T:
            giau[al[0]].append((al, seq))

    hang, khong_do = [], []
    for al, r in kinh.iterrows():
        if ba.get(al, 0) >= T:          # khong ngheo o nguong nay
            continue
        if al not in ps.index:
            khong_do.append(al)
            continue
        seq = ps.loc[al, "seq"]
        ung = giau[r.locus]
        if not ung:
            khong_do.append(al)
            continue
        kc = min(hamming(seq, s) for _, s in ung)
        hang.append({"allele2": al, "locus": r.locus, "af": r.af,
                     "BA": int(ba.get(al, 0)), "kc": kc, "nguong": T})
    return pd.DataFrame(hang), khong_do


def main():
    print("=== NGUON DU LIEU ===")
    ps = doc_pseudoseq()
    ba = doc_counts()
    kinh = doc_kinh()

    tat_ca, tom_tat = [], []
    for T in NGUONG:
        df, khong_do = chay_mot_nguong(T, kinh, ba, ps)
        tat_ca.append(df)
        nhan = "  [DOI CHIEU voi script 14]" if T == DOI_CHIEU else ""
        print("\n" + "=" * 70)
        print(f"NGUONG BA = {T}{nhan}")
        print(f"  allele Kinh ngheo (BA < {T}) do duoc: {len(df)}"
              f"   |  khong do duoc: {len(khong_do)}")
        if khong_do:
            print(f"    (khong co pseudosequence): {', '.join(khong_do[:8])}"
                  + (" ..." if len(khong_do) > 8 else ""))
        print("=" * 70)
        for loc in LOCI:
            s = df[df.locus == loc]
            if s.empty:
                print(f"  HLA-{loc}: khong co allele ngheo")
                continue
            tong = s.af.sum()
            print(f"  HLA-{loc}  (n={len(s):2d}, AF ngheo {tong:.3f})")
            for lo, hi, ten in BIN:
                w = s[(s.kc >= lo) & (s.kc <= hi)].af.sum()
                pct = 100 * w / tong if tong else 0
                if lo == 0:
                    tom_tat.append({"nguong": T, "locus": loc, "n": len(s),
                                    "af_ngheo": round(tong, 4),
                                    "pct_trong_2_thay_the": round(pct, 1)})
                print(f"     kc {lo}-{hi if hi < 99 else '+'}  {ten:8s}: {pct:5.1f}%")

    # ------------------------------------------------------------- tom tat
    tt = pd.DataFrame(tom_tat)
    print("\n" + "=" * 70)
    print("TOM TAT -- % (trong so AF) allele ngheo nam TRONG VONG 2 THAY THE")
    print("=" * 70)
    bang = tt.pivot(index="nguong", columns="locus",
                    values="pct_trong_2_thay_the")
    n_bang = tt.pivot(index="nguong", columns="locus", values="n")
    print(f"{'nguong':>8s}  " + "  ".join(f"{'HLA-'+l:>14s}" for l in LOCI))
    for T in NGUONG:
        o = []
        for l in LOCI:
            v = bang.loc[T, l] if l in bang.columns else float("nan")
            n = n_bang.loc[T, l] if l in n_bang.columns else 0
            o.append(f"{v:8.1f}% (n={int(n):2d})")
        sao = "  <-- script 14" if T == DOI_CHIEU else ""
        print(f"{T:>8d}  " + "  ".join(o) + sao)

    print("\nDOI CHIEU: dong nguong=50 phai khop output cua script 14")
    print("  script 14 bao: HLA-A 100.0%, HLA-B 87.5%, HLA-C 73.6%")
    try:
        got = {l: round(float(bang.loc[DOI_CHIEU, l]), 1) for l in LOCI}
        ky_vong = {"A": 100.0, "B": 87.5, "C": 73.6}
        lech = {l: (got[l], ky_vong[l]) for l in LOCI
                if abs(got[l] - ky_vong[l]) > 0.15}
        if lech:
            print(f"  >>> LECH: {lech}  -- DUNG LAI, mot trong hai script sai")
        else:
            print(f"  -> KHOP ({got}). Logic script nay trung script 14,")
            print(f"     nen cac nguong khac dang tin.")
    except Exception as e:
        print(f"  khong doi chieu duoc: {e}")

    dest = RES / "pseudoseq_sensitivity.csv"
    pd.concat(tat_ca, ignore_index=True).to_csv(dest, index=False)
    print(f"\n-> {dest.relative_to(ROOT)}")

    print("\nDOC KET QUA:")
    print("  Nguong cao hon = phep thu khac nghiet hon (tap ngheo rong ra VA")
    print("  tap da dac trung hep lai cung luc).")
    print("  Neu % o T=200 van cao => ket luan 'phep ngoai suy la ngan' khong")
    print("  phu thuoc vao con so 50, va README co the noi vay.")
    print("  Neu % sut manh => ket luan PHU THUOC nguong, phai viet ro trong")
    print("  README thay vi chi ghi 'chua kiem do nhay'.")


if __name__ == "__main__":
    main()
