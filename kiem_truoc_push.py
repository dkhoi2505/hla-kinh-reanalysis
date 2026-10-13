"""Kiem README truoc khi push. Chi bao cao, khong sua gi.
Chi giu cac phep kiem ma sua_lan8.py CHUA lam:
  1. Moi duong dan README dan toi co TON TAI va duoc GIT THEO DOI khong?
     (loi da mac that: mot commit truoc day README dan toi 3 file chua add)
  2. Script trong repo ma README khong nhac ten -- da bo loc backtick gay nhieu
  3. Doan van bi lap y nguyen (rui ro that sau 11 lan chen chuoi)
  4. Bang markdown lech so cot
"""
import pathlib, re, subprocess, sys
from collections import Counter

t = pathlib.Path("README.md").read_text(encoding="utf-8")
dong = t.split("\n")
loi = []
try:
    theo_doi = set(subprocess.run(["git", "ls-files"], capture_output=True,
                                  text=True, check=True).stdout.split())
except Exception as e:
    theo_doi = None; print(f"(khong chay duoc git ls-files: {e})")

print("--- 1. duong dan README dan toi ---")
dd = sorted(set(re.compile(
    r"`((?:data|scripts|results|docs|public)/[A-Za-z0-9_.\-/]+)`").findall(t)))
for d in dd:
    if not pathlib.Path(d).exists():
        loi.append(f"{d}: FILE KHONG TON TAI"); print(f"  [KHONG CO] {d}")
    elif theo_doi is not None and d not in theo_doi:
        loi.append(f"{d}: GIT CHUA THEO DOI"); print(f"  [CHUA ADD] {d}")
print(f"  {len(dd)} duong dan trong backtick, {len(loi)} co van de")

print("--- 2. script trong repo ma README khong nhac ---")
sot = [f.name for f in sorted(pathlib.Path("scripts").iterdir())
       if f"scripts/{f.name}" not in t] if pathlib.Path("scripts").exists() else []
if sot:
    loi.append(f"README khong nhac: {', '.join(sot)}")
    print(f"  [THIEU] {', '.join(sot)}")
else:
    print("  ok -- moi script deu duoc nhac ten")

print("--- 3. doan bi lap y nguyen ---")
doan = [p.strip() for p in t.split("\n\n")
        if len(p.strip()) > 120 and not p.strip().startswith("|")]
lap = [(v, k[:60]) for k, v in Counter(doan).items() if v > 1]
for v, k in lap:
    loi.append(f"doan xuat hien {v} lan: {k}"); print(f"  [LAP x{v}] {k}...")
if not lap: print(f"  ok -- {len(doan)} doan dai, khong doan nao trung")

print("--- 4. bang markdown ---")
i = nb = 0
while i < len(dong):
    if dong[i].strip().startswith("|") and i + 1 < len(dong) \
            and re.match(r"^\s*\|[\s:|-]+\|\s*$", dong[i + 1]):
        nb += 1; dau = i; cot = []
        while i < len(dong) and dong[i].strip().startswith("|"):
            cot.append(dong[i].strip().strip("|").count("|") + 1); i += 1
        c = Counter(cot)
        if len(c) > 1:
            loi.append(f"bang dong {dau+1} lech cot {dict(c)}")
            print(f"  [LECH COT] dong {dau+1}: {dict(c)}")
    else: i += 1
print(f"  {nb} bang")

print("\n" + "=" * 56)
if loi:
    print(f"CO {len(loi)} VAN DE -- KHONG NEN PUSH:")
    for x in loi: print(f"  - {x}")
    sys.exit(1)
print("TAT CA DAT -- san sang push")
