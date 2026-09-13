# HLA Kinh reanalysis — nhat ky

## Moi truong
- WSL2, Ubuntu, 8 core, RAM 3.7GiB (da nang .wslconfig: 5GB + swap 8GB)
- conda env `hla`: optitype 1.5.0, razers3 3.5.12, glpk 5.0, samtools, seqtk, fastp 1.3.6
- Python 3.13.15 — chay duoc, khong xung dot

## Du lieu
- PRJNA609593, 101 run, 1 run/nguoi, library_strategy=AMPLICON
- Tong 10.64 GB nen
- Read goc 151bp, da trim truoc khi nop SRA (80-86% con >=150bp)
- Q40 = 0% -> quality da bi binning, khong phai gia tri may goc
- Insert size 211-266bp — KHONG khop mo ta "~2kb" trong paper
  -> suy doan: co buoc tagmentation trung gian (CHUA xac minh)

## Pilot: 3 mau theo min/median/max read count
- SRR11212911: 96,030 read (min)
- SRR11212875: 603,077 read (median)
- SRR11212960: 2,872,506 read (max)

## Ket qua SRR11212911
Genotype GIONG NHAU o ca 2 muc do phu:
  A*24:03, A*34:01 | B*15:25, B*15:35 | C*04:03, C*07:02

| muc      | read dung | thoi gian | RAM dinh |
|----------|-----------|-----------|----------|
| 20k cap  | 1,253     | 1m03s     | 1.43 GB  |
| toan bo  | 5,643     | 15m25s    | 4.44 GB  |

-> Chi phi PHI TUYEN: read x4.6 nhung thoi gian x14.7
-> Subsample la BAT BUOC cho 101 mau

## Phat hien: vi sao chi ~3% read dung duoc
- BAM co 7,059,824 alignment / 9,088 read duy nhat
- TB 777 allele/read, min 1, max 2823 (tren 11,179 allele)
- Da hinh HLA tap trung o exon 2-3; du lieu DNA phu ca intron
  -> phan lon read roi vao vung bao ton, khong phan biet duoc allele
- => Tang read co loi ich giam dan nhanh. Day la thuoc tinh bai toan,
     khong phai loi cong cu.

## Con lai
- [ ] Chay 2 mau pilot con lai
- [ ] Thi nghiem downsampling day du (Chang 4)
- [ ] Xac minh tagmentation qua tai lieu TruSight HLA

## Chang 4 — thi nghiem on dinh

### Thang do phu (3 mau x 5 muc, seed 100)
- SRR11212911, SRR11212960: on dinh tu 2,000 cap (~116-136 read dung)
- SRR11212875: on dinh tu 5,000 cap (152 read dung)
- gap_pct KHONG tang theo do phu -> loai bo lam thang do tin cay

### Thu seed (20k cap, seed 100/200/300)
- SRR11212911, SRR11212960: 3 seed GIONG NHAU
- SRR11212875: PHAN KY o C*06:02 vs C*06:03
- Tang len 40k, 60k cap: VAN phan ky, con DAO CHIEU
  (60k: seed 100 -> C*06:02; seed 200,300 -> C*06:03)
- Chenh diem giua 2 dap an: chi 0.10-0.32%

### Nguyen nhan (da xac minh tren file tham chieu)
- Tham chieu OptiType: 6,889 allele dich -> 11,179 muc (no 1.62x)
- 10,779/11,179 muc (96.4%) dung intron MUON tu allele ho hang
- Chi ~400 allele co trinh tu gen THAT
- C*06:02 co 2 muc voi intron that; C*06:03 KHONG co muc nao
  -> intron cua C*06:03 trong tham chieu CHINH LA intron cua C*06:02
  -> read roi vao intron khop deu ca hai, khong phan biet duoc
- Du lieu amplicon DNA phu toan gene -> phan lon read roi vao intron
- => Mot so cap allele khong tach duoc, BAT KE do phu

### He qua
- Khong the dat 1 nguong subsample "an toan" cho moi mau
- Phai chay MOI mau voi nhieu seed va danh dau call khong on dinh

## Buoc 4.5 — kiem tra allele Kinh trong tham chieu

Gia thuyet cua Claude: allele pho bien o nguoi Kinh thieu trinh tu gen
that hon allele chau Au -> thien lech quan the o tang tham chieu.

KET QUA: BAC BO.
- 19/19 allele class I pho bien o nguoi Kinh deu CO trinh tu gen that
- KHONG allele Kinh nao di muon intron tu allele 2-field khac
- Nguoc lai, tat ca deu la DONOR: A*24:02 cho 197 allele muon,
  C*07:02 cho 135, C*04:01 cho 129, B*40:01 cho 120, A*11:01 cho 116

Dien giai: day la hieu ung PHO BIEN, khong phai hieu ung dia ly.
Allele pho bien -> duoc giai trinh tu gen som -> thanh khung intron
cho cac allele hiem cung nhom.

[SUY LUAN - can kiem chung o Chang 5]
Rui ro gọi sai tap trung o ca the mang ALLELE HIEM, khong phai allele
pho bien. Vi du: SRR11212911 mang A*34:01 (hiem o DNA);
SRR11212875 mang C*06:03 (khong co gen that) -> chinh la ca phan ky.
Kiem chung: sau khi chay 101 mau x 3 seed, xem call khong hoi tu co
tap trung o allele hiem khong.

## Chang 5 — 101 mau x 3 seed (303 lan chay, 0 loi)

### On dinh theo locus
HLA-A: 81/101 (80.2%)  <- kem nhat
HLA-B: 93/101 (92.1%)
HLA-C: 93/101 (92.1%)

Phan loai mau: cao 74 | trung binh 19 | thap 8

### Allele gay bat on (do dung: allele co o seed nay, vang o seed khac)
allele hiem (<=3 mau)   : 39.1% (43/110)
allele pho bien (>3 mau):  8.4% (43/512)
-> chenh 4.7 lan

Gay bat on nhieu nhat: A*24:02 (8 lan), A*11:01 (7), A*02:01 (4)
= dung 2 donor intron lon nhat (197 va 116 allele muon)

Allele dao dong 100%: hau het la bien the hiem cung nhom voi donor
  A*11:06/10/19 (nhom A*11:01)
  A*24:08/20/28/29/63 (nhom A*24:02)
  A*02:11/12/44/45 (nhom A*02:01)
  B*40:08/49, B*48:02/08

### KET LUAN MUC TIEU 1
Chuoi nhan qua da khep kin:
tham chieu 96.4% intron muon -> allele hiem muon intron cua allele
pho bien -> du lieu amplicon phu toan gene, read chu yeu o intron ->
allele hiem khong tach duoc -> dao dong 4.7 lan cao hon

Locus A kem on dinh nhat vi mang 2 donor lon nhat.

### LUU Y
"Tin cay cao" = PRECISION (lap lai duoc), KHONG phai ACCURACY.
Accuracy chi do duoc khi so voi ket qua CMB (Chang 6-7).

### Loi phuong phap da mac va sua
Lan dau do "allele co trong mau bat on khong" -> sai, vi allele
cang pho bien cang chac chan True. Phai do "allele co DAO DONG
giua cac seed khong".

## Chang 6 — so tan suat voi Bang 1 (Do et al. 2020)

### Tuong quan (81 allele co o ca hai, muc 2-field)
HLA-A: n=24  Spearman 0.990  Pearson 0.995
HLA-B: n=40  Spearman 0.969  Pearson 0.996
HLA-C: n=17  Spearman 0.984  Pearson 0.999
Gop  : n=81  Spearman 0.983  Pearson 0.997

### Sai khac
Toi co 87 allele, Bang 1 co 85, trung 81.

6 allele CHI toi co (moi cai dung 1 nguoi):
  A*02:11, A*02:12, A*02:45, B*39:02, B*40:49, B*48:08
  -> 6/6 nam trong danh sach allele DAO DONG 100% o Chang 5
  -> DOI CHUNG DOC LAP: day la artefact, khong phai phat hien

4 allele CHI Bang 1 co:
  C*04:82 -> 0 muc trong tham chieu OptiType (khong the goi duoc)
  A*33:01 (20 muc), C*03:17 (11 muc), B*55:18 (2 muc) -> co ma bi sot

### Cap lech doi xung
A*11:01: toi 43 allele, bai 46 (-3)
A*11:02: toi  8 allele, bai  5 (+3)
-> dung 3 nguoi bi chuyen tu allele pho bien sang bien the hiem
A*33:03: toi 24, bai 22 (+2) / A*33:01: toi 0, bai 2 (-2)
-> lech NGUOC chieu -> sai so khong co huong co dinh

### KET LUAN MUC TIEU 1
OptiType tai lap duoc phan bo tan suat allele cua phan mem thuong mai
o muc Spearman 0.98-0.99, du chay tren du lieu khong dung loai thiet ke
va chi dung ~1.4% so read goc.
Toan bo sai khac tap trung o allele hiem, va da duoc du bao truoc boi
thi nghiem on dinh 3-seed.

## Buoc 6.2 — khao sat tham so --beta

### Chieu 1: 3 nguoi goi la A*11:02 dong hop tu
beta 0.001 -> A*11:01 A*11:02  (di hop tu)
beta 0.009 -> A*11:02 x2       (MAC DINH)
beta 0.05  -> A*11:02 x2
=> 3/3 doi, cung huong. Call mac dinh la FALSE HOMOZYGOTE.

Khop voi 3 phep do doc lap deu ra con so 3:
- Bang tan suat: A*11:01 thieu 3, A*11:02 thua 3
- HWE locus A: thua dung 3 ca dong hop tu (14 quan sat vs 11 ky vong)
- Thi nghiem beta: dung 3 nguoi do chuyen ve di hop tu

### HWE dong hop tu (quan sat / ky vong)
HLA-A: 14 / 11.0 = 1.27   <- thua, chi o locus A
HLA-B:  6 /  5.8 = 1.04
HLA-C:  9 / 11.3 = 0.79
Chenh tuyet doi chi 3 ca -> chua du y nghia thong ke voi n=101.
Phat bieu dung: "co xu huong thua dong hop tu o locus A, phu hop voi
co che da mo ta, nhung chua du manh de ket luan."

### Chieu 2 (kiem tra nguoc): 5 mau dong hop tu o allele PHO BIEN
beta 0.001 tao ra 5 FALSE HETEROZYGOTE trong 5 mau:
  C*07:02 x2 -> C*07:02 / C*07:123
  A*24:02 x2 -> A*24:02 / A*24:07
  A*02:07 x2 -> A*02:07 / A*31:01
  A*02:07 x2 -> A*02:07 / A*30:02
  C*01:02 x2 -> C*01:02 / C*02:19
Allele moi deu la allele HIEM -> cung loai false positive o Chang 6.

### KET LUAN
beta cao  -> false homozygote
beta thap -> false heterozygote
KHONG co gia tri beta nao dung cho moi mau.
Ly do: read intron khong phan biet duoc 2 allele -> ti le read giua
chung mat y nghia -> beta la nguong ap len mot dai luong da vo nghia.
=> Van de nam o THAM CHIEU, khong sua duoc bang tham so.

### LUU Y PHUONG PHAP
3 mau A*11:02 duoc gan nhan "TIN CAY CAO" o Chang 5 (3 seed giong nhau)
nhung van SAI. Thi nghiem 3-seed do PRECISION, mu hoan toan truoc sai
lech do tham so.
=> Day la chieu bat on THU HAI, 3-seed khong phat hien duoc.

### Loi phuong phap da mac va sua (ghi de nho)
1. Do "allele co trong mau bat on khong" -> sai, vi allele cang pho bien
   cang chac chan True. Phai do "allele co DAO DONG giua cac seed khong".
2. Dem "nguoi" vs dem "allele" - nguoi dong hop tu mang 2 ban nhung
   chi la 1 nguoi. Luon noi ro don vi dem.

## Chang 9 — viet lai script kiem toan tham chieu

Script 07_reference_audit.py thay cho lenh ad-hoc truoc do.
Khac biet: CHI tinh locus A/B/C (loai E/F/G/H/J/K/L/V - class I khong
co dien va gia gen). Day la pham vi dung cho du an.

SO LIEU CHINH THUC (A/B/C):
  11,047 muc / 6,797 allele dich -> he so no 1.63
  10,713 muc dung intron muon = 97.0%
  Chi 334 muc (3.0%) co trinh tu gen that
  Theo locus: A 97.4% | B 95.8% | C 97.7%

(So cu 11,179 / 6,889 / 96.4% la tinh ca pseudogene -> khong dung nua)

DONOR LON NHAT (dem so allele 2-field di muon):
  A*02:01   263
  A*02:43N  202   <- allele NULL van la donor lon thu 2
  A*24:02   198
  B*07:02   151
  A*03:01   145
  C*07:01   145 / C*07:06 139 / C*07:02 136  <- nhom C*07 tong 420

Nhan xet: tham chieu xay theo tieu chi "co du lieu intron", KHONG theo
tieu chi quan trong sinh hoc. Mot allele null (A*02:43N) van lam khung
xuong cho 202 allele khac.
