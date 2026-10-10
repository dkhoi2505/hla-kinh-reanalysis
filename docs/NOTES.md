# HLA Kinh reanalysis - nhat ky

## Moi truong
- WSL2, Ubuntu, 8 core, RAM 3.7GiB (da nang .wslconfig: 5GB + swap 8GB)
- conda env `hla`: optitype 1.5.0, razers3 3.5.12, glpk 5.0, samtools, seqtk, fastp 1.3.6
- Python 3.13.15 - chay duoc, khong xung dot

## Du lieu
- PRJNA609593, 101 run, 1 run/nguoi, library_strategy=AMPLICON
- Tong 10.64 GB nen
- Read goc 151bp, da trim truoc khi nop SRA (80-86% con >=150bp)
- Q40 = 0% -> quality da bi binning, khong phai gia tri may goc
- Insert size 211-266bp - KHONG khop mo ta "~2kb" trong paper
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

## Chang 4 - thi nghiem on dinh

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

## Buoc 4.5 - kiem tra allele Kinh trong tham chieu

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

## Chang 5 - 101 mau x 3 seed (303 lan chay, 0 loi)

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

## Chang 6 - so tan suat voi Bang 1 (Do et al. 2020)

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

## Buoc 6.2 - khao sat tham so --beta

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

## Chang 9 - viet lai script kiem toan tham chieu

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

## Chang 8 — du lieu huan luyen MHCflurry cho allele nguoi Kinh

### 8.1 Dem truc tiep (trung vi co trong so AF)
         BA      tong peptide
HLA-A  4,515     8,746
HLA-B    163     4,271
HLA-C     74     2,858
-> thieu hut tap trung o BA, khong phai thieu du lieu noi chung.
   Locus C co 2,858 peptide nhung chi 74 phep do ai luc (gan het la EL).

### 8.2 Khoang cach pseudosequence (34 aa) toi allele gan nhat co >=50 BA
KHONG allele nao bi co lap. Xa nhat la B*13:02 (kc=6).
   HLA-A: 100%  AF ngheo o kc 0-2
   HLA-B: 87.5% o kc 0-2
   HLA-C: 73.6% o kc 0-2, 26.4% o kc 3-5

Cac ca tuong nghiem trong o 8.1 hoa ra khong:
   B*07:05 (AF .069, BA=0) -> B*07:02, kc=1, BA 5,614
     (khac dung 1 vi tri: D->N o vi tri 22)
   A*29:01 (AF .084, BA=0) -> A*29:02, kc=0, BA 2,665
   A*24:20 (AF .015, BA=0) -> A*24:02, kc=0, BA 4,549
   C*04:82 -> C*04:01, kc=0 (pseudoseq GIONG HET)

=> SUA LAI ket luan 8.1: "BA=0" KHONG dong nghia "mo hinh chua tung thay".
   Phai noi: khong co du lieu truc tiep, nhung pseudosequence gan mot
   allele giau du lieu nen mo hinh khai quat hoa duoc.

### 8.3 BA "chong lung" = max(BA cua no, BA cua allele gan nhat)
         truc tiep -> chong lung
HLA-A      4,515  ->  4,737   (x1.05)
HLA-B        163  ->  1,775   (x10.9)
HLA-C         74  ->    200   (x2.7)

Ganh nang AF co >= 2,000 BA chong lung:
   HLA-A 75.7% | HLA-B 43.6% | HLA-C 0.0%

P75 cua locus C (291) < P25 cua locus B (486) - hai phan bo gan nhu
khong chong lan.

=> Khai quat hoa pan-allele CUU duoc locus B (x10.9) nhung KHONG cuu
   duoc locus C, vi nguon tham chieu cua C cung ngheo.
   C*07:02 (AF .218, hang 2), C*01:02 (.134, hang 4), C*08:01 (.129,
   hang 5) deu nam o vung du lieu mong.

### CHUA KET LUAN DUOC
[SUY LUAN] Mau hinh A >> B >> C rat co the la thien lech theo LOCUS,
khong phai theo QUAN THE. Trong toan bo du lieu HLA cua MHCflurry,
locus C chi chiem 11.6% (A 39%, B 49%).
Can quan the chau Au tham chieu moi phan biet duoc.

### Loi phuong phap da mac
- Nguong "BA < 100" o ban dau la tu dat, khong co nguon. Bo, thay bang
  phan vi co trong so.
- Nguong 50 BA de goi "giau du lieu" qua thap: C*14:02 (244 BA) duoc
  xep ngang A*11:01 (6,931 BA). Day la ly do phai dung BA chong lung
  thay vi chi dung khoang cach.
## Chang 8.4 -- Quan the doi chieu (AFND Germany DKMS)

### Nguon
AFND (allelefrequencies.net), population 3767 "Germany DKMS - German donors"
n = 3,456,066 | phuong phap NGS | nam 2021 | standard=g (Gold only)
URL: hla6006a.asp?hla_locus=<A|B|C>&hla_locus_type=Classical
     &hla_population=3767&standard=g&hla_order=order_2&page=<N>
Luu y: hla_order=order_2 la sap theo TEN allele, khong phai theo tan suat.

Vi sao can quan the doi chieu: trong du lieu Kinh, hai gia thuyet long vao nhau.
Mot allele Kinh dong thoi (i) thuoc mot locus nao do va (ii) hiem o chau Au,
nen khong tach duoc "locus bias" khoi "population inequity" bang du lieu Kinh.

Nguong 0.005 = 1/202 = gioi han phat hien cua n=101 (202 nhiem sac the).
=> Hai quan the duoc cat o cung mot san, khong phai nguong tuy chon.

### Ket qua quet
locus  ban_ghi   AF>0   tong_AF   giu(>=0.005)  tong_AF_giu
A        835      337   0.99543        17         0.96940
B      1,071      474   0.99564        31         0.95670
C        902      332   0.99518        17         0.97880
                                    tong 65 allele | file day du 1,143 dong

Bien nguong (de biet nguong cat o dau co "sat" khong):
  A: thap nhat giu A*30:02 0.0051  | cao nhat bi loai A*66:01 0.0044
  B: thap nhat giu B*44:05 0.0050  | cao nhat bi loai B*45:01 0.0044
  C: thap nhat giu C*12:02 0.0085  | cao nhat bi loai C*16:02 0.0040

### Nam lop tu kiem (quan trong -- day la ly do tin duoc so lieu)
1. So ban ghi khop dung con so site tu khai (835 / 1,071 / 902).
2. 0 allele trung lap trong ca ba locus.
3. Locus A khop 100% voi 17 dong minh trich tay truoc do (tong 0.9694).
   Day la chot kiem soat: quet lai locus da lam tay, khop thi tin phan con lai.
4. germany_af.tsv la tap con dung cua germany_af_full.tsv:
   0 dong thieu, 0 gia tri lech, dung 65 allele >= 0.005.
5. Chay script 15 voi file day du 78 dong va 1,143 dong cho Spearman GIONG HET.
   Dung nhu ky vong: 1,065 dong thua deu ngoai tap Kinh nen khong duoc phep
   doi ket qua. Neu doi => join sai.

### Van de con MO
Tong AF cac dong khac 0 chi dat ~0.995 o ca ba locus, khong phai 1.0
(thieu ~0.45%). Chua ro nguyen nhan. Lam tron cua AFND khong du giai thich
muc nay. KHONG bia ly do. Can hoi hoac tra lai AFND neu dung so nay de
chuan hoa (renormalise).

### Hai bai hoc ky thuat
1. TOOL TIMEOUT KHONG DUNG SCRIPT DANG CHAY.
   Quet locus C bi timeout o 45s; script van chay nen va nap trung trang 5-6.
   Phat hien bang so hoc: store co 400, them 3 trang, thanh 900 (phai la 700).
   Kiem 200 cap trung: 0 cap co AF lech nhau -> dedupe theo ten allele vo hai.
   QUY TAC: kiem trang thai truoc khi chay lai bat ky vong lap GHI du lieu.

2. KHONG NUONG NGUONG VAO FILE DU LIEU.
   Ban dau script gan af_duc = 0 cho allele duoi nguong (vi file chi co
   >= 0.005). Hau qua: A*33:03 co AF Duc THAT la 0.0037 bi ghi thanh 0,
   ty le thanh "inf", va mot khoi ties gia bop meo Spearman. Do la bia so lieu.
   Sua: hai file, nguong cat trong code. Script ghi "<0.005" + "?" va BO QUA
   Phan 3 neu thieu file day du, thay vi tra ve mot con so sai.

### Phat hien
[XAC LAP] Khoang cach du lieu A >> C la hieu ung LOCUS, khong phai quan the.
  Giu tan suat o Duc gan nhu khong doi, khoang cach van con:
    C*07:02  AF Duc 0.1393 (allele C pho bien thu 2 o Duc) ->   194 BA
    A*24:02  AF Duc 0.0911 (IT pho bien hon)               -> 4,549 BA  (23x)
    C*03:04  AF Duc 0.0816                                 ->   316 BA
    C*01:02  AF Duc 0.0351                                 ->    25 BA
    A*11:01  AF Duc 0.0528                                 -> 6,931 BA

[XAC LAP] "Hiem o chau Au" KHONG keo theo "ngheo du lieu" -- doi chung nguoc:
    A*02:03  AF Duc 0.000054 (gan nhu khong ton tai)       -> 6,382 BA
    B*46:01  AF Duc 0.0002                                 -> 1,831 BA
    B*58:01  AF Duc 0.0079                                 -> 3,770 BA

[XAC LAP] A*33:03 (AF Kinh 0.109, hang 7) co AF Duc 0.0037 -- duoi san 0.005.
  Khong phai bo sot dong khi trich tay. Mot allele top-10 cua Kinh gan nhu
  vang mat o quan the Duc.

[SUY LUAN] Thieu hut nam o BA, khong phai o tong du lieu:
    C*01:02 BA 25 / tong 2,789 | C*03:02 BA 0 / tong 1,338
    C*08:01 BA 44 / tong 2,053 | B*38:02 BA 0 / tong 3,170
  Khop lich su: thoi ky assay IC50 tap trung HLA-A va -B; cac nghien cuu
  mono-allelic EL hien dai moi phu deu HLA-C.
  => Thu hep tiep moi lo o 8.1, bang bang chung DOC LAP voi 8.3 (khoang cach
     pseudosequence). Day la lan tu sua thu hai cua chang nay.
  CAN XAC MINH: dinh nghia cot "tong" trong 11_mhcflurry_counts.py.
  Dang gia dinh tong = BA + EL. Neu khong phai thi suy luan nay phai sua.

### File
data/germany_af.tsv           65 dong, AF >= 0.005. Du cho Phan 1 + Phan 2.
data/germany_af_full.tsv   1,143 dong, moi allele AF > 0. Can cho Phan 3.
scripts/15_locus_vs_population.py
results/locus_vs_population.csv   (script tu sinh)
DA XOA: data/germany_af_top20kinh.tsv -- du thua, la tap con cua file day du
        (da kiem 20/20 khop). Hai nguon so cho cung mot thu se gay lech.

### Thiet ke script 15
Phan 1  Khoa quan the, doc theo locus: chi allele pho bien o Duc.
        Day la allele ma cac quan the sinh ra du lieu huan luyen dang mang,
        nen chenh lech con lai giua ba locus khong the quy cho quan the.
Phan 2  Khoa locus, doc theo quan the: trong tung locus, pho bien o Duc
        vs pho bien o Kinh ma hiem o Duc.
Phan 3  Tuong quan Spearman AF Duc <-> so phep do BA, trong tung locus.
        Spearman tinh bang numpy rank (env hla thieu scipy).

## Chang 8.4 -- SUA LOI (ghi sau, de nguyen muc 8.4 o tren lam ho so loi)

Ba muc trong 8.4 o tren phai coi la BI THAY THE boi muc nay.

### SUA 1 -- Ket luan "[XAC LAP] la hieu ung LOCUS, khong phai quan the" la SAI
Phan "locus" dung. Phan "khong phai quan the" SAI.

Loi suy luan: nhay tu "co hai phan vi du (A*02:03, B*46:01 hiem o Duc ma nhieu
du lieu)" sang "vay quan the khong anh huong". Phan vi du chi chung minh mot
yeu to KHONG PHAI nguyen nhan duy nhat; no khong chung minh yeu to do KHONG CO
tac dong trung binh. Phep do tac dong trung binh la Spearman, va no noi nguoc lai.

### KET LUAN DUNG -- ca hai hieu ung deu that va CONG DON
[XAC LAP] Co hieu ung LOCUS (Phan 1, khoa quan the o 65 allele pho bien o Duc):
  locus  n   trung vi BA   BA=0    ty trong BA
  A     17       4,300     0/17       62.3%
  B     31         290     3/31       35.3%
  C     17         177     3/17        2.4%
  Khoang cach trung vi A:C = 24x, tren chinh quan the sinh ra du lieu huan luyen.

[XAC LAP] Co hieu ung QUAN THE (Phan 2 + 3, khoa locus):
  Trong locus A: pho bien o Duc trung vi 4,300 BA  vs  pho bien Kinh/hiem Duc
  trung vi 20 BA  -> 215x.  Cung locus, chi khac quan the.
  Spearman AF_duc <-> BA duong o ca ba locus: A +0.666, B +0.689, C +0.770
  (chung +0.593, n=63).

  Truong hop xau nhat = giao cua hai hieu ung: allele locus C, pho bien o Kinh,
  hiem o Duc -> trung vi BA = 0.

### SUA 2 -- Thieu hut nam o BA, nay da do duoc (truoc chi la [SUY LUAN])
[XAC LAP] Tren 65 allele pho bien o Duc:
  locus  tong BA   tong EL   EL/BA
  A       79,825   158,134     2.0
  B       45,235   245,829     5.4
  C        3,120    87,891    28.2
  HLA-C khong thieu du lieu -- no thieu BA. 87,891 EL ma chi 3,120 BA.
  Khop lich su: thoi ky assay IC50 tap trung A/B; mono-allelic EL moi phu deu C.

### SUA 3 -- Dong "CAN XAC MINH dinh nghia cot tong" -- DA GIAI
[XAC LAP] tong = BA + EL. Kiem 3/3 voi bang 8.3:
  C*01:02  25 + 2,764 = 2,789 (bang 8.3 ghi 2789)
  C*03:02   0 + 1,338 = 1,338 (ghi 1338)
  B*38:02   0 + 3,170 = 3,170 (ghi 3170)

### SUA 4 -- LOI THIET KE NGUONG (quan trong, anh huong moi so o tren)
AF cua Kinh bi luong tu hoa: allele thay k lan co AF = k/202.
1/202 = 0.004950495..., NHUNG bai bao lam tron thanh 0.00495 -- NHO HON 1/202.
=> Nguong >= 0.005 LOAI SACH moi singleton (k=1). San thuc te cua Kinh thanh
   2/202 = 0.0099, trong khi san cua Duc la 0.005.
   Hai quan the KHONG duoc cat o cung mot san -- pha vo chinh ly le bien minh
   cho nguong 0.005. Day la ly do "63 allele" thay vi 85.

Sua trong script ban 2: nguong chinh 0.0049 (ngay duoi 1/202) cho ca hai ben,
nen Kinh nhan ca singleton. Phia Duc khong doi (khong allele Duc nao nam trong
khoang 0.0049-0.0050). Them Phan 4 chay lai o nguong Kinh 0.005 va 0.0099 de
kiem ket luan co phu thuoc lua chon nguong khong.

### SUA 5 -- "90 dong -> 85 allele" KHONG phai loi du lieu
kinh_af.tsv la ten 3 truong; 5 cap gop thanh cung mot ten 2 truong, AF duoc
cong lai. Bang chung dung: tong AF = 2.99992 ~ 3.0 (ba locus x 1.0).
Checkpoint "phai ra 90 allele" la checkpoint dat sai -- checkpoint dung la
tong AF ~ 3.0.

### Bai hoc chung
Phan vi du bac bo tinh TAT YEU, khong bac bo tinh TRUNG BINH. Muon ket luan
"X khong anh huong" thi phai do tuong quan tren toan tap, khong duoc dua vao
vai truong hop nguoc.

## Chang 8.4 -- BAN CHOT (doc cai nay; hai khoi 8.4 o tren la ho so loi)

Hai khoi "Chang 8.4" va "Chang 8.4 -- SUA LOI" o tren duoc giu lai lam vet
loi, KHONG dung de viet README. Khoi nay tu chua du moi thu can thiet.

=============================================================================
## A. KET LUAN
=============================================================================

Cau hoi: du lieu huan luyen peptide-MHC mong di theo huong A >> B >> C.
Vi LOCUS (HLA-C it duoc nghien cuu) hay vi QUAN THE (allele ngoai chau Au
bi bo roi)? Hai thu nay long vao nhau trong du lieu Kinh -- mot allele Kinh
dong thoi thuoc mot locus VA hiem o chau Au.

Cach tach: dung quan the Duc (DKMS, n=3,456,066) lam don bay, khoa lan luot
tung bien.

[XAC LAP] CA HAI hieu ung deu that va CONG DON.
  Khoa quan the -> van con khoang cach giua locus  => co hieu ung LOCUS
  Khoa locus    -> van con khoang cach giua quan the => co hieu ung QUAN THE
  Truong hop xau nhat = giao cua hai: allele locus C, pho bien o Kinh,
  hiem o Duc -> trung vi BA = 0.

[XAC LAP] Thieu hut khu tru o BA, KHONG phai o toan bo du lieu.
  HLA-C chiem 17.9% du lieu EL nhung chi 2.4% du lieu BA (lech 7.5 lan).

[XAC LAP] 15/46 (33%) allele pho bien o Kinh ma hiem o Duc vang HOAN TOAN
  khoi bang huan luyen MHCflurry. Phia Duc: 0/65.
  Moi allele trong tap lien quan ma MHCflurry chua tung thay du lieu nao
  deu la allele cua nguoi Kinh.
  LUU Y: MHCflurry la mo hinh pan-allele -- "vang khoi bang huan luyen"
  KHONG dong nghia "khong du doan duoc". No co nghia du doan cho nhung
  allele do la ngoai suy thuan tu pseudosequence. Chang 8.3 da cho thay
  chung gan voi ho hang duoc dac trung tot, nen rui ro thuc te duoc giam bot.

=============================================================================
## B. SO LIEU
=============================================================================

### Phan 1 -- khoa quan the (chi 65 allele pho bien o Duc, AF >= 0.005)
locus  n   trung vi BA  p25   p75   BA=0    tong BA   tong EL   EL/BA
A     17      4,300   2,546  5,601  0/17     79,825   158,134    2.0
B     31        290      24  2,921  3/31     45,235   245,829    5.4
C     17        177      25    244  3/17      3,120    87,891   28.2

Ty trong (cung allele, cung quan the, chi khac cong nghe do):
        BA      EL
  A   62.3%   32.2%
  B   35.3%   50.0%
  C    2.4%   17.9%

### Phan 2 -- khoa locus (trung vi BA: Duc-pho-bien vs Kinh-rieng)
nguong      HLA-A              HLA-B              HLA-C
k>=1    4300 vs 7  (614x)  290 vs 0 (vo cuc)  177 vs 0 (vo cuc)
k>=2    4300 vs 20 (215x)  290 vs 0 (vo cuc)  177 vs 0 (vo cuc)
k>=3    4300 vs 14 (307x)  290 vs 4 (82.9x)   177 vs 0 (vo cuc)

### Phan 3/4 -- Spearman AF_duc <-> BA
nguong   n    chung    HLA-A    HLA-B    HLA-C
k>=1    85   +0.712   +0.752   +0.762   +0.815
k>=2    63   +0.593   +0.666   +0.689   +0.770
k>=3    50   +0.573   +0.605   +0.578   +0.775
=> 9/9 o deu duong. Hieu ung quan the KHONG phu thuoc lua chon nguong.

### Phan 5 -- chan doan BA=0
pho bien o Duc : 6 BA=0 | vang khoi bang 0  | trong bang nhung BA=0 6
pho bien o Kinh: 28 BA=0 | vang khoi bang 15 | trong bang nhung BA=0 13
=> Khong allele Duc nao bi thieu khoi bang. Nen trung vi o Phan 1 KHONG phai
   tao tac cua join thieu -- chung la thieu hut that. Phan 1 dung vung.
   Vi du: C*16:01 co 5,788 EL va 0 BA. B*49:01: 5,959 EL, 0 BA.

### Do phu cua tap Duc da cat (mau so = 1.0 = toan bo nhiem sac the)
  HLA-A: giu 96.9%  | AFND bao tong 0.99543, hut 0.46% | cat nguong bo 2.6%
  HLA-B: giu 95.7%  | AFND bao tong 0.99564, hut 0.44% | cat nguong bo 3.9%
  HLA-C: giu 97.9%  | AFND bao tong 0.99518, hut 0.48% | cat nguong bo 1.6%

### CHIEU vung, DO LON khong vung
[XAC LAP]      Chieu: Duc-pho-bien >> Kinh-rieng, moi locus, moi nguong.
[KHONG XAC LAP] Do lon. Ty le locus A nhay 614x -> 215x -> 307x KHONG don
  dieu theo nguong. Do la nhieu trung vi o mau nho (n = 9-14), khong phai
  tin hieu. README KHONG duoc trich mot con so ty le duy nhat -- phai bao
  trung vi kem n cua tung nhom.

=============================================================================
## C. BON KHAI NIEM (ghi lai de 3 thang sau con doc duoc)
=============================================================================

### C1. Tong AF ~ 3.0 nghia la gi
Khong phai quy uoc, ma la phep DEM. Moi nguoi co 2 ban sao moi locus.
101 nguoi => 202 nhiem sac the tai HLA-A. Moi mot trong 202 cai do mang
DUNG MOT allele HLA-A. Nen tong tan suat tat ca allele HLA-A = 202/202 = 1.0.
Tan suat allele la mot PHAN HOACH tap nhiem sac the. Ba locus => tong 3.0.

Mot con so bat hai loai loi nguoc nhau:
  < 3.0  thieu dong (go sot, hoac typing khong goi duoc)
  > 3.0  trung dong (dem hai lan)

Ap dung:
  kinh     2.99992  lech 0.00008 -- nhieu lam tron thuan (5 c.s. x 90 dong)
  duc_full 2.98625  lech 0.01375 = 0.458% -- CHUA RO NGUYEN NHAN (xem muc D)
  duc      2.90490  PHAI nho hon 3.0, do thiet ke (da cat duoi 0.005)

Voi tap DA LOC, tong AF khong con la phep kiem dung/sai -- no tro thanh
thuoc do DO PHU. Mau so phai la 1.0, khong phai tong AFND bao (~0.995);
chia cho 0.995 se am tham hap thu phan hut cua AFND vao phan "giu lai".

### C2. Ty trong BA / EL theo locus nghia la gi
Mot minh ty trong BA gan nhu vo nghia -- khong co moc so sanh hien nhien
(chia deu 33%? Khong, vi locus B co 31 allele con A va C moi ben 17).

Suc manh nam o chO dat HAI ty trong canh nhau. Hai dong tinh tren CUNG 65
allele, CUNG 3 locus, CUNG mot quan the. Moi thu y het -- chi khac CONG NGHE
DO. Nen moi chenh lech giua hai dong KHONG the giai thich bang:
  - "locus B co nhieu allele hon"     -> giong nhau o ca hai dong
  - "quan the khac nhau"              -> cung mot quan the
  - "HLA-C it quan trong sinh hoc"    -> neu vay EL cung phai thap, nhung
                                          EL = 17.9%
Chi con mot loi giai thich: ban than assay BA da bo qua HLA-C.

Vi sao bao CA trung vi LAN ty trong -- hai cau hoi khac nhau:
  trung vi -> "mot allele dien hinh co bao nhieu du lieu?" Ben voi ngoai lai.
  ty trong -> "tong cong suc nghien cuu da chay ve dau?" BI VAI ALLELE
              KHONG LO CHI PHOI (A*02:01 mot minh chiem phan lon locus A).
Dung mot cai ma bo cai kia la tu bit mot mat.

### C3. Spearman -- cach tinh va y nghia
Cach tinh: (1) thay moi gia tri bang HANG trong cot cua no (dong hang nhan
hang trung binh); (2) tinh Pearson thong thuong tren HANG. Script tu tinh
bang numpy vi env hla thieu scipy.

Vi sao dung hang: af_duc trai 5 bac do lon (1e-6 den 0.28), BA trai 0 den
~12,000, ca hai lech rat manh. Pearson tren gia tri goc se bi A*02:01 va
A*11:01 keo di gan het. Spearman chi hoi: "allele xep hang cao hon ve tan
suat o Duc co xu huong xep hang cao hon ve luong du lieu khong?" -- bat bien
voi moi bien doi don dieu.

Y nghia trong bai: trong so cac allele pho bien o Kinh, tan suat o Duc giai
thich phan lon THU TU ve luong du lieu BA. Tinh trong tung locus rieng, nen
locus da bi giu co dinh -- do la cai bien no tu phat bieu ve locus thanh
phat bieu ve quan the.

Nhung NO KHONG NOI:
  - Khong phai nhan qua. Tan suat o Duc khong GAY RA du lieu. Ca hai cung
    phan anh mot nguyen nhan thu ba: quan the nao da duoc nghien cuu.
  - Khong loai tru ngoai le. A*02:03 co AF Duc 0.000054 ma 6,382 BA.
    Spearman noi ve THU TU TRUNG BINH.
    => Day dung la loi da mac: thay 2 ngoai le -> ket luan "quan the khong
       anh huong". Ngoai le bac bo tinh TAT YEU, khong bac bo tinh TRUNG
       BINH. Muon noi "X khong anh huong" phai do tuong quan toan tap.
  - CHUA co khoang tin cay / p-value. n chi 13-41 moi locus. Lam duoc bang
    permutation test (xao tron ngau nhien mot cot hang nghin lan) khong can
    scipy. CHUA LAM -- ghi vao han che.
  - Dong hang (DA DO; du doan ban dau SAI):
      k>=1: af_duc 4 (5% mau)  | BA 28 (33% mau)
      k>=2: af_duc 4 (6% mau)  | BA 18 (29% mau)
      k>=3: af_duc 4 (8% mau)  | BA 14 (28% mau)
    Du doan sai: "22 allele singleton co cung AF nen tao khoi dong hang 26%".
    Sai vi Spearman tinh tren af_duc <-> BA; AF cua KINH KHONG tham gia phep
    tinh -- no chi dung de CHON tap allele. Cac singleton Kinh co af_duc khac
    nhau (A*02:03 = 0.000054, A*33:03 = 0.0037), nen khong dong hang.
    Khoi dong hang THAT la BA = 0 (28 allele): nhieu allele cung khong co
    phep do BA nao -> cung hang.
    Huong anh huong khong doi: dong hang lam giam phuong sai cot hang va KEO
    |Spearman| VE 0, nen +0.712 co le la DANH GIA THAP, khong phai phong dai.
    Chua dinh luong duoc muc anh huong.
    BAI HOC: "bien nao dang thuc su nam trong phep tinh" -- dung gia dinh,
    phai do.

### C4. Nguong 0.005 lay o dau
Tra loi thang: TU CHON, roi dung ly le bien minh SAU.
  Ban 1 viet trong note: "Nguong 0.005 = 1/202 = gioi han phat hien".
  Nhung 1/202 = 0.004950495, KHONG phai 0.005. Phuong trinh do SAI.
  No song sot nhieu luot vi GAN DUNG -- lech 1%, du nho de khong ai thay.

Cong bang: 0.5% LA mot quy uoc pho bien trong di truyen quan the HLA de
dinh nghia "allele thuong gap". Neu da viet "0.005 la quy uoc 0.5%" thi dung.
Nhung da viet "0.005 = gioi han phat hien" -- va do la phat bieu sai.

Phan biet quan trong: nguong QUY UOC (chon, phai neu ro la chon) khac nguong
DAN XUAT (tinh ra tu thiet ke nghien cuu). Da trinh bay cai thu nhat nhu
cai thu hai.

=============================================================================
## D. LUONG TU HOA -- bai hoc chuyen giao duoc
=============================================================================

AF = k/202 voi k la mot phep DEM (nguyen). Mau so BAT DONG.
=> AF chi co 203 gia tri kha di. AF = 0.007 la BAT KHA THI trong bo nay.
=> Truc AF khong phai duong lien, ma la CAI THANG, bac cach nhau 0.00495.
   Giua hai bac la chan khong.

   AF
 0.01980 |====  k=4   11 allele
 0.01485 |====  k=3   11 allele   } giu
 0.00990 |====  k=2   13 allele   }
         | - - - - - - - - - - -  nguong 0.005
 0.00495 |====  k=1   22 allele   <-- BI LOAI SACH
         | - - - - - - - - - - -  nguong 0.0049
 0       |      k=0

Cai thang that (tu du lieu Kinh):
  k=1: 22 allele | k=2: 13 | k=3: 11 | k=4: 11 | k=5: 1 | k=6: 4 | k>=7: 23
  cat o k>=1: 85 | k>=2: 63 | k>=3: 50      (tong 85)

### HE QUA 1 -- nguong tro thanh HAM BAC THANG
Ca mot DAI nguong cho ket qua y het, roi nhay bac:
  nguong 0.0001 - 0.00495 -> k>=1 -> 85 allele
  nguong 0.0050 - 0.00990 -> k>=2 -> 63 allele
  nguong 0.0099 - 0.01485 -> k>=3 -> 50 allele

### HE QUA 2 -- GAN BANG VE SO khac GAN BANG VE TAC DUNG
  y dinh: cat o 1/202 = 0.004950 | da viet: 0.005
  sai lech ve so:      1%        | sai lech ve mau: 26% (85 -> 63)
Tren dai luong lien tuc, sai 1% thi he qua cung co 1%. Tren dai luong luong
tu hoa, sai 1% co the nhay CA MOT BAC.
=> Truc giac "gan dung thi gan du" CHI hop le tren dai luong lien tuc.

### HE QUA 3 -- cung mot con so, HAI PHEP TOAN khac nhau
            Kinh              Duc
  n         101               3,456,066
  nhiem sac the  202          6,912,132
  bac cach nhau  0.00495      0.00000014
O Duc, bac sat nhau den muc THUC TE LA LIEN TUC (dich nguong tu 0.0049 len
0.005 doi 0 allele -- da kiem). Nen "AF >= 0.005" la:
  voi Duc : phep cat lien tuc, vo hai
  voi Kinh: phep cat roi rac, loai ca mot bac
Cung ky hieu, hai ban chat. DO MOI LA LOI THAT -- khong phai chuyen 0.0049
hay 0.005, ma la dung mot nguong chung ma khong nhan ra no hanh xu khac nhau.

### BAY DA SAP HAI LAN
Da kiem bang Decimal: gia tri trong bai bao bi LAM TRON XUONG so voi phan so
that, voi dung nhung k can quan tam:
  k=1: 1/202 = 0.004950495, bai bao ghi 0.00495  -> 0.00495 < 1/202
  k=2: 2/202 = 0.009900990, bai bao ghi 0.00990  -> 0.00990 < 2/202
  k=3: 3/202 = 0.014851485, bai bao ghi 0.01485  -> 0.01485 < 3/202
(Luu y: khong phai moi k deu lam tron xuong -- k=11,12,16 lam tron LEN.
 Nhung dung ba k dung lam nguong thi deu xuong.)

=> NGUONG DAT DUNG BANG k/202 SE LOAI CHINH BAC k DO. Moi lan.

Lan 2: ban 2 cua script dat nguong do nhay = 2/202 va dan nhan
"(loai singleton)", nghi rang cho k>=2. THUC TE cho k>=3 (50 allele, khong
phai 63). Nhan SAI.
  Tin tot: do nhay thuc te da phu BA muc cat that (k>=1, 2, 3) -- rong hon
           du dinh. Ket luan vung qua ca ba => manh hon tuong.
  Tin xau: bay sap ngay trong code viet ra de tranh no.
  Y nghia: bay khong sap vi thieu hieu biet. No sap vi VAN VIET NGUONG DUOI
           DANG AF. Hieu roi ma dung cong cu cu thi van sap.

### CACH SUA (ban 3 cua script)
  k = round(af_kinh * N_CHROM);  giu = k >= 1
Ba loi the:
  1. Khong co epsilon -- so nguyen so voi so nguyen.
  2. Mien nhiem voi lam tron cua bai bao: round(0.00495 * 202) = 1 bat ke
     bai bao ghi 0.00495 hay 0.004950495.
  3. Tu giai thich -- "allele quan sat duoc it nhat 1 lan trong 202 nhiem
     sac the" la cau nguoi doc kiem duoc; "AF >= 0.0049" thi khong.
Phia Duc giu nguyen AF >= 0.005 nhung GOI DUNG TEN: quy uoc 0.5%.

### NGUYEN TAC
Khi mau so cua mot ty le la CO DINH va NHO, ty le do bi luong tu hoa. Moi
phep loc / so sanh / xep hang tren no phai phat bieu bang TU SO (phep dem),
khong bang ty le. Ranh gioi "nho": bac cach nhau co dang ke so voi nguong
ta muon dat hay khong.

=============================================================================
## E. CON MO / CHUA LAM
=============================================================================
1. Tong AF cac dong khac 0 cua AFND chi ~0.995, hut ~0.45% o ca ba locus.
   Chua ro nguyen nhan; lam tron khong du giai thich. KHONG bia ly do.
   Can hoi AFND neu dung so nay de chuan hoa.
2. Chua co khoang tin cay cho Spearman (permutation test -- lam duoc, chua lam).
3. Chua dinh luong anh huong cua khoi dong hang len Spearman.
4. Co che cua hieu ung locus: gia thuyet la thoi ky assay IC50 can MHC tinh
   sach, va HLA-C bieu hien be mat thap hon A/B. CHUA KIEM phan "kho tinh
   sach" trong y van. Neu dua vao README phai gan nhan gia thuyet.
5. Can chay lai script ban 3 de lay: do phu theo mau so 1.0, kich thuoc khoi
   dong hang, va nhan nguong dung. CAC SO LIEU O MUC B VAN DUNG -- chi nhan
   nguong la sai (0.0049/0.005/0.009901 thuc chat la k>=1/2/3).


## Chang 8.5 -- Tra y van hai cau chua kiem (2026-10-10)

### 1. Nguong 50 trong 14_pseudoseq_distance.py -- CAN CU KHONG DUNG DUOC
Docstring cu viet: nguong 50 "theo tien le MHCflurry 2.0 (O'Donnell et al.
2020), khi ho thay du doan cua allele hiem bang allele ho hang gan nhat co
>= 50 mau". Da tra:
  - MHCflurry 1.x (Cell Syst 2018): nguong la 200 phep do moi allele, KHONG
    phai 50. Bai do noi ro "no information is shared between alleles" va
    khuyen dung NetMHCpan cho allele it du lieu -- nguoc lai voi y "du doan
    bang allele ho hang".
  - MHCflurry 2.0 la pan-allele; nguong du lieu toi thieu moi allele khong
    ap dung theo kieu do.
  - Khong truy cap duoc toan van 2.0 (ScienceDirect robots, PMC captcha,
    EuropePMC robots). => "khong xac minh duoc", khong phai "da bac bo".
[XAC LAP] Con so 50 khong co can cu duoc neu. Da sua docstring va da them
  mot dong han che vao README.
LOAI LOI: lan thu BA cung kieu -- so tron + ly le gan vao sau
  (0.005 -> 2/202 -> 50). Lan nay nang hon hai lan truoc vi no QUY KET cho
  mot bai bao cu the. Quy tac moi: khong viet ten tac gia canh mot con so
  tru khi da mo bai bao do ra va doc dung cau.

### 2. Lich su assay va HLA-C -- CO NGUON, da sua cach dien dat
Cau cu: "in-vitro IC50 ... concentrated on HLA-A and -B, while mono-allelic
mass-spectrometry studies cover HLA-C more evenly" -- [SUY LUAN], chua nguon.
Da tra:
  - Abelin et al. 2017 (mono-allelic MS dau tien, 16 allele): toan HLA-A va
    HLA-B, KHONG co HLA-C nao duoc nhac ten. => cau cu SAI neu hieu la
    "mono-allelic MS noi chung phu HLA-C".
  - Sarkizova et al. 2020 (Nat Biotechnol 38:199-209), 95 dong te bao
    mono-allelic: 31 HLA-A, 40 HLA-B, 21 HLA-C, 3 HLA-G. Chon allele theo
    cho THIEU du lieu gan ket trong co so du lieu cong khai. Bao cao muc
    cai thien lon nhat o HLA-C (+20% PPV) va HLA-G (+38%), va goi do la
    "poorly characterized alleles".
[XAC LAP] Phu HLA-C bang MS den tu cac nghien cuu mono-allelic QUY MO LON
  (Sarkizova 2020), khong phai tu mono-allelic MS noi chung.
[KHONG XAC LAP] Nguyen nhan vi sao HLA-C it du lieu BA tu dau. Sarkizova
  khong noi, va README gio ghi ro la chua xac lap.
