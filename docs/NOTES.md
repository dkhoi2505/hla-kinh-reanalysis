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
Rui ro goi sai tap trung o ca the mang ALLELE HIEM, khong phai allele
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

## Chang 8 -- du lieu huan luyen MHCflurry cho allele nguoi Kinh

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


## Chang 8.6 -- Doi chieu Sarkizova 2020 va hai sua chu (2026-10-12)

Khoi da tu mo bai Sarkizova et al. 2020 (Nat Biotechnol 38:199-209) va doi
chieu ba diem. CA BA DEU DUNG:
  1. "the collection of 95 cell lines (31 HLA-A, 40 HLA-B, 21 HLA-C and
     3 HLA-G)" -- Results, tieu muc "Systematic LC-MS/MS profiling of HLA
     class I ligands from mono-allelic cell lines".
  2. Chon allele theo cho thieu du lieu gan ket cong khai -- co trong bai.
  3. "+20% PPV for HLA-C and +38% for HLA-G" -- Results, tieu muc
     "Generation and performance of allele-and-length-specific and
     pan-allele-pan-length models".

### SUA 1 -- 'predicts' -> 'is associated with' (README dong 9)
Spearman do TUONG QUAN, khong do nhan qua. Tu 'predicts' trong tieng Anh de
bi doc thanh nhan qua. README khong co bang chung nhan qua nao, nen dung tu
yeu hon cho dung. Chi mot cho trong README dung tu nay theo nghia do; cac
cho khac la 'predictor'/'prediction' chi cong cu, giu nguyen.

### SUA 2 -- them 'among the classical loci' (README dong 266)
Diem 3 o tren chinh la cho cau cu YEU: bai bao bao HLA-G duoc +38%, CAO HON
HLA-C (+20%). Cau "reported its largest prediction gains at HLA-C" vi vay
bat be duoc. HLA-G la class I khong co dien, OptiType khong goi, nen no
ngoai pham vi bai nay -- gioi han vao 'classical loci' (A, B, C) vua dung
voi bai bao vua dung voi pham vi cua repo.
BAI HOC: khi trich mot so sanh ("lon nhat", "cao nhat") tu bai bao khac,
phai kiem xem bai do con hang muc nao cao hon khong, ke ca hang muc nam
ngoai pham vi cua minh.

### Them scripts/16_pseudoseq_sensitivity.py
Chay lai phan tich khoang cach pseudosequence o nguong BA 25 / 50 / 100 / 200
de tra loi: ba phan tram o Ket qua 4 (100 / 87.5 / 73.6%) co phu thuoc con so
50 tuy tien khong?
Thiet ke: khi nguong T tang thi tap "ngheo" rong ra VA tap "da dac trung" hep
lai cung luc, ca hai deu lam khoang cach tang => T=200 la phep thu khac
nghiet nhat.
Tu kiem: script tinh lai tu dau, khong dung lai logic script 14, nhung o T=50
PHAI tai tao dung 100 / 87.5 / 73.6. Lech thi no bao ">>> LECH" va phai dung
lai. (Da thu voi du lieu gia lap: canh bao nay no dung nhu thiet ke.)
KET QUA: chua chay tren du lieu that.


## Chang 8.7 -- Do nhay nguong NGUONG_BA: ket qua (2026-10-12)

Chay scripts/16_pseudoseq_sensitivity.py tren du lieu that.

DOI CHIEU DAT: o T=50 script 16 tai tao DUNG 100.0 / 87.5 / 73.6 cua script 14,
du no tinh lai tu dau va khong dung lai logic script 14. => hai script dong y,
cac nguong khac dang tin.

% (trong so AF) allele Kinh ngheo nam trong vong 2 thay the:
  nguong   HLA-A          HLA-B          HLA-C
    25     100.0% (n= 8)   86.3% (n=19)   96.3% (n= 8)
    50     100.0% (n= 9)   87.5% (n=21)   73.6% (n=11)   <- script 14
   100     100.0% (n= 9)   82.8% (n=23)   73.6% (n=11)
   200      96.3% (n=11)   65.8% (n=24)   53.6% (n=14)

[XAC LAP] Ket luan "phep ngoai suy la ngan" KHONG vung deu theo locus.
  HLA-A: vung o moi nguong (100 -> 96.3%).
  HLA-B: sut vua phai (87.5 -> 65.8%).
  HLA-C: sut manh (73.6 -> 53.6%). O T=200 chi hon mot nua chut.
  => Khong duoc viet "the extrapolation is short" khong kem dieu kien.
     Da sua README: neu ro ket qua theo tung locus va tung nguong.

[XAC LAP] PHAT HIEN MOI, va no lam cau chuyen locus MANH HON:
  Locus thieu du lieu BA nhat (HLA-C: 2.4% ty trong BA) cung la locus ma
  phuong an du phong pseudosequence YEU NHAT. Hai han che CONG DON cung
  chieu, khong bu tru nhau.
  Truoc do README ngam y phuong an du phong bu lai duoc cho thieu du lieu.
  Dieu do dung voi HLA-A, KHONG dung voi HLA-C.

[SUY LUAN] Vi sao HLA-C sut manh: chinh vi it allele HLA-C dat nguong, nen
  tap "da dac trung" o locus C mong, va allele ngheo phai voi xa hon de tim
  ho hang. Day la cung mot nguyen nhan voi thieu hut BA, nhin tu goc khac.

GHI CHU: o T=25 HLA-C len 96.3% roi tut xuong 73.6% o T=50 -- khong don dieu.
  Co che: C*01:02 co dung BA=25, nen o T=25 no CHUYEN TU nhom ngheo SANG nhom
  da dac trung va thanh diem tham chieu cho cac allele C khac. Tap ngheo co
  lai tu 11 xuong 8. Day la hanh vi dung cua thiet ke, khong phai loi.

BAI HOC: mot nguong tuy tien khong chi can duoc THU NHAN, no can duoc THU.
  Neu chi ghi "nguong nay tuy tien" roi dung lai, ta se bo lo mot phat hien
  that (hai han che cong don). Phep thu do nhay re hon nhieu so voi gia tri
  no mang lai.

--------------------------------------------------------------------
## Chang 9.1 -- MUC 1: loi lam tron (7.5 -> 7.3)

### Van de
README ghi ti trong EL cua HLA-C gap "7.5 lan" ti trong BA. Con so do
chi ra duoc khi chia HAI PHAN TRAM DA LAM TRON: 17.9 / 2.4 = 7.4583.

### Gia tri that
So tho lay tu chinh bang Result 1 cua README (Total BA, Total EL):
    HLA-A  79,825 BA   158,134 EL
    HLA-B  45,235 BA   245,829 EL
    HLA-C   3,120 BA    87,891 EL
    => tong BA = 128,180 ; tong EL = 491,854
    ti trong HLA-C: BA = 3120/128180 = 2.4341% ; EL = 87891/491854 = 17.8693%
    ti so that = 7.3413  ->  7.3

Ba duong tinh doc lap deu cho 7.3413:
  1. %EL / %BA (hai ti trong chinh xac)
  2. (EL_C/BA_C) x (tBA/tEL)  -- dai so tuong duong, khong qua ti trong
  3. phan so chinh xac tuyet doi 14443421/1967416
Duong 3 dung so nguyen nen loai han kha nang sai so dau phay dong.
Sai 1.6% o gia tri, khong phai sai o chu so cuoi.

### Kiem tra so 0
Tai tao ca sau phan tram cua bang Result 3 tu sau so tho: 6/6 khop
(62.3 / 35.3 / 2.4 va 32.2 / 50.0 / 17.9). Nen so tho va bang nhat quan
voi nhau; loi CHI nam o ti so.

### Da sua
  a) Key findings muc 6 : "a 7.5-fold" -> "a 7.3-fold"
  b) Result 3 cau ket   : "is 7.5 times" -> "is 7.3 times", VA them cau
     noi ro ti so tinh tu tong goc.
Ly do them cau o (b): bang Result 3 chi in 1 chu so thap phan, nen nguoi
doc tu chia 17.9/2.4 se duoc 7.5. Neu khong noi ro thi SAU khi sua README
se trong nhu co loi noi tai -- te hon truoc khi sua. Cau them bien mot cho
de bi bat loi thanh mot cho the hien phan biet duoc gia tri goc va gia tri
hien thi.

### Ra soat toan README: con ti so nao tinh tu gia tri da lam tron?
Kiem moi so sanh dang ti so. Bay ti so tinh tu SO THO, tat ca dung:
  EL:BA tai A/B/C  1.9810 / 5.4345 / 28.1702  -> 2.0 / 5.4 / 28.2
  24-fold   = 4300/177        = 24.2938
  29.9-fold = 2872506/96030   = 29.9126
  14.7x     = 925s/63s        = 14.6825
  4.6x reads                  = 4.6000

Tim them MOT truong hop cung ho loi: "about 4.7 times as often" (Key
findings) lay tu 39.1% / 8.4%, ca hai da lam tron. Chan tren chan duoi
cua ti so that la [4.5976, 4.6886) -- khoang nay CAT NGANG nguong 4.65,
nen TRUOC khi do thi khong the biet 4.6 hay 4.7.

Da do tu so dem goc (results/stability_by_allele.csv):
  hiem     68 allele,  110 xuat_hien,  43 bat_on  ->  39.0909%
  pho bien 40 allele,  512 xuat_hien,  43 bat_on  ->   8.3984%
  ti so chinh xac = 256/55 = 4.654545  -> 4.7
  ti so tu % lam tron      = 4.654762  -> 4.7
  chenh lech hai duong     = 0.000216

=> "4.7" DUNG, KHONG SUA.
Ly do cau truc khien no ben: gay_bat_on cua hai nhom BANG NHAU (43 va 43),
nen hai so 43 triet tieu va ti so rut con 512/110 = 256/55. Mot ti so rut
gon thanh phan so hai so nguyen nho thi rat ben voi lam tron.
Day la truong hop PHUONG PHAP LONG nhung KET QUA DUNG -- dung do may, khong
do can than. Van dung khi kiem, va dung khi khong khang dinh truoc khi do.

Kiem nhat quan phu, dat: dinh nghia "rare = 3 or fewer occurrences" doi hoi
hiem <= 68x3 = 204 (co 110, dat) va pho bien >= 40x4 = 160 (co 512, dat).

### Bai hoc ghi lai
1. Ti so phai tinh tu SO THO, khong bao gio tu gia tri da in ra. Ba lan
   mac loi nay trong du an (0.005 = 1/202, 2/202, va 7.5) deu cung mot goc:
   lay gia tri HIEN THI lam gia tri GOC.
2. Khi phat hien mot ti so tinh tu % da lam tron, DUNG ket luan no sai.
   Tinh chan tren chan duoi truoc. Neu khoang khong cat nguong lam tron thi
   con so van dung, chi phuong phap la long.
3. Sua mot con so trong van ban ma de nguyen bang so lam tron ben tren co
   the TAO RA mot ve mau thuan moi. Luon doc lai doan da sua bang mat nguoi
   doc, khong chi kiem chuoi.
4. Loi phuong phap: trong script va, dat dieu kien kiem
   README.count("87.5%") == 2 ma CHUA HE DO. Do ra la 1, script bao LECH
   ngay lan chay dau. Da doi thanh bat bien dung: so lan "87.5%" KHONG DOI
   truoc va sau khi va -- ve nay khong can biet con so la bao nhieu.
--------------------------------------------------------------------

--------------------------------------------------------------------
## Chang 9.2 -- MUC 2: phan tram pseudosequence co trong so tan suat

### Van de
Sau con so 100 / 87.5 / 73.6 / 96.3 / 65.8 / 53.6 deu CO TRONG SO TAN SUAT
(frequency-weighted) nhung README khong noi, va cau chu doc nhu ti le dem
allele.

### Tu kiem (scripts tam, khong commit)
Tinh lai CA HAI cach tu results/pseudoseq_sensitivity.csv, moi nguong:
                co trong so        dem thuan
  T= 25  A 100.0  B 86.3  C 96.3   A 100.0 (8/8)   B 78.9 (15/19)  C 87.5 (7/8)
  T= 50  A 100.0  B 87.5  C 73.6   A 100.0 (9/9)   B 81.0 (17/21)  C 81.8 (9/11)
  T=100  A 100.0  B 82.8  C 73.6   A 100.0 (9/9)   B 73.9 (17/23)  C 81.8 (9/11)
  T=200  A  96.3  B 65.8  C 53.6   A  90.9 (10/11) B 70.8 (17/24)  C 50.0 (7/14)

=> Cach CO TRONG SO tai tao dung 6/6 con so README. Vay README dang dung
   trong so ma khong khai bao.

### Phat hien them 1: dao chieu o HAI nguong, khong phai mot
Chi hai nguong README trich (T=50, T=200) tung duoc kiem. Chay ca bon:
  T= 25  cung chieu  (C tot hon B o ca hai cach)
  T= 50  DAO CHIEU   (trong so: C te hon B | dem thuan: C tot hon B)
  T=100  DAO CHIEU   (nhu tren)
  T=200  cung chieu  (C te hon B o ca hai cach)
T=100 cung dao. Van de khong chi o mot nguong.

### Phat hien them 2 (nang hon): cau ket luan cu qua manh
Cau cu: "The locus with the largest training-data deficit is also the locus
where the pan-allele fallback has to reach furthest."
Dieu kien kiem duoc: HLA-C phai la locus THAP NHAT trong ba.
  locus thap nhat, 8 o (4 nguong x 2 cach):
    T= 25  trong so -> B   | dem thuan -> B
    T= 50  trong so -> C   | dem thuan -> B
    T=100  trong so -> C   | dem thuan -> B
    T=200  trong so -> C   | dem thuan -> C
  => "C thap nhat trong ba locus"  dung  4/8 o
  => "C thap hon A"                dung  8/8 o

Cau cu SAI o T=25 NGAY CA KHI co trong so (B 86.3 thap hon C 96.3). README
chua bao gio in hai con so do, du T=25 nam trong dai da chay.

Da thu hep tuyen bo xuong cap A so C (ben 8/8). Day cung la dung cap ma
Result 1 dung de do thieu hut (trung vi BA 4.300 so 177, khoang cach 24 lan),
nen thu hep khong lam yeu lap luan "compound" ma lam no khop cap do hon.

### Co che cua viec dao chieu
Huong lech giua hai cach noi len ALLELE NAO nam xa:
  HLA-B T=50: trong so 87.5 CAO hon dem thuan 81.0
              => allele nam xa la allele IT pho bien, trong so lam B dep hon
  HLA-C T=50: trong so 73.6 THAP hon dem thuan 81.8
              => allele nam xa la allele PHO BIEN hon trung binh, trong so
                 lam C xau hon
Huong nay ben: chenh 6.5 va 8.2 diem phan tram, con lam tron chi gay sai so
0.05 diem. DO LON thi KHONG tinh duoc tu day -- muon biet hai allele HLA-C
nam xa co AF bao nhieu thi can AF tung allele, khong duoc suy tu 73.6% da
lam tron (dung la loi cua MUC 1).

Day la ly do thuc chat de chon trong so lam cach chinh: o HLA-C, allele ma
mo hinh phai ngoai suy xa nhat lai la allele TUONG DOI PHO BIEN trong quan
the nay. Dem thuan che mat dieu do; trong so phoi ra. Cau hoi co y nghia la
"bao nhieu phan ganh nang HLA cua quan the nam trong tam voi cua phuong an
du phong", va moi nhiem sac the mang dung mot allele moi locus, nen trong so
tan suat la phep tinh tra loi dung cau do.

### Da sua
  a) Result 4 cau dau  : doi khung tu "100% cua cac allele" sang "100% cua
     TAN SUAT allele thieu du lieu", de con so khong con doc nhu ti le dem.
  d) Result 4, DOAN MOI tach rieng: khai bao frequency-weighted, noi ro vi
     sao chon trong so, in ba con so dem thuan 100 (9/9) / 81.0 (17/21) /
     81.8 (9/11), neu co che (allele xa o B la allele hiem, o C la allele
     pho bien hon), va chot rang cap A so C ben 8/8 o.
     Ly do tach doan: gop het vao mot doan thi doan do dai 290 tu, qua day
     cho README se nam trong CV. Doan 1 giu phat hien + ket luan; doan 2 lo
     phuong phap.
  b) Result 4 cau ket  : thu hep tu "C thap nhat trong ba locus" (4/8) xuong
     "C thap hon A o ca bon nguong va ca hai cach dem" (8/8), va noi thang
     rang thu tu B so C khong on dinh nen khong tuyen bo nao dua vao no.
  c) Limitations nguong 50: them cau khai bao trong so va neu ro cap A so C
     la cap dung duoi ca hai cach.

### Bai hoc ghi lai
1. Mot con so co trong so ma khong khai bao thi nguoi doc se hieu thanh so
   dem. Moi ti le phai noi ro trong so theo cai gi.
2. Khi hai cach do cho thu tu KHAC NHAU, do chua chac la nhieu. Huong lech
   mang thong tin: no cho biet nhom nam ngoai la nhom nang hay nhe ve trong
   so. Doc huong lech truoc khi goi no la nhieu.
3. Kiem do nhay phai kiem MOI gia tri da chay, khong chi gia tri da in ra.
   Lan nay T=100 bi bo sot vi README khong in no.
4. Mot tuyen bo dang "X la cai nhat trong ba" yeu hon han dang "X thap hon
   Y". Truoc khi viet dang superlative, dem xem no dung bao nhieu o.
--------------------------------------------------------------------

--------------------------------------------------------------------
## Chang 9.3 -- MUC 3: khoi Reproducing

### Tu doc ma, khong doc mo ta
Doc toan van 5 script .sh (01, 02, 03, 04, 05) thay vi tin bang liet ke.
Ly do: MUC 3 la mo ta HANH VI, khong phai con so, nen khong co cach kiem nao
ngoai mo ma ra doc. Day cung la cho da sai mot lan trong du an (gan nguong 50
cho O'Donnell et al. 2020 ma chua doc bai).

### Khoi Reproducing cu KHONG dut o cho duoc neu ban dau
Mot nhan dinh can bac bo: "02_filter.sh CHI chay fastp", va "lam theo khoi
nay se chay fastp roi nhay toi 06_gather.py, khong co gi de gop".
KHONG DUNG. 02_filter.sh dong 16 la:  bash scripts/05_run_all.sh
Script do goi 05 o dong cuoi, nen chay 02 thi THAT SU chay fastp roi 303 lan
OptiType -- dung nhu README mo ta. Khoi Reproducing cu KHONG dut o cho do.

### Nhung moi lo do dung, vi ly do khac
README khong sai su that, no CHE MOT PHU THUOC. Ba he qua that:
  1. N=20000 va SEEDS="100 200 300" dat o 05_run_all.sh dong 5-6. README noi
     nhieu ve viec chon 20.000 cap read nhung nguoi doc muon doi thi khong
     biet mo file nao.
  2. Muon chay lai chi phan typing thi khong co duong. 05_run_all.sh co san
     resume (dong 14-16) ma README khong chi.
  3. 01_fetch.sh dong 22-31 CUNG chay fastp, dung cung vong nhu 02_filter.sh.
     Nen "01 tai, 02 loc" la sai: 01 = tai + kiem toan ven + loc,
     02 = loc lai (bo qua cai da xong) + typing.

### XEP LAI muc nghiem trong
Sau script thieu la dung, nhung cai nang nhat KHONG phai 05_run_all.sh:
  12_anonymize.py   -- sinh ra public/, thu ma Repository layout mo ta la san
                       pham cong bo. Lam theo khoi cu thi public/ KHONG BAO GIO
                       ton tai. => NANG NHAT, va la cho de bi bo sot nhat.
  09b_beta_summary  -- 09_beta_experiment.sh CO trong khoi, script tong hop cua
                       no KHONG. Chay thi nghiem roi bo do.
  05_run_all.sh     -- chay duoc qua 02, nhung tham so bi che.
  03_ladder + 03b   -- sinh cac con so chi phi README dang trich.
  04_seedtest.sh    -- grep README: khong dong nao trich results/seedtest/.

### Thi nghiem beta KHONG phai phu
grep README ra: dong 27, ca muc "### The --beta trade-off" (dong 137-167, hai
bang), dong 288 va 302. Nen 09 va 09b phai o DUONG CHINH. Du tinh ban dau
la xep chung vao nhom tuy chon -- sai, da sua truoc khi de xuat.

### Hai con so trong Methods
1. "4.6x the reads" -- DA DO THAT. Mau nho nhat sau fastp (SRR11212911) co
   91.801 cap read; 91801/20000 = 4.59 -> 4.6. DUNG.
   Nhung README in 96.030 (read THO) o dong 80, nen nguoi doc tu chia se ra
   4.80. Hai con so dung hai don vi ma README khong noi. Da them menh de
   "(91,801 cleaned pairs against 20,000)".
   LOI PHUONG PHAP o MUC 1: da bao "4.6x reads = 4.6000 ok" voi phep kiem
   (20000 x 4.6)/20000 -- MOT HANG DANG THUC, khong kiem gi ca, luon ra 4.6
   bat ke su that. Phep kiem vo nghia ma bao ok.
2. "using all reads took 15m25s and 4.44 GB peak RAM" -- KHONG ra tu script
   nao trong repo. 03_ladder.sh co LEVELS="2000 5000 10000 20000 40000", va
   ls results/ladder/ chi co dung 5 muc do x 3 mau = 15 lan chay. Toan bo
   read la 91.801 cap, gap hon doi muc cao nhat cua ladder.
   Ma day lai la con so README dung de bien minh cho ca quyet dinh subsample.
   Da them menh de noi thang: do trong mot lan chay rieng, khong script nao
   tai tao duoc.

### Tong ke: ca 18 script deu co cho
  duong chinh Phan I : 01, 05, 06, 07, 08, 09, 09b, 10, 12   (9)
  Phan II            : 11, 13, 14, 15, 16                    (5)
  do chi phi         : 03, 03b                               (2)
  khong o duong nao  : 02, 04  -- co giai thich tai sao       (2)
Khong con script nao trong repo ma README khong nhac ten.

### Bai hoc ghi lai
1. Mo ta hanh vi chi kiem duoc bang cach doc ma. Khong co phep tinh nao thay
   the duoc, khac han MUC 1 va MUC 2.
2. Cong cu quet tu khoa tu viet cho ket qua sai hai lan trong mot lan chay:
   - DUONG TINH GIA: 01_fetch.sh va 02_filter.sh bao "goi fastp" vi chuoi
     "fastp" nam trong TEN THU MUC logs/fastp, khong phai lenh goi.
   - AM TINH GIA: 05_run_all.sh bao khong goi OptiType, vi no goi
     "optitype run" CHU THUONG dang lenh con, con regex tim chu hoa.
   Phat hien duoc nho mot mau thuan so hoc (06_gather.py noi "303 ket qua" ma
   khong script nao co ve sinh ra chung), khong nho su can than.
   => Quet ma bang tu khoa thi dung -i, va KET QUA RONG nghia la "regex khong
      khop", KHONG phai "script khong goi".
3. Mot phep kiem co the "dat" ma khong kiem gi ca neu hai ve cua no cung xuat
   phat tu mot con so. Phep kiem phai lay so tu NGUON DOC LAP voi con so dang
   kiem. (4.6x la vi du: kiem bang 20000x4.6/20000 thi vo nghia; kiem bang
   cach dem read trong file .fq.gz thi moi co nghia.)
4. Mot con so trich trong README ma khong script nao sinh ra duoc la mot lo
   hong tai tao, du con so do dung. Noi thang con hon de nguoi doc di tim
   script khong ton tai.
--------------------------------------------------------------------

--------------------------------------------------------------------
## Chang 9.4 -- MUC 4: dinh vi so voi Atkins 2024 va Glynn 2025

### Da TU DOC toan van ca hai bai, khong qua tom tat cua ai
Quy tac Phan 3A cua bao cao Chang 8: khong trich khi chua doc. Lan nay lam
dung. Hai bai doc tu ban full text:
  Atkins TK, Solanki A, Vasmatzis G, Cornette J, Riedel M. Front Immunol.
    2024;14:1288105. doi:10.3389/fimmu.2023.1288105
  Glynn E, Ghersi D, Singh M. PNAS 2025;122(8):e2405106122.
    doi:10.1073/pnas.2405106122

### Atkins 2024 -- hai diem can tinh chinh
Xac nhan tung diem:
  nguoi Viet          : 5.1% (A), 27.7% (B), 12.1% (C), 33.6% (DRB1)
  chau Au da trang    : 0.4%, 0.9%, 0.6%, 2.6%
  8 allele moi        : 7 class I + DRB1*12:02, KHONG co HLA-B nao
  7 allele class I    : A*24:07, A*34:01, A*34:02, A*36:01,
                        C*03:02, C*04:03, C*14:03
  danh gia tren       : ELUTED LIGAND. Ho noi thang "did not perform any
                        experiments to characterize the model's ability to
                        predict binding affinity"
TINH CHINH 1: ket luan cua ho KHONG phai "khong co
  hieu ung". Class I moi hieu nang CAO HON trung binh, class II moi THAP HON,
  va CA HAI khac biet deu CO Y NGHIA THONG KE. Ho danh gia do lon qua nho de
  dang ke. "Hieu ung co nhung nho", khong phai "khong co".
TINH CHINH 2: du lieu eluted ligand ho dung lay tu
  Sarkizova et al. -- dung bai README nay dang dan o Result 3.

### Glynn 2025 -- day moi la bai nang
Trung lap DANG KE voi Phan II, va do tot hon (6.59 trieu nguoi):
  ti le pseudosequence CO du lieu binding: A 5.6%, B 5.1%, C 2.4%
  ~87% pseudosequence class I khong co du lieu binding nao
  benh nhan TCGA ghi la White co du lieu NHIEU HON CO Y NGHIA so voi Black,
    Asian, Native Hawaiian/Pacific Islander
  train set: 43.626 BA so 378.039 MA + 518.573 SA (EL)
  expected PPV thap nhat o HLA-C
Doi chieu voi repo:
  hieu ung locus     -> Result 1 va Result 3 cua repo do cung dieu nay
  hieu ung quan the  -> Result 2 (Spearman +0.57 den +0.82) do cung dieu nay
  87% ngoai suy      -> Result 4 va phan tich khoang cach pseudosequence
Glynn KHONG co nguoi Viet (gop vao "Asian and Pacific Islander").
Glynn DEM rieng BA va EL nhung KHONG bao cao hieu nang tach theo loai du lieu.

### HAI BAY SO TRUNG -- quan trong nhat cua muc nay
1. 2.4% cua Glynn = ti le PSEUDOSEQUENCE HLA-C co du lieu binding.
   2.4% cua repo  = ti le PHEP DO BA thuoc HLA-C.
   HAI DAI LUONG KHAC NHAU TINH CO TRUNG SO. Da viet mot cau trong README
   noi thang dieu nay, de nguoi doc khong tu ghep.
2. 33.6% cua Atkins = ti le CA NHAN Viet mang allele DRB1 thieu du lieu.
   33% cua repo     = ti le ALLELE class I (15/46) vang mat khoi bang
                      MHCflurry.
   Cung la hai dai luong khac nhau trung so. Khac predictor, khac mau so,
   khac locus. KHONG duoc dat canh nhau.

### Giao diem allele -- lien he truc tiep giua hai cong trinh
4 trong 7 allele class I Atkins danh gia CO trong cohort Kinh:
   A*24:07  AF 0.02970  k=6
   A*34:01  AF 0.01485  k=3
   C*03:02  AF 0.08911  k=18
   C*04:03  AF 0.07426  k=15
Hai allele C cong lai = 33/202 nhiem sac the (16.3%) cua locus C. Khong phai
allele ben le.
Khong co trong cohort: A*34:02, A*36:01, C*14:03.
LUU Y: C*14:02 CO trong cohort (AF 0.00495) nhung Atkins test C*14:03 --
allele khac, de lan.

### Ra soat tuyen bo ve chat luong du doan: README SACH, khong can sua
Quet 7 tu khoa (degrad, less accurate, worse, unreliab, quality of, accuracy,
practical concern) tren toan README: chi mot dong trung, la "Stability is not
accuracy" o Limitations -- cau do noi ve do chinh xac TYPING, khong phai chat
luong du doan peptide-MHC.
=> README khong chua tuyen bo "gap lam giam chat luong du doan". No da duoc
   viet dung pham vi tu dau: Result 4 noi thang "The statement is about
   direct evidence, not about whether output is produced."
=> Cau hoi "lam diu hay giu nguyen tuyen bo" phan lon khong con. Khong co tuyen bo nao can lam diu. Cai thieu la hai viec khac:
   thua nhan y van di truoc, va noi ro "dem du lieu khong phai do hieu nang".

### Da sua
  A) dau Phan II: doan dinh vi. Noi thang "What follows is not a new
     disparity", kem con so cua ca hai bai, kem cau phan biet hai dai luong
     2.4%.
  B) cuoi Result 3: Atkins do tren EL, dung assay ma phan tich nay thay
     KHONG thieu; 4 allele giao diem kem AF; va cau tu gioi han pham vi
     "this analysis counts how much data exists, not how well a model
     performs on it".
  C) References: them Atkins 2024 va Glynn 2025 sau Sarkizova.
  D) tom tat dau README: them "Both are documented elsewhere (Glynn et al.
     2025; Atkins et al. 2024); what is added here is a split by assay type."

Ly do them (D) du khong ai yeu cau: doan tom tat dau README la
thu Thay Minh doc truoc tien. Neu no tuyen bo hai hieu ung nhu phat hien moi
ma khong nhac y van, thi do dung la cho phoi nhat, va sua doan Phan II o duoi
khong che duoc.

### Danh gia tinh moi -- da ha bot, co y
Phan II KHONG moi o hai ket qua lon nhat. Hieu ung locus va hieu ung quan the
deu da co trong Glynn 2025, do tot hon.
Con lai la cua repo:
  1. Phan biet BA so EL tren CUNG tap allele, CUNG quan the (Result 3).
     Glynn dem rieng nhung khong bao cao hieu nang tach theo loai. Atkins
     danh gia CHI tren EL va tu noi khong kiem BA. Nen cau hoi "vi sao mot
     danh gia dua tren EL lai khong thay hieu ung, trong khi thieu hut nam o
     BA?" khong bai nao tra loi.
  2. Quan the Kinh cu the, typing do phan giai cao. Glynn khong co nguoi
     Viet; Atkins co nhung tu du lieu hien tuy My.
  3. Khoang cach pseudosequence kem phan tich do nhay 4 nguong. Glynn uoc
     luong hieu nang bang MHCPerf; repo do khoang cach Hamming. Bo sung
     nhau, khong trung.
Ha tuyen bo tinh moi LA viec nen lam: hoi dong hoac nguoi huong dan biet hai
bai nay se nghi toi chung ngay. Dinh vi dung gay an tuong hon han tuyen bo
moi roi bi hoi lai.

### Bai hoc ghi lai
1. Phep kiem "tong AF ~ 3.0" KHONG xac thuc duoc nguon. Mot ban kinh_af.tsv
   cu chi co 74 allele 2-truong (that: 85) van qua sach: moi AF la boi so
   chinh xac cua 1/202 va moi locus cong dung 202. Phep kiem do bat duoc hang
   THIEU hoac NHAN DOI trong cung mot lan trich xuat, KHONG phan biet duoc
   hai lan trich xuat khac nhau.
   => Them mot phep kiem nua: DEM SO ALLELE, khong chi cong AF.
2. Suyt mac loi cung ho MUC 1 o dang khac: dinh viet "C*03:02 o tan suat
   0.089, allele HLA-C pho bien thu tu" -- tinh tren dung cai file vua chung
   minh khong tin duoc. MUC 1 la tinh tren GIA TRI DA LAM TRON; lan nay la
   tinh tren NGUON CHUA XAC THUC. Cung mot goc: dung so ma khong truy nguon.
3. Khi hai nghien cuu bao cung mot con so, kiem DAI LUONG truoc khi vui. Hai
   lan trong muc nay (2.4% va 33.6%) deu la trung so chu khong phai trung
   ket qua. Mot README dat chung canh nhau se bi bat ngay.
4. Mot bai bao noi "khong tim thay bang chung" khong co nghia la "khong co
   hieu ung". Atkins tim thay hieu ung CO Y NGHIA THONG KE o ca hai huong,
   roi danh gia do lon qua nho. Trich "ho khong tim thay gi" la trich sai.
--------------------------------------------------------------------

--------------------------------------------------------------------
## Chang 9.5 -- MUC 5: kiem dinh hoan vi cho cac tuong quan Spearman

### Tu viet script, tinh lai tu dau
Viet scripts/17_permutation_test.py, tinh lai rho TU DAU (khong dung lai
logic script 15), 10.000 hoan vi moi o, hat giong 100.

### Tai tao: 12/12 o KHOP bang Spearman cua README
Script 17 va script 15 doc lap xac nhan nhau. Script 17 co co che tu dung:
neu lech bat ky o nao thi in ro va KHONG ghi file ket qua. Da test co che do
tren du lieu mock -- no phat hien lech ca 12 o va tu choi ghi file.

### Xac thuc phan hoan vi o BON muc (truoc khi chay tren du lieu that)
  1. xep_hang va spearman doi chieu dap an tinh tay:
     [1..5] vs [1..5]      -> +1.000
     [1..5] vs [5..1]      -> -1.000
     [1..5] vs [2,1,4,3,5] -> +0.800  (tinh tay: 8/10)
     [1,2,3,4] vs [1,1,2,2]-> +0.894  (tinh tay: 4/sqrt(20))
     bat bien qua bien doi don dieu: rho(x,y) = rho(x,log y) = rho(x,y^3)
  2. tuong quan hoan toan -> p = gia tri san 1/(N+1)
  3. khong co lien he      -> p lon, rho nam TRONG khoang rong
  4. do lech chuan phan phoi rong = 0.1850 so voi cong thuc giai tich
     1/sqrt(n-1) = 0.1857 cho n=30. KHOP toi ba chu so.
Muc 4 la xac thuc manh nhat: neu phep lay mau hoan vi sai thi gia tri do
khong khop cong thuc giai tich duoc.

### Loi phuong phap: hai gia tri ky vong trong phep kiem don vi
Viet ky vong xep_hang([3,1,2]) = [2,1,3] -- SAI, dung la [3,1,2].
Viet ky vong spearman([1,2,3,4],[1,1,2,2]) = 0.9487 -- SAI, dung la 0.8944
(0.9487 = 3/sqrt(10), mot phep tinh khac nho lan).
Ham DUNG o ca hai ca; gia tri ky vong sai. Da tinh tay bang phan so chinh xac
de xac dinh ben nao sai.
Day la lan thu NAM trong dot nay khang dinh truoc khi do.
CAI CUU: dung quy uoc nguoi dung dat ra -- IN GIA TRI THU DUOC CANH GIA TRI
KY VONG de doi chieu bang mat, KHONG dung assert. Neu dung assert thi chi
thay "that bai" va rat co the da di "sua" cai ham dang dung.

### KET QUA -- p thay doi hai bac do lon theo nguong
Neu chi chay k >= 1 (nguong de nhat) thi ca bon o deu cho "p < 0,0001".
Chay ca 12 o:
            k>=1     k>=2     k>=3
  gop     0.0001   0.0001   0.0001
  HLA-A   0.0001   0.0023   0.0190   <-- o mong nhat
  HLA-B   0.0001   0.0001   0.0060
  HLA-C   0.0001   0.0009   0.0034
O k>=3 HLA-A cho p = 0.019 -- van duoi 0.05 nhung LON HON HAI BAC DO LON so
voi con so o k >= 1. Neu README chi viet "p < 0.0001"
thi se la noi qua o 5 trong 9 o.

Khoang 95% cua phan phoi rong rong ra dung nhu phai:
  n=85 [-0.216,+0.219]   n=41 [-0.311,+0.311]   n=25 [-0.403,+0.396]
  n=19 [-0.461,+0.453]   n=13 [-0.545,+0.559]
O n=13, ngau nhien thuan du suc tao ra |rho| toi 0.56 -- chinh la ly do muc
han che do ton tai. Ca 12 o quan sat deu nam ngoai khoang rong, nhung o mong
nhat la HLA-A k>=3: quan sat +0.605 so voi bien tren +0.519, du 0.086.

Khoi dong hang: 23.1% den 36.8%. O k>=1 gop la 32.9% = dung 28/85, xac nhan
can cu cua muc han che "Tied ranks" dang co.

### "Group sizes are 13 to 41 per locus" -- README DUNG
  k>=1  A 25  B 41  C 19
  k>=2  A 19  B 28  C 16
  k>=3  A 15  B 22  C 13
min = 13 (HLA-C o k>=3), max = 41 (HLA-B o k>=1). Dem chi o k>=1 thi thay
min la 19; con so 13 den tu nguong chat nhat. KHONG CAN SUA GI.

### LOI KHAI NIEM trong chinh README, da sua
Cau cu: "A permutation test would give confidence intervals without
additional dependencies".
SAI VE PHUONG PHAP. Hoan vi kiem GIA THUYET RONG; no KHONG cho khoang tin
cay. Muon khoang tin cay cho rho thi phai bootstrap. Da sua trong e2.

### KHONG dong han muc han che, chi thu hep
Mot de nghi la "dong muc han che No uncertainty estimates". Da chi THU HEP.
Ly do: hoan vi tra loi "co phai ngau nhien khong", KHONG tra loi "rho that
nam trong khoang nao". Hai cau hoi khac nhau; cau thu hai can bootstrap.
Neu README viet "da co uoc luong bat dinh" thi trong phong van se bi hoi
"khoang tin cay cua rho la bao nhieu" va khong co cau tra loi. Viet dung cai
da lam -- da kiem gia thuyet rong, chua chan khoang -- thi tra loi duoc, va
con neu duoc buoc tiep theo.

### Da sua
  e1) Result 2: them doan ket qua hoan vi, kem khoang p 0.0001-0.019, o mong
      nhat, luu y cac phep kiem long nhau nen khong hieu chinh boi so duoc,
      va viec hoan vi xu ly dung dong hang.
  e2) Limitations: thay "No uncertainty estimates on the correlations" bang
      "The correlations are tested but not bounded" -- va sua luon loi khai
      niem ve confidence interval.
  e3) Limitations muc dong hang: them cau noi ro dong hang KHONG anh huong
      p-value hoan vi, chi anh huong GIA TRI cua rho.
  e4) Reproducing: them script 17.

### Bai hoc ghi lai
1. Uoc luong p kieu cong mot: p = (1 + so lan vuot)/(1 + N). Voi N=10.000 va
   0 lan vuot, p = 0.0001, KHONG phai 0. Mot kiem dinh N lan khong the chung
   minh p nho hon 1/(N+1). Viet "p < 0.0001" khi chi chay 10.000 lan la noi
   qua.
2. Kiem do nhay phai chay MOI nguong, khong chi nguong de nhat. Tab chien
   luoc chay k>=1 va ket luan "tat ca xa ngoai khoang rong"; o k>=3 thi p
   lon hon hai bac do lon. Lan thu hai trong dot nay bo sot cac gia tri khong
   duoc in ra (lan truoc la T=100 o MUC 2).
3. Hoan vi xu ly dong hang dung vi no xao tron chinh vecto hang. Cong thuc
   tiem can cua Spearman gia dinh khong dong hang, nen voi 23-37% dong hang
   no cho p sai. Chon phep kiem theo cau truc du lieu, khong theo tien le.
4. Mot de nghi "dong muc han che nay" can kiem xem phep do moi co THAT SU
   tra loi dung cau hoi cua muc han che do khong. O day khong: no tra loi
   mot cau hoi lan can. Thu hep dung hon la dong.
--------------------------------------------------------------------
