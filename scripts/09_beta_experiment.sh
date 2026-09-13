#!/usr/bin/env bash
#
# Thi nghiem tham so --beta (homozygosity detection) cua OptiType.
#
# Cau hoi: co gia tri beta nao cho ket qua dung tren MOI mau khong?
#
# Hai chieu kiem tra:
#   Chieu 1 - mau goi la DONG HOP TU o allele HIEM (nghi ngo false homozygote)
#             ha beta -> co tach thanh di hop tu khong?
#   Chieu 2 - mau goi la DONG HOP TU o allele PHO BIEN (tin la dung)
#             ha beta -> co bi tach NHAM thanh di hop tu khong?
#
# Neu chieu 1 sua duoc ma chieu 2 khong hong -> beta thap tot hon.
# Neu ca hai deu doi -> beta la danh doi, khong phai loi chinh duoc.

set -uo pipefail

N=20000
SEED=100
OUT=results/beta_experiment
LOG=logs/beta_experiment

# Chieu 1: nghi la false homozygote (allele hiem)
CHIEU1="SRR11212885 SRR11212891 SRR11212973"
BETA1="0.001 0.009 0.05"

# Chieu 2: dong hop tu o allele pho bien (kiem tra nguoc)
CHIEU2="SRR11212897 SRR11212967 SRR11212938 SRR11212904 SRR11212971"
BETA2="0.001"

mkdir -p "$OUT" "$LOG"

chay() {
  local s=$1 b=$2
  local tag="${s}_beta${b}"
  [ -f "${OUT}/${tag}_result.tsv" ] && { echo "[bo qua] $tag"; return; }

  seqtk sample -s${SEED} data/clean/${s}_1.fq.gz $N | gzip > /tmp/${tag}_1.fq.gz
  seqtk sample -s${SEED} data/clean/${s}_2.fq.gz $N | gzip > /tmp/${tag}_2.fq.gz

  echo "[chay] $tag"
  if optitype run -i /tmp/${tag}_1.fq.gz -i /tmp/${tag}_2.fq.gz \
       --dna -o "$OUT" -p "$tag" --beta "$b" --enumerate 3 --threads 4 \
       > "${LOG}/${tag}.log" 2>&1; then
    echo "   OK"
  else
    echo "   LOI"
  fi
  rm -f /tmp/${tag}_*.fq.gz
}

echo "=== CHIEU 1: nghi false homozygote ==="
for s in $CHIEU1; do for b in $BETA1; do chay "$s" "$b"; done; done

echo ""
echo "=== CHIEU 2: kiem tra nguoc ==="
for s in $CHIEU2; do for b in $BETA2; do chay "$s" "$b"; done; done

echo ""
echo "=== TONG HOP ==="
python scripts/09b_beta_summary.py
