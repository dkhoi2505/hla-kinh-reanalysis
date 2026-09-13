#!/usr/bin/env bash
set -uo pipefail

SAMPLES="SRR11212911 SRR11212875 SRR11212960"
LEVELS="2000 5000 10000 20000 40000"
SEED=100

mkdir -p data/sub results/ladder logs/ladder

for s in $SAMPLES; do
  for n in $LEVELS; do
    tag="${s}_${n}_s${SEED}"

    # da co ket qua thi bo qua -> chay lai script khong mat cong
    if [ -f "results/ladder/${tag}_result.tsv" ]; then
      echo "[bo qua] $tag"
      continue
    fi

    # subsample neu chua co
    if [ ! -f "data/sub/${tag}_2.fq.gz" ]; then
      echo "[subsample] $tag"
      seqtk sample -s${SEED} data/clean/${s}_1.fq.gz $n | gzip > data/sub/${tag}_1.fq.gz
      seqtk sample -s${SEED} data/clean/${s}_2.fq.gz $n | gzip > data/sub/${tag}_2.fq.gz
    fi

    echo "[chay] $tag"
    if /usr/bin/time -v optitype run \
        -i data/sub/${tag}_1.fq.gz \
        -i data/sub/${tag}_2.fq.gz \
        --dna -o results/ladder -p "$tag" \
        --enumerate 3 --threads 4 \
        > logs/ladder/${tag}.log 2>&1; then
      echo "   -> OK"
    else
      echo "   -> LOI (xem logs/ladder/${tag}.log)"
    fi
  done
done

echo "XONG"
