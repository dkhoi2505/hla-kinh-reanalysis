#!/usr/bin/env bash
set -uo pipefail

mkdir -p data/sub results/main logs/main
N=20000
SEEDS="100 200 300"

tong=0; xong=0; loi=0
for s in $(cut -f1 data/runs_sorted.tsv); do
  for seed in $SEEDS; do
    tong=$((tong+1))
    tag="${s}_${N}_s${seed}"

    if [ -f "results/main/${tag}_result.tsv" ]; then
      xong=$((xong+1)); continue
    fi

    if [ ! -f "data/sub/${tag}_2.fq.gz" ]; then
      seqtk sample -s${seed} data/clean/${s}_1.fq.gz $N | gzip > data/sub/${tag}_1.fq.gz
      seqtk sample -s${seed} data/clean/${s}_2.fq.gz $N | gzip > data/sub/${tag}_2.fq.gz
    fi

    if optitype run -i data/sub/${tag}_1.fq.gz -i data/sub/${tag}_2.fq.gz \
         --dna -o results/main -p "$tag" --enumerate 3 --threads 4 \
         > logs/main/${tag}.log 2>&1; then
      xong=$((xong+1))
      echo "[$xong/$tong] OK  $tag"
    else
      loi=$((loi+1))
      echo "[$xong/$tong] LOI $tag"
    fi

    # xoa file sub sau khi chay xong de tiet kiem dia
    rm -f data/sub/${tag}_1.fq.gz data/sub/${tag}_2.fq.gz
  done
done

echo "=== XONG: $xong thanh cong, $loi loi / $tong ==="
