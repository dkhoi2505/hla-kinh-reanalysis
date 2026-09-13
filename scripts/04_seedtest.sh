#!/usr/bin/env bash
set -uo pipefail
mkdir -p results/seedtest logs/seedtest

for s in SRR11212911 SRR11212875 SRR11212960; do
  for seed in 200 300; do
    tag="${s}_20000_s${seed}"
    [ -f "results/seedtest/${tag}_result.tsv" ] && { echo "[bo qua] $tag"; continue; }

    if [ ! -f "data/sub/${tag}_2.fq.gz" ]; then
      seqtk sample -s${seed} data/clean/${s}_1.fq.gz 20000 | gzip > data/sub/${tag}_1.fq.gz
      seqtk sample -s${seed} data/clean/${s}_2.fq.gz 20000 | gzip > data/sub/${tag}_2.fq.gz
    fi

    echo "[chay] $tag"
    optitype run -i data/sub/${tag}_1.fq.gz -i data/sub/${tag}_2.fq.gz \
      --dna -o results/seedtest -p "$tag" --enumerate 3 --threads 4 \
      > logs/seedtest/${tag}.log 2>&1 && echo "   OK" || echo "   LOI"
  done
done
