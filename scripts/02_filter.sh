#!/usr/bin/env bash
set -uo pipefail
mkdir -p data/clean logs/fastp

for s in $(cut -f1 data/runs_sorted.tsv); do
  [ -f "data/clean/${s}_2.fq.gz" ] && continue
  echo "[fastp] $s"
  fastp -i data/raw/${s}_1.fastq.gz -I data/raw/${s}_2.fastq.gz \
        -o data/clean/${s}_1.fq.gz  -O data/clean/${s}_2.fq.gz \
        --length_required 100 --detect_adapter_for_pe --thread 4 \
        --json logs/fastp/${s}.json --html logs/fastp/${s}.html \
        2> logs/fastp/${s}.log
done

echo "=== FASTP XONG: $(ls data/clean/*_1.fq.gz | wc -l)/101 ==="
bash scripts/05_run_all.sh
