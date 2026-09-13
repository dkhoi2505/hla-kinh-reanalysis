#!/usr/bin/env bash
set -uo pipefail

mkdir -p data/raw data/clean logs/fastp

# tao danh sach URL cho tat ca 101 run (bo qua mau pilot da co)
tail -n +2 data/ena_full.tsv | cut -f2 | tr ';' '\n' | sed 's|^|https://|' > data/all_urls.txt
echo "So file can tai: $(wc -l < data/all_urls.txt)"

cd data/raw
wget -c -nv -i ../all_urls.txt
cd ../..

echo "=== TAI XONG, kiem tra toan ven ==="
loi=0
for f in data/raw/*.fastq.gz; do
  gzip -t "$f" 2>/dev/null || { echo "HONG: $f"; loi=$((loi+1)); }
done
echo "So file hong: $loi"
[ $loi -gt 0 ] && { echo "DUNG LAI — tai lai file hong truoc"; exit 1; }

echo "=== LOC BANG FASTP ==="
for s in $(cut -f1 data/runs_sorted.tsv); do
  [ -f "data/clean/${s}_2.fq.gz" ] && continue
  echo "[fastp] $s"
  fastp -i data/raw/${s}_1.fastq.gz -I data/raw/${s}_2.fastq.gz \
        -o data/clean/${s}_1.fq.gz  -O data/clean/${s}_2.fq.gz \
        --length_required 100 --detect_adapter_for_pe --thread 4 \
        --json logs/fastp/${s}.json --html logs/fastp/${s}.html \
        2> logs/fastp/${s}.log
done

echo "XONG: $(ls data/clean/*_1.fq.gz | wc -l) mau da loc"
