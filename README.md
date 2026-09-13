# Intron-borrowed reference sequences limit OptiType HLA class I typing on targeted amplicon data

A reanalysis of PRJNA609593 (Kinh Vietnamese, n = 101)

OptiType reproduces published population-level HLA class I allele frequencies from targeted amplicon DNA data (Spearman 0.983), but one class of calls — rare alleles sharing an allele group with a common one — is unstable in a way that parameter tuning cannot fix. This repository documents the mechanism, quantifies it, and shows why the usual remedy fails.

**Scope:** HLA class I (A, B, C) only. The source study also typed DRB1 and DQB1; OptiType does not support class II.

---

## Key findings

**1. The reference is almost entirely composite.**

Of 11,047 class I entries in OptiType's genomic reference, 10,713 (97.0%) pair real exon sequence with intron sequence borrowed from a related allele, because IMGT/HLA has full genomic sequence for only a minority of alleles. Only 334 entries (3.0%) carry the allele's own genomic sequence.

| Locus | Entries | Target alleles | Intron-borrowed |
|---|---|---|---|
| HLA-A | 3,823 | 2,192 | 97.4% |
| HLA-B | 3,374 | 2,852 | 95.8% |
| HLA-C | 3,850 | 1,753 | 97.7% |

**2. This matters specifically for amplicon DNA.**

Long-range PCR amplicon data covers the whole gene, so most reads fall in introns — that is, in sequence that does not reflect the target allele. In one sample, 9,088 distinct reads produced 7,059,824 alignments: a mean of 777 reference entries per read, ranging from 1 to 2,823.

**3. Instability is concentrated, and parameter tuning trades one error for another.**

Running each sample three times with independent read subsets, rare alleles (≤3 occurrences) changed between runs in 39.1% of appearances, versus 8.4% for common alleles. Lowering `--beta` fixed three suspected false homozygotes but created five false heterozygotes in five test samples.

---

## Data

| | |
|---|---|
| Source | Do MD *et al.*, Front Genet 2020;11:383 |
| Accession | [PRJNA609593](https://www.ebi.ac.uk/ena/browser/view/PRJNA609593) |
| Samples | 101 runs, one per individual |
| Library | `AMPLICON` — long-range PCR of HLA-A/-B/-C/-DRB1/-DQB1, Illumina MiniSeq |
| Size | 10.64 GB compressed |
| Reference typing | Assign TruSight HLA v2.0 (commercial), 3-field |

Read properties below were measured from the data, not taken from the publication:

- Read length mode 151 bp with a trimmed tail; 80–86% of reads ≥150 bp
- Insert size peak 211–266 bp
- No Q40 bases, indicating quality-score binning upstream
- Read count per sample varies 29.9-fold (96,030 to 2,872,506)

---

## Methods

Reads were downloaded from ENA, filtered with `fastp` (`--length_required 100`, adapter auto-detection), subsampled to **20,000 read pairs** with `seqtk`, and typed with **OptiType 1.5.0** (`--dna`, razers3 3.5.12, GLPK 5.0).

Each sample was run with **three independent seeds** (100, 200, 300) to separate reproducible calls from calls that depend on which reads happened to be drawn. 303 runs completed, 0 failures.

Subsampling was necessary rather than convenient. On the smallest sample, using all reads took 15m25s and 4.44 GB peak RAM, versus 1m03s and 1.43 GB at 20,000 pairs — for an identical genotype. Cost scales non-linearly: 4.6× the reads produced 14.7× the runtime.

Comparison against the published Table 1 was done at **2-field resolution**. OptiType cannot resolve 3-field differences (synonymous and intronic variants), and the 2020 allele names predate several IMGT/HLA renamings. The study's 3-field claims are therefore **not** independently checked here.

---

## Results

### Population allele frequency

| Locus | n alleles | Spearman | Pearson |
|---|---|---|---|
| HLA-A | 24 | 0.990 | 0.995 |
| HLA-B | 40 | 0.969 | 0.996 |
| HLA-C | 17 | 0.984 | 0.999 |
| Pooled | 81 | **0.983** | 0.997 |

87 alleles were called here, 85 appear in Table 1, and 81 are shared.

**Six alleles called only here** — `A*02:11`, `A*02:12`, `A*02:45`, `B*39:02`, `B*40:49`, `B*48:08`. Each occurs in exactly one individual, and all six appear in the set of alleles that fluctuated across seeds. This is an independent check that they are artifacts.

**Four alleles only in Table 1.** `C*04:82` has no entry in OptiType's reference and therefore cannot be called. `A*33:01` (20 entries), `C*03:17` (11 entries) and `B*55:18` (2 entries) are present in the reference but were missed.

### Call stability

| Locus | Stable across 3 seeds |
|---|---|
| HLA-A | 81/101 (80.2%) |
| HLA-B | 93/101 (92.1%) |
| HLA-C | 93/101 (92.1%) |

74 samples were stable at all three loci, 19 at two, 7 at one, and 1 at none.

HLA-A is the least stable locus. It also carries the two largest intron donors in the reference — `A*02:01` lends intron sequence to 263 other 2-field alleles and `A*24:02` to 198 — and the alleles that fluctuate most are rare members of those same groups: `A*11:06`, `A*11:10`, `A*11:19`, `A*24:08`, `A*24:20`, `A*02:11`, `A*02:12`, `A*02:45`.

The alleles most often involved in an unstable call are `A*24:02` (8 occasions), `A*11:01` (7) and `A*02:01` (4) — the group anchors themselves.

### Homozygote excess

| Locus | Observed | Expected (HWE) | Ratio |
|---|---|---|---|
| HLA-A | 14 | 10.7 | 1.31 |
| HLA-B | 6 | 5.7 | 1.05 |
| HLA-C | 9 | 11.3 | 0.79 |

The absolute excess at HLA-A is 3.3 individuals. This is consistent with the mechanism described above but **not statistically significant at n = 101**. The source study reported no HWE deviation, and these data do not contradict that.

### The `--beta` trade-off

`--beta` is OptiType's homozygosity detection threshold. Three samples called `A*11:02` homozygous were re-run across beta values:

| beta | Call |
|---|---|
| 0.001 | `A*11:01` / `A*11:02` |
| 0.009 (default) | `A*11:02` / `A*11:02` |
| 0.05 | `A*11:02` / `A*11:02` |

All three samples changed in the same direction, consistent with false homozygosity at the default setting.

However, five samples homozygous for *common* alleles were then re-run at beta 0.001, and five loci were split that should not have been:

| Sample locus | Default (0.009) | beta 0.001 |
|---|---|---|
| C | `C*07:02` / `C*07:02` | `C*07:02` / `C*07:123` |
| A | `A*02:07` / `A*02:07` | `A*02:07` / `A*31:01` |
| A | `A*24:02` / `A*24:02` | `A*24:02` / `A*24:07` |
| A | `A*02:07` / `A*02:07` | `A*02:07` / `A*30:02` |
| C | `C*01:02` / `C*01:02` | `C*01:02` / `C*02:19` |

Every newly introduced allele is rare.

High beta produces false homozygotes; low beta produces false heterozygotes. No single value is correct for all samples, because the read ratio that beta thresholds is itself uninformative when two alleles share intron sequence in the reference. The problem is in the reference, not in the parameter.

---

## Limitations

**No per-sample ground truth.** The source publication reports haplotype and allele frequency tables, not individual genotypes. Comparison is therefore population-level only. Where the two methods disagree, this analysis **cannot determine which is correct** — Assign TruSight draws on the same IMGT/HLA database and faces the same intron problem. Disagreements are reported as disagreements, not as errors by either side.

**Stability is not accuracy.** The three-seed design measures precision (reproducibility), not correctness. Three samples labelled high-confidence by that criterion were shown by the beta experiment to be wrong. A tool with a systematic bias would be perfectly stable and still wrong.

**Single tool, single reference version.** Only OptiType was tested. HLA-HD, Kourami and HISAT-genotype construct their references differently and may behave differently. The finding about intron-borrowed sequences applies to the reference shipped with OptiType 1.5.0 via bioconda.

**Class II not addressed.** DRB1 and DQB1 carry the most distinctive alleles in this population (`DRB1*12:02`, `DQB1*03:01`) but fall outside OptiType's scope.

**Library preparation inference is unverified.** Measured insert size (211–266 bp) is roughly ten-fold shorter than the ~2 kb fragmentation described in the source methods. A tagmentation step between the two is the most plausible explanation, but this has **not** been checked against the TruSight HLA technical documentation and is stated here as a hypothesis, not a finding.

**Statistical power.** Several observations rest on small counts: 3 excess homozygotes, 6 discordant alleles, 5 samples in the reverse beta test.

---

## Reproducing

```
conda env create -f environment-hla.yml
conda activate hla

bash   scripts/01_fetch.sh           # download from ENA (~10.6 GB)
bash   scripts/02_filter.sh          # fastp, then 101 x 3 OptiType runs (~5 h)
python scripts/06_gather.py          # 303 results -> consensus table
python scripts/07_reference_audit.py
python scripts/08_stability.py
bash   scripts/09_beta_experiment.sh
python scripts/10_compare_af.py
```

Scripts are idempotent — completed work is skipped on re-run.

Hardware used: WSL2, 8 cores, 3.7 GiB RAM (configured with 5 GB memory and 8 GB swap).

---

## Repository layout

| Path | Contents |
|---|---|
| `data/` | Accession lists, published Table 1 frequencies |
| `scripts/` | Numbered in execution order |
| `results/` | Summary tables (CSV) |
| `public/` | Genotype tables, sample IDs anonymised |
| `docs/NOTES.md` | Working log, in Vietnamese |

Per-sample genotypes are published with randomised sample codes rather than SRA accessions. The underlying reads are public, but linking a clinically meaningful genotype to a specific accession adds nothing to the analysis. The accession mapping is not distributed.

Raw reads, intermediate files and per-run OptiType outputs are not committed; `scripts/01_fetch.sh` retrieves them.

---

## Status

Objective 1 — typing and comparison — is complete.

A second objective is in progress: auditing per-allele training data coverage in peptide–MHC binding predictors for alleles common in this population. `results/mhcflurry_training_counts.csv` is a partial result and its interpretation is not yet settled.

**Related work:** [luad-deg-analysis](https://github.com/dkhoi2505/luad-deg-analysis) — differential expression in lung adenocarcinoma, including an MHC class I antigen presentation panel. Both repositories approach antigen presentation from opposite ends: one in tumour tissue, one at the population level.

---

## Citation

Do MD, Le LGH, Nguyen VT, Dang TN, Nguyen NH, Vu HA, Mai TP. High-Resolution HLA Typing of HLA-A, -B, -C, -DRB1, and -DQB1 in Kinh Vietnamese by Using Next-Generation Sequencing. *Front Genet.* 2020;11:383. doi:[10.3389/fgene.2020.00383](https://doi.org/10.3389/fgene.2020.00383)

This repository is an independent reanalysis of public data. It is not affiliated with, endorsed by, or produced in collaboration with the authors of that study.
