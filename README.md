# OptiType HLA class I typing on targeted amplicon data: high population-level concordance, unstable rare-allele calls

A reanalysis of PRJNA609593 (Kinh Vietnamese, n = 101)

OptiType reproduces the published HLA class I allele frequencies of this cohort from long-range PCR amplicon data, with Spearman 0.983 against the commercial pipeline used in the source study. Individual calls are less reliable. Across three independent read subsamples per individual, rare alleles changed between runs about 4.7 times as often as common ones, and every allele called here but absent from the published table came from an unstable call.

These patterns are consistent with limitations the OptiType authors described for closely related alleles (Szolek *et al.* 2014). The contribution of this repository is to measure them on a data type the tool was not benchmarked on, and to show that the obvious parameter fix trades one error for another.

**Scope:** HLA class I (A, B, C) only. The source study also typed DRB1 and DQB1; OptiType does not support class II.

---

## Key findings

**1. Population-level frequencies are reproduced.**

At 2-field resolution, allele frequencies from 101 individuals correlate with the published Table 1 at Spearman 0.983 (HLA-A 0.990, HLA-B 0.969, HLA-C 0.984).

**2. Instability concentrates in rare members of common allele groups.**

Each sample was typed three times from independent subsets of 20,000 read pairs. Rare alleles (≤3 occurrences) changed between runs in 39.1% of their appearances; common alleles in 8.4%. The most unstable alleles are rare relatives of common ones — `A*11:06`, `A*11:10`, `A*24:08`, `A*24:20`, `A*02:11` — and all six alleles called here but absent from Table 1 are among them.

**3. Tuning `--beta` does not fix it.**

Lowering OptiType's homozygosity threshold turned three suspected false homozygotes heterozygous. In five other samples, homozygous for common alleles, the same setting split five loci — each time into a rare allele. The trade-off is not specific to OptiType: other HLA callers rely on similarly empirical thresholds.

---

## Background: OptiType's reference and amplicon data

Three features of how OptiType builds its reference explain these results. All are described in the original methods (Szolek *et al.* 2014).

- **Only exons 2 and 3, with flanking intron, are used.** This window encodes the peptide-binding groove and is the only region sequenced for nearly all known alleles; restricting the reference to it gives every allele an equal chance of being called. Entries in the shipped reference are about 1.5 kb long.
- **Missing intron sequence is reconstructed from the nearest fully sequenced relative.** The authors report 10,779 reconstructed sequences. The reference shipped with OptiType 1.5.0 contains exactly 10,779 reconstructed entries — 10,713 of 11,047 (97.0%) at HLA-A, -B and -C. The authors validated the reconstruction at 99.89% intron similarity, and found that including intron sequence reduced typing error 2.7- to 3.9-fold on exome data.
- **Alleles never reported in population databases are excluded** before optimisation (allelefrequencies.net, dbMHC).

Two consequences matter for amplicon data that covers whole genes:

- **Most reads have no target.** Reads from exon 1, exons 4–8, distal intron, and all class II loci fall outside the reference window. In one sample, about 10% of R1 reads mapped.
- **Close relatives can be indistinguishable over most of their sequence.** When an allele's intron is reconstructed from a close relative — `C*06:03` from `C*06:02`, for example — the two are identical across that intron in the reference and differ only in exon sequence. The call then depends on the few reads that happen to cover a differing position. In one sample, 9,088 R1 reads produced 7,059,824 alignments: a mean of 777 reference entries per read, ranging from 1 to 2,823.

The authors named both resulting failure modes: ambiguity between alleles that differ only in poorly covered segments, and homozygous calls when two highly similar alleles form a heterozygous locus. This analysis observes both on amplicon data.

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

Reads were downloaded from ENA, filtered with `fastp` (`--length_required 100`, adapter auto-detection), subsampled to **20,000 read pairs** with `seqtk`, and typed with **OptiType 1.5.0** (`--dna`, razers3 3.5.12, GLPK 5.0, default `--beta 0.009`).

Each sample was run with **three independent seeds** (100, 200, 300) to separate reproducible calls from calls that depend on which reads happened to be drawn. 303 runs completed, 0 failures.

Subsampling was necessary rather than convenient. On the smallest sample, using all reads took 15m25s and 4.44 GB peak RAM, versus 1m03s and 1.43 GB at 20,000 pairs — for an identical genotype. Cost scales non-linearly: 4.6× the reads produced 14.7× the runtime.

Comparison against the published Table 1 was done at **2-field resolution**. OptiType does not resolve 3-field differences, and allele names have changed across IMGT/HLA releases. The source study's 3-field calls are therefore **not** independently checked here.

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

**Six alleles called only here** — `A*02:11`, `A*02:12`, `A*02:45`, `B*39:02`, `B*40:49`, `B*48:08`. Each occurs in exactly one individual, and all six are among the alleles that changed between seeds. Two independent measurements point at the same calls.

**Four alleles only in Table 1.** `C*04:82` has no entry in OptiType's reference and cannot be called. `A*33:01` (20 entries), `C*03:17` (11) and `B*55:18` (2) are present in the reference but were not called.

### Call stability

| Locus | Stable across 3 seeds |
|---|---|
| HLA-A | 81/101 (80.2%) |
| HLA-B | 93/101 (92.1%) |
| HLA-C | 93/101 (92.1%) |

74 samples were stable at all three loci, 19 at two, 7 at one, and 1 at none.

HLA-A is the least stable locus. It carries the two largest intron donors in the reference — `A*02:01` supplies reconstructed intron to 263 other 2-field alleles, `A*24:02` to 198 — and the alleles that fluctuate most are rare members of those groups: `A*11:06`, `A*11:10`, `A*11:19`, `A*24:08`, `A*24:20`, `A*02:11`, `A*02:12`, `A*02:45`.

The alleles most often involved in an unstable call are `A*24:02` (8 occasions), `A*11:01` (7) and `A*02:01` (4) — the group anchors themselves.

### Homozygote excess

| Locus | Observed | Expected (HWE) | Ratio |
|---|---|---|---|
| HLA-A | 14 | 10.7 | 1.31 |
| HLA-B | 6 | 5.7 | 1.05 |
| HLA-C | 9 | 11.3 | 0.79 |

The absolute excess at HLA-A is 3.3 individuals. This is consistent with the zygosity failure mode described above but **not statistically significant at n = 101**. The source study reported no HWE deviation, and these data do not contradict that.

### The `--beta` trade-off

`--beta` sets how many additional reads a heterozygous solution must explain before OptiType prefers it over a homozygous one. The authors chose 0.009 by cross-validation on exome data from the 1000 Genomes Project.

Three samples called `A*11:02` homozygous were re-run across beta values:

| beta | Call |
|---|---|
| 0.001 | `A*11:01` / `A*11:02` |
| 0.009 (default) | `A*11:02` / `A*11:02` |
| 0.05 | `A*11:02` / `A*11:02` |

All three changed in the same direction — the pattern expected when two highly similar alleles form a heterozygous locus. The population data point the same way: the default calls give 3 fewer `A*11:01` and 3 more `A*11:02` than Table 1, and re-calling these three individuals as heterozygous would bring both alleles to exactly their Table 1 counts. This suggests, but does not prove, that the default calls were wrong.

Five samples homozygous for *common* alleles were then re-run at beta 0.001. Five loci were split:

| Locus | Default (0.009) | beta 0.001 |
|---|---|---|
| C | `C*07:02` / `C*07:02` | `C*07:02` / `C*07:123` |
| A | `A*02:07` / `A*02:07` | `A*02:07` / `A*31:01` |
| A | `A*24:02` / `A*24:02` | `A*24:02` / `A*24:07` |
| A | `A*02:07` / `A*02:07` | `A*02:07` / `A*30:02` |
| C | `C*01:02` / `C*01:02` | `C*01:02` / `C*02:19` |

Every newly introduced allele is rare, which points to spurious reads rather than a genuine second allele.

Why a threshold is needed at all: maximising the number of explained reads can never be hurt by adding a second allele, so without a penalty the optimum is always heterozygous. The OptiType authors note that the unpenalised formulation favours heterozygous solutions because of spurious hits such as sequencing errors. Other sources of such reads include PCR misincorporation and cross-mapping from related loci.

High beta risks false homozygotes; low beta risks false heterozygotes. This is the sensitivity–specificity trade-off of any cut-off: when a true second allele is supported by only a few discriminating reads — as for close relatives — its signal overlaps with spurious reads, and no threshold separates them cleanly. On this data no single value is right for every sample.

The trade-off is not specific to OptiType. arcasHLA calls a locus homozygous when the minor allele's non-shared reads fall below 15% of the major allele's, and culls alleles below 10% of the top abundance, a cut-off it adopted from HISAT-genotype (Orenbuch *et al.* 2020). Like OptiType's β, these values were chosen empirically.

---

## Limitations

**The mechanisms are not new.** OptiType's reference construction, intron reconstruction and both failure modes observed here are described in Szolek *et al.* 2014. This analysis quantifies them on long-range PCR amplicon data; it does not identify a new mechanism.

**No per-sample ground truth.** The source publication reports allele and haplotype frequency tables, not individual genotypes. Comparison is therefore population-level only. Where the two methods disagree, this analysis **cannot determine which is correct**. Disagreements are reported as disagreements, not as errors by either side.

**Stability is not accuracy.** The three-seed design measures precision, not correctness. Three samples labelled high-confidence by that criterion are, on the evidence of the beta experiment and the population counts, probably wrong. A tool with a systematic bias would be perfectly stable and still wrong.

**Reference age.** The number of reconstructed entries in the shipped reference matches the figure published in 2014 exactly, which suggests it is the original build from IMGT/HLA Release 3.14.0 (July 2013). This is an inference from the count, not a confirmed fact. If correct, alleles named after 2013 cannot be called — `C*04:82` has no entry — and the comparison with a 2019-era commercial pipeline spans six years of nomenclature.

**Reference window.** Alleles that differ only outside exons 2–3 cannot be separated by OptiType, by design. Reference data for the remaining exons would help resolve some rare alleles.

**Single tool.** Only OptiType was tested. HLA-HD, Kourami and HISAT-genotype construct their references differently and may behave differently.

**Class II not addressed.** DRB1 and DQB1 carry the most distinctive alleles in this population (`DRB1*12:02`, `DQB1*03:01`) but fall outside OptiType's scope.

**Library preparation inference is unverified.** Measured insert size (211–266 bp) is roughly ten-fold shorter than the ~2 kb fragmentation described in the source methods. A tagmentation step is the most plausible explanation but has **not** been checked against the TruSight HLA technical documentation.

**Laboratory artefacts are not separated from computational ones.** Allelic dropout during PCR — preferential amplification of one allele — also produces false homozygotes, and polymerase errors produce spurious reads. These data cannot distinguish such effects from threshold effects in the typing algorithm.

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

A second objective is in progress: how much training data do peptide–MHC binding predictors actually have for the alleles common in this population? `results/mhcflurry_training_counts.csv` holds raw per-allele counts; the population-weighted analysis is not finished.

**Related work:** [luad-deg-analysis](https://github.com/dkhoi2505/luad-deg-analysis) — differential expression in lung adenocarcinoma, including an MHC class I antigen presentation panel.

---

## References

Do MD, Le LGH, Nguyen VT, Dang TN, Nguyen NH, Vu HA, Mai TP. High-Resolution HLA Typing of HLA-A, -B, -C, -DRB1, and -DQB1 in Kinh Vietnamese by Using Next-Generation Sequencing. *Front Genet.* 2020;11:383. doi:[10.3389/fgene.2020.00383](https://doi.org/10.3389/fgene.2020.00383)

Szolek A, Schubert B, Mohr C, Sturm M, Feldhahn M, Kohlbacher O. OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics.* 2014;30(23):3310–3316. doi:[10.1093/bioinformatics/btu548](https://doi.org/10.1093/bioinformatics/btu548)

Orenbuch R, Filip I, Comito D, Shaman J, Pe'er I, Rabadan R. arcasHLA: high-resolution HLA typing from RNAseq. *Bioinformatics.* 2020;36(1):33–40. doi:[10.1093/bioinformatics/btz474](https://doi.org/10.1093/bioinformatics/btz474)

This repository is an independent reanalysis of public data. It is not affiliated with, endorsed by, or produced in collaboration with the authors of any cited work.
