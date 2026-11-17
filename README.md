# Long-read shotgun metagenomics of bovine faecal microbiota and MAP

Analysis scripts supporting a study of bovine faecal microbiota and *Mycobacterium avium* subsp. *paratuberculosis* (MAP) infection using Oxford Nanopore long-read shotgun metagenomics.

## Scripts

| Script | Analysis |
|---|---|
| `01_Pilot_Analysis.R` | Paired 0/24 h pilot taxonomic composition, alpha and beta diversity |
| `02_Taxonomic_Analysis.R` | Taxonomic composition, differential-abundance analyses, five-classifier comparison |
| `03_Alpha_Diversity.R` | Alpha diversity, MAP-status, farm and classifier comparisons |
| `03b_Farm_Adjusted_Alpha_Diversity.R` | Farm-adjusted sensitivity analyses for alpha diversity (farm-level permutation, farm means, optional mixed model) |
| `04_Beta_Diversity.R` | Bray–Curtis, PCoA, PERMANOVA, dispersion, farms and classifier comparisons |
| `04b_Farm_Adjusted_Beta_IndVal.R` | Farm-adjusted sensitivity analyses for beta diversity and indicator taxa (farm-level exact permutation) |
| `05_IndVal.R` | Indicator taxa with `indicspecies::multipatt` |
| `06_KMA_MAP_Detection.R` | IS900/F57 marker detection and association with diagnostics/sequencing yield |

## Input data and folder layout

**The study input tables are deposited on Zenodo:** 10.5281/zenodo.24200566

Download the Zenodo `data.zip` and extract it at the **GitHub repository root** so that there is a `data/` directory alongside the `.R` files. The scripts use the original, descriptive filenames in the ZIP; **do not rename these files**. See the separate Zenodo `README_DATA.md` for formats, lineage naming and the metadata requirements.

```
repository/
├── 01_Pilot_Analysis.R
├── 02_Taxonomic_Analysis.R
├── 03_Alpha_Diversity.R
├── 03b_Farm_Adjusted_Alpha_Diversity.R
├── 04_Beta_Diversity.R
├── 04b_Farm_Adjusted_Beta_IndVal.R
├── 05_IndVal.R
├── 06_KMA_MAP_Detection.R
├── README.md
└── data/                                
    ├── metadata.txt
    ├── pilot/
    │   ├── Pilot_meta.txt
    │   └── Pilot_metaphlan_gtdb.txt
    ├── taxonomy/
    │   ├── MetaPhlAn_GTDB/
    │   ├── MetaPhlAn_default/
    │   ├── Kraken2/
    │   ├── Kraken2_conf10/
    │   └── ganon2/
    └── map_detection/
        ├── kma_IS900_F57_neu.xlsx          (Zenodo)
        └── S1_Sample_Master.xlsx           (Supplementary Table S1, see below)
```

Each taxonomy folder contains prefixed abundance and normalised profiles, for example `MetaPhlAn_GTDB_Genus_abundance.tsv`, `MetaPhlAn_GTDB_Genus_normalized.tsv`, `MetaPhlAn_Genus_normalized.tsv` under MetaPhlAn_default, or `Kraken2_conf10_Species_normalized.tsv`. The ganon2 folder additionally contains three `*_counts.tsv` files.

## Running the analyses

From the repository root in R, with the required packages installed:

```r
source("01_Pilot_Analysis.R")
source("02_Taxonomic_Analysis.R")
source("03_Alpha_Diversity.R")
source("03b_Farm_Adjusted_Alpha_Diversity.R")
source("04_Beta_Diversity.R")
source("04b_Farm_Adjusted_Beta_IndVal.R")
source("05_IndVal.R")
source("06_KMA_MAP_Detection.R")
```

Scripts generate their outputs in separate analysis-specific folders. For standalone execution, e.g. `Rscript 03_Alpha_Diversity.R`, launch from the repository root. Run alpha before beta if using its exported sample-inclusion QC file. The beta script can also determine the same inclusion set from metadata and the explicit exclusions.

Required packages across the scripts include `tidyverse`, `readxl`, `readr`, `openxlsx`, `writexl`, `vegan`, `ape`, `patchwork`, `indicspecies`, `permute`, and `broom`; `lme4` is optional for script 03b. Record `sessionInfo()` after a successful reproduction.

To run script 06, download Supplementary Table S1 from the article's supplementary files, save the unmodified workbook as `data/map_detection/S1_Sample_Master.xlsx`, and run the script from the repository root.

Script 06 calls a marker detected when KMA covers at least 90% of the marker reference sequence (`breadth_threshold`). Comparator status is taken at sample level from S1: a sample is reference-positive if its study qPCR or culture was positive, and reference-negative if it had a negative qPCR or culture result and no positive one. Reference-negative samples are reported in two sets: animals from herds without known MAP history (primary) and samples from previously test-positive animals (secondary). Samples without a study qPCR or culture result are not used for agreement. Besides the Wilcoxon, Spearman and logistic-regression results, the script writes `kma_agreement.csv` (PPA/NPA with exact 95% CIs), `detection_by_group.csv`, `detection_by_ct_bin_descriptive.csv` (equal-width Ct bins, descriptive only), `kma_reported_hits_qc.csv` (breadth, depth and identity of every reported hit) and `kma_threshold_sensitivity.csv` (any reported hit, ≥50% and ≥90% breadth).

## SNP position check

`snp_verification/check_snp_positions.py` (Python) checks the positions of a SNP table against the reads of the corresponding BAM files and flags positions inside deletions or with mixed base signals. See `snp_verification/README.md`.

For the published SNP distances, positions classified as deletions (`-`) or ambiguous (`N`) were excluded, and the SNP distance between two samples is the number of positions at which both carry an unambiguous base (A/C/G/T) and the bases differ; insertions and deletions are not counted. Positions within the structurally complex region CP033688.1:1,790,000–1,890,710 (split alignments) were masked and excluded from the distances. These two steps are not performed by `check_snp_positions.py` and were applied to its corrected output table.

## Authorship and provenance

These analysis scripts were written with assistance from ChatGPT (OpenAI) and Claude (Anthropic), reviewed by the first author, and consolidated for publication. Figure formatting may differ from manuscript illustrations.

## Citations

Will be inserted when available.
