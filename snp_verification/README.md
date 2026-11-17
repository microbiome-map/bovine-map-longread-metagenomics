# SNP position check (isolate vs metagenome comparison)

`check_snp_positions.py` verifies the variable positions of a SNP table against the reads. For each position and sample it reports the base counts (A/C/G/T), the number of reads with a deletion and a class:

| class | meaning (default thresholds) |
|---|---|
| NOCOV | fewer than 10 reads |
| DEL | at least 50 % of reads carry a deletion (position lies inside a deletion, not a SNP) |
| MIX | no single base makes up at least 80 % of base calls (ambiguous) |
| OK | clear single-base call |

It was used to check a SNP matrix built from ONT reads (minimap2 alignment to CP033688.1, bcftools variant calling, vSNP matrix script) for positions that are deletions or mixed signals rather than SNPs. The counts reproduce those shown in the IGV coverage pop-up.

## Usage

```
pip install pysam openpyxl
python check_snp_positions.py --table SNP_table.xlsx \
  --bam isolate=isolate_nodup.bam --bam metagenome=metagenome_nodup.bam \
  --contig CP033688.1 --out per_position_counts.csv
```

BAM files must be sorted and indexed (`.bai`). Thresholds can be changed with `--min-depth`, `--del-frac` and `--mixed-frac`. Duplicate reads were removed before this step; pileup uses pysam's default read filters.

## Provenance

This script was written with assistance from Claude (Anthropic) and reviewed by the author.
