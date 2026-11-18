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

## Genome-wide comparison

`genome_wide_compare.py` compares two samples (e.g. an isolate and its matched metagenome) at every reference position, not only at the positions of a SNP table. Only reads that clearly come from the target organism are used: primary alignments with MAPQ >= 20 and minimap2 per-read divergence (`de` tag) < 0.10 (adjustable). A consensus base is called where at least 10 filtered reads are present and one base reaches at least 80 %; consensus bases are compared wherever both samples have one. Positions where the second sample carries an additional base (>= 10 %, >= 3 reads) absent from the first are listed for inspection in IGV. Indels are not assessed.

```
pip install pysam numpy
python genome_wide_compare.py --bam isolate=isolate_nodup.bam --bam metagenome=metagenome_nodup.bam \
  --contig CP033688.1 --outdir gw_isolate_vs_metagenome      # add --max-de 0.03 for a stricter run
```
**Input order matters.** Give `--bam` exactly twice, each as `NAME=PATH` (no spaces around `=`):
first the **isolate** (sample A), then the **metagenome** (sample B). The `--contig` value is the
sequence name in the BAM header (e.g. `CP033688.1`), not the FASTA file. The consensus comparison
gives the same result in either order, but `minor_alleles.tsv` and the last line of `summary.txt`
only report extra bases in sample B, so with the metagenome second they show bases present in the
metagenome but absent from the isolate.

Outputs: `summary.txt`, `differences.tsv`, `minor_alleles.tsv`, and cached per-chunk counts in `counts/` (an interrupted run resumes where it stopped). Options: `--min-mapq`, `--max-de`, `--min-aln-len`, `--min-depth`, `--min-frac`, `--minor-frac`, `--minor-reads`.

## High-identity reads at ambiguous positions

`count_high_identity_reads.py` counts, at selected positions, how many reads carry an expected base (e.g. the isolate base), separately for all reads and for reads with low per-read divergence (`de` < 0.10 by default), with and without restriction to primary alignments and MAPQ >= 20. It was used for the 22 positions with mixed base signals in the metagenome Dia5 (`mixed_positions_Dia5.tsv`, isolate 26-127 bases).

```
python count_high_identity_reads.py --bam metagenome_nodup.bam \
  --positions mixed_positions_Dia5.tsv --contig CP033688.1 --out high_identity_counts.tsv
```

## Per-read divergence histogram

`extract_read_divergence.py` writes the per-read divergence (minimap2 `de` tag) of all mapped, non-duplicate reads of a BAM file; `plot_read_divergence.R` (R, ggplot2) plots one histogram per sample (primary alignments, log-scaled counts, dashed lines at 3 % and 10 %) and writes a summary table with the share of reads below each threshold. The plot shows whether a divergence threshold separates target-organism reads from reads of other taxa.

```
python extract_read_divergence.py isolate_nodup.bam isolate isolate_de.tsv
python extract_read_divergence.py metagenome_nodup.bam metagenome metagenome_de.tsv
cat isolate_de.tsv > read_divergence.tsv; tail -n +2 metagenome_de.tsv >> read_divergence.tsv
Rscript plot_read_divergence.R read_divergence.tsv read_divergence_histogram \
  "isolate=Isolate;metagenome=Metagenome"
```

## Provenance

These scripts were written with assistance from Claude (Anthropic) and reviewed by the author.
