#!/usr/bin/env python3
"""Write the per-read divergence (minimap2 `de` tag) of all mapped reads in a BAM file.

Usage
  python extract_read_divergence.py sample.bam SAMPLE_NAME out.tsv

Duplicate and secondary alignments are skipped. Output columns:
  sample, read, primary (1 = primary, 0 = supplementary), mapq, aln_len (reference span, bp), de
Concatenate the outputs of several samples (keeping one header) as input for plot_read_divergence.R.
"""
import sys
import pysam

if len(sys.argv) != 4:
    sys.exit(__doc__)
bam_path, name, out = sys.argv[1:4]
with pysam.AlignmentFile(bam_path) as bam, open(out, "w") as fh:
    fh.write("sample\tread\tprimary\tmapq\taln_len\tde\n")
    for r in bam.fetch():
        if r.is_unmapped or r.is_duplicate or r.is_secondary or not r.has_tag("de"):
            continue
        fh.write(f"{name}\t{r.query_name}\t{0 if r.is_supplementary else 1}\t{r.mapping_quality}\t"
                 f"{r.reference_end - r.reference_start}\t{r.get_tag('de'):.4f}\n")
