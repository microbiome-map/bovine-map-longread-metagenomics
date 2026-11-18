#!/usr/bin/env python3
"""Count reads carrying a given base at selected positions, split by read filter.

Used to check positions that are ambiguous (mixed base signals) in a metagenome:
do the reads that clearly derive from the target organism (low per-read divergence,
minimap2 `de` tag) carry the same base as the matched isolate?

Input
  --bam        metagenome BAM (sorted, indexed, duplicates removed or flagged)
  --positions  TSV with columns `position` (1-based) and `expected_base`
               (e.g. the isolate base at each position)
  --contig     reference sequence name in the BAM header (e.g. CP033688.1)
  --max-de     divergence threshold for "high-identity" reads (default 0.10)

Output (TSV, one row per position plus a TOTAL row); for each read filter the number of
reads carrying the expected base (`*_expected`) and all reads with a base call (`*_n`):
  all      all non-duplicate, non-secondary alignments
  hi       alignments with de < --max-de (primary and supplementary)
  hi_prim  as `hi`, primary alignments only
  hi_mapq  as `hi_prim`, MAPQ >= --min-mapq

Example
  python count_high_identity_reads.py --bam metagenome_nodup.bam \
      --positions mixed_positions.tsv --contig CP033688.1 --out high_identity_counts.tsv
"""
import argparse, csv, sys
import pysam

FILTERS = ["all", "hi", "hi_prim", "hi_mapq"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--max-de", type=float, default=0.10)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--out", default="high_identity_counts.tsv")
    a = p.parse_args()

    with open(a.positions) as fh:
        pos = [(int(r["position"]), r["expected_base"].strip().upper()) for r in csv.DictReader(fh, delimiter="\t")]

    rows, total = [], [0] * (2 * len(FILTERS))
    with pysam.AlignmentFile(a.bam) as bam:
        for position, base in pos:
            c = [0] * (2 * len(FILTERS))
            for col in bam.pileup(a.contig, position - 1, position, truncate=True, stepper="nofilter",
                                  min_base_quality=0, min_mapping_quality=0, max_depth=10**7,
                                  ignore_overlaps=False, ignore_orphans=False):
                for pr in col.pileups:
                    if pr.is_del or pr.is_refskip or pr.query_position is None:
                        continue
                    r = pr.alignment
                    if r.is_unmapped or r.is_duplicate or r.is_secondary:
                        continue
                    b = r.query_sequence[pr.query_position].upper()
                    if b not in "ACGT":
                        continue
                    de = r.get_tag("de") if r.has_tag("de") else 1.0
                    hi = de < a.max_de
                    prim = not r.is_supplementary
                    keep = [True, hi, hi and prim, hi and prim and r.mapping_quality >= a.min_mapq]
                    for k, ok in enumerate(keep):
                        if ok:
                            c[2 * k + 1] += 1
                            c[2 * k] += (b == base)
            total = [x + y for x, y in zip(total, c)]
            rows.append([position, base] + c)

    header = ["position", "expected_base"] + [f"{f}_{s}" for f in FILTERS for s in ("expected", "n")]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
        w.writerow(["TOTAL", ""] + total)
    print(f"{len(rows)} positions; high-identity reads (de < {a.max_de}) with expected base: "
          f"{total[2]} of {total[3]}", file=sys.stderr)


if __name__ == "__main__":
    main()
