#!/usr/bin/env python3
"""
Per-position read-support check for a SNP table (e.g. a vSNP-style SNP matrix).

For every variable position listed in the header row of the table, and for every
BAM file, the script counts the bases (A/C/G/T) and the deletions (reads with a
deletion at that position) and classifies the position:

    NOCOV  fewer than --min-depth reads (bases + deletions)
    DEL    >= --del-frac of the reads carry a deletion
    MIX    no single base accounts for >= --mixed-frac of the base calls
    OK     otherwise

Usage
-----
python check_snp_positions.py \
    --table SNP_table.xlsx \
    --bam isolate=isolate_nodup.bam --bam metagenome=metagenome_nodup.bam \
    --contig CP033688.1 --out per_position_counts.csv

Requirements: python >= 3.8, pysam, openpyxl. BAM files must be sorted and indexed.
The table must have positions such as "CP033688.1:1251875" in the first row (cells B1 onward).
"""
import argparse
import csv
import sys

import openpyxl
import pysam


def read_positions(table_path):
    ws = openpyxl.load_workbook(table_path, read_only=True).active
    header = next(ws.iter_rows(min_row=1, max_row=1, values_only=True))
    return [int(str(h).split(":")[1]) for h in header[1:] if h]


def count_position(bam, contig, pos):
    bases = {"A": 0, "C": 0, "G": 0, "T": 0}
    deletions = 0
    for column in bam.pileup(contig, pos - 1, pos, truncate=True, stepper="all",
                             min_base_quality=0, min_mapping_quality=0,
                             max_depth=1_000_000, ignore_overlaps=False,
                             ignore_orphans=False):
        for read in column.pileups:
            if read.is_del:
                deletions += 1
            elif not read.is_refskip and read.query_position is not None:
                base = read.alignment.query_sequence[read.query_position].upper()
                if base in bases:
                    bases[base] += 1
    return bases, deletions


def classify(bases, deletions, min_depth, del_frac, mixed_frac):
    n_bases = sum(bases.values())
    total = n_bases + deletions
    if total < min_depth:
        return "NOCOV"
    if deletions / total >= del_frac:
        return "DEL"
    if max(bases.values()) / n_bases < mixed_frac:
        return "MIX"
    return "OK"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="SNP table (.xlsx) with positions in row 1")
    ap.add_argument("--bam", action="append", required=True, metavar="NAME=PATH",
                    help="sample name and BAM path; repeat for each sample")
    ap.add_argument("--contig", default="CP033688.1")
    ap.add_argument("--out", default="per_position_counts.csv")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--del-frac", type=float, default=0.5)
    ap.add_argument("--mixed-frac", type=float, default=0.8)
    args = ap.parse_args()

    positions = read_positions(args.table)
    rows = []
    for item in args.bam:
        name, path = item.split("=", 1)
        with pysam.AlignmentFile(path, "rb") as bam:
            for pos in positions:
                bases, dele = count_position(bam, args.contig, pos)
                n = sum(bases.values())
                rows.append([name, pos, n, bases["A"], bases["C"], bases["G"], bases["T"], dele,
                             round(100 * dele / (n + dele), 1) if n + dele else 100.0,
                             classify(bases, dele, args.min_depth, args.del_frac, args.mixed_frac)])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "pos", "bases", "A", "C", "G", "T", "DEL", "pct_reads_with_deletion", "class"])
        w.writerows(rows)
    print(f"{len(rows)} rows written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
