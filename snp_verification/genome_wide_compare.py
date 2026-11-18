#!/usr/bin/env python3
"""Genome-wide consensus comparison of two samples aligned to the same reference.

Purpose: check whether an isolate and its matched metagenome differ anywhere in the
genome, not only at the positions of a SNP table. Only reads that clearly belong to
the target organism are used (primary alignments, MAPQ >= --min-mapq, minimap2
per-read divergence `de` < --max-de, optional minimum aligned length).

For every reference position and sample, a consensus base is called if the position
has >= --min-depth filtered reads and one base reaches >= --min-frac of them.
Positions with a consensus base in both samples are compared.

Outputs (in --outdir):
  summary.txt             settings, number of comparable positions, number of differences
  differences.tsv         positions where the two consensus bases differ (with base counts)
  minor_alleles.tsv       positions where sample B carries a base at >= --minor-frac
                          (and >= --minor-reads reads) that is (almost) absent in sample A
  counts/                 cached per-chunk base counts (lets an interrupted run resume)

Indels are not assessed. BAM files must be sorted and indexed.

Example:
  python genome_wide_compare.py --bam isolate=isolate_nodup.bam --bam metagenome=metagenome_nodup.bam \
      --contig CP033688.1 --outdir gw_isolate_vs_metagenome
"""
import argparse, glob, os
import numpy as np
import pysam

BASES = "ACGT"


def keep_read(r, min_mapq, max_de, min_len):
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
        return False
    if r.mapping_quality < min_mapq:
        return False
    if not r.has_tag("de") or r.get_tag("de") >= max_de:
        return False
    if r.reference_end - r.reference_start < min_len:
        return False
    return True


def count_sample(name, bam_path, contig, outdir, chunk, a):
    """Per-position A/C/G/T counts (shape 4 x L) from filtered reads; cached per chunk."""
    bam = pysam.AlignmentFile(bam_path)
    L = bam.get_reference_length(contig)
    cdir = os.path.join(outdir, "counts")
    os.makedirs(cdir, exist_ok=True)
    cb = lambda r: keep_read(r, a.min_mapq, a.max_de, a.min_aln_len)
    for i, start in enumerate(range(0, L, chunk)):
        f = os.path.join(cdir, f"{name}_{i:04d}.npy")
        if os.path.exists(f):
            continue
        end = min(start + chunk, L)
        c = bam.count_coverage(contig, start, end, quality_threshold=0, read_callback=cb)
        np.save(f, np.array(c, dtype=np.int32))
        print(f"{name}: {end:,}/{L:,} bp", flush=True)
    parts = sorted(glob.glob(os.path.join(cdir, f"{name}_*.npy")))
    return np.concatenate([np.load(p) for p in parts], axis=1)


def consensus(c, min_depth, min_frac):
    tot = c.sum(0)
    ok = (tot >= min_depth) & (c.max(0) >= min_frac * np.maximum(tot, 1))
    return c.argmax(0), ok, tot


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", action="append", required=True, help="NAME=PATH, given exactly twice (A then B)")
    p.add_argument("--contig", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--max-de", type=float, default=0.10, help="max minimap2 per-read divergence (de tag)")
    p.add_argument("--min-aln-len", type=int, default=0, help="min aligned length on the reference (bp)")
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--min-frac", type=float, default=0.8)
    p.add_argument("--minor-frac", type=float, default=0.10)
    p.add_argument("--minor-reads", type=int, default=3)
    p.add_argument("--chunk", type=int, default=250000)
    a = p.parse_args()
    if len(a.bam) != 2:
        p.error("give --bam exactly twice")
    (na, pa), (nb, pb) = [x.split("=", 1) for x in a.bam]
    os.makedirs(a.outdir, exist_ok=True)

    A = count_sample(na, pa, a.contig, a.outdir, a.chunk, a)
    B = count_sample(nb, pb, a.contig, a.outdir, a.chunk, a)
    ca, oka, ta = consensus(A, a.min_depth, a.min_frac)
    cb, okb, tb = consensus(B, a.min_depth, a.min_frac)
    both = oka & okb
    diff = np.where(both & (ca != cb))[0]

    with open(os.path.join(a.outdir, "differences.tsv"), "w") as fh:
        fh.write(f"position\t{na}_base\t{nb}_base\t{na}_ACGT\t{nb}_ACGT\n")
        for i in diff:
            fh.write(f"{i+1}\t{BASES[ca[i]]}\t{BASES[cb[i]]}\t{','.join(map(str, A[:, i]))}\t{','.join(map(str, B[:, i]))}\n")

    fa = A / np.maximum(ta, 1)
    fb = B / np.maximum(tb, 1)
    minor = (fb >= a.minor_frac) & (B >= a.minor_reads) & (fa < 0.02)
    mpos = np.where((ta >= a.min_depth) & (tb >= a.min_depth) & minor.any(0))[0]
    with open(os.path.join(a.outdir, "minor_alleles.tsv"), "w") as fh:
        fh.write(f"position\t{nb}_minor_base\t{nb}_minor_frac\t{na}_ACGT\t{nb}_ACGT\n")
        for i in mpos:
            j = int(np.argmax(np.where(minor[:, i], fb[:, i], -1)))
            fh.write(f"{i+1}\t{BASES[j]}\t{fb[j, i]:.2f}\t{','.join(map(str, A[:, i]))}\t{','.join(map(str, B[:, i]))}\n")

    L = A.shape[1]
    summary = (
        f"Samples: A={na} ({pa}), B={nb} ({pb}); contig {a.contig} ({L:,} bp)\n"
        f"Read filters: primary, MAPQ >= {a.min_mapq}, de < {a.max_de}, aligned length >= {a.min_aln_len} bp\n"
        f"Consensus: >= {a.min_depth} reads and one base >= {a.min_frac:.0%}\n"
        f"Mean filtered depth: {na} {ta.mean():.1f}x, {nb} {tb.mean():.1f}x\n"
        f"Comparable positions (consensus in both): {both.sum():,} ({both.sum()/L:.1%} of reference)\n"
        f"Consensus differences: {len(diff)}\n"
        f"Positions with a minor base only in {nb} (>= {a.minor_frac:.0%}, >= {a.minor_reads} reads): {len(mpos)}\n"
    )
    with open(os.path.join(a.outdir, "summary.txt"), "w") as fh:
        fh.write(summary)
    print(summary)


if __name__ == "__main__":
    main()
