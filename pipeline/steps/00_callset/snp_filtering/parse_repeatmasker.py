#!/usr/bin/env python3
"""
Parse RepeatMasker .out or .gff into a classified 4-column BED.

Output columns: chrom  start  end  repeat_class
Repeat classes: Simple, TE, rRNA, NUMT, Other

Usage:
    python3 parse_repeatmasker.py \
        --input genome.fa.out \
        --output repeats_classified.bed

    python3 parse_repeatmasker.py \
        --input genome.fa.out.gff \
        --output repeats_classified.bed
"""

import argparse
import re
import sys
from collections import Counter


# ─────────────────────────────────────────────
# Classification patterns
# ─────────────────────────────────────────────

TE_PATTERNS = re.compile(
    r'LINE|SINE|LTR|DNA|RC|Retroposon|Transposon|Helitron|TIR|MITE|PLE|Maverick',
    re.IGNORECASE
)
SIMPLE_PATTERNS = re.compile(
    r'Simple_repeat|Low_complexity|Satellite|Tandem|Microsatellite',
    re.IGNORECASE
)
RRNA_PATTERNS = re.compile(
    r'\brRNA\b|snRNA|scRNA|srpRNA|tRNA|Small_RNA',
    re.IGNORECASE
)
NUMT_PATTERNS = re.compile(
    r'NUMT|mitochondri|^MT|chrMT',
    re.IGNORECASE
)


def classify(repeat_class: str, repeat_family: str) -> str:
    combined = f"{repeat_class}/{repeat_family}"
    if NUMT_PATTERNS.search(combined):
        return "NUMT"
    if RRNA_PATTERNS.search(combined):
        return "rRNA"
    if SIMPLE_PATTERNS.search(combined):
        return "Simple"
    if TE_PATTERNS.search(combined):
        return "TE"
    return "Other"


# ─────────────────────────────────────────────
# GFF parser
# ─────────────────────────────────────────────

def parse_gff(fh, out):
    n = 0
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip().split("\t")
        if len(fields) < 9:
            continue
        chrom = fields[0]
        start = int(fields[3]) - 1   # GFF 1-based → BED 0-based
        end   = int(fields[4])
        attrs = fields[8]

        repeat_class  = fields[2]
        repeat_family = ""

        # RepeatMasker GFF: Target "Motif:AluSx" or "Repeat:LINE/L1"
        m = re.search(r'Target[= ]"?([^"\s]+)"?', attrs)
        if m:
            parts = m.group(1).replace("Motif:", "").replace("Repeat:", "").split("/")
            if len(parts) >= 2:
                repeat_class  = parts[0]
                repeat_family = parts[1]
            else:
                repeat_class = parts[0]

        cls = classify(repeat_class, repeat_family)
        out.write(f"{chrom}\t{start}\t{end}\t{cls}\n")
        n += 1
    return n


# ─────────────────────────────────────────────
# .out parser
# ─────────────────────────────────────────────

def parse_out(fh, out):
    """
    RepeatMasker .out fixed-width columns (1-indexed):
      1:  SW score
      2:  % divergence
      3:  % deletions
      4:  % insertions
      5:  query sequence (chrom)
      6:  query start
      7:  query end
      8:  query left (in parens)
      9:  strand (+/C)
      10: repeat name
      11: repeat class/family
      12-14: repeat position info
      15: ID
    """
    n = 0
    for i, line in enumerate(fh):
        if i < 3:   # skip header lines
            continue
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) < 11:
            continue

        chrom  = fields[4]
        start  = int(fields[5]) - 1   # 1-based → 0-based
        end    = int(fields[6])

        repeat_class_family = fields[10]   # e.g. "DNA/TcMar-Tc1" or "Simple_repeat"
        parts         = repeat_class_family.split("/")
        repeat_class  = parts[0]
        repeat_family = parts[1] if len(parts) > 1 else ""

        cls = classify(repeat_class, repeat_family)
        out.write(f"{chrom}\t{start}\t{end}\t{cls}\n")
        n += 1
    return n


# ─────────────────────────────────────────────
# Format detection
# ─────────────────────────────────────────────

def detect_format(path: str) -> str:
    """Return 'gff' or 'out' based on file content."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                # RepeatMasker GFF typically starts with ##gff-version
                if "gff" in line.lower():
                    return "gff"
                continue
            # First non-comment line: GFF has 9 tab-separated fields
            fields = line.rstrip().split("\t")
            if len(fields) == 9:
                return "gff"
            # .out has space-separated fields, first real data line ~15 fields
            return "out"
    return "out"


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Parse RepeatMasker output into a classified BED file"
    )
    parser.add_argument("--input", required=True,
                        help="RepeatMasker .out or .gff file")
    parser.add_argument("--output", required=True,
                        help="Output classified BED (chrom start end class)")
    parser.add_argument("--format", choices=["gff", "out", "auto"], default="auto",
                        help="Input format (default: auto-detect)")
    args = parser.parse_args()

    fmt = args.format
    if fmt == "auto":
        fmt = detect_format(args.input)
        print(f"  Detected format: {fmt}")

    counts: Counter = Counter()
    with open(args.input) as fh, open(args.output, "w") as out:
        if fmt == "gff":
            n = parse_gff(fh, out)
        else:
            n = parse_out(fh, out)

    # Count classes in output for QC
    with open(args.output) as fh:
        for line in fh:
            cls = line.rstrip().split("\t")[-1]
            counts[cls] += 1

    print(f"  Repeat intervals written: {n}")
    for cls, count in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"    {cls}: {count}")
    print(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
