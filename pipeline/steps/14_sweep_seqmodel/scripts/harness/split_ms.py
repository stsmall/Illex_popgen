#!/usr/bin/env python3
"""Split a gzipped ms file into chunks of R blocks each, streaming (low memory).

Two passes: count `//` blocks, then stream blocks into chunk files with a correct
`ms <nhap> <count>` header on each. Prints each chunk path it writes.

Usage: split_ms.py <in.msOut.gz> <outbase> <R> <nhap>
  -> <outbase>.chunk0001.msOut.gz, <outbase>.chunk0002.msOut.gz, ...
"""
import gzip
import sys


def main():
    if len(sys.argv) != 5:
        sys.exit("usage: split_ms.py <in.msOut.gz> <outbase> <R> <nhap>")
    ms_path, outbase, R, nhap = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

    total = 0
    with gzip.open(ms_path, "rt") as fh:
        for line in fh:
            if line.startswith("//"):
                total += 1
    if total == 0:
        sys.exit(f"split_ms.py: no blocks in {ms_path}")

    def chunk_count(c):
        return min(R, total - c * R)

    blk = -1
    cur = -1
    o = None
    with gzip.open(ms_path, "rt") as fh:
        for line in fh:
            if line.startswith("//"):
                blk += 1
                c = blk // R
                if c != cur:
                    if o:
                        o.close()
                    cur = c
                    path = f"{outbase}.chunk{c + 1:04d}.msOut.gz"
                    o = gzip.open(path, "wt", compresslevel=1)
                    o.write(f"ms {nhap} {chunk_count(c)}\n0\n\n")
                    print(path)
                o.write(line)
            elif o is not None:
                o.write(line)
    if o:
        o.close()


if __name__ == "__main__":
    main()
