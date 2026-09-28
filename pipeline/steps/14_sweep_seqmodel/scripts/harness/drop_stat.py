#!/usr/bin/env python3
"""drop_stat.py -- remove all <stat>_win* columns from a diploSHIC fvec, preserving order.

Used to drop diplo_ZnS (the LD statistic that is computed ~40x differently between the
fvecSim training path and the fvecVcf scan path -- see diploshic-zns-train-scan-mismatch).
Applied CONSISTENTLY to training fvecs (before makeTrainingSets/train) and to scan/calib
fvecs (before predict) so the model and inputs share the identical reduced feature layout.

Usage: drop_stat.py <in.fvec> <out.fvec> <stat1>[,<stat2>,...]
"""
import sys


def main():
    inp, out, stats = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    with open(inp) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        keep = [i for i, c in enumerate(header)
                if not any(c == s or c.startswith(s + "_win") for s in stats)]
        dropped = len(header) - len(keep)
        with open(out, "w") as o:
            o.write("\t".join(header[i] for i in keep) + "\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                o.write("\t".join(f[i] for i in keep) + "\n")
    sys.stderr.write(f"[drop] {inp}: dropped {dropped} cols ({stats}) -> "
                     f"{len(keep)} cols -> {out}\n")


if __name__ == "__main__":
    main()
