#!/usr/bin/env bash
# combine_ms.sh OUT N_HAP FILE [FILE ...]
# Merge per-sim gzipped-ms files into one: a single `ms N_HAP TOTAL` header followed
# by every `//` block from every input (each input's own header is stripped).
# Uses zcat/gzip (C) -- Python's gzip made this take ~1h/position on 240MB inputs.
# Prints TOTAL block count.
set -euo pipefail
out=$1; nhap=$2; shift 2
[ "$#" -eq 0 ] && { printf 'ms %d 0\n0\n\n' "$nhap" | gzip > "$out"; echo 0; exit 0; }
total=0
for f in "$@"; do total=$((total + $(zcat "$f" | grep -c '^//' || true))); done
{ printf 'ms %d %d\n0\n\n' "$nhap" "$total"
  for f in "$@"; do zcat "$f" | awk '/^\/\//{p=1} p'; done
} | ${COMPRESS:-gzip -1} > "$out"
echo "$total"
