#!/usr/bin/env bash
# Re-run only the corrupted .addTheta (missing 'final model' = truncated) estimations, single-writer,
# then mv into rand dirs, verify, and aggregate with Stairpainter (full classpath, headless).
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography
BP=$D/illex_baker633.blueprint.sh
SWE=/home/ssmall/programs/stairway/stairway_plot_v2.1.1/stairway_plot_es
CP="$SWE:$SWE/gral-core-0.11.jar:$SWE/swarmops.jar:$SWE/VectorGraphics2D-0.9.3.jar"
cd "$D"

# 1) collect commands for corrupted (N,nrand) pairs
: > "$D/repair_cmds.txt"
for r in 316 632 948 1264; do
  for f in $D/stairway_out/rand$r/*.addTheta; do
    grep -q 'final model' "$f" && continue
    base=$(basename "$f"); N=$(echo "$base" | sed 's/illex_baker633-//; s/\..*//')
    rm -f "$f" "$D/stairway_out/input/${base}"        # remove truncated copy (both locations)
    grep -F "input/illex_baker633-${N} ${r} 0.67" "$BP" >> "$D/repair_cmds.txt"
  done
done
NREP=$(wc -l < "$D/repair_cmds.txt")
echo "[$(date)] re-running $NREP corrupted estimations (-P24)"

# 2) re-run (single xargs, P24)
xargs -d '\n' -P 24 -I CMDLINE bash -c CMDLINE < "$D/repair_cmds.txt"

# 3) mv fresh outputs into rand dirs
for r in 316 632 948 1264; do mv -f $D/stairway_out/input/*.${r}_0.67.addTheta $D/stairway_out/rand$r/ 2>/dev/null; done

# 4) verify no more truncated files
bad=0; for r in 316 632 948 1264; do for f in $D/stairway_out/rand$r/*.addTheta; do grep -q 'final model' "$f" || bad=$((bad+1)); done; done
echo "[$(date)] still missing final model: $bad"
[[ "$bad" -eq 0 ]] || { echo "STILL CORRUPTED — abort"; exit 1; }

# 5) aggregate + plot
mkdir -p $D/stairway_out/final
java -Djava.awt.headless=true -cp "$CP" Stairpainter $D/illex_baker633.blueprint > $D/logs/stairpainter.log 2>&1
echo "[$(date)] Stairpainter exit $?"
tail -3 $D/logs/stairpainter.log
wc -l $D/stairway_out/*.final.summary 2>/dev/null
