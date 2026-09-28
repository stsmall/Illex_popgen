#!/usr/bin/env bash
# Stairway estimation at higher parallelism (P48), skipping already-done .addTheta, then Step-2 tail.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography
BP=$D/illex_baker633.blueprint.sh
cd "$D"
PAR=48

# 1) kill the v1 orchestrator + its xargs + all running estimation java (by PID)
for p in $(pgrep -f run_stairway_parallel.sh) $(pgrep -f Stairway_fold_training_testing7); do kill "$p" 2>/dev/null; done
sleep 4

# 2) build list of estimation cmds still needing to run (no .addTheta yet)
grep 'Stairway_fold_training_testing7' "$BP" > "$D/stairway_est_all.txt"
: > "$D/stairway_est_todo.txt"
while IFS= read -r line; do
  inp=$(awk '{for(i=1;i<=NF;i++) if($i ~ /input\/illex_baker633-/){print $i; break}}' <<<"$line")
  nr=$(awk '{for(i=1;i<=NF;i++) if($i ~ /input\/illex_baker633-/){print $(i+1); break}}' <<<"$line")
  [[ -s "${inp}.${nr}_0.67.addTheta" ]] || printf '%s\n' "$line" >> "$D/stairway_est_todo.txt"
done < "$D/stairway_est_all.txt"
NTODO=$(wc -l < "$D/stairway_est_todo.txt")
NDONE=$(ls $D/stairway_out/input/*.addTheta 2>/dev/null | wc -l)
echo "[$(date)] $NDONE done, $NTODO todo -> running -P$PAR"

# 3) run remaining estimations in parallel
xargs -d '\n' -P "$PAR" -I CMDLINE bash -c CMDLINE < "$D/stairway_est_todo.txt"
echo "[$(date)] estimation done: $(ls $D/stairway_out/input/*.addTheta 2>/dev/null|wc -l)/800 addTheta"

# 4) Step-2 tail (mv into rand*/, aggregate, summary_plot2, Stairpainter)
STEP2=$(grep -n '# Step 2' "$BP" | head -1 | cut -d: -f1)
tail -n +"$STEP2" "$BP" > "$D/stairway_tail.sh"
bash "$D/stairway_tail.sh"
echo "[$(date)] STAIRWAY COMPLETE"
ls -la $D/stairway_out/*.final.summary $D/stairway_out/*.png $D/stairway_out/*.pdf 2>/dev/null
