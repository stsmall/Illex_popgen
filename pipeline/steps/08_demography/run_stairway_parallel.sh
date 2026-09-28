#!/usr/bin/env bash
# Parallelize Stairway Plot 2 estimation: the 800 independent Stairway_fold_training_testing7
# calls run at -P16 (16 cores, -Xmx1g each), then run the Step-2 tail (mv + aggregate + plot) serially.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography
BP=$D/illex_baker633.blueprint.sh
cd "$D"

# 1) kill any running serial blueprint run + its java (by PID; pattern lives only in this file, no self-match)
for p in $(pgrep -f illex_baker633.blueprint.sh) $(pgrep -f Stairway_fold_training_testing7); do kill "$p" 2>/dev/null; done
sleep 3

# 2) split: estimation commands (Step 1) vs the Step-2 tail (mv/aggregate/plot)
STEP2=$(grep -n '# Step 2' "$BP" | head -1 | cut -d: -f1)
grep 'Stairway_fold_training_testing7' "$BP" > "$D/stairway_est_cmds.txt"
tail -n +"$STEP2" "$BP" > "$D/stairway_tail.sh"
NEST=$(wc -l < "$D/stairway_est_cmds.txt")
echo "[$(date)] $NEST estimation cmds -> running -P16; tail from line $STEP2"

# 3) run estimations in parallel (16 concurrent)
xargs -d '\n' -P 16 -I CMDLINE bash -c CMDLINE < "$D/stairway_est_cmds.txt"
echo "[$(date)] estimation done: $(ls $D/stairway_out/input/*.addTheta 2>/dev/null|wc -l)/$NEST addTheta"

# 4) run the Step-2 tail (mv into rand*/, aggregate, Stairway_output_summary_plot2, Stairpainter)
bash "$D/stairway_tail.sh"
echo "[$(date)] STAIRWAY COMPLETE"
ls -la $D/stairway_out/*.final.summary $D/stairway_out/*.png $D/stairway_out/*.pdf 2>/dev/null
