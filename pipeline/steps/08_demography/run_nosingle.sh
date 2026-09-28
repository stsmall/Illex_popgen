#!/usr/bin/env bash
# No-singleton Stairway run (smallest_size_of_SFS_bin=2), clean single-writer pipeline.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography
BP=$D/illex_baker633_nosingle.blueprint.sh
SWE=/home/ssmall/programs/stairway/stairway_plot_v2.1.1/stairway_plot_es
CP="$SWE:$SWE/gral-core-0.11.jar:$SWE/swarmops.jar:$SWE/VectorGraphics2D-0.9.3.jar"
export _JAVA_OPTIONS="-Djava.awt.headless=true"
cd "$D"
mkdir -p $D/stairway_out_nosingle/final
STEP2=$(grep -n '# Step 2' "$BP" | head -1 | cut -d: -f1)
grep 'Stairway_fold_training_testing7' "$BP" > "$D/nosingle_est.txt"
NEST=$(wc -l < "$D/nosingle_est.txt")
echo "[$(date)] nosingle: $NEST estimations -P48 (single xargs)"
xargs -d '\n' -P 48 -I CMDLINE bash -c CMDLINE < "$D/nosingle_est.txt"
# verify none truncated
bad=$(for f in $D/stairway_out_nosingle/input/*.addTheta; do grep -q 'final model' "$f" || echo x; done | wc -l)
echo "[$(date)] estimations done; truncated=$bad"
# Step-2 tail (mv into rand dirs) then Stairpainter + plot.sh
tail -n +"$STEP2" "$BP" > "$D/nosingle_tail.sh"
bash "$D/nosingle_tail.sh" > "$D/logs/nosingle_tail.log" 2>&1
bash "$D/illex_baker633_nosingle.blueprint.plot.sh" > "$D/logs/nosingle_plot.log" 2>&1
echo "[$(date)] NOSINGLE COMPLETE"
wc -l $D/stairway_out_nosingle/Illex_baker633_nosingleton.final.summary 2>/dev/null
