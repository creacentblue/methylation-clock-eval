#!/usr/bin/env bash
# One-command rerun of the whole pipeline. Usage: ./run_all.sh   (set SKIP_DOWNLOAD=1 to reuse data/raw)
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d .venv ]; then
  python3 -m venv .venv
  .venv/bin/pip install -r requirements.txt
fi
PY=.venv/bin/python
[ "${SKIP_DOWNLOAD:-0}" = "1" ] || scripts/01_download.sh
$PY scripts/02_clock_features.py
$PY scripts/02b_cell_reference.py       # 450K cell-type reference from Reinius 2012 (GSE35069)
$PY scripts/03_parse.py
$PY scripts/04_run_clocks.py
$PY scripts/05_evaluate.py
$PY scripts/06_train_elasticnet.py   # optional stretch
$PY scripts/07_sample_identity.py    # GSE87571 732 vs 421: repeat-individual check, Sample_Group, dedup flag
$PY scripts/08a_estimate_cells.py    # Houseman-style cell proportions
$PY scripts/08_cell_composition.py   # acceleration vs cells, cell-adjusted comparison, rankings
$PY scripts/09_subset_metrics.py     # GSE87571 dedup / batch subset metrics
echo "Done. See results/metrics.csv, results/own_model_metrics.csv, results/plots/"
