PY=.venv/bin/python
.PHONY: all download features cellref parse clocks evaluate train identity cells cellanalysis subsets env clean-results
all: env download features cellref parse clocks evaluate train identity cells cellanalysis subsets
env:
	test -d .venv || (python3 -m venv .venv && .venv/bin/pip install -r requirements.txt)
download:
	scripts/01_download.sh
features:
	$(PY) scripts/02_clock_features.py
cellref:
	$(PY) scripts/02b_cell_reference.py
parse:
	$(PY) scripts/03_parse.py
clocks:
	$(PY) scripts/04_run_clocks.py
evaluate:
	$(PY) scripts/05_evaluate.py
train:
	$(PY) scripts/06_train_elasticnet.py
identity:
	$(PY) scripts/07_sample_identity.py
cells:
	$(PY) scripts/08a_estimate_cells.py
cellanalysis:
	$(PY) scripts/08_cell_composition.py
subsets:
	$(PY) scripts/09_subset_metrics.py
clean-results:
	rm -f results/*.csv results/plots/*.png
