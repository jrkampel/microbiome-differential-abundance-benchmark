# Simulated route, end to end. Real-data route is in the README.
PY := python3
export PYTHONPATH := src

.PHONY: all simulate decontam ordination da random sim test clean

all: simulate decontam ordination da random sim

simulate:
	$(PY) scripts/01_simulate_dataset.py

decontam:
	$(PY) scripts/02_decontam.py

ordination:
	$(PY) scripts/03_ordination.py

da:
	$(PY) scripts/04_differential_abundance.py

random:
	$(PY) scripts/05_random_effects.py

sim:
	$(PY) scripts/06_simulation_study.py --replicates 4 --donors 30

test:
	$(PY) -m pytest tests -q

clean:
	rm -rf results/figures/* results/tables/* data/processed/*.tsv
