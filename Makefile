# kyfan-tucker: test tiers and the paper-cited computations.  Uses .venv/bin/python if present (see README.md "Install").
PY ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
export PYTHONPATH := .

.PHONY: test-fast test-slow test-veryslow test-all reproduce-paper gf2solve

gf2solve:            ## build the sparse GF(2) solver (needed by test_sparse, test_realize*, test_rule_vs_solver)
	cargo build --release --manifest-path gf2solve/Cargo.toml

test-fast:           ## ~1.5 min: every proof check and paper number that fits (random-labeling certificate checks)
	$(PY) -m pytest -m "not slow and not veryslow" -q --durations=10

test-slow:           ## ~6 min: exhaustive certificate checks, CP-SAT optimality of 304, S^3 sparse dual, n = 7 tree shapes
	$(PY) -m pytest -m "slow" -q --durations=10

test-veryslow:       ## ~20 min: S^3 degree-3 dense GF(2) lower bound (z9), S^2 Q degree-5 certificate, realization m = 10, m = 8 sphere verify
	$(PY) -m pytest -m "veryslow" -q --durations=10

test-all: test-fast test-slow test-veryslow

reproduce-paper:     ## the three paper-cited computations + the five proof checks (fast + slow tiers of those files)
	$(PY) -m pytest -m "not veryslow" -q --durations=0 \
	    tests/test_paper_numbers.py tests/test_min_size.py tests/test_flag_hitting.py tests/test_spread_mcmc.py \
	    tests/test_extension_lemma.py tests/test_rule_vs_solver.py tests/test_realize_explicit.py tests/test_tower.py
