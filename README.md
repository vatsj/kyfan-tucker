# kyfan-tucker

Nullstellensatz / Sherali–Adams degree of Tucker's lemma and Ky Fan's lemma on the subdivided cross-polytope, over F_2 and Q — as a toy instance of "lemma-shaped advice" for heuristic estimators.

- `CLAUDE.md` — full brief: definitions, settled results, theorems, tasks. Read first.
- `results.md` — every number with the test that produces it and its runtime.
- `review.md` — Phase-0 review of the original scripts (findings, conventions check).
- `kyfan/` — the package (complex parametrized by m, labelings, group actions, SA-dual builder, primal, solvers).
- `tests/` — the settled table as assertions. `./reproduce.sh` runs the fast subset (~2 min), `./reproduce.sh all` everything (~15 min).

Setup: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt` (numpy, pytest). Python ≥ 3.10.
