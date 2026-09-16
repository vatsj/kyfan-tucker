# kyfan-tucker

Nullstellensatz / Sherali–Adams degree of Tucker's lemma and Ky Fan's lemma on the subdivided cross-polytope, over F_2 and Q — as a toy instance of "lemma-shaped advice" for heuristic estimators.

- `CLAUDE.md` — full brief: definitions, settled results, theorem, tasks. Read first.
- `results.md` — every number with the command that produced it.
- `scripts/` — the computations. Run from inside `scripts/` (they `exec()` each other by relative path; refactoring into a package is Task 0 in the brief).
- `reproduce.sh` — runs the fast subset and prints the lines to compare against `results.md`.

Requirements: Python ≥ 3.10, numpy. Nothing else.
