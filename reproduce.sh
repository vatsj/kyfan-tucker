#!/usr/bin/env bash
# Reproduce every settled number in results.md. The tests ARE the table (tests/*.py, one assertion per number).
#   fast subset (~2 min):   ./reproduce.sh
#   everything (~10 min):   ./reproduce.sh all      (adds the S^3 degree-3 lower bound, z9: dense GF(2), ~7 min)
set -e; cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}; [ -x "$PY" ] || PY=python3
if [ "$1" = all ]; then "$PY" -m pytest -m "" -v; else "$PY" -m pytest -v; fi
