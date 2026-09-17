#!/usr/bin/env bash
# Run the test tiers (the tests are the settled table: HISTORY.md, docs/README.md).
#   fast + slow (~5 min):   ./reproduce.sh
#   everything (~20 min):   ./reproduce.sh all      (adds the veryslow tier: S^3 degree-3 dense GF(2) lower bound ~7 min, ...)
# Equivalent make targets: make test-fast / test-slow / test-veryslow / reproduce-paper.
set -e; cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}; [ -x "$PY" ] || PY=python3
if [ "$1" = all ]; then PYTHONPATH=. "$PY" -m pytest -m "" -v; else PYTHONPATH=. "$PY" -m pytest -v; fi
