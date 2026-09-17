#!/usr/bin/env bash
# Run the test tiers (the tests are the settled table: HISTORY.md, docs/README.md).
#   fast + slow (~8 min):   ./reproduce.sh
#   everything (~50 min):   ./reproduce.sh all      (adds the veryslow tier, ~40 min: S^7 chain gadget, z9 dense GF(2), Q degree 5, ...)
# Equivalent make targets: make test-fast / test-slow / test-veryslow / reproduce-paper.
set -e; cd "$(dirname "$0")"
PY=${PY:-.venv/bin/python}; [ -x "$PY" ] || PY=python3
if [ "$1" = all ]; then PYTHONPATH=. "$PY" -m pytest -m "" -v; else PYTHONPATH=. "$PY" -m pytest -v; fi
