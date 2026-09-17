# kyfan-tucker

Companion repository for the paper on the F_2 Nullstellensatz / Sherali–Adams degree of Tucker and Ky Fan. Start with
`README.md`, then `docs/README.md` (theorem → proof file → test → runtime). `HISTORY.md` is the development log; the original
working brief is `docs/AGENT_BRIEF.md`.

Conventions for working in this repository: use the virtual environment (`python3 -m venv .venv && .venv/bin/pip install -r
requirements.txt`); `make test-fast` before committing; `gf2solve/` is built with `cargo` on first use. Do not alter settled
numbers — every number in `docs/` and `README.md` is produced by a test in `tests/`.
