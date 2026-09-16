#!/usr/bin/env bash
# Fast reproduction (~5 min). Compare each block against results.md.
set -e; cd "$(dirname "$0")/scripts"
echo "== build.py (expect: Euler=2; octahedron A+ histogram {1: 192}; Tucker solutions 0; Ky Fan all odd)"; python3 build.py
echo; echo "== nsdeg.py (expect: d=2 no refutation both fields; d=3 F_2 True, F_p False)"; python3 nsdeg.py
echo; echo "== pseudo2b.py (expect: unrestricted 597 free; distinct-only NO; label-sym NO; hemisphere NO)"; python3 pseudo2b.py 2>&1 | grep -v '^complex\|^free'
echo; echo "== tower.py (expect: 0 violations of local lemma; identity 400/400 for m=3,4,5)"; python3 tower.py
echo; echo "Slow (run separately): primal.py (~1 min), z9a.py + z9b.py (~15 min, S^3 lower bound)"
