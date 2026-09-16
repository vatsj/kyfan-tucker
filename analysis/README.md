# analysis/ — exploratory drivers (Tasks 2, 3 and the restriction gadgets)

Run from the repo root with `PYTHONPATH=. .venv/bin/python analysis/<script> ...`. Each `.out` file is the output of the
command that produced it (`.pkl` files are ignored by git; regenerate by rerunning). Everything here is exploratory — the
settled numbers live in `tests/`; results.md points at the relevant `.out` files.

| script | command | output | what it shows |
|---|---|---|---|
| `s3_degree3_types.py` | `analysis/s3_degree3_types.py` (16 min) | `s3_degree3_types.out` | 808 full-group types of size-3 partial labelings on S^3; knockout: 4 necessary types |
| `s3_degree3_marginals.py` | (3 min) | `s3_degree3_marginals.out` | marginal knockout on size-1, size-2 types |
| `s3_degree3_restrictions.py` | | `s3_degree3_restrictions.out` | label-sym / hemisphere / distinct-magnitude restrictions (all NO) |
| `s3_degree3_greedy.py` | (20 min, needs `s3_degree3_types.pkl`) | `s3_degree3_greedy.out` | greedy minimal type family (379 of 808) |
| `class_knockout.py` | `analysis/class_knockout.py 3` / `4` | `class_knockout_m3.out`, `class_knockout_m4.out`, `class_knockout_m4_invariant.out`, `class_necessary_suffice_m4.out` | class-level (simplex?, #edges, AAA/AAB/ABC) knockouts and minimal families |
| inline (see results.md) | | `simplex_support.out`, `magnitude_support.out`, `s2_degree2_greedy.out` | further support restrictions at n=2,3 |
| `base_labeling.py`, `s3_degree3_base.py` | `analysis/base_labeling.py 3 6` / `4 4` | `base_labeling_m3.out`, `base_labeling_m4.out`, `base_labeling_m4_random.out`, `s3_degree3_base.out` | "base labeling + corrections" supports (NO) |
| inline | | `ball_m3.out`, `ball_m4.out` | Task 3 ball degrees (m=3: all 80 L_eq; m=4: 45 L_eq) |
| inline | | `residual_m3.out`, `residual_m3_small.out`, `residual_m4_chains.out`, `residual_m4_Wv.out`, `residual_m5_Wv.out` | restriction experiments leading to the chain gadget |
| `gadget.py` | `analysis/gadget.py m ntries` | `gadget_m5.out` | first chain-gadget construction (random equatorial labelings; fails at m=5) |
| `gadget2.py` | `analysis/gadget2.py m ntries "1:1"` | `gadget2_m6.out` | chain gadget with a forbidden value on the star of sigma's rank-1 vertex (found the m=5 gadget) |
| `gadget3.py` | `analysis/gadget3.py 6 3 "(1;(2;x,x),(3;x,(4;x,x)))" "5,4,3,2,1"` | `gadget3_m6.out`, `gadget3_m7.out`, `gadget3_m8.out` | **realizes a target tree gadget at any m; the S^5, S^6, S^7 lower bounds** |
| `verify_gadget.py` | `analysis/verify_gadget.py 5 0` (needs `gadget2_m5.pkl`) | `verify_gadget_m5_0.out` | the first independent verification of the S^4 gadget (superseded by `kyfan.gadget.verify`) |
| inline | | `abstract_trees.out`, `mu_structure.out` | abstract tree gadgets (degree n+1, unique pseudo-solution, support 3^n); a global mu of support 3^(n-1) generating it |
