# Estimation version of Ky Fan: what a low-degree estimator outputs without advice

Question (the note's "next step"): draw an antipodal labeling of S^n with k ≥ n+1 label magnitudes uniformly at random. Ky Fan's
lemma says E[(−1)^A | valid] = −1 exactly, where A counts positively alternating top simplices. Equivalently, with
v_e = 1[edge e is not complementary] and s_σ = (−1)^[σ alternates, either sign],

    N := E[ Π_e v_e · Π_σ s_σ ] = − E[ Π_e v_e ] =: −Z .

Both are expectations of products of local factors of independent uniform variables (one label per antipodal class).

## What was run

- `enum_gen.py`, `enum.c`: exact enumeration on S^2 with k = 3 (13 variables, 6 labels). 22,193,664 of 6^13 labelings are valid
  (Z = 1.699271e-3); A is odd on every one of them (histogram A = 1, 3, …, 15), mean 2.21, variance 2.30.
- `lowdeg.py`: a generic low-degree deduction–projection estimator. g is kept as a sum of ANOVA components on at most d
  variables; each step multiplies in the next local factor and projects back (orthogonally in L^2) onto degree ≤ d; the estimate
  is the constant term. It is exact at d = number of variables (checked on S^1 by brute force). This is my instantiation of the
  deduction–projection scheme, with "low degree" as the sketch; it is not ARC's code, and it uses no advice.
- `run_s2.py`, `run_s3.py`: the runs below.

## Results

S^2, k = 3 (truth: N/Z = −1):

| d | Z estimate / exact Z | N estimate / Z estimate (edges first) | simplices first | shuffled order |
|---|---|---|---|---|
| 0, 1 | 0.83 | +0.059 | +0.059 | +0.059 |
| 2 | 1.08 | +0.116 | +0.117 | +0.098 |
| 3 | 1.03 | +0.035 | +0.035 | +0.034 |
| 4 | 1.00 | −0.072 | −0.077 | −0.019 |
| 5 | 1.00 | −0.184 | −0.184 | −0.092 |
| 6 | 1.00 | −0.314 | −0.313 | −0.216 |

S^3, k = 4 (40 variables; truth −1): ratio +0.0105 at d = 0, 1 and +0.051 at d = 2 (`s3.out`).

Reading: the probability of validity is right to 0.3% by d = 4, while the conditional parity is essentially 0 and has the wrong
sign through d = 3; it moves toward −1 only slowly as d approaches the number of variables (−0.31 at d = 6 of 13; that run takes 12 minutes). Full independence gives
Π_σ E[s_σ] = (1 − 12/216)^24 = +0.059. Cumulants of A do not help either: with the exact mean and variance of A among valid
labelings a Gaussian treatment gives magnitude e^{−π²·2.30/2} ≈ 1e-5, and higher truncations diverge.

Caveats: one small sphere with exact truth, one more with the truth known from the lemma; one estimator family; and the
unconditioned quantity N is itself exponentially small, so this is a question about relative accuracy in the low-probability
regime (an estimator that outputs 0 for N already matches sampling at any feasible budget).
