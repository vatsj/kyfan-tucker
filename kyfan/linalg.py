"""Linear algebra over F_2 and F_p.

Row format for the sparse builders: list of (coeffs: dict col -> int, rhs: int).
- gf2_dense: bit-packed dense Gauss-Jordan in numpy (rows x ncols/64 words). Fine to ~40k unknowns; a real sparse
  GF(2) solver is Task 1 and should plug in behind the same `GF2Result` interface.
- fp_sparse: dict-based sparse elimination mod p (triangular, no fill control).
- gf2_bitset_solve / fp_dense_solve: for the random-point primal (rows are labelings, dense-ish).
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class GF2Result:
    consistent: bool
    rank: int
    pivots: list            # pivot columns, in order
    ncols: int
    solution: object = None  # np.uint8 array (particular solution, free variables = 0) or None if inconsistent

    @property
    def n_free(self):
        """Dimension of the solution space (number of free columns) when consistent."""
        return self.ncols - self.rank


def _pack_rows(rows, ncols):
    R = len(rows)
    W = (ncols + 1 + 63) // 64
    M = np.zeros((R, W), dtype=np.uint64)
    one = np.uint64(1)
    for r, (f, rhs) in enumerate(rows):
        for c, x in f.items():
            if x % 2:
                M[r, c // 64] ^= one << np.uint64(c % 64)
        if rhs % 2:
            M[r, ncols // 64] ^= one << np.uint64(ncols % 64)
    return M


def pack_bit_rows(bit_rows, ncols):
    """Pack 0/1 rows (each of length ncols+1, last entry = rhs) into the uint64 word layout used by gf2_dense_packed:
    bit c of a row lives in word c//64 at position c%64 (little-endian host assumed, as in z9a.py)."""
    W = (ncols + 1 + 63) // 64
    M = np.zeros((len(bit_rows), W), dtype=np.uint64)
    for r, bits in enumerate(bit_rows):
        packed = np.packbits(np.asarray(bits, dtype=np.uint8), bitorder='little').tobytes().ljust(W * 8, b'\0')
        M[r] = np.frombuffer(packed, dtype=np.uint64)
    return M


def gf2_dense(rows, ncols, want_solution=False):
    """Gauss-Jordan over F_2 on rows (dict, rhs). rhs is column `ncols`. Returns GF2Result."""
    return gf2_dense_packed(_pack_rows(rows, ncols), ncols, want_solution)


def gf2_dense_packed(M, ncols, want_solution=False):
    """Gauss-Jordan on a packed uint64 matrix (modified in place); rhs is bit `ncols`."""
    R = M.shape[0]
    one = np.uint64(1)
    r = 0
    pivots = []
    for c in range(ncols):
        if r == R:
            break
        w, b = c // 64, np.uint64(c % 64)
        colbits = (M[r:, w] >> b) & one
        nz = np.nonzero(colbits)[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        mask = ((M[:, w] >> b) & one).astype(bool)
        mask[r] = False
        idx = np.nonzero(mask)[0]
        if len(idx):
            M[idx] ^= M[r]
        pivots.append(c)
        r += 1
    rhs_rest = (M[r:, ncols // 64] >> np.uint64(ncols % 64)) & one
    consistent = not bool(rhs_rest.any())
    sol = None
    if consistent and want_solution:
        sol = np.zeros(ncols, dtype=np.uint8)
        for k, c in enumerate(pivots):
            sol[c] = int((M[k, ncols // 64] >> np.uint64(ncols % 64)) & one)
    return GF2Result(consistent, r, pivots, ncols, sol)


def fp_sparse(rows, ncols, p=1000003):
    """Sparse elimination mod p. Returns (consistent, rank)."""
    piv = {}   # leading col -> (row dict normalized, rhs)
    rank = 0
    for f, rhs in rows:
        f = {c: x % p for c, x in f.items() if x % p}
        rhs %= p
        while f:
            c = min(f)
            if c in piv:
                pf, prhs = piv[c]
                x = f[c]
                for cc, y in pf.items():
                    v = (f.get(cc, 0) - x * y) % p
                    if v:
                        f[cc] = v
                    else:
                        f.pop(cc, None)
                rhs = (rhs - x * prhs) % p
            else:
                inv = pow(f[c], p - 2, p)
                piv[c] = ({cc: (y * inv) % p for cc, y in f.items()}, (rhs * inv) % p)
                rank += 1
                break
        else:
            if rhs:
                return False, rank
    return True, rank


def gf2_bitset_solve(rows_bits, rhs, ncols=None):
    """rows_bits: python-int bitsets. Solve A x = rhs over F_2. Returns x as an int bitset (free vars 0) or None."""
    piv = {}
    for r, b in zip(rows_bits, rhs):
        b &= 1
        while r:
            c = r.bit_length() - 1
            if c in piv:
                r ^= piv[c][0]
                b ^= piv[c][1]
            else:
                piv[c] = (r, b)
                break
        else:
            if b:
                return None
    x = 0
    for c in sorted(piv):
        r, b = piv[c]
        val = b ^ (bin(r & x & ~(1 << c)).count('1') & 1)
        if val:
            x |= 1 << c
    return x


def fp_dense_solve(A, b, p=1000003):
    """Dense Gauss-Jordan mod p. A: (R,C) int64 with entries < p. Returns particular solution or None."""
    A = A.copy() % p
    b = b.copy() % p
    R, C = A.shape
    r = 0
    pivcols = []
    for c in range(C):
        if r == R:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        pi = r + nz[0]
        if pi != r:
            A[[r, pi]] = A[[pi, r]]
            b[[r, pi]] = b[[pi, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        b[r] = (b[r] * inv) % p
        f = A[:, c].copy()
        f[r] = 0
        A = (A - np.outer(f, A[r])) % p
        b = (b - f * b[r]) % p
        pivcols.append(c)
        r += 1
    if np.any(b[r:] != 0):
        return None
    x = np.zeros(C, dtype=np.int64)
    for k, c in enumerate(pivcols):
        x[c] = b[k]
    return x


def real_lstsq_residual(A, b):
    """Relative residual of least squares over R: < 1e-9 counts as a certificate, > 1e-3 as none (CLAUDE.md §2)."""
    x, *_ = np.linalg.lstsq(A.astype(float), b.astype(float), rcond=None)
    return float(np.linalg.norm(A @ x - b) / np.sqrt(len(b)))
