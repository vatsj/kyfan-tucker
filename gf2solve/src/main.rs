//! Sparse Gaussian elimination over GF(2) with Markowitz pivoting.
//!
//! Input (little-endian binary):  u64 nrows, u64 ncols, u64 nnz, i64 indptr[nrows+1], u32 indices[nnz], u8 rhs[nrows].
//! Column indices within a row need not be sorted or unique (duplicates cancel mod 2).
//! Output (text): `consistent 0|1`, `rank r`, `nfree f`, then optionally
//!   `solution c1 c2 ...`   (support of a particular solution, free variables = 0)
//!   `pivots c1 c2 ...`     (pivot columns in elimination order)
//!   `null k: c1 c2 ...`    (one line per null-space basis vector; free column k set to 1)
//!
//! usage: gf2solve IN OUT [--solution] [--pivots] [--nullspace] [--stop-early] [--quiet]
//!
//! Rows are sorted Vec<u32>; a column -> rows index is kept with lazy invalidation; the pivot column is the
//! active column of minimum count (heap with stale entries), and the pivot row is the shortest row in it.
//! This is plain structured Gaussian elimination: fill-in is the enemy, so low-count columns go first.

use std::collections::BinaryHeap;
use std::cmp::Reverse;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::time::Instant;

fn read_u64(r: &mut impl Read) -> u64 {
    let mut b = [0u8; 8];
    r.read_exact(&mut b).unwrap();
    u64::from_le_bytes(b)
}

fn read_vec<T: Copy, const N: usize>(r: &mut impl Read, n: usize, f: fn([u8; N]) -> T) -> Vec<T> {
    let mut buf = vec![0u8; n * N];
    r.read_exact(&mut buf).unwrap();
    buf.chunks_exact(N).map(|c| f(c.try_into().unwrap())).collect()
}

/// symmetric difference of two sorted vectors
fn xor_rows(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            out.push(a[i]);
            i += 1;
        } else if a[i] > b[j] {
            out.push(b[j]);
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: gf2solve IN OUT [--solution] [--pivots] [--nullspace] [--stop-early] [--quiet]");
        std::process::exit(2);
    }
    let want_solution = args.iter().any(|a| a == "--solution");
    let want_pivots = args.iter().any(|a| a == "--pivots");
    let want_null = args.iter().any(|a| a == "--nullspace");
    let stop_early = args.iter().any(|a| a == "--stop-early");
    let quiet = args.iter().any(|a| a == "--quiet");
    let t0 = Instant::now();

    // ---- read ----
    let mut f = BufReader::new(File::open(&args[1]).unwrap());
    let nrows = read_u64(&mut f) as usize;
    let ncols = read_u64(&mut f) as usize;
    let nnz = read_u64(&mut f) as usize;
    let indptr: Vec<i64> = read_vec(&mut f, nrows + 1, i64::from_le_bytes);
    let indices: Vec<u32> = read_vec(&mut f, nnz, u32::from_le_bytes);
    let mut rhs: Vec<u8> = read_vec(&mut f, nrows, |b: [u8; 1]| b[0]);
    drop(f);

    let mut rows: Vec<Vec<u32>> = Vec::with_capacity(nrows);
    for r in 0..nrows {
        let mut v: Vec<u32> = indices[indptr[r] as usize..indptr[r + 1] as usize].to_vec();
        v.sort_unstable();
        // cancel duplicates mod 2
        let mut w: Vec<u32> = Vec::with_capacity(v.len());
        let mut k = 0;
        while k < v.len() {
            let mut cnt = 1;
            while k + cnt < v.len() && v[k + cnt] == v[k] {
                cnt += 1;
            }
            if cnt % 2 == 1 {
                w.push(v[k]);
            }
            k += cnt;
        }
        rows.push(w);
    }
    if !quiet {
        eprintln!("read {} x {} nnz {} [{:.1}s]", nrows, ncols, nnz, t0.elapsed().as_secs_f64());
    }

    // ---- column index ----
    let mut col_rows: Vec<Vec<u32>> = vec![Vec::new(); ncols];
    let mut col_count: Vec<u32> = vec![0; ncols];
    for (r, row) in rows.iter().enumerate() {
        for &c in row {
            col_rows[c as usize].push(r as u32);
            col_count[c as usize] += 1;
        }
    }
    let mut row_active = vec![true; nrows];
    let mut col_active = vec![true; ncols];
    let mut heap: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::with_capacity(ncols);
    for c in 0..ncols {
        if col_count[c] > 0 {
            heap.push(Reverse((col_count[c], c as u32)));
        }
    }
    // rows that are already empty
    let mut inconsistent = false;
    for r in 0..nrows {
        if rows[r].is_empty() {
            row_active[r] = false;
            if rhs[r] & 1 == 1 {
                inconsistent = true;
            }
        }
    }

    // ---- elimination ----
    let mut pivots: Vec<(u32, u32)> = Vec::new(); // (row, col)
    let mut steps = 0usize;
    let mut fill_max = 0usize;
    let mut last_report = Instant::now();
    while let Some(Reverse((cnt, c))) = heap.pop() {
        let c = c as usize;
        if !col_active[c] || col_count[c] != cnt {
            continue; // stale
        }
        // compact the column's row list and pick the shortest active row containing c
        let mut live: Vec<u32> = Vec::with_capacity(col_rows[c].len());
        let mut best: Option<(usize, usize)> = None; // (len, row)
        let mut cand = std::mem::take(&mut col_rows[c]);
        cand.sort_unstable();
        cand.dedup();
        for &r in &cand {
            let r = r as usize;
            if row_active[r] && rows[r].binary_search(&(c as u32)).is_ok() {
                live.push(r as u32);
                let len = rows[r].len();
                if best.map_or(true, |(bl, _)| len < bl) {
                    best = Some((len, r));
                }
            }
        }
        col_rows[c] = live;
        let real_count = col_rows[c].len() as u32;
        if real_count == 0 {
            col_count[c] = 0;
            col_active[c] = false;
            continue;
        }
        if real_count != cnt {
            col_count[c] = real_count;
            heap.push(Reverse((real_count, c as u32)));
            continue;
        }
        let (_, prow) = best.unwrap();
        // eliminate c from all other rows in the column
        let pivot = std::mem::take(&mut rows[prow]);
        let prhs = rhs[prow] & 1;
        row_active[prow] = false;
        col_active[c] = false;
        let targets: Vec<u32> = col_rows[c].iter().copied().filter(|&r| r as usize != prow).collect();
        for &r in &targets {
            let r = r as usize;
            let newrow = xor_rows(&rows[r], &pivot);
            // update column counts: columns of pivot (other than c) toggle membership in r
            for &pc in &pivot {
                if pc as usize == c {
                    continue;
                }
                let pcu = pc as usize;
                if rows[r].binary_search(&pc).is_ok() {
                    col_count[pcu] -= 1;
                } else {
                    col_count[pcu] += 1;
                    col_rows[pcu].push(r as u32);
                }
                heap.push(Reverse((col_count[pcu], pc)));
            }
            rows[r] = newrow;
            rhs[r] ^= prhs;
            if rows[r].is_empty() {
                row_active[r] = false;
                if rhs[r] & 1 == 1 {
                    inconsistent = true;
                }
            } else if rows[r].len() > fill_max {
                fill_max = rows[r].len();
            }
        }
        // the pivot row leaves the active set: its other columns lose one active row
        for &pc in &pivot {
            if pc as usize == c {
                continue;
            }
            col_count[pc as usize] -= 1;
            if col_count[pc as usize] == 0 {
                col_active[pc as usize] = false;
            } else {
                heap.push(Reverse((col_count[pc as usize], pc)));
            }
        }
        col_count[c] = 0;
        col_rows[c].clear();
        rows[prow] = pivot;
        pivots.push((prow as u32, c as u32));
        steps += 1;
        if !quiet && last_report.elapsed().as_secs_f64() > 5.0 {
            let active_nnz: usize = rows.iter().enumerate().filter(|(r, _)| row_active[*r]).map(|(_, v)| v.len()).sum();
            eprintln!(
                "  pivots {} active nnz {} longest row {} inconsistent {} [{:.0}s]",
                steps, active_nnz, fill_max, inconsistent, t0.elapsed().as_secs_f64()
            );
            last_report = Instant::now();
        }
        if inconsistent && stop_early {
            break;
        }
    }
    let rank = pivots.len();
    let consistent = !inconsistent;
    if !quiet {
        eprintln!(
            "done: rank {} consistent {} longest row {} [{:.1}s]",
            rank, consistent, fill_max, t0.elapsed().as_secs_f64()
        );
    }

    // ---- output ----
    let mut out = BufWriter::new(File::create(&args[2]).unwrap());
    writeln!(out, "consistent {}", consistent as u8).unwrap();
    writeln!(out, "rank {}", rank).unwrap();
    writeln!(out, "nfree {}", ncols - rank).unwrap();
    if want_pivots {
        write!(out, "pivots").unwrap();
        for &(_, c) in &pivots {
            write!(out, " {}", c).unwrap();
        }
        writeln!(out).unwrap();
    }
    let mut is_pivot_col = vec![false; ncols];
    for &(_, c) in &pivots {
        is_pivot_col[c as usize] = true;
    }
    // back substitution: x[c_k] = rhs_k ^ sum_{j in row_k, j != c_k} x_j, in reverse pivot order.
    // Row k contains no earlier pivot columns; later pivot columns are already solved; free columns are given.
    let back_sub = |free_one: Option<u32>, use_rhs: bool| -> Vec<u32> {
        let mut x = vec![0u8; ncols];
        if let Some(fc) = free_one {
            x[fc as usize] = 1;
        }
        for &(r, c) in pivots.iter().rev() {
            let mut v = if use_rhs { rhs[r as usize] & 1 } else { 0 };
            for &j in &rows[r as usize] {
                if j != c {
                    v ^= x[j as usize];
                }
            }
            x[c as usize] = v;
        }
        (0..ncols as u32).filter(|&i| x[i as usize] == 1).collect()
    };
    if want_solution && consistent {
        let sol = back_sub(None, true);
        write!(out, "solution").unwrap();
        for c in sol {
            write!(out, " {}", c).unwrap();
        }
        writeln!(out).unwrap();
    }
    if want_null {
        for fc in 0..ncols as u32 {
            if is_pivot_col[fc as usize] {
                continue;
            }
            let v = back_sub(Some(fc), false);
            write!(out, "null {}:", fc).unwrap();
            for c in v {
                write!(out, " {}", c).unwrap();
            }
            writeln!(out).unwrap();
        }
    }
}
