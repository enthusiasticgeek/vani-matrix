# vani-matrix — TODO

All matrices are flat row-major `Vec<f64>`: `M[i * cols + j]` for row i, col j.
Explicit dimension args everywhere — no metadata hidden in the Vec.

---

## Compiler dependencies

| Ticket | Needed for |
|--------|-----------|
| None currently | All operations use existing arithmetic builtins |

---

## v0.1.0 — Implementation plan (complete)

All functions, tests, examples, and safety annotations below are implemented
and passing (`vanic run` exit 0 on every file in `tests/` and `examples/`).

### Section 1: Constructors and accessors

- [x] `mat_zeros(rows, cols)` — push 0.0 `rows*cols` times; WCET = rows×cols × 8 cycles
- [x] `mat_ones(rows, cols)` — push 1.0 `rows*cols` times
- [x] `mat_identity(n)` — push 1.0 at diagonal positions, else 0.0; WCET = n² × 12 cycles
- [x] `mat_from_diag(diag, n)` — n×n matrix with `diag[i]` on diagonal, 0 elsewhere
- [x] `mat_copy(A, rows, cols)` — deep copy via push loop
- [x] `mat_get(A, cols, i, j)` — return `A[i * cols + j]`; O(1)
- [x] `mat_set(A mut ref, cols, i, j, v)` — `set(mut ref A, i*cols+j, v)`; O(1); returns i64 0
- [x] `mat_row(A, cols, i, n_cols)` — extract row i into new Vec; WCET = n_cols × 12 cycles
- [x] `mat_col(A, rows, cols, j)` — extract column j into new Vec; WCET = rows × 12 cycles

### Section 2: Arithmetic (element-wise)

- [x] `mat_add(A, B, rows, cols)` — C[k] = A[k] + B[k]; WCET = rows×cols × 10 cycles
- [x] `mat_sub(A, B, rows, cols)` — C[k] = A[k] - B[k]
- [x] `mat_scale(A, rows, cols, s)` — C[k] = A[k] * s
- [x] `mat_hadamard(A, B, rows, cols)` — C[k] = A[k] * B[k] (element-wise product)
- [x] `mat_transpose(A, rows, cols)` — result is cols×rows; T[j*rows+i] = A[i*cols+j]; WCET = rows×cols × 12 cycles

### Section 3: Vector ops and norms

- [x] `dot_n(x, y, n)` — Σ x[i]*y[i]; WCET = n × 10 cycles (named `dot_n`, not `vec_dot`, to avoid conflicting with the compiler builtin)
- [x] `vec_norm2(x, n)` — sqrt(vec_dot(x, x, n)); WCET = n × 10 + sqrt
- [x] `vec_outer(x, y, m, n)` — m×n matrix R[i*n+j] = x[i]*y[j]; WCET = m×n × 10 cycles
- [x] `mat_norm_fro(A, rows, cols)` — sqrt(Σ A[k]²); WCET = rows×cols × 12 cycles
- [x] `mat_norm_max(A, rows, cols)` — max |A[k]|; WCET = rows×cols × 8 cycles

### Section 4: Multiplication

- [x] `mat_mul(A, B, n)` — square n×n × n×n; triple loop O(n³); WCET = n³ × 10 cycles
- [x] `mat_mul_rect(A, B, m, k, n)` — m×k × k×n → m×n; C[i*n+j] = Σₗ A[i*k+l]*B[l*n+j]
- [x] `mat_vec_mul(A, x, m, n)` — m×n matrix × n-Vec → m-Vec; y[i] = Σⱼ A[i*n+j]*x[j]
- [x] `mat_pow_n(A, n, p)` — integer matrix power via repeated squaring; p=0 → identity

### Section 5: 2×2 and 3×3 (closed-form, no loops)

- [x] `mat_det_2x2(A)` — A[0]*A[3] - A[1]*A[2]; WCET = 27 cycles (measured)
- [x] `mat_inv_2x2(A)` — (1/det) × adjugate; returns zeros if |det| < ε
- [x] `mat_det_3x3(A)` — cofactor expansion along row 0; WCET = 102 cycles (measured)
- [x] `mat_inv_3x3(A)` — (1/det) × adjugate matrix; returns zeros if |det| < ε

### Section 6: General n×n (Gauss-Jordan with partial pivoting)

- [x] `mat_inv_n(A, n)` — augmented [A | I] → row reduce → [I | A⁻¹]; O(n³); returns zeros if singular
- [x] `mat_det_n(A, n)` — via lu_factor: product of U diagonal × pivot sign
- [x] `mat_solve(A, b, n)` — solve Ax=b via augmented Gauss-Jordan [A | b] → [I | x]

### Section 7: LU decomposition (Doolittle, partial pivoting)

Packed format: `lu_factor` returns `Vec<f64>` of length `n*n + n`.
- Elements `0..n²−1`: combined LU (L below diagonal with implicit 1s; U on/above diagonal)
- Elements `n²..n²+n−1`: pivot row indices as f64

- [x] `lu_factor(A, n)` — Doolittle with partial pivoting; O(n³); packed output
- [x] `lu_det(lu_packed, n)` — product of U diagonal × (−1)^{swap_count from pivots}
- [x] `lu_solve(lu_packed, b, n)` — apply pivots, forward sub (L), back sub (U); O(n²)

### Section 8: Cholesky decomposition (symmetric positive-definite)

- [x] `cholesky_factor(A, n)` — lower triangular L such that A = LLᵀ;
      L[i*n+j] = (A[i*n+j] − Σₖ<ⱼ L[i*n+k]×L[j*n+k]) / L[j*n+j] for j<i;
      L[j*n+j] = sqrt(A[j*n+j] − Σₖ<ⱼ L[j*n+k]²);
      returns zeros if diagonal goes negative (not SPD); O(n³/3)
- [x] `cholesky_solve(L, b, n)` — solve Ax=b given L: forward sub Ly=b, back sub Lᵀx=y; O(n²)

---

## Safety / WCET annotations

- [x] `#[bounded_stack(bytes=N)]` on all v0.1.0 functions after `vanic stack-depth` run
- [x] `#[bounded_stack(bytes=N)]` on all v0.2.0 functions (budgets set to `vanic check`'s
      exact reported worst-case, not hand estimates -- the Householder QR and SVD
      bidiagonalization functions needed the largest budgets in the library so far,
      up to 1072 bytes for `mat_svd_bidiag`)
- [x] `#[wcet(cycles=N)]` on closed-form leaf functions (mat_get, mat_set, mat_det_2x2,
      mat_det_3x3). `mat_inv_2x2`/`mat_inv_3x3` do NOT get `#[wcet]`: they call
      `mat_zeros`, whose loop bound is a runtime parameter, so the compiler's static
      estimator can't prove a finite bound even though the actual call sites only
      ever pass 2 or 3.

---

## Tests

- [x] `tests/test_construct.vani` — zeros/ones/identity/from_diag/copy/get/set/row/col
- [x] `tests/test_arithmetic.vani` — add/sub/scale/hadamard/transpose; verify element-wise
- [x] `tests/test_multiply.vani` — mat_mul identity, mat_mul_rect, mat_vec_mul, mat_pow_n (p=0,1,2,3)
- [x] `tests/test_inverse.vani` — det_2x2, inv_2x2 × original = I, det_3x3, inv_3x3, mat_det_n,
      mat_inv_n, mat_solve (including singular-matrix cases)
- [x] `tests/test_decomp.vani` — lu_factor/lu_det/lu_solve residual checks; cholesky_factor:
      L×Lᵀ = A, cholesky_solve, non-SPD detection

## Examples

- [x] `examples/linear_system.vani` — 3×3 Ax=b via mat_solve and mat_inv_n
- [x] `examples/least_squares.vani` — overdetermined 5×2 system via normal equations (XᵀX)⁻¹Xᵀy
- [x] `examples/cholesky_demo.vani` — 3×3 SPD matrix → Cholesky factor → solve + verify

---

## v0.2.0 -- Implemented ✓

No dependency on vani-probability: checked both repos' TODOs and found no
technical reason for v0.2.0 to wait. If anything the dependency runs the
other way — vani-probability's v0.4.0 (`mlr_fit`, `cov_matrix`,
`correlation_matrix`, `pca_power_iter`) needs `mat_inv_n`, `mat_mul_rect`,
`mat_transpose`, and Cholesky from *this* repo, all of which now ship in
v0.1.0, and its planned `pca_power_iter` implements its own power iteration
rather than calling into vani-matrix's eigenvalue functions. So vani-matrix
v0.2.0 can start independently whenever it's prioritized.

- [x] `mat_eig_power(A, n, max_iter, tol)` — dominant eigenvalue/vector via power iteration
- [x] `mat_eig_deflate(A, n, lambda, v)` — deflate to get next eigenpair
- [x] `mat_qr_householder(A, m, n)` — QR decomposition (Householder reflectors)
- [x] `mat_svd_bidiag(A, m, n)` — bidiagonalisation step for SVD (not a full SVD solver --
      diagonalizing the resulting bidiagonal form needs an iterative implicit-shift QR
      sweep, intentionally out of scope here)
- [x] `mat_cond(A, n)` — condition number (max eigenvalue / min eigenvalue via power + inverse power)
- [x] `tests/test_eig.vani` — positive/negative/symmetric eigenvalues, deflation, cond
- [x] `tests/test_qr_svd.vani` — QR against the canonical Wikipedia 3×3 example (validated
      via orthogonality + reconstruction, not exact sign match), non-square QR, 4×3 SVD
      bidiagonalization (orthogonality + reconstruction)
- [x] `examples/eigen_qr_svd_demo.vani`

**Bug found and fixed during validation**: `mat_eig_power`'s original uniform starting
vector (`[1/√n, ..., 1/√n]`) can land exactly on a *non-dominant* eigenvector of a
symmetric matrix by construction — this reliably happened after `mat_eig_deflate` on a
symmetric input (deflating `[[2,1],[1,2]]`'s dominant eigenpair leaves a matrix whose
zero eigenvalue's eigenvector is exactly the old uniform start), silently converging to
the wrong eigenvalue, and the same degeneracy independently corrupted `mat_cond` (which
calls `mat_eig_power` on `A⁻¹` internally). Fixed by starting from a normalized
`[1, 2, ..., n]` vector instead, which breaks the symmetry that caused it. Caught by
testing the eigen/QR/SVD example end-to-end rather than only running the unit tests
that happened not to exercise this exact case — now locked in as a dedicated
regression test in `test_eig.vani`.
