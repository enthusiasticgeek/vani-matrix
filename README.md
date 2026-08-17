# vani-matrix

Dense linear algebra library for the [vāṇī compiler](https://github.com/enthusiasticgeek/vani-compiler).

Matrices are stored as flat row-major `Vec<f64>`: element (i, j) of an m×n matrix M is `M[i*n + j]`.
All functions take explicit dimension arguments so no metadata is hidden inside the Vec.

**API reference / tutorial:** <https://enthusiasticgeek.github.io/vani-matrix/>

## Add to your project

```toml
# vani.toml
[deps]
matrix = { registry = "kosh", version = "^0.2" }
```

```sh
vanic add matrix
vanic build
```

## What's included (v0.2.0 — complete; see TODO.md)

| Module | Functions |
|---|---|
| Construct | `mat_zeros`, `mat_ones`, `mat_identity`, `mat_from_diag`, `mat_copy` |
| Access | `mat_get`, `mat_set`, `mat_row`, `mat_col` |
| Arithmetic | `mat_add`, `mat_sub`, `mat_scale`, `mat_hadamard`, `mat_transpose` |
| Vector ops | `dot_n`, `vec_norm2`, `vec_outer`, `mat_norm_fro`, `mat_norm_max` |
| Multiply | `mat_mul`, `mat_mul_rect`, `mat_vec_mul`, `mat_pow_n` |
| 2×2 / 3×3 | `mat_det_2x2`, `mat_inv_2x2`, `mat_det_3x3`, `mat_inv_3x3` |
| General n×n | `mat_det_n`, `mat_inv_n`, `mat_solve` |
| LU decomp | `lu_factor`, `lu_det`, `lu_solve` |
| Cholesky | `cholesky_factor`, `cholesky_solve` |
| Eigenvalues | `mat_eig_power`, `mat_eig_deflate`, `mat_cond` |
| QR / SVD | `mat_qr_householder`, `mat_svd_bidiag` |

## Matrix encoding

```
// 3×3 matrix A:
//   | a b c |
//   | d e f |   →  Vec<f64> of length 9
//   | g h i |      A[0*3+0]=a  A[0*3+1]=b  ...  A[2*3+2]=i
```

All functions index as `A[i * cols + j]` for row i, column j.

## LU packing convention

`lu_factor(A, n)` returns a `Vec<f64>` of length `n*n + n`:
- Elements `0 .. n*n − 1`: packed LU (L below diagonal, implicit 1s; U on/above diagonal)
- Elements `n*n .. n*n+n − 1`: pivot row indices stored as `f64`

`lu_solve` and `lu_det` expect this packed format directly.

## Eigenvalue / QR / SVD packing conventions

`mat_eig_power(A, n, max_iter, tol)` returns a `Vec<f64>` of length `n+1`: element 0 is
the eigenvalue, elements `1..n+1` are the unit eigenvector.

`mat_qr_householder(A, m, n)` returns `Q` (`m*m` elements) followed by `R` (`m*n`
elements), concatenated.

`mat_svd_bidiag(A, m, n)` returns `U` (`m*m`) followed by `B` (`m*n`, bidiagonal)
followed by `V` (`n*n`), concatenated. This is the bidiagonalization *step* toward a
full SVD (`A = U B Vᵀ`), not a full SVD solver -- diagonalizing `B` further needs an
iterative implicit-shift QR sweep, which is out of scope for v0.2.0.

## What this library does NOT provide

These are already vāṇī compiler builtins — call them directly, no import needed:

`abs` `sqrt` `pow` `exp` `log` `sin` `cos` `f64_tgamma` `f64_lgamma`

## License

MIT
