# Knowledge Graph ¡ª mini-linear-algebra-control

## L1: Definitions (Complete)
| Entry | Implementation |
|-------|---------------|
| Matrix data type (dense, row-major) | matrix.h: Matrix struct |
| Vector data type | vector.h: Vector struct |
| Determinant | matrix.c: matrix_determinant() |
| Trace | matrix.c: matrix_trace() |
| Rank | matrix.c: matrix_rank() |
| Matrix inverse | matrix.c: matrix_inverse() |
| Dot product | vector.c: vector_dot() |
| Cross product | vector.c: vector_cross() |
| Vector norms (L1, L2, Linf, Lp) | vector.c: vector_norm_*() |
| Matrix norms (Frobenius, induced) | norms.c: matrix_norm_*() |
| Controllability matrix | control_linalg.c: control_controllability_matrix() |
| Observability matrix | control_linalg.c: control_observability_matrix() |
| Gramians (controllability, observability) | control_linalg.c: *_gramian() |
| Eigenvalues, eigenvectors | eigen.h, eigen.c |
| Characteristic polynomial | eigen.c: charpoly_faddeev_leverrier() |
| Spectral radius / abscissa | eigen.c: eigen_spectral_*() |
| Matrix exponential | matrix_exp.c: matrix_exp() |

## L2: Core Concepts (Complete)
| Entry | Implementation |
|-------|---------------|
| Linear transformations | matrix.c: matrix_multiply() |
| Span and basis | vector.c: vector_span_dimension(), vector_gram_schmidt() |
| Linear independence | vector.c: vector_is_linearly_independent() |
| Orthogonality | matrix.c: matrix_is_orthogonal(), vector_gram_schmidt() |
| Controllability concept | control_linalg.c: control_is_controllable() |
| Observability concept | control_linalg.c: control_is_observable() |
| Duality (controllability/observability) | control_linalg.c: observer_gain_luenberger() |
| State transition | matrix_exp.c: matrix_exp_t() |
| Pseudoinverse (Moore-Penrose) | matrix.c: matrix_pseudoinverse() |
| Kronecker product | matrix.c: matrix_kronecker() |
| Subspace iteration | eigen.c: eigen_subspace_iteration() |

## L3: Engineering Quantities (Complete)
| Entry | Implementation |
|-------|---------------|
| Condition number (L1, L2, Linf, Frobenius) | norms.c: matrix_condition_*() |
| Spectral radius rho(A) | eigen.c: eigen_spectral_radius() |
| Spectral abscissa alpha(A) | eigen.c: eigen_spectral_abscissa() |
| Hankel singular values | control_linalg.c: control_hankel_singular_values() |
| LU growth factor | norms.c: lu_growth_factor() |
| Orthogonality error | norms.c: orthogonality_error() |
| Linear system residual | norms.c: linear_system_residual() |
| Lyapunov residual | norms.c: lyapunov_residual() |

## L4: Conservation Laws / Theorems (Complete)
| Entry | Implementation |
|-------|---------------|
| Cayley-Hamilton theorem | matrix.c: matrix_polynomial() |
| Spectral theorem (symmetric eig) | eigen.c: eigen_symmetric_decompose() |
| Schur decomposition theorem | decompositions.c: schur_decompose() |
| SVD existence theorem | decompositions.c: svd_jacobi() |
| Cholesky for SPD matrices | decompositions.c: cholesky_decompose() |
| LU existence for invertible matrices | decompositions.c: lu_decompose() |
| Gram-Schmidt constructive basis proof | vector.c: vector_gram_schmidt() |
| PBH controllability/observability test | control_linalg.c: control_pbh_*_test() |
| Lyapunov stability theorem | control_linalg.c: lyapunov_solve_continuous() |
| Matrix exponential properties | matrix_exp.c |

## L5: Algorithms / Engineering Methods (Complete)
| Entry | Implementation |
|-------|---------------|
| Gaussian elimination + partial pivoting | decompositions.c: lu_decompose() |
| QR via Modified Gram-Schmidt | decompositions.c: qr_decompose_mgs() |
| QR via Householder reflections | decompositions.c: qr_decompose_householder() |
| Givens rotations | decompositions.c: qr_givens_rotate() |
| Cholesky decomposition | decompositions.c: cholesky_decompose() |
| LDL^T decomposition | decompositions.c: ldl_decompose() |
| SVD via one-sided Jacobi | decompositions.c: svd_jacobi() |
| Hessenberg reduction | decompositions.c: hessenberg_reduce() |
| QR algorithm for eigenvalues | eigen.c: eigen_qr_symmetric/general() |
| Power iteration | eigen.c: eigen_power_iteration() |
| Inverse iteration | eigen.c: eigen_inverse_iteration() |
| Rayleigh quotient iteration | eigen.c: eigen_rayleigh_quotient_iteration() |
| Faddeev-Leverrier (char poly) | eigen.c: charpoly_faddeev_leverrier() |
| Scaling+squaring (matrix exp) | matrix_exp.c: matrix_exp() |
| Bartels-Stewart (Lyapunov) | control_linalg.c: lyapunov_solve_continuous() |
| Hamiltonian/Schur (Riccati) | control_linalg.c: riccati_solve_care() |
| Normal equations (least squares) | solvers.c: least_squares_normal_equations() |
| QR least squares | solvers.c: least_squares_qr() |
| SVD least squares | solvers.c: least_squares_svd() |
| Conjugate gradient | solvers.c: solve_conjugate_gradient() |
| Jacobi preconditioned CG | solvers.c: solve_cg_jacobi_preconditioned() |
| GMRES | solvers.c: solve_gmres() |

## L6: Canonical Problems (Complete)
| Entry | Implementation |
|-------|---------------|
| Controllability analysis | control_linalg.c + examples/controllability_check.c |
| Observability analysis | control_linalg.c + examples/controllability_check.c |
| Stabilizability / Detectability | control_linalg.c: control_is_stabilizable/detectable() |
| Kalman decomposition | control_linalg.c: control_kalman_decomposition() |
| Pole placement (Ackermann) | control_linalg.c + examples/pole_placement_demo.c |
| Pole placement (Bass-Gura) | control_linalg.c |
| ZOH discretization | matrix_exp.c: discretize_zoh() |
| Stability analysis (Hurwitz/Schur) | eigen.c: eigen_is_hurwitz/schur() |
| Minimum-norm solution | solvers.c: solve_minimum_norm() |

## L7: Applications (Partial+)
| Entry | Implementation |
|-------|---------------|
| LQR optimal control | control_linalg.c: lqr_gain() + examples/lqr_design_demo.c |
| Luenberger observer | control_linalg.c: observer_gain_luenberger() |
| Kalman filter gain | control_linalg.c: kalman_gain() |
| Ridge regression | solvers.c: least_squares_ridge() |
| Non-negative least squares | solvers.c: least_squares_nonnegative() |
| Iterative solvers (CG, GMRES) | solvers.c |

## L8: Advanced Methods (Partial+)
| Entry | Implementation |
|-------|---------------|
| Balanced realization | control_linalg.c: balanced_realization() |
| Balanced truncation (model reduction) | control_linalg.c: balanced_truncation() |
| Total least squares | solvers.c: least_squares_total() |
| Matrix logarithm | matrix_exp.c: matrix_log() |
| Matrix cosine/sine | matrix_exp.c: matrix_cos/sin() |
| Frechet derivative (matrix exp) | matrix_exp.c: matrix_exp_frechet() |
| Matrix sign function | eigen.c: matrix_sign_function() |
| Matrix square root | eigen.c: matrix_sqrt() |
| L2,1 norm (structured sparsity) | norms.c: matrix_norm_l21() |

## L9: Research Frontiers (Partial)
| Entry | Implementation |
|-------|---------------|
| Fractional matrix powers | matrix_exp.c: matrix_pow_real() |
| Nuclear norm (low-rank recovery) | norms.c: matrix_norm_nuclear() |
| Truncated SVD | decompositions.c: svd_truncate() |
| Generalized eigenvalue problem | eigen.c: eigen_generalized() |
