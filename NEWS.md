# rclsp 1.0.0

## Changes
* Preserved rcond in stored results.
* Added custom `.cond()` function with a relative singular-value cutoff.
* Added `cond_tolerance` argument to `clsp()`.
* Used Euclidean/Frobenius norms for corr() changes in zhat, z, and x.
* Fixed a bootstrap bug in ttest().

## Bug fixes
* Improved `alpha` handling in `clsp()` to correctly support numeric scalars,
  numeric vectors, and `NULL`.
* Fixed candidate `alpha` selection when non-finite values are supplied.
* Restricted automatic `alpha` selection to cases where the final convex
  correction step is enabled.

# rclsp 0.5.0

## Bug fixes
* Updated CVXR integration for CVXR 1.8.x compatibility.

# rclsp 0.4.0

## Bug fixes
* Fixed potential projector Z dimension mismatch under repeated solve() calls.

# rclsp 0.3.0

## Bug fixes
* Achieved full numerical parity between R and Python implementations.
* Replaced MASS::ginv() with explicit SVD-based pseudoinverse logic.
* Aligned Bott–Duffin tolerance handling with NumPy semantics.
* Corrected row-wise reconstruction of X from the solution vector z.
* Fixed zero-diagonal constraint handling for allocation (TM/AP) problems.

# rclsp 0.2.0

## Bug fixes
* Removed duplicated warnings.
