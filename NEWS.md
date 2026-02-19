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
