#' Convex Least Squares Programming (CLSP) estimator
#'
#' @param problem character scalar, optional
#'   Structural template for matrix construction. One of:
#'   - 'ap'   or 'tm' : allocation or tabular matrix problem.
#'   - 'cmls' or 'rp' : constrained modular least squares or RP-type.
#'   - ''     or other: general CLSP problems (user-defined \code{C} and/or
#'   \code{M}).
#'
#' @param C,S,M numeric matrix or NULL
#'   Blocks of the constraint matrix
#'   \eqn{A = \begin{bmatrix} \mathbf{C} & \mathbf{S} \\ \mathbf{M} & \mathbf{Q}
#'            \end{bmatrix}}.}. If \code{C} and/or \code{M} are provided, the
#'   matrix \code{A} is constructed accordingly. If both are NULL and \code{A}
#'   is not yet defined, an error is raised.
#'
#' @param b numeric vector or NULL
#'   Right-hand side vector. Must have as many rows as \code{A}. Required.
#'
#' @param m,p integer scalar or NULL
#'   Dimensions of \eqn{\mathbf{X} \in \mathbb{R}^{m \times p}}, relevant for
#'   allocation problems ('ap').
#'
#' @param i,j integer scalar, default = 1
#'   Grouping sizes for row and column sum constraints in AP problems.
#'
#' @param zero_diagonal logical scalar, default = FALSE
#'   If TRUE, enforces structural zero diagonals via identity truncation.
#'
#' @param r integer scalar, default = 1
#'   Number of refinement iterations for the pseudoinverse-based estimator.
#'   When \eqn{r > 1}, the slack block \code{Q} is updated iteratively to
#'   improve feasibility in underdetermined or ill-posed systems.
#'
#' @param Z numeric matrix or NULL
#'   A symmetric idempotent matrix (projector) defining the subspace for
#'   Bott–Duffin pseudoinversion. If NULL, the identity matrix is used, reducing
#'   to the Moore–Penrose case.
#'
#' @param rcond numeric scalar or logical scalar, default = FALSE
#'   Regularization parameter for the Moore–Penrose and Bott–Duffin inverses,
#'   providing numerically stable inversion and ensuring convergence of singular
#'   values. If TRUE, an automatic tolerance equal to \code{tolerance} is
#'   applied. If set to a numeric value, it specifies the relative cutoff below
#'   which small singular values are treated as zero.
#'
#' @param tolerance numeric scalar or NULL, default = NULL
#'   Convergence tolerance for NRMSE change between iterations.
#'
#' @param iteration_limit integer scalar or NULL, default = NULL
#'   Maximum number of iterations allowed in the refinement loop.
#'
#' @param final logical scalar, default = TRUE
#'   If TRUE, a convex programming problem is solved to refine \code{zhat}.
#'   The resulting solution \code{z} minimizes a weighted \eqn{\ell_1/\ell_2}
#'   norm around \code{zhat} subject to \eqn{\mathbf{A}\mathbf{z} = \mathbf{b}}.
#'
#' @param alpha numeric scalar, numeric vector, or NULL, default = NULL
#'   Regularization parameter controlling the convex correction type:
#'   - \eqn{\alpha = 0}: Lasso (\eqn{\ell_1} norm)
#'   - \eqn{\alpha = 1}: Ridge (\ell_1} norm)
#'   - \eqn{0 < \alpha < 1}: Elastic Net
#'   If a numeric scalar is provided, that value is used after clipping to
#'   [0, 1]. If a numeric vector is provided, each candidate is evaluated via a
#'   full solve, and the \eqn{\alpha} with the smallest NRMSE is selected. If
#'   NULL, \eqn{\alpha is chosen, based on an error rule:
#'   \eqn{\alpha = \min(1, \text{NRMSE}_{\alpha=0} / (\text{NRMSE}_{\alpha=0} +
#'        \text{NRMSE}_{\alpha=1} + \text{tolerance}))}
#'
#' @param ... optional
#'   Additional arguments passed to the CVXR solver backend.
#'
#' @return
#'   An object of class \code{"clsp"} representing the fitted
#'   Convex Least Squares Programming (CLSP) model.
#'   The object is a list containing all initialized fields and
#'   solver results.
#'   Class-specific methods such as \code{summary.clsp()},
#'   \code{corr.clsp()}, and \code{ttest.clsp()} can be used to
#'   extract, analyze, and summarize the results.
#'
#' @section Initialized Fields:
#'   The returned \code{"clsp"} object is a named list containing
#'   the following elements:
#'
#' \describe{
#'   \item{A}{numeric matrix. Canonical design matrix with block structure
#'     \eqn{A = \begin{bmatrix} C & S \\ M & Q \end{bmatrix}}.}
#'
#'   \item{C_idx}{integer vector of length 2. Indices defining the row and
#'     column ranges of the \code{C} block inside \code{A}. Used for matrix
#'     slicing and partitioning.}
#'
#'   \item{b}{numeric vector. Right-hand side vector for the system
#'     \eqn{\mathbf{A}\mathbf{z} = \mathbf{b}}.}
#'
#'   \item{Z}{numeric matrix or NULL. Projection matrix for Bott–Duffin
#'     inversion; must be symmetric and idempotent. Defaults to identity
#'     (Moore–Penrose).}
#'
#'   \item{tolerance}{numeric scalar. Convergence tolerance for NRMSE change
#'     between iterations.}
#'
#'   \item{iteration_limit}{integer scalar. Maximum number of iterations
#'     allowed in the refinement loop.}
#'
#'   \item{zhat}{numeric vector. Unregularized pseudoinverse estimate of
#'     \code{z} from Step 1.}
#'
#'   \item{final}{logical scalar. Whether to run the convex refinement step.
#'     If TRUE, the estimate is regularized via convex programming.}
#'
#'   \item{alpha}{numeric scalar or NULL. Regularization parameter:
#'     \eqn{\alpha=0} → Lasso, \eqn{\alpha=1} → Ridge,
#'     \(0<\alpha<1\) → Elastic Net.}
#'
#'   \item{z}{numeric vector. Final estimate after regularization. If
#'     \code{final = FALSE}, equals \code{zhat}.}
#'
#'   \item{x}{numeric matrix. Variable component extracted from \code{z},
#'     reshaped to \eqn{m \times p}.}
#'
#'   \item{y}{numeric vector. Slack component of \code{z} representing
#'     inequality residuals.}
#'
#'   \item{r}{integer scalar. Number of refinement iterations performed during
#'     Step 1. Iteration stops when NRMSE stabilizes or the limit is reached.}
#'
#'   \item{kappaC,kappaB,kappaA}{numeric scalars. Condition numbers of the
#'     constraint block C, projected estimator B\eqn{^{(r)}}, and canonical
#'     matrix A\eqn{^{(r)}} respectively.}
#'
#'   \item{rmsa}{numeric scalar. Total RMSA (Root Mean Square Adjustment)
#'     over all rows.}
#'
#'   \item{rmsa_i,rmsa_dkappaC,rmsa_dkappaB,rmsa_dkappaA,
#'         rmsa_dnrmse,rmsa_dzhat,rmsa_dz,rmsa_dx}{numeric vectors.
#'     Diagnostics showing changes in RMSA, condition numbers, NRMSE, and
#'     estimates when rows of \code{[C | S]} are removed and the model
#'     re-estimated.}
#'
#'   \item{r2_partial}{numeric scalar. Partial R² computed over the M block.}
#'
#'   \item{nrmse,nrmse_partial}{numeric scalars. Normalized RMSE for the full
#'     system and for the M block respectively.}
#'
#'   \item{nrmse_ttest}{numeric vector. NRMSE samples generated via Monte
#'     Carlo simulation for empirical t-testing.}
#'
#'   \item{z_lower,z_upper,x_lower,x_upper,y_lower,y_upper}{numeric vectors.
#'     Lower and upper bounds of condition-weighted confidence bands derived
#'     from \eqn{\kappa(A)} and residual norms.}
#'
#'   \item{seed}{integer scalar. Random seed used for reproducible Monte
#'     Carlo diagnostics.}
#'
#'   \item{distribution}{function. Generates random samples for simulation,
#'     typically \code{rnorm(n)}.}
#' }
#'
#' @export
clsp <- function(problem = "",
                 C = NULL,
                 S = NULL,
                 M = NULL,
                 b = NULL,
                 m = NULL,
                 p = NULL,
                 i = 1L,
                 j = 1L,
                 zero_diagonal = FALSE,
                 r = 1L,
                 Z = NULL,
                 rcond = FALSE,
                 tolerance = NULL,
                 iteration_limit = NULL,
                 final = TRUE,
                 alpha = NULL,
                 ...) {
  
  # ---- default normalization --------------------------------------------
  if (isFALSE(rcond)) rcond <- sqrt(.Machine$double.eps)
  if (is.null(tolerance)) tolerance <- sqrt(.Machine$double.eps)
  if (is.null(iteration_limit)) iteration_limit <- 50L
  
  # ---- RNG ---------------------------------------------------------------
  seed <- 123456789L
  set.seed(seed)
  distribution <- function(n) stats::rnorm(n, mean = 0, sd = 1)
  
  # ---- initialize object -------------------------------------------------
  object <- list(
    problem         = problem,
    C               = C,
    S               = S,
    M               = M,
    b               = b,
    m               = m,
    p               = p,
    i               = i,
    j               = j,
    zero_diagonal   = zero_diagonal,
    r               = r,
    Z               = Z,
    rcond           = rcond,
    tolerance       = tolerance,
    iteration_limit = iteration_limit,
    final           = final,
    alpha           = alpha,
    # placeholders and derived fields
    A              = NULL,
    C_idx          = c(NA_integer_, NA_integer_),
    zhat           = NULL,
    z              = NULL,
    x              = NULL,
    y              = NULL,
    kappaC         = NULL,
    kappaB         = NULL,
    kappaA         = NULL,
    rmsa           = NULL,
    rmsa_i         = numeric(),
    rmsa_dkappaC   = numeric(),
    rmsa_dkappaB   = numeric(),
    rmsa_dkappaA   = numeric(),
    rmsa_dnrmse    = numeric(),
    rmsa_dzhat     = numeric(),
    rmsa_dz        = numeric(),
    rmsa_dx        = numeric(),
    r2_partial     = NULL,
    nrmse          = NULL,
    nrmse_partial  = NULL,
    nrmse_ttest    = numeric(),
    z_lower        = NULL,
    z_upper        = NULL,
    x_lower        = NULL,
    x_upper        = NULL,
    y_lower        = NULL,
    y_upper        = NULL,
    seed           = seed,
    distribution   = distribution,
    solver         = "auto",
    status         = NULL,
    value          = NULL
  )
  
  # ---- delegate to solver -----------------------------------------------
  if (!is.null(M) && !is.null(b)) {
    object <- .clsp_solver(object)
  }
  
  class(object) <- "clsp"
  object
}