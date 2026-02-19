#' Convex Least Squares Programming (CLSP) estimator.
#'
#' @description
#' The Convex Least Squares Programming (CLSP) estimator solves
#' underdetermined, ill-posed, or structurally constrained least-squares
#' problems using a modular two-step approach. The first step computes a
#' pseudoinverse-based estimate, and the second step applies a convex
#' correction (Lasso, Ridge, or Elastic Net) to ensure numerical stability,
#' constraint enforcement, and interpretability.
#'
#' @details
#' This estimator unifies pseudoinverse-based least squares with convex
#' programming correction. The pseudoinverse step computes an initial solution
#' \eqn{\mathbf{z}^{(r)}} iteratively via the Moore–Penrose or Bott–Duffin
#' inverse. The convex step then refines \eqn{\boldsymbol{z}} by minimizing a
#' mixed \eqn{\ell_1/\ell_2} norm under equality constraints
#' \eqn{\mathbf{A}\mathbf{z} = \mathbf{b}}. The method supports allocation
#' problems (AP), constrained modular least squares (CMLS), and general CLSP
#' formulations.
#'
#' @param problem character scalar, optional  
#'   Structural template for matrix construction. One of:
#'   - `'ap'`   or `'tm'`: allocation or tabular matrix problem.
#'   - `'cmls'` or `'rp'`: constrained modular least squares or RP-type.
#'   - `''`     or other:  general CLSP problems (user-defined
#'     \eqn{\boldsymbol{C}} and/or \eqn{\boldsymbol{M}}).
#'
#' @param C,S,M numeric matrix or \code{NULL}  
#'   Blocks of the constraint matrix
#'   \eqn{\mathbf{A} = \begin{bmatrix} \mathbf{C} & \mathbf{S} \\
#'                                     \mathbf{M} & \mathbf{Q}
#'                     \end{bmatrix}}.
#'   If \eqn{\boldsymbol{C}} and/or \eqn{\boldsymbol{M}} are provided, the
#'   matrix \eqn{\boldsymbol{A}} is constructed accordingly. If both are
#'   \code{NULL} and \eqn{\boldsymbol{A}} is not yet defined, an error is
#'   raised.
#'
#' @param b numeric vector or \code{NULL}  
#'   Right-hand-side vector. Must have as many rows as \eqn{\boldsymbol{A}}.
#'   Required.
#'
#' @param m,p integer scalar or \code{NULL}  
#'   Dimensions of \eqn{\mathbf{X} \in \mathbb{R}^{m \times p}}, relevant for
#'   allocation problems ('ap').
#'
#' @param i,j integer scalar, default = \code{1}  
#'   Grouping sizes for row and column-sum constraints in AP problems.
#'
#' @param zero_diagonal logical scalar, default = \code{FALSE}  
#'   If \code{TRUE}, enforces structural zero diagonals via identity truncation.
#'
#' @param r integer scalar, default = \code{1}  
#'   Number of refinement iterations for the pseudoinverse-based estimator.
#'   When \eqn{r > 1}, the slack block \eqn{\boldsymbol{Q}} is updated
#'   iteratively to improve feasibility in underdetermined or ill-posed systems.
#'
#' @param Z numeric matrix or \code{NULL}  
#'   A symmetric idempotent matrix (projector) defining the subspace for
#'   Bott–Duffin pseudoinversion. If \code{NULL}, the identity matrix is used,
#'   reducing to the Moore–Penrose case.
#'
#' @param rcond numeric scalar or logical scalar, default = \code{FALSE}  
#'   Regularization parameter for the Moore–Penrose and Bott–Duffin inverses,
#'   providing numerically stable inversion and ensuring convergence of singular
#'   values. If \code{TRUE}, an automatic tolerance equal to \code{tolerance} is
#'   applied. If set to a numeric value, it specifies the relative cutoff below
#'   which small singular values are treated as zero.
#'
#' @param tolerance numeric scalar or \code{NULL}, default = \code{NULL}  
#'   Convergence tolerance for NRMSE change between iterations.
#'
#' @param iteration_limit integer scalar or \code{NULL}, default = \code{NULL}  
#'   Maximum number of iterations allowed in the refinement loop.
#'
#' @param final logical scalar, default = \code{TRUE}  
#'   If \code{TRUE}, a convex programming problem is solved to refine
#'   \code{zhat}. The resulting solution \eqn{\boldsymbol{z}} minimizes a
#'   weighted \eqn{\ell_1/\ell_2} norm around \eqn{\widehat{\mathbf{z}}}
#'   subject to \eqn{\mathbf{A}\mathbf{z} = \mathbf{b}}.
#'
#' @param alpha numeric scalar, numeric vector, or \code{NULL},
#'   default = \code{NULL}  
#'   Regularization parameter:  
#'   - \eqn{\alpha = 0}: Lasso (\eqn{\ell_1} norm)  
#'   - \eqn{\alpha = 1}: Ridge (\eqn{\ell_2} norm)  
#'   - \eqn{0 < \alpha < 1}: Elastic Net.  
#'   If a numeric scalar is provided, that value is used after clipping to
#'   \eqn{[0,1]}. If a numeric vector is provided, each candidate is evaluated
#'   via a full solve, and the \eqn{\alpha} with the smallest NRMSE is selected.
#'   If \code{NULL}, \eqn{\alpha} is chosen automatically according to
#'   \deqn{\alpha =
#'          \min\left(1,\,
#'          \frac{\mathrm{NRMSE}_{\alpha=0}}
#'               {\mathrm{NRMSE}_{\alpha=0} +
#'                \mathrm{NRMSE}_{\alpha=1} +
#'                \mathrm{tolerance}}\right)}.
#'
#' @param ... Optional.  
#'   Additional arguments passed to the \pkg{CVXR} solver backend.
#'
#' @return
#' An object of class \code{"clsp"} representing the fitted
#' Convex Least Squares Programming (CLSP) model.
#' The object is a named list containing all initialized fields and
#' solver results.  
#' Class-specific methods such as \code{summary.clsp()},
#' \code{corr.clsp()}, and \code{ttest.clsp()} can be used to
#' extract, analyze, and summarize the results.
#'
#' @seealso \link[CVXR]{CVXR}
#'
#' @examples
#' \dontrun{
#'   ## Example: CMLS (RP) estimation with stationary-point constraints
#'
#'   set.seed(123456789)
#'
#'   # sample (dataset)
#'   k  <- 500L                                        # number of observations
#'   p  <- 6L                                          # number of regressors
#'   c0 <- 1                                           # sum of coefficients
#'
#'   D        <- matrix(NA_real_, nrow = k, ncol = p)
#'   D[, 1]   <- 1.0                                   # constant
#'   D[, 2:p] <- matrix(rnorm(k * (p - 1)), k, p - 1)
#'
#'   b_true   <- rnorm(p)
#'   b_true   <- (b_true / sum(b_true)) * c0           # normalize to sum = c
#'
#'   e        <- matrix(rnorm(k), ncol = 1)
#'   y        <- D %*% b_true + e
#'
#'   # build blocks for CLSP (CMLS)
#'   b <- rbind(
#'       matrix(c0, ncol = 1),                         # sum of coefficients
#'       matrix(0,  nrow = k - 2, ncol = 1),
#'       matrix(0,  nrow = k - 1, ncol = 1),
#'       matrix(y,  ncol = 1)
#'   )
#'
#'   C <- rbind(
#'       matrix(1, nrow = 1, ncol = p),                # row of ones
#'       diff(D, differences = 2),                     # 2nd differences
#'       diff(D, differences = 1)                      # 1st differences
#'   )
#'
#'   # diagonal sign-matrix for 2nd differences
#'   S <- rbind(
#'       matrix(0, nrow = 1,      ncol = k - 2),
#'       diag(sign(diff(as.numeric(y), differences = 2))),
#'       matrix(0, nrow = k - 1,  ncol = k - 2)
#'   )
#'
#'   # model
#'   model <- rclsp::clsp(
#'       problem = "cmls",
#'       b       = b,
#'       C       = C,
#'       S       = S,
#'       M       = D,
#'       r       = 1L,                                 # no refinement
#'       alpha   = 1.0                                 # MNBLUE solution
#'   )
#'
#'   # results
#'   print("true beta (x_M):")
#'   print(round(b_true, 4))
#'
#'   print("beta hat (x_M hat):")
#'   print(round(model$x, 4))
#'
#'   print(model)
#'
#'   # bootstrap t-test
#'   tt <- rclsp::ttest(
#'       model,
#'       sample_size  = 30L,
#'       seed         = 123456789L,
#'       distribution = rnorm,
#'       partial      = TRUE
#'   )
#'
#'   print("Bootstrap t-test:")
#'   print(tt)
#' }
#'
#' @importFrom methods as
#' @importFrom stats rnorm sd setNames
#'
#' @export
clsp <- function(problem="", C=NULL, S=NULL, M=NULL, b=NULL, m=NULL, p=NULL,
                 i=1L, j=1L, zero_diagonal=FALSE, r=1L, Z=NULL, rcond=FALSE,
                 tolerance=NULL, iteration_limit=NULL, final=TRUE, alpha=NULL,
                 ...) {
  object <- list(
    A               = NULL,
    C_idx           = c(NA_integer_, NA_integer_),
    b               = NULL,
    Z               = NULL,
    tolerance       = if (is.null(tolerance)) sqrt(.Machine$double.eps)
    else as.numeric(tolerance),
    iteration_limit = if (is.null(iteration_limit)) 50L
    else as.integer(iteration_limit),
    r               = 0L,
    zhat            = NULL,
    final           = isTRUE(final),
    alpha           = if (is.null(alpha)) NULL else as.numeric(alpha),
    z               = NULL,
    x               = NULL,
    y               = NULL,
    kappaC          = NA_real_,
    kappaB          = NA_real_,
    kappaA          = NA_real_,
    rmsa            = NA_real_,
    rmsa_i          = NULL,
    rmsa_dkappaC    = NULL,
    rmsa_dkappaB    = NULL,
    rmsa_dkappaA    = NULL,
    rmsa_dnrmse     = NULL,
    rmsa_dzhat      = NULL,
    rmsa_dz         = NULL,
    rmsa_dx         = NULL,
    r2_partial      = NA_real_,
    nrmse           = NA_real_,
    nrmse_partial   = NA_real_,
    nrmse_ttest     = NULL,
    z_lower         = NULL,
    z_upper         = NULL,
    x_lower         = NULL,
    x_upper         = NULL,
    y_lower         = NULL,
    y_upper         = NULL,
    seed            = 123456789L,
    distribution    = rnorm
  )
  
  # solve the CLSP problem
  dots <- list(...)
  object <- do.call(.solve, c(list(object, problem=problem,
                                   C=C, S=S, M=M, b=b, m=m, p=p, i=i, j=j,
                                   zero_diagonal=zero_diagonal,
                                   r=r, Z=Z, rcond=rcond, tolerance=tolerance,
                                   iteration_limit=iteration_limit,
                                   final=final, alpha=alpha),
                              dots))
  
  class(object) <- "clsp"
  object
}
################################################################################
# Ancillary functions
################################################################################
.solve.instance <- function(object, problem="", C=NULL, S=NULL, M=NULL, b=NULL,
                            m=NULL, p=NULL, i=1L, j=1L, zero_diagonal=FALSE,
                            r=1L, Z=NULL, rcond=FALSE, tolerance=NULL,
                            iteration_limit=NULL, final=NULL, alpha=NULL, ...) {
  # (A), (b) Construct a conformable canonical form for the CLSP estimator
  if (!is.null(C) || !is.null(M) || (!is.null(m) && !is.null(p))) {
    object <- canonize.clsp(object, problem, C, S, M, NULL, b,
                            m, p, i, j, zero_diagonal)
  } else if (is.null(object$A))
    stop("At least one of C, M, m, or p must be provided.")
  if (nrow(object$A) != nrow(object$b))
    stop(sprintf(paste0("The matrix A and vector b must have the same number ",
                        "of rows: A has %d, b has %d"),
                 nrow(object$A), nrow(object$b)))
  
  # (zhat) (Iterated if r > 1) first-step estimate
  if      (r < 1L)   stop("Number of refinement iterations r must be \u2265 1.")
  if      (!is.null(Z))               object$Z <- Z
  else if (is.null(object$Z))         object$Z <- diag(ncol(object$A))
  else                                object$Z <- object$Z[1:ncol(object$A),
                                                           1:ncol(object$A)]
  if      (!is.null(tolerance))       object$tolerance       <- tolerance
  if      (!is.null(iteration_limit)) object$iteration_limit <- iteration_limit
  if      (!isTRUE(all.equal(object$Z,                    t(object$Z),
                             tolerance=object$tolerance)) ||
           !isTRUE(all.equal(object$Z %*% object$Z,         object$Z,
                             tolerance=object$tolerance)) ||
           nrow(object$Z) != ncol(object$A))
    stop(sprintf(paste0("Matrix Z must be symmetric, idempotent and match ",
                        "the number of columns in A: expected (%d,%d), ",
                        "got (%d,%d)"), ncol(object$A), ncol(object$A),
                 nrow(object$Z), ncol(object$Z)))
  for (n_iter in seq_len(if (nrow(object$A) > object$C_idx[1]) r else 1L)) {
    # save NRMSE from the previous step, construct Q and Z
    if (n_iter > 1L) {
      res        <- object$b - object$A %*% object$zhat
      nrmse_prev <- .nrmse.r2(object, res=res)
      Q          <- diag(as.numeric(-sign(res[(object$C_idx[1]    +
                                                 1):nrow(object$A), ,
                                              drop=FALSE])))
      object     <- canonize.clsp(object, problem, C, S, M, Q, b, m, p, i, j,
                                  zero_diagonal)
      Z_delta <- ncol(object$A) - nrow(object$Z)
      if (Z_delta > 0L) {                              # augment Z by I
        object$Z <- rbind(cbind(object$Z,
                                matrix(0, nrow(object$Z), Z_delta)),
                          cbind(matrix(0, Z_delta, ncol(object$Z)),
                                diag(Z_delta)))
      }
    }
    # solve via the Bott–Duffin inverse
    object$zhat <- with(svd(M <- object$Z %*% crossprod(object$A) %*% object$Z),
                        v %*% diag(ifelse(d > ((     if (isFALSE(rcond))
                          max(dim(M))      *
                            .Machine$double.eps
                          else if (isTRUE(rcond))
                            object$tolerance
                          else rcond) * max(d)), 1/d, 0),
                          length(d)) %*% t(u)            %*%
                          object$Z %*% t(object$A)
    ) %*% object$b
    object$nrmse <- .nrmse.r2(object, res=object$b - object$A %*% object$zhat)
    # break on convergence
    object$r <- n_iter
    if (n_iter > 1L && (abs(object$nrmse - nrmse_prev) < object$tolerance ||
                        n_iter > object$iteration_limit))              break
  }    
  if (!all(is.finite(object$zhat))) {
    object$zhat <- NA_real_
    stop("Pseudoinverse estimate zhat failed")
  }
  
  # (z) Final solution (if available), or set object$z = object$zhat
  if (!is.null(final)) object$final <- isTRUE(final)
  if (!is.null(alpha)) object$alpha <- max(0, min(1, as.numeric(alpha)))
  if (isTRUE(object$final)) {
    # build a convex problem (p_cvx) and its solver (c_cvx)
    A_csc <- as(Matrix::Matrix(object$A, sparse=TRUE), "dgCMatrix")
    z_cvx <- CVXR::Variable(ncol(A_csc))
    d_cvx <- z_cvx - as.numeric(object$zhat)
    if         (isTRUE(all.equal(object$alpha, 0))) {  # Lasso
      f_obj <- CVXR::norm1(d_cvx)
      s_cvx <- "ECOS"
    } else  if (isTRUE(all.equal(object$alpha, 1))) {  # Ridge
      f_obj <- CVXR::sum_squares(d_cvx)
      s_cvx <- "OSQP"
    } else  {                                          # Elastic Net
      f_obj <- (1 - object$alpha) * CVXR::norm1(d_cvx) +
        object$alpha * CVXR::sum_squares(d_cvx)
      s_cvx <- "SCS"
    }
    c_cvx <- list(A_csc %*% z_cvx == as.numeric(object$b))
    p_cvx <- CVXR::Problem(CVXR::Minimize(f_obj), c_cvx)
    # solve
    dots       <- list(...)                            # pass arguments
    dots$rcond <- NULL
    solution   <- try(do.call(p_cvx$solve, c(list(solver=s_cvx, verbose=FALSE),
                                             dots)),              silent=TRUE)
    if (inherits(solution, "try-error") || is.null(z_cvx$value)) {
      warning(sprintf("Step 2 infeasible (%s); falling back", p_cvx@status),
              call. = FALSE)
      object$z     <- object$zhat
    } else {
      object$z     <- as.numeric(z_cvx$value)
      object$nrmse <- .nrmse.r2(object)
    }
  } else {
    object$z   <- object$zhat
  }
  
  # (x), (y) Variable and slack components of z
  object$x <- matrix(object$z[1:object$C_idx[2], ,        drop=FALSE],
                     nrow=if (!is.null(m)) m else object$C_idx[2],
                     ncol=if (!is.null(p)) p else 1,      byrow=TRUE)
  object$y <- (if (ncol(object$A) > object$C_idx[2])
    object$z[(object$C_idx[2] + 1):length(object$z), ,    drop=FALSE]
    else matrix(numeric(0), ncol=1))
  
  # (kappaC), (kappaB), (kappaA) Condition numbers
  object$kappaC <- kappa(object$A[1:object$C_idx[1], ,    drop=FALSE])
  object$kappaB <- kappa(object$A %*% MASS::ginv(object$A[1:object$C_idx[1], ,
                                                          drop=FALSE]))
  object$kappaA <- kappa(object$A)
  
  # (r2_partial), (nrmse_partial) M-block-based statistics
  if (nrow(object$A) > object$C_idx[1]) {
    object$r2_partial    <- .nrmse.r2(object, r2=TRUE, partial=TRUE)
    object$nrmse_partial <- .nrmse.r2(object,          partial=TRUE)
  }
  
  # (z_lower), (z_upper) Condition-weighted confidence band
  b_norm <- sqrt(sum(object$b^2))
  dz     <- if (isTRUE(all.equal(b_norm, 0))) Inf                       else
    object$kappaA                 *
    sqrt(sum((object$b - object$A %*%
                matrix(object$z, ncol=1))^2)) /
    b_norm
  object$z_lower <- object$z * (1 - dz)
  object$z_upper <- object$z * (1 + dz)
  
  # (x_lower), (x_upper), (y_lower), (y_upper)
  object$x_lower <- matrix(object$z_lower[1:object$C_idx[2], , drop=FALSE],
                           nrow=if (!is.null(m)) m else object$C_idx[2],
                           ncol=if (!is.null(p)) p else 1,     byrow=TRUE)
  object$x_upper <- matrix(object$z_upper[1:object$C_idx[2], , drop=FALSE],
                           nrow=if (!is.null(m)) m else object$C_idx[2],
                           ncol=if (!is.null(p)) p else 1,     byrow=TRUE)
  object$y_lower <- (if (ncol(object$A) > object$C_idx[2])
    object$z_lower[(object$C_idx[2] + 1):length(object$z), ,
                   drop=FALSE]
    else matrix(numeric(0), ncol=1))
  object$y_upper <- (if (ncol(object$A) > object$C_idx[2])
    object$z_upper[(object$C_idx[2] + 1):length(object$z), ,
                   drop=FALSE]
    else matrix(numeric(0), ncol=1))
  
  object
}
.solve <- function(object, tolerance=NULL, alpha=NULL, ...) {
  dots           <- list(...)
  dots$tolerance <- NULL
  dots$final     <- NULL
  dots$alpha     <- NULL
  to.alpha       <- function(a) max(0, min(1, as.numeric(a)))
  to.nrmse       <- function(n) if (is.finite(n)) n else Inf
  
  # process alpha
  if      (!is.null(tolerance))
    object$tolerance <- tolerance
  if      (!is.null(alpha) && (is.finite(alpha) && is.atomic(alpha) &&
                               length(alpha) == 1L))
    object$alpha     <- to.alpha(alpha)
  else if (!is.null(alpha) && (is.list(alpha)   || is.atomic(alpha) &&
                               length(alpha) >  1L)) {
    alpha  <- sapply(alpha, to.alpha)
    result <- numeric(length(alpha))
    idx    <- 1L
    for (a in alpha) {
      result[idx] <- suppressWarnings(to.nrmse(do.call(.solve.instance,
                                     c(list(object, tolerance=object$tolerance,
                                            final=NULL, alpha=a), dots))$nrmse))
      idx         <- idx + 1L
    }
    object$alpha  <- if (length(result) > 0L) alpha[which.min(result)]  else
      NULL
  }
  if (is.null(object$alpha)) {                         # error rule
    nrmse.alpha0  <- suppressWarnings(to.nrmse(do.call(.solve.instance,
                                     c(list(object, tolerance=object$tolerance,
                                            final=NULL, alpha=0), dots))$nrmse))
    nrmse.alpha1  <- suppressWarnings(to.nrmse(do.call(.solve.instance,
                                     c(list(object, tolerance=object$tolerance,
                                            final=NULL, alpha=1), dots))$nrmse))
    denominator   <- nrmse.alpha0 + nrmse.alpha1 + object$tolerance
    object$alpha  <- if (is.finite(denominator) && denominator > 0)
                     to.alpha(nrmse.alpha0 / denominator)               else
        0.5
  }
  
  do.call(.solve.instance, c(list(object, tolerance=object$tolerance,
                                  final=NULL, alpha=object$alpha), dots))
}
