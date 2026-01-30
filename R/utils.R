#' Construct the canonical design matrix \eqn{A = [C | S; M | Q]} for CLSP.
#'
#' This method assembles the constraint matrix A from user-supplied or
#' internally generated components — C, S, M, and Q — and assigns the
#' corresponding right-hand side vector b. It is a required pre-step before
#' solving a Convex Least Squares Programming (CLSP) problem.
#'
#' Depending on the specified problem type, it can generate allocation,
#' tabular matrix, or modular constraints and enforce optional diagonal
#' exclusions. All missing blocks are padded to ensure conformability.
#'
#' @param object An object of class \code{"clsp"}.
#'
#' @param problem Character, optional. Structural template for matrix
#'   construction. One of:
#'   - `'ap'`   or `'tm'`: allocation or tabular matrix problem.
#'   - `'cmls'` or `'rp'`: constrained modular least squares or RP-type.
#'   - `''`     or other: General CLSP problems (user-defined C and/or M).
#'
#' @param C,S,M Numeric matrix or \code{NULL}. Blocks of the constraint
#'   matrix \eqn{A = [C | S; M | Q]}. If \code{C} and/or \code{M} are
#'   provided, the matrix A is constructed accordingly. If both are
#'   \code{NULL} and A is not yet defined, an error is raised.
#'
#' @param Q Numeric matrix or \code{NULL}. Externally supplied residual slack
#'   matrix used to adjust inequality constraints in M. Required only when
#'   \eqn{r > 1}. Encodes the sign pattern of residuals from the previous
#'   iteration and is used to construct the \eqn{[C | S; M | Q]} canonical
#'   form. Defaults to a conformable zero matrix on the first iteration.
#'
#' @param b Numeric vector or \code{NULL}. Right-hand side vector. Must have
#'   as many rows as \eqn{A}. Required.
#'
#' @param m,p Integer or \code{NULL}. Dimensions of
#'   \eqn{X \in \mathbb{R}^{m \times p}}, relevant for allocation problems
#'   ('ap').
#'
#' @param i,j Integer, default = \code{1}. Grouping sizes for row and column
#'   sum constraints in AP problems.
#'
#' @param zero_diagonal Logical, default = \code{FALSE}. If \code{TRUE},
#'   enforces structural zero diagonals via identity truncation.
#'
#' @section Attributes Set:
#' \describe{
#'   \item{\code{A}}{Numeric matrix. Canonical design matrix constructed from
#'     (C, S, M, Q).}
#'   \item{\code{C_idx}}{Integer vector of length 2 indicating the size of
#'     the C block.}
#'   \item{\code{b}}{Numeric vector. Conformable right-hand side vector.}
#' }
#'
#' @return An updated object of class \code{"clsp"}.
#' @export
canonize      <- function(object, problem="", C=NULL, S=NULL, M=NULL, Q=NULL,
                          b=NULL, m=NULL, p=NULL, i=1L, j=1L,
                          zero_diagonal=FALSE) UseMethod("canonize")
#' @export
canonize.clsp <- function(object, problem="", C=NULL, S=NULL, M=NULL, Q=NULL,
                          b=NULL, m=NULL, p=NULL, i=1L, j=1L,
                          zero_diagonal=FALSE) {
  # (b) Ensure the right-hand side is defined and set `object$b`
  if (is.null(b)) stop("Right-hand side vector b must be provided.")
  object$b <- matrix(as.numeric(b), ncol=1, byrow=TRUE)
  
  if (!is.null(C)) storage.mode(C) <- "double"
  if (!is.null(S)) storage.mode(S) <- "double"
  if (!is.null(M)) storage.mode(M) <- "double"
  if (!is.null(Q)) storage.mode(Q) <- "double"
  # (A) Option 1. AP (TM) problems with an optional zero diagonal
  if (grepl("ap|tm", tolower(problem))) {
    m <- as.integer(m)
    p <- as.integer(p)
    i <- as.integer(i)
    j <- as.integer(j)
    if (is.null(m) || is.null(p)) stop("Both m and p must be specified.")
    if ((m %% i) != 0) stop(sprintf("m  %d must be divisible by i  %d", m, i))
    if ((p %% j) != 0) stop(sprintf("p  %d must be divisible by j  %d", p, j))
    # Construct the C block using Kronecker product
    row_groups <- kronecker(kronecker(diag(m %/% i),   matrix(1, 1, i)),
                            matrix(1, 1, p))
    col_groups <- kronecker(matrix(1, 1, m), kronecker(diag(p %/% j),
                                                       matrix(1, 1, j)))
    if (!is.null(C) && ncol(C) != ncol(row_groups))
      stop(sprintf("C must have %d columns", ncol(row_groups)))
    C <- if (!is.null(C)) rbind(row_groups, col_groups, C)       else
      rbind(row_groups, col_groups)
    # append an optional identity matrix to M, remove duplicates
    if (isTRUE(zero_diagonal)) {
      M_flag <- is.null(M) || length(M) == 0
      M_diag <- matrix(0, nrow=min(m, p), ncol=m * p)
      for (k in seq_len(min(m, p))) M_diag[k, (k - 1) * p + k] <- 1
      tmp    <- if (M_flag) M_diag                               else
                   rbind(M, M_diag)
      ord <- do.call(order, as.data.frame(tmp))
      idx <- ord[!duplicated(tmp[ord, , drop=FALSE])]
      M   <- tmp[idx, ,                 drop=FALSE]
      object$b <- rbind(
        object$b[seq_len(nrow(C)), ,    drop=FALSE],
        (if (M_flag)
          matrix(0, min(m, p), 1)
         else
           rbind(object$b[(nrow(C) + 1):nrow(object$b), , drop=FALSE],
                 matrix(0, min(m, p), 1))
        )[idx, , drop = FALSE]
      )
    }
  }
  
  # (A) Option 2. CMLS and RP problems
  if (grepl("cmls|rp", tolower(problem)) && (is.null(C) || is.null(M)))
    stop("Both C and M must be provided.")
  
  # (A) Option 3. General problems
  if (is.null(C) && is.null(M)) stop("At least one of C or M must be provided.")
  
  # (A) Convert missing blocks to conformable zero matrices
  n_col <- if (!is.null(C)) ncol(C) else ncol(M)
  C     <- if (!is.null(C)) C else matrix(0, nrow=0,       ncol=n_col)
  M     <- if (!is.null(M)) M else matrix(0, nrow=0,       ncol=n_col)
  S     <- if (!is.null(S)) S else matrix(0, nrow=nrow(C), ncol=0)
  Q     <- if (!is.null(Q)) Q else matrix(0, nrow=nrow(M), ncol=0)
  
  if (nrow(C) != nrow(S))
    stop(sprintf("C and S must have the same number of rows: %d vs %d",
                 nrow(C), nrow(S)))
  if (ncol(C) != ncol(M))
    stop(sprintf("C and M must have the same number of columns: %d vs %d",
                 ncol(C), ncol(M)))
  
  # (A) Pad C and Q with zeros and set object$.A and object$C_idx
  object$A     <- rbind(cbind(C, S, matrix(0, nrow=nrow(C), ncol=ncol(Q))),
                        cbind(M,    matrix(0, nrow=nrow(M), ncol=ncol(S)), Q))
  object$C_idx <- c(nrow(C), ncol(C))
  
  object
}
#' Compute the structural correlogram of the CLSP constraint system.
#'
#' This method performs a row-deletion sensitivity analysis on the canonical
#' constraint matrix \eqn{[C | S]}, denoted as \eqn{C_{\text{canon}}}, and
#' evaluates the marginal effect of each constraint row on numerical stability,
#' angular alignment, and estimator sensitivity.
#'
#' For each row \eqn{i} in \eqn{C_{\text{canon}}}, it computes:
#' \itemize{
#'   \item The Root Mean Square Alignment (\eqn{\mathrm{RMSA}_i}) with all
#'         other rows \eqn{j \ne i}.
#'   \item The change in condition numbers \eqn{\kappa(C)}, \eqn{\kappa(B)}, and
#'         \eqn{\kappa(A)} when row \eqn{i} is deleted.
#'   \item The effect on estimation quality: changes in NRMSE, \eqn{\hat{z}},
#'         \eqn{z}, and \eqn{x}.
#' }
#'
#' Additionally, it computes the total RMSA statistic across all rows,
#' summarizing the overall angular alignment of the constraint block.
#'
#' @param object An object of class \code{"clsp"}.
#'
#' @param reset Logical, default = \code{FALSE}.
#'   If \code{TRUE}, forces recomputation of all diagnostic values.
#'
#' @param threshold Numeric, default = \code{0}.
#'   If positive, limits the output to constraints with
#'   \eqn{\mathrm{RMSA}_i \ge \text{threshold}}.
#'
#' @return
#' A named list containing per-row diagnostic values:
#' \describe{
#'   \item{constraint}{Vector of constraint indices (1-based).}
#'   \item{rmsa_i}{List of \eqn{\mathrm{RMSA}_i} values.}
#'   \item{rmsa_dkappaC}{List of \eqn{\Delta\kappa(C)} after deleting row i.}
#'   \item{rmsa_dkappaB}{List of \eqn{\Delta\kappa(B)} after deleting row i.}
#'   \item{rmsa_dkappaA}{List of \eqn{\Delta\kappa(A)} after deleting row i.}
#'   \item{rmsa_dnrmse}{List of \eqn{\Delta\mathrm{NRMSE}} after deleting row i.}
#'   \item{rmsa_dzhat}{List of \eqn{\Delta\hat{z}} after deleting row i.}
#'   \item{rmsa_dz}{List of \eqn{\Delta z} after deleting row i.}
#'   \item{rmsa_dx}{List of \eqn{\Delta x} after deleting row i.}
#' }
#'
#' @export
corr      <- function(object, reset=FALSE, threshold=0) UseMethod("corr")
#' @export
corr.clsp <- function(object, reset=FALSE, threshold=0) {
  # (RMSA) Total RMSA
  if (is.null(object$rmsa)                     || isTRUE(reset)) {
    k           <- object$C_idx[1]
    p           <- object$C_idx[2]
    C_canon     <- object$A[1:k, , drop=FALSE]
    norms       <- sqrt(rowSums(C_canon^2))
    object$rmsa <- {
      if (k < 2) {
        NA_real_
      } else {
        cos <- c()
        for (i in 1:(k - 1))
          for (j in (i + 1):k) cos <- c(cos, sum(C_canon[i, ] * C_canon[j, ]) /
                                          (norms[i] * norms[j]))
        sqrt(2 / (k * (k - 1)) * sum(cos^2))
      }
    }
  }
  
  # (RMSA) Constraint-wise RMSA, changes in condition numbers, and GoF
  if (length(object$rmsa_i) != object$C_idx[1] || isTRUE(reset)) {
    tmp                 <- unserialize(serialize(object, NULL))
    k                   <- object$C_idx[1]
    p                   <- object$C_idx[2]
    C_canon             <- object$A[1:k, , drop = FALSE]
    norms               <- sqrt(rowSums(C_canon^2))
    object$rmsa_i       <- rep(NA_real_, k)
    object$rmsa_dkappaC <- rep(NA_real_, k)
    object$rmsa_dkappaB <- rep(NA_real_, k)
    object$rmsa_dkappaA <- rep(NA_real_, k)
    object$rmsa_dnrmse  <- rep(NA_real_, k)
    object$rmsa_dzhat   <- vector("list", k)
    object$rmsa_dz      <- vector("list", k)
    object$rmsa_dx      <- vector("list", k)
    for (i in 1:k) {
      tmp$A                  <- object$A[-i, , drop = FALSE]
      tmp$b                  <- object$b[-i, , drop = FALSE]
      tmp$C_idx              <- c(k - 1, p)
      suppressWarnings(tmp   <- .solve(tmp))
      object$rmsa_i[i]       <- {
        cos_i <- c()
        for (j in 1:k)
          if (j != i) cos_i <- c(cos_i, sum(C_canon[i, ] * C_canon[j, ]) /
                                   (norms[i] * norms[j]))
        sqrt(1 / (k - 1) * sum(cos_i^2))
      }
      object$rmsa_dkappaC[i] <- tmp$kappaC - object$kappaC
      object$rmsa_dkappaB[i] <- tmp$kappaB - object$kappaB
      object$rmsa_dkappaA[i] <- tmp$kappaA - object$kappaA
      object$rmsa_dnrmse[i]  <- tmp$nrmse  - object$nrmse
      object$rmsa_dzhat[[i]] <- tmp$zhat   - object$zhat
      object$rmsa_dz[[i]]    <- tmp$z      - object$z
      object$rmsa_dx[[i]]    <- matrix(tmp$x,    ncol = 1, byrow=TRUE) -
        matrix(object$x, ncol = 1, byrow=TRUE)
    }
  }
  
  # Return the correlogram
  indices <- which(!is.na(object$rmsa_i) & object$rmsa_i >= threshold)
  list(
    constraint   = indices,
    rmsa_i       = object$rmsa_i[indices],
    rmsa_dkappaC = object$rmsa_dkappaC[indices],
    rmsa_dkappaB = object$rmsa_dkappaB[indices],
    rmsa_dkappaA = object$rmsa_dkappaA[indices],
    rmsa_dnrmse  = object$rmsa_dnrmse[indices],
    rmsa_dzhat   = object$rmsa_dzhat[indices],
    rmsa_dz      = object$rmsa_dz[indices],
    rmsa_dx      = object$rmsa_dx[indices]
  )
}
#' Perform bootstrap or Monte Carlo t-tests on the NRMSE statistic from
#' the CLSP estimator.
#'
#' This function either (a) resamples residuals via a nonparametric bootstrap
#' to generate an empirical NRMSE sample, or (b) produces synthetic right-hand
#' side vectors \code{b} from a user-defined or default distribution and
#' re-estimates the model. It tests whether the observed NRMSE significantly
#' deviates from the null distribution of resampled or simulated NRMSE values.
#'
#' @param object An object of class \code{"clsp"}.
#'
#' @param reset Logical, default = \code{FALSE}.
#'   If \code{TRUE}, forces recomputation of the NRMSE null distribution.
#'
#' @param sample_size Integer, default = \code{50}.
#'   Size of the Monte Carlo simulated sample under H0.
#'
#' @param seed Integer or \code{NULL}, default = \code{NULL}.
#'   Optional random seed to override the default.
#'
#' @param distribution Function or \code{NULL}, default = \code{NULL}.
#'   Distribution for generating synthetic \code{b} vectors. One of:
#'   \code{rnorm}, \code{runif}, or a custom RNG function. Defaults to
#'   standard normal.
#'
#' @param partial Logical, default = \code{FALSE}.
#'   If \code{TRUE}, runs the t-test on the partial NRMSE: during simulation,
#'   the C-block entries are preserved and the M-block entries are simulated.
#'
#' @param simulate Logical, default = \code{FALSE}.
#'   If \code{TRUE}, performs a parametric Monte Carlo simulation by generating
#'   synthetic right-hand side vectors \code{b}. If \code{FALSE} (default),
#'   executes a nonparametric bootstrap procedure on residuals without
#'   re-estimation.
#'
#' @return
#' A named list containing test results and null distribution statistics:
#' \describe{
#'   \item{p_one_left}{P(nrmse \eqn{\le} null mean)}
#'   \item{p_one_right}{P(nrmse \eqn{\ge} null mean)}
#'   \item{p_two_sided}{2-sided t-test p-value}
#'   \item{nrmse}{Observed value}
#'   \item{mean_null}{Mean of null distribution}
#'   \item{std_null}{Standard deviation of null distribution}
#' }
#'
#' @export
ttest      <- function(object, reset=FALSE, sample_size=50L,
                       seed=NULL, distribution=NULL,
                       partial=FALSE, simulate=FALSE) UseMethod("ttest")
#' @export
ttest.clsp <- function(object, reset=FALSE, sample_size=50L,
                       seed=NULL, distribution=NULL,
                       partial=FALSE, simulate=FALSE) {
  # Set the seed, RNG configuration, and distribution
  if (!is.null(seed))        object$seed  <- seed
  if (is.null(distribution)) distribution <- object$distribution
  else if (!is.function(distribution))
    stop("distribution must be a function accepting one integer argument n.")
  set.seed(object$seed)
  
  # (t-test) Bootstrap-resampled or simulated NRMSE distribution under H0
  if (length(object$nrmse_ttest) != sample_size || isTRUE(reset)) {
    if (isTRUE(partial) && nrow(object$A) == object$C_idx[1]) {
      warning("No M-block present in A; falling back to full NRMSE t-test.")
      partial <- FALSE
    }
    object$nrmse_ttest <- rep(NA_real_, sample_size)
    # (re)generate a nonparametric bootstrap sample
    if (!isTRUE(simulate)) {
      res <- object$b - object$A %*% object$zhat
      for (i in seq_len(sample_size)) {
        res_bs                <- matrix(sample(as.vector(res), size=length(res),
                                               replace=TRUE), ncol=1,
                                        byrow=TRUE)
        object$nrmse_ttest[i] <-as.numeric(.nrmse.r2(object, res=res_bs,
                                                     partial=isTRUE(partial)))
      }
      # (re)generate a parametric Monte Carlo sample
    } else {
      tmp <- unserialize(serialize(object, NULL))
      for (i in seq_len(sample_size)) {
        tmp$b                 <- if (!isTRUE(partial))
          matrix(distribution(length(object$b)), ncol=1, byrow=TRUE)  else
            rbind(matrix(object$b[1:object$C_idx[1]],
                         ncol=1, byrow=TRUE),
                  matrix(distribution(length(object$b) - object$C_idx[1]),
                         ncol=1, byrow=TRUE))            # simulate b_M only
        suppressWarnings(tmp  <- .solve(tmp))
        object$nrmse_ttest[i] <- if (!isTRUE(partial)) tmp$nrmse      else
          tmp$nrmse_partial
      }
    }
  }
  
  # Return the t-test
  nrmse_null <- as.numeric(object$nrmse_ttest)
  mean_null  <- mean(nrmse_null,      na.rm=TRUE)
  std_null   <- stats::sd(nrmse_null, na.rm=TRUE)
  if (!is.finite(std_null) || std_null == 0) {
    p_left   <- p_right <- p_two <- NA_real_
  } else {
    t_stat   <- ((if (!isTRUE(partial)) object$nrmse             else
      object$nrmse_partial) - mean_null) / (std_null / sqrt(sample_size))
    p_left   <- stats::pt(t_stat, df=sample_size-1)
    p_right  <- 1 - p_left
    p_two    <- 2 * stats::pt(-abs(t_stat), df=sample_size-1)
  }
  list(
    p_one_left  = p_left,
    p_one_right = p_right,
    p_two_sided = p_two,
    nrmse       = if (!isTRUE(partial)) object$nrmse             else
      object$nrmse_partial,
    mean_null   = mean_null,
    std_null    = std_null
  )
}
#' @export
print.clsp <- function(x, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else cat("clsp(...)\n")
}
#' @export
summary.clsp <- function(object, ...) {
  # Declare "out"
  out <- list(
    # Metadata
    call            = if (!is.null(object$call)) object$call     else NULL,
    tolerance       = object$tolerance,
    iteration_limit = object$iteration_limit,
    inverse         = ifelse(!is.null(object$Z) &&
                             !isTRUE(all.equal(object$Z, diag(nrow(object$Z)))),
                             "Bott-Duffin", "Moore-Penrose"),
    r               = object$r,
    final           = object$final,
    alpha           = object$alpha,
    # Estimates
    zhat            = object$zhat,
    z               = object$z,
    x               = object$x,
    y               = object$y,
    # Numerical stability
    kappaC          = object$kappaC,
    kappaB          = object$kappaB,
    kappaA          = object$kappaA,
    # Goodness of fit
    nrmse           = object$nrmse,
    nrmse_partial   = object$nrmse_partial,
    r2_partial      = object$r2_partial,
    # Confidence bands
    z_lower         = .summary(object$z_lower),
    x_lower         = .summary(object$x_lower),
    y_lower         = .summary(object$y_lower),
    z_upper         = .summary(object$z_upper),
    x_upper         = .summary(object$x_upper),
    y_upper         = .summary(object$y_upper),
    # Simulation
    seed            = object$seed,
    distribution    = object$distribution
  )
  
  # Expand "out" (if possible)
  if (!is.null(object$rmsa) && !is.na(object$rmsa)) {
    out$rmsa         <- object$rmsa
    out$rmsa_i       <- .summary(object$rmsa_i)
    out$rmsa_dkappaC <- .summary(object$rmsa_dkappaC)
    out$rmsa_dkappaB <- .summary(object$rmsa_dkappaB)
    out$rmsa_dkappaA <- .summary(object$rmsa_dkappaA)
    out$rmsa_dnrmse  <- .summary(object$rmsa_dnrmse)
    out$rmsa_dzhat   <- .summary(object$rmsa_dzhat)
    out$rmsa_dz      <- .summary(object$rmsa_dz)
    out$rmsa_dx      <- .summary(object$rmsa_dx)
  }
  
  # Return "out"
  class(out) <- "summary.clsp"
  out
}
#' @export
print.summary.clsp <- function(x, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else cat("clsp(...)\n")
  
  cat("\nEstimator Configuration:\n")
  .print(x, "inverse",   label="Generalized inverse")
  .print(x, "r",         label="Iterations (r)")
  .print(x, "tolerance", label="Tolerance")
  .print(x, "final",     label="Final correction")
  .print(x, "alpha",     label="Regularization (\u03B1)")
  cat("\nNumerical Stability:\n")
  .print(x, c("kappaC", "kappaB", "kappaA"))
  if (!is.null(x$rmsa) && !is.na(x$rmsa)) {
    .print(x, "rmsa")
    .print(x,  c("i", "dkappaC", "dkappaB", "dkappaA",
                 "dnrmse", "dzhat", "dz", "dx"),  "rmsa_")
  }
  cat("\nGoodness of Fit:\n")
  .print(x, c("r2_partial", "nrmse", "nrmse_partial"))
  .print(x, c("z_lower", "z_upper", "x_lower", "x_upper", "y_lower", "y_upper"))
  invisible(x)
}
################################################################################
# Ancillary functions
################################################################################
.pop.sd <- function(x){                                # population SD
  n  <- length(x)
  if (n == 0) return(NA_real_) else if (n == 1) return(0) else {
    sd <- stats::sd(x, na.rm=TRUE) * sqrt((n - 1) / n)
    if (is.nan(sd))  NA_real_  else sd
  }
}
.nrmse.r2 <- function(object, res=NULL, sdb=NULL, r2=FALSE, partial=FALSE) {
  if (isTRUE(partial) && nrow(object$A) == object$C_idx[1]) return(NA_real_)
  res <- if (is.null(res)) object$b - object$A %*% object$z else res
  res <- if (!isTRUE(partial))      res                     else
    res[(object$C_idx[1] + 1):nrow(object$A), ,
        drop=FALSE]
  b   <- if (!isTRUE(partial)) object$b else
    object$b[(object$C_idx[1] + 1):nrow(object$b), ,
             drop=FALSE]
  sdb <- if (is.null(sdb)) .pop.sd(as.numeric(b))           else sdb
  if (isTRUE(r2)) {
    tss <- sum((as.numeric(b) - mean(as.numeric(b)))^2)
    if      (isTRUE(all.equal(tss, 0))) NA_real_            else
      1 -  sum(res^2)  /  tss
  } else if (isTRUE(all.equal(sdb, 0))) NA_real_            else
    sqrt(sum(res^2)) /  sqrt(nrow(b)) / sdb
}
.summary <- function(v) {
  if (is.null(v) || all(is.na(v)))
    return(setNames(rep(NA_real_, 4), c("min", "max", "mean", "sd")))
  v <- as.numeric(v)
  c(min = min(v, na.rm=TRUE), max=max(v, na.rm=TRUE),
    mean=mean(v, na.rm=TRUE), sd = sd(v, na.rm=TRUE))
}
.print <- function(x, attributes, base="", label=NULL) {
  attributes <- as.vector(attributes)
  if (!is.null(label)) label <- as.vector(label)
  for (idx in seq_along(attributes)) {
    attr <- attributes[idx]
    key  <- paste0(base, attr)
    v    <- x[[key]]
    lab  <- if (!is.null(label)) label[idx] else attr
    lab  <- paste0(lab, ":", strrep(" ", 20 - nchar(lab, type="width") - 1))
    if (length(v) == 1) {                              # print value
      width <- 15
      cat(sprintf("  %-20s %s\n", lab, .format(x, v, width)))
    } else {                                           # print list
      width <- 10
      fmt <- paste0("  %-14s  min=%",  width, "s", "  max=%",  width, "s",
                    "  mean=%", width, "s", "  sd=%",   width, "s\n")
      cat(sprintf(fmt, lab, .format(x, v[["min"]],  width),
                  .format(x, v[["max"]],  width),
                  .format(x, v[["mean"]], width),
                  .format(x, v[["sd"]],   width)))
    }
  }
}
.format <- function(object, x, width=10) {
  fmt <- paste0("%", width, "s")
  if (is.null(x) || is.na(x)) return(sprintf(fmt, "NA"))
  if (!is.numeric(x))         return(sprintf(fmt, as.character(x)))
  if (abs(x) <= 1e10 && x == as.integer(x))
    return(sprintf(fmt, as.character(as.integer(x))))
  dp <- max(3, floor(width / 2) - 1)
  ep <- max(2, dp - 1)
  fixed <- formatC(x, digits=dp, format="f", drop0trailing=FALSE)
  if (!(abs(x) >= 10^(-(dp + 1)) && abs(x) <  10^(width - dp) &&
        nchar(fixed) <= width))
    return(sprintf(fmt, formatC(x, digits=ep, format="e")))
  sprintf(fmt, fixed)
}
