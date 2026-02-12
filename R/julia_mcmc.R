#' Bayesian Spatial Factor Analysis
#'
#' Provides methods for Bayesian spatial factor analysis using
#' latent spatial factors. The implementation is written in Julia and uses
#' nearest-neighbor Gaussian process approximations, sparse matrix computations,
#' and iterative solvers to enable scalable Markov chain Monte Carlo sampling.
#' It integrates R-based spatial utilities for spatial ordering and neighbor
#' selection.
#'
#' @param coords Numeric matrix of spatial coordinates with dimension
#' \eqn{n \times 2}.
#' @param X Numeric design matrix of observed covariates with
#' dimension \eqn{n \times p}. Rows correspond to spatial locations
#' and columns correspond to predictors.
#' @param Y Numeric response matrix with dimension \eqn{n \times q},
#' where each row corresponds to a spatial location and each column
#' to a response variable.
#' @param K Integer specifying the number of latent spatial factors.
#' @param phi_K Numeric vector of length \eqn{K} giving the
#' spatial range parameter for each latent factor.
#' @param m Integer specifying the number of neighbors for nearest neighbor
#' calculation.
#' @param N_sam Integer specifying the number of MCMC samples.
#'
#' @return A list containing MCMC samples for regression and factor
#' loadings (\code{gamma_samples}), response variances
#' (\code{sigma_samples}), and latent spatial factors
#' (\code{F_samples}). Rows correspond to MCMC iterations and columns to
#' model parameters.
#' @importFrom JuliaCall julia_setup julia_source julia_call
#' @import spNNGP
#' @import GPvecchia

#' @export
#'
#' @examples
#' set.seed(1)
#'
#' n <- 50
#' p <- 2
#' q <- 5
#' K <- 2
#'
#' # Spatial coordinates
#' coords <- cbind(runif(n), runif(n))
#'
#' # Covariates and responses
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- matrix(rnorm(n * q), n, q)
#'
#' # Spatial range parameters (one per factor)
#' phi_K <- rep(1.0, K)
#'
#' # Run MCMC
#' samples <- julia_mcmc(
#'   coords = coords,
#'   X = X,
#'   Y = Y,
#'   phi_K = phi_K,
#'   K = K,
#'   m = 10,
#'   N_sam = 50
#' )
#'
#' # Inspect output
#' str(samples)
#'
julia_mcmc<-function(coords,X,Y,phi_K,K,m=10, N_sam = 5) { ## returns F, gamma (composed of lambda and beta), and sigma


  # Locate Julia backend shipped with the package
  jl_file <- system.file("julia", "PBSF.jl", package = "PBSF")

  ## --------------------------------------------------
  ## 2. Input validation (USER-FACING CHECKS)
  ## --------------------------------------------------

  # ensures the PBSF.jl file is found
  if (!nzchar(jl_file) || !file.exists(jl_file)) {
    stop("Could not find PBSF.jl in installed PBSF package.")
  }
  ## -------------------------
  ## 1. Type checks first
  ## -------------------------
  if (!is.matrix(coords)) {
    stop("coords must be a matrix")
  }
  if (!is.matrix(X)) {
    stop("X must be a matrix")
  }

  if (!is.matrix(Y)) {
    stop("Y must be a matrix")
  }

  if (!is.numeric(coords)) {
    stop("coords must be numeric")
  }

  if (!is.numeric(X)) {
    stop("X must be numeric")
  }

  if (!is.numeric(Y)) {
    stop("Y must be numeric")
  }

## -------------------------
## 2. Dimension checks
## -------------------------
    n <- nrow(coords)

    if (ncol(coords) != 2) {
      stop("coords must have exactly 2 columns for spNNGP.")
    }

    if (n == 0) {
      stop("coords must have at least one row")
    }

    if (nrow(X) == 0 || ncol(X) == 0) {
      stop("X must have at least one row and one column")
    }

    if (nrow(Y) == 0 || ncol(Y) == 0) {
      stop("Y must have at least one row and one column")
    }

    if (nrow(coords) != nrow(X) || nrow(coords) != nrow(Y)) {
      stop("coords, X, and Y must have the same number of rows")
    }


## -------------------------
## 3. Parameter checks
## -------------------------
  if (!is.numeric(K) || length(K) != 1 || K %% 1 != 0 || K <= 0) {
    stop("K must be a single integer")
  }

  if (!is.numeric(phi_K) || !is.atomic(phi_K) || is.matrix(phi_K) || length(phi_K) != K) {
    stop("phi_K must be a numeric vector of length K")
  }

  if (!is.numeric(m) || length(m) != 1 || m %% 1 != 0 || m <= 0) {
    stop("m must be a single positive integer")
  }


  if (m >= n) {
    stop("m must be less than the number of spatial locations.")
  }

  if (!is.numeric(N_sam) || length(N_sam) != 1 || N_sam %% 1 != 0 || N_sam <= 0) {
    stop("N_sam must be a single positive integer")
  }


## -------------------------
## 4. Missing data check
## -------------------------
    if (any(!is.finite(coords)) ||
        any(!is.finite(X)) ||
        any(!is.finite(Y)) ||
        any(!is.finite(phi_K))) {
      stop("All inputs must be non missing (no NA, NaN, or Inf values)") # ?? specify which one has missing
    }



  # Source Julia code (safe to call multiple times)
  JuliaCall::julia_setup()
  JuliaCall::julia_source(jl_file)
  samples<-JuliaCall::julia_call("wrapper",coords,X,Y,phi_K,as.integer(K),m=as.integer(m), N_sam=as.integer(N_sam))

  return(samples)

}

