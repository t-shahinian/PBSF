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
#' \eqn{n \times d}, where \eqn{d} represents the number of spatial dimensions.
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

### ensures coords, X, and Y are entered as matrices ###
  if (!is.matrix(coords)) {
    stop("coords must be a matrix")
  }

  if (!is.numeric(coords)) {
    stop("coords must be numeric")
  }


### ensures no missing data ###
  if (any(!is.finite(coords)) ||
      any(!is.finite(X)) ||
      any(!is.finite(Y)) ||
      any(!is.finite(phi_K))) {
    stop("All inputs must be non missing (no NA, NaN, or Inf values)") # ?? specify which one has missing
  }

### ensures k is a single positive integer ###
  # numeric is important to ensure numbers as characters aren't entered
  if (!is.numeric(K) || length(K) != 1 || K %% 1 != 0 || K <= 0) {
    stop("K must be a single integer")
  }

  ### making sure coords, X, and Y have at least one row and one column
  if (nrow(coords) == 0 || ncol(coords) == 0) {
    stop("coords must have at least one row and one column")
  }

  if (nrow(X) == 0 || ncol(X) == 0) {
    stop("X must have at least one row and one column")
  }

  if (nrow(Y) == 0 || ncol(Y) == 0) {
    stop("Y must have at least one row and one column")
  }


### row count consistency across inputs ###
  # ensure coords, X, Y have consistent row counts
  if (nrow(coords) != nrow(X) || nrow(coords) != nrow(Y)) {
    stop("coords, X, and Y must have the same number of rows")
  }

### consistency with phi_K and K ###
  if (!is.numeric(phi_K) || !is.atomic(phi_K) || is.matrix(phi_K) || length(phi_K) != K) {
    stop("phi_K must be a numeric vector of length K")
  }

### m and N_sam are entered correctly ###
  if (!is.numeric(m) || length(m) != 1 || m %% 1 != 0 || m <= 0) {
    stop("m must be a single positive integer")
  }

  if (!is.numeric(N_sam) || length(N_sam) != 1 || N_sam %% 1 != 0 || N_sam <= 0) {
    stop("N_sam must be a single positive integer")
  }

  # Source Julia code (safe to call multiple times)
  JuliaCall::julia_setup()
  JuliaCall::julia_source(jl_file)
  samples<-JuliaCall::julia_call("wrapper",coords,X,Y,phi_K,as.integer(K),m=as.integer(m), N_sam=as.integer(N_sam))

  samples <- list(
    gamma_samples = samples$γ_samples,
    sigma_samples = samples$Σ_samples,
    F_samples     = samples$F_samples
  )

  return(samples)

}

