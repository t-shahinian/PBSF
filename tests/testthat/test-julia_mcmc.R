## make sure everything works properly and output looks the way it should
## for underlying and wrapper functions
## can put a bunch into one test_that, or can separate
## i.e. dimension checking etc (ask gpt what it thinks)


## check that you get informative messages for misuse cases (i.e. they give
## decimal for integer, etc)
## just for the wrapper function


library(testthat)

# ---- Helper: create a minimal valid input set ----
  # use this to create correct input, then in the tests we purposely
  # change the value we are testing
make_valid_inputs <- function(n = 15, p = 2, q = 3, K = 2) {
  coords <- matrix(runif(n * 2), nrow = n, ncol = 2)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * q), nrow = n, ncol = q)
  phi_K <- rep(1.0, K)

  list(coords = coords, X = X, Y = Y, phi_K = phi_K, K = K)
}


# m < n (number of neighbors is less than sample size)
test_that("m must be less than n", {
  good <- make_valid_inputs(n = 5)

  expect_error(
    julia_mcmc(
      good$coords,
      good$X,
      good$Y,
      good$phi_K,
      good$K,
      m = 5
    ),
    "m must be less than"
  )
})

# ---- Test: coords must be a numeric 2-column matrix ----
test_that("coords must be a numeric 2-column matrix", {
  good <- make_valid_inputs()

  # Case 1: not a matrix (vector)
  bad1 <- good
  bad1$coords <- runif(5)
  expect_error(
    julia_mcmc(bad1$coords, bad1$X, bad1$Y, bad1$phi_K, bad1$K),
    "coords must be a matrix"
  )

  # Case 2: data.frame instead of matrix
  bad2 <- good
  bad2$coords <- data.frame(x = runif(5))
  expect_error(
    julia_mcmc(bad2$coords, bad2$X, bad2$Y, bad2$phi_K, bad2$K),
    "coords must be a matrix"
  )

  # Case 3: wrong number of columns (1 column)
  bad3 <- good
  bad3$coords <- matrix(runif(nrow(good$coords)), ncol = 1)
  expect_error(
    julia_mcmc(bad3$coords, bad3$X, bad3$Y, bad3$phi_K, bad3$K),
    "coords must have exactly 2 columns"
  )

  # Case 4: wrong number of columns (3 columns)
  bad4 <- good
  bad4$coords <- matrix(runif(nrow(good$coords) * 3),
                        nrow = nrow(good$coords),
                        ncol = 3)
  expect_error(
    julia_mcmc(bad4$coords, bad4$X, bad4$Y, bad4$phi_K, bad4$K),
    "coords must have exactly 2 columns"
  )
})


# ---- Test: inputs must not contain NA/NaN/Inf ----
test_that("all inputs must be fully observed", {
  good <- make_valid_inputs()

  # Case 1: NA in phi_K
  bad1 <- good
  bad1$phi_K[1] <- NA
  expect_error(
    julia_mcmc(bad1$coords, bad1$X, bad1$Y, bad1$phi_K, bad1$K),
    "All inputs must be non missing"
  )

  # Case 2: NA in coords
  bad2 <- good
  bad2$coords[1, 1] <- NA
  expect_error(
    julia_mcmc(bad2$coords, bad2$X, bad2$Y, bad2$phi_K, bad2$K),
    "All inputs must be non missing"
  )

  # Case 3: Inf in X
  bad3 <- good
  bad3$X[1, 1] <- Inf
  expect_error(
    julia_mcmc(bad3$coords, bad3$X, bad3$Y, bad3$phi_K, bad3$K),
    "All inputs must be non missing"
  )

  # Case 4: NaN in Y
  bad4 <- good
  bad4$Y[1, 1] <- NaN
  expect_error(
    julia_mcmc(bad4$coords, bad4$X, bad4$Y, bad4$phi_K, bad4$K),
    "All inputs must be non missing"
  )
})


# ---- Test 5: K must be a single positive integer ----
test_that("K must be a single positive integer", {
  good <- make_valid_inputs()

  # non-numeric
  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = "3"
    ),
    "K must be a single integer"
  )

  # length > 1
  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = c(1, 2)
    ),
    "K must be a single integer"
  )

  # non-integer numeric
  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = 1.5
    ),
    "K must be a single integer"
  )

  # non-positive
  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = 0
    ),
    "K must be a single integer"
  )
})

# ---- Test: X and Y must be non-empty matrices ----
test_that("X and Y must be non-empty matrices", {
  good <- make_valid_inputs()

  # X: zero columns
  bad1 <- good
  bad1$X <- matrix(numeric(0), nrow = nrow(good$coords), ncol = 0)
  expect_error(
    julia_mcmc(bad1$coords, bad1$X, bad1$Y, bad1$phi_K, bad1$K),
    "X must have at least one row and one column"
  )

  # X: zero rows
  bad2 <- good
  bad2$X <- matrix(numeric(0), nrow = 0, ncol = 1)
  expect_error(
    julia_mcmc(bad2$coords, bad2$X, bad2$Y, bad2$phi_K, bad2$K),
    "X must have at least one row and one column"
  )

  # Y: zero columns
  bad3 <- good
  bad3$Y <- matrix(numeric(0), nrow = nrow(good$coords), ncol = 0)
  expect_error(
    julia_mcmc(bad3$coords, bad3$X, bad3$Y, bad3$phi_K, bad3$K),
    "Y must have at least one row and one column"
  )

  # Y: zero rows
  bad4 <- good
  bad4$Y <- matrix(numeric(0), nrow = 0, ncol = 1)
  expect_error(
    julia_mcmc(bad4$coords, bad4$X, bad4$Y, bad4$phi_K, bad4$K),
    "Y must have at least one row and one column"
  )
})

# ---- Test 7: row counts must be consistent ----
test_that("coords, X, and Y must have the same number of rows", {
  good <- make_valid_inputs()

  # Make X have a different number of rows
  bad <- good
  bad$X <- matrix(
    rnorm((nrow(good$coords) - 1) * ncol(good$X)),
    nrow = nrow(good$coords) - 1,
    ncol = ncol(good$X)
  )

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
    "coords, X, and Y must have the same number of rows"
  )
})

# ---- Test: phi_K must be a numeric vector of length K ----
test_that("phi_K must be a numeric vector of length K", {
  good <- make_valid_inputs()

  # Case 1: wrong length
  bad1 <- good
  bad1$phi_K <- rep(1, good$K - 1)

  expect_error(
    julia_mcmc(bad1$coords, bad1$X, bad1$Y, bad1$phi_K, bad1$K),
    "phi_K must be a numeric vector of length K"
  )

  # Case 2: non-numeric
  bad2 <- good
  bad2$phi_K <- "a"

  expect_error(
    julia_mcmc(bad2$coords, bad2$X, bad2$Y, bad2$phi_K, bad2$K),
    "phi_K must be a numeric vector of length K"
  )

  # Case 3: matrix instead of vector
  bad3 <- good
  bad3$phi_K <- matrix(1, nrow = good$K, ncol = 1)

  expect_error(
    julia_mcmc(bad3$coords, bad3$X, bad3$Y, bad3$phi_K, bad3$K),
    "phi_K must be a numeric vector of length K"
  )
})

# ---- Test: m and N_sam must be single positive integers ----
test_that("m and N_sam must be single positive integers", {
  good <- make_valid_inputs()

  # m negative
  expect_error(
    julia_mcmc(good$coords, good$X, good$Y, good$phi_K, good$K, m = -1),
    "m must be a single positive integer"
  )

  # m non-integer
  expect_error(
    julia_mcmc(good$coords, good$X, good$Y, good$phi_K, good$K, m = 1.5),
    "m must be a single positive integer"
  )

  # N_sam zero
  expect_error(
    julia_mcmc(good$coords, good$X, good$Y, good$phi_K, good$K, N_sam = 0),
    "N_sam must be a single positive integer"
  )

  # N_sam non-integer
  expect_error(
    julia_mcmc(good$coords, good$X, good$Y, good$phi_K, good$K, N_sam = 2.3),
    "N_sam must be a single positive integer"
  )
})



# ---- Test: julia_mcmc output structure and dimensions ----
test_that("julia_mcmc returns correctly structured MCMC output", {
  skip_on_cran()

  good <- make_valid_inputs(n = 5, p = 2, q = 3, K = 2)
  N_sam <- 4

  expect_no_error(
    res <- julia_mcmc(
      good$coords,
      good$X,
      good$Y,
      good$phi_K,
      good$K,
      m = 2,
      N_sam = N_sam
    )
  )

  # Structure
  expect_type(res, "list")
  expect_named(res, c("gamma_samples", "sigma_samples", "F_samples"))
  expect_true(is.matrix(res$gamma_samples))
  expect_true(is.matrix(res$sigma_samples))
  expect_true(is.matrix(res$F_samples))

  # Dimensions
  n <- nrow(good$coords)
  p <- ncol(good$X)
  q <- ncol(good$Y)
  K <- good$K

  expect_equal(nrow(res$gamma_samples), N_sam)
  expect_equal(nrow(res$sigma_samples), N_sam)
  expect_equal(nrow(res$F_samples), N_sam)

  expect_equal(ncol(res$gamma_samples), (p + K) * q)
  expect_equal(ncol(res$sigma_samples), q)
  expect_equal(ncol(res$F_samples), n * K)
})



