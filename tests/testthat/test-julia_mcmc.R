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
# ---- Test 1: coords must be a matrix ----
test_that("coords must be a matrix", {
  bad <- make_valid_inputs()
  bad$coords <- data.frame(x = runif(5))   # not a matrix

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
    "coords must be a matrix"
  )
})

# ---- Test 2: coords must not be a vector ----
test_that("coords cannot be a vector", {
  bad <- make_valid_inputs()
  bad$coords <- runif(5)   # numeric vector

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
    "coords must be a matrix"
  )
})

# ---- Test 3: coords must be numeric ----
test_that("coords must be numeric", {
  bad <- make_valid_inputs()
  bad$coords <- matrix(letters[1:5], ncol = 1)

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
    "coords must be numeric"
  )
})

# ---- coords must be 2 dimensional ---- #
test_that("coords must have exactly 2 columns", {
  good <- make_valid_inputs()

  # Case 1: only 1 column
  bad1 <- good
  bad1$coords <- matrix(runif(nrow(good$coords)), ncol = 1)

  expect_error(
    julia_mcmc(
      coords = bad1$coords,
      X = bad1$X,
      Y = bad1$Y,
      phi_K = bad1$phi_K,
      K = bad1$K
    ),
    "coords must have exactly 2 columns"
  )

  # Case 2: 3 columns
  bad2 <- good
  bad2$coords <- matrix(runif(nrow(good$coords) * 3),
                        nrow = nrow(good$coords),
                        ncol = 3)

  expect_error(
    julia_mcmc(
      coords = bad2$coords,
      X = bad2$X,
      Y = bad2$Y,
      phi_K = bad2$phi_K,
      K = bad2$K
    ),
    "coords must have exactly 2 columns"
  )
})

# ---- Test 4: missing or non-finite values are rejected ----
test_that("all inputs must be fully observed", {
  bad <- make_valid_inputs()
  bad$phi_K[1] <- NA   # break exactly one rule

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
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

# ---- Test 6: X and Y must have at least one row and one column ----
test_that("X, and Y must be non-empty matrices", {
  good <- make_valid_inputs()

  # X: zero columns
  bad2 <- good
  bad2$X <- matrix(numeric(0), nrow = nrow(good$coords), ncol = 0)
  expect_error(
    julia_mcmc(
      coords = bad2$coords,
      X = bad2$X,
      Y = bad2$Y,
      phi_K = bad2$phi_K,
      K = bad2$K
    ),
    "X must have at least one row and one column"
  )

  # Y: zero rows
  bad3 <- good
  bad3$Y <- matrix(numeric(0), nrow = 0, ncol = 1)
  expect_error(
    julia_mcmc(
      coords = bad3$coords,
      X = bad3$X,
      Y = bad3$Y,
      phi_K = bad3$phi_K,
      K = bad3$K
    ),
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

# ---- Test 8: phi_K length must equal K ----
test_that("phi_K must be a numeric vector of length K", {
  good <- make_valid_inputs()

  # wrong length
  bad <- good
  bad$phi_K <- rep(1, bad$K - 1)

  expect_error(
    julia_mcmc(
      coords = bad$coords,
      X = bad$X,
      Y = bad$Y,
      phi_K = bad$phi_K,
      K = bad$K
    ),
    "phi_K must be a numeric vector of length K"
  )
})

# ---- Test 8: m must be a single positive integer ---
test_that("m must be a single positive integer", {
  good <- make_valid_inputs()

  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = good$K,
      m = -1
    ),
    "m must be a single positive integer"
  )
})

# ---- Test 10: N_sam must be a single positive integer ----
test_that("N_sam must be a single positive integer", {
  good <- make_valid_inputs()

  expect_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = good$K,
      N_sam = 0
    ),
    "N_sam must be a single positive integer"
  )
})

######## ???? everything works except these two tests. fix later

# ---- Test 11: julia_mcmc returns correctly structured output ----
test_that("julia_mcmc returns structured MCMC output", {
  skip_on_cran()

  good <- make_valid_inputs()
  res <- julia_mcmc(
    coords = good$coords,
    X = good$X,
    Y = good$Y,
    phi_K = good$phi_K,
    K = good$K,
    m = 2,
    N_sam = 3
  )

  expect_type(res, "list")
  expect_named(res, c("gamma_samples", "sigma_samples", "F_samples"))

  expect_true(is.matrix(res$gamma_samples))
  expect_true(is.matrix(res$sigma_samples))
  expect_true(is.matrix(res$F_samples))
})


### ---- Test 12: MCMC output dimensions are correct ----
test_that("MCMC output dimensions match model specification", {
  skip_on_cran()

  good <- make_valid_inputs(n = 5, p = 2, q = 3, K = 2)
  N_sam <- 4

  res <- julia_mcmc(
    coords = good$coords,
    X = good$X,
    Y = good$Y,
    phi_K = good$phi_K,
    K = good$K,
    m = 2,
    N_sam = N_sam
  )

  # extract dimensions from inputs
  n <- nrow(good$coords)
  p <- ncol(good$X)
  q <- ncol(good$Y)
  K <- good$K

  # rows = number of MCMC samples
  expect_equal(nrow(res$gamma_samples), N_sam)
  expect_equal(nrow(res$sigma_samples), N_sam)
  expect_equal(nrow(res$F_samples), N_sam)

  # columns follow model algebra
  expect_equal(ncol(res$gamma_samples), (p + K) * q)
  expect_equal(ncol(res$sigma_samples), q)
  expect_equal(ncol(res$F_samples), n * K)
})

# ---- Test 13: minimal valid call runs without error ----
test_that("minimal valid example runs without error", {
  skip_on_cran()

  good <- make_valid_inputs()

  expect_no_error(
    julia_mcmc(
      coords = good$coords,
      X = good$X,
      Y = good$Y,
      phi_K = good$phi_K,
      K = good$K,
      m = 2,
      N_sam = 2
    )
  )
})

# ---- Test 14: returned object has exactly expected names ----
test_that("returned object has exactly expected components", {
  skip_on_cran()

  good <- make_valid_inputs()

  res <- julia_mcmc(
    coords = good$coords,
    X = good$X,
    Y = good$Y,
    phi_K = good$phi_K,
    K = good$K,
    m = 2,
    N_sam = 2
  )

  expect_setequal(
    names(res),
    c("gamma_samples", "sigma_samples", "F_samples")
  )
})
