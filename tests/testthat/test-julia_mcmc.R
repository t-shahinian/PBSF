## make sure everything works properly and output looks the way it should
## for underlying and wrapper functions
## can put a bunch into one test_that, or can separate
## i.e. dimension checking etc (ask gpt what it thinks)
test_that("multiplication works", {


  expect_equal(2 * 2, 4)




})





## check thaty ou get informative messages for misuse cases (i.e. they give
## decimal for integer, etc)
## just for the wrapper function
