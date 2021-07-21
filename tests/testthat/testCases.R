# Tests
library(testthat)

test_that("bit masking is working", {
  e = create.embryo(100, 0.1, 0.1)

  # Check masking of chr 1
  e = set.aneuploid(e, 1, 1)
  expect_equal(bitwAnd(e$isAneuploid[1], 1), 1)

  # Check masking of chr 2
  e = set.aneuploid(e, 1, 2)
  expect_equal(bitwAnd(e$isAneuploid[1], 2), 2)


  expect_equal(is.aneuploid(e, 1, 1), T)
  expect_equal(is.aneuploid(e, 1, 2), T)
  expect_equal(is.aneuploid(e, 1, 3), F)
  expect_equal(is.aneuploid(e, 1, 4), F)
})
