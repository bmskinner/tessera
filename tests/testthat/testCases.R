# Tests
library(testthat)

test_that("bit masking is working", {
  e = create.embryo(100, 0, 0.1)

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


test_that("counting aneuploid chrs", {
  e = create.embryo(100, 0, 0.1)
  e = set.aneuploid(e, 1, 1)
  e = set.aneuploid(e, 2, 1)

  expect_equal(count.aneuploid(e, 1), 2)
  expect_equal(count.aneuploid(e, 2), 0)

  e = set.aneuploid(e, 3, 1)
  expect_equal(count.aneuploid(e, 1), 3)
})

test_that("taking one biopsy", {
  e = create.embryo(100, 0, 0.1)
  e = set.aneuploid(e, 1, 1)
  e = set.aneuploid(e, 2, 1)

  # The biopsy here should contain only 2 aneuploid cells for chr1
  # and no other aneuploidies
  c = take.one.biopsy(e, 5, 1, 1)
  expect_equal(c, 2)

  for(i in 2:31){
    expect_equal(take.one.biopsy(e,
                                 n.sampled.cells=5,
                                 index.cell=1,
                                 chromosome=i),
                 0)
  }

})
