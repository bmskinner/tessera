# Tests
library(testthat)

test_that("isAneuploid is working for all chromosomes", {

  for(cell in 1:100){ # test each cell in turn
    e = create.embryo(100, 0, 0) # no aneuploidy

    # chrs should not be aneuploid before setting
    for(i in 1:31){
      expect_equal(is.aneuploid(e, chromosome=i, cell.index=cell), F)
    }

    # set all chrs aneuploid
    for(i in 1:31){
      e = set.aneuploid(e, chromosome=i, cell.index=cell)
    }

    # all chrs should now be aneuploid
    for(i in 1:31){
      expect_equal(is.aneuploid(e, chromosome=i, cell.index=cell), T)
    }
  }
  # print(e$isAneuploid)
})

test_that("isAneuploid bit-masking works for all chrs", {
  e = create.embryo(100, 0, 0) # no aneuploidy


  # First check if the bit masking is working directly
  # Check setting aneuploidy of all chrs in cell 1
  for(i in 1:31){
    e = set.aneuploid(e, chromosome=i, cell.index=1)
  }

  for(i in 1:31){
    expect_equal(bitwAnd(e$isAneuploid[1], i), i)
  }

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
  e = create.embryo(100, c(0), c(0.1))
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

test_that("embryo created with correct chromosomes", {
  props = c(0, 0.1, 0.5, 0.4, 0.2)
  disps = c(0, 0, 0, 0, 0)

  e = create.embryo(100, props, disps)

  for(chr in 1:length(props)){
    exp = props[chr] * 100
    expect_equal(count.aneuploid(e, chr), exp)
  }
})
