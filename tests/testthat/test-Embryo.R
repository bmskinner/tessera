# Tests

test_that("Embryo does not have aneuploid chromosomes when not specified", {
  
  # Assume 31 chromosome pairs. None should be aneuuploid
  e <- Embryo(n.cells = 100, n.chrs = 31, prop.aneuploid = 0, dispersal = 0)
  for (i in 1:31) {
    expect_equal(tessera::countAneuploidCells(e, i), 0)
  }
  expect_equal(tessera::countAneuploidCells(e), 0)
})

test_that("Embryo has aneuploid specified chromosomes", {
  
  # Set each chromosome in turn
  for (i in 1:31) {
    if (i == 1) {
      chrs <- c(1, rep(0, 30))
    }
    if (i > 1 && i < 31) {
      chrs <- c(rep(0, i - 1), 1, rep(0, 31 - i))
    }
    if (i == 31) {
      chrs <- c(rep(0, 30), 1)
    }
    
    e <- Embryo(n.cells = 100, n.chrs = 31, prop.aneuploid = chrs, dispersal = 0)
    expect_equal(tessera::countAneuploidCells(e, i), 100)
  }
})

test_that("Number of aneuploid cells matches input aneuploidy level", {
  for (i in seq(0, 1, 0.05)) {
    e <- Embryo(n.cells = 100, n.chrs = 1, prop.aneuploid = i, dispersal = 0)
    expect_equal(tessera::countAneuploidCells(e), i * 100)
  }
})

test_that("RNG seed makes reproducible embryos", {
  
  # Embryos with the same seed should be strictly identical
  e <- Embryo(n.cells = 100, n.chrs = 1, prop.aneuploid = 0.5, dispersal = 0, rng.seed = 42)
  e2 <- Embryo(n.cells = 100, n.chrs = 1, prop.aneuploid = 0.5, dispersal = 0, rng.seed = 42)
  expect_identical(e, e2)
  
  # Embryos with different seeds should be different
  e3 <- Embryo(n.cells = 100, n.chrs = 1, prop.aneuploid = 0.5, dispersal = 0, rng.seed = 43)
  expect_false(isTRUE(all.equal(e, e3)))
})

test_that("Clustered embryos have all aneuploid cells adjacent to another", {
  has.adjacent.aneuploid <- function(embryo, cell.index, chromosome) {
    adj.list <- embryo@neighbours[[paste0("n", cell.index)]]
    isAenu <- embryo@ploidy[, chromosome] != embryo@euploidy
    return(any(adj.list & isAenu))
  }
  
  for (aneu in seq(0.05, 0.95, 0.05)) {
    # clustered embryo
    e <- Embryo(n.cells = 100, n.chrs = 1, prop.aneuploid = aneu, dispersal = 0)
    
    # Check any aneuploid cells have an aneuploid neighbour
    for (i in 1:100) {
      if (e@ploidy$chr1[i] != e@euploidy) {
        expect_true(has.adjacent.aneuploid(e, i, 1))
      }
    }
  }
})
