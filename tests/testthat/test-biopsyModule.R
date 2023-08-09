#' #' @import testthat
#' test_that("reactives and output updates", {
#'   shiny::testServer(tesseraServer, {
#'     session$setInputs(
#'       n.cells = 100, proportion = 0.1, dispersal = 0.5, n.samples = 4, biopsy.size = 4,
#'       n.transfer = 3, n.pool = 7, new.embryo = 1
#'     )
#'     print(calculateRanks())
#'     print(output$pool.data)
#'     
#'     # expect_equal(output$pool.data, 100)
#'   })
#' })
