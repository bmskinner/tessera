#' Run the Tessera Shiny app
#'
#' Launches a webview containing the embryo model
#'
#' @export
runTessera<- function(...) {
  shiny::shinyApp(ui = tesseraUI, server = tesseraServer)
}

