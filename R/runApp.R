#' Run the Tessera Shiny app
#'
#' Launches a webview containing the embryo model
#'
#' @export
runTessera<- function(...) {
  appDir <- system.file("shiny-examples", "tessera", package = "tessera")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `tessera`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", ...)
}
