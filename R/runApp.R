# Run the Shiny app

#' @export
runModel<- function() {
  appDir <- system.file("shiny-examples", "tessera", package = "tessera")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `tessera`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
