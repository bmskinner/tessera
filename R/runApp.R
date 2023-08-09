#' Run the Tessera Shiny app
#'
#' Launches a webview containing the embryo model
#' @import shinytheme
#' @export
runTessera <- function(...) {
  aboutUi <- mainPanel(
    p("This model visualises the possible tropectoderm biopsies obtained from a blastocyst in Preimplantation Genetic Testing for Aneuploidy (PGT-A) during IVF."),
    p("It is intended to help demonstrate the impact that the proportion and distribution of aneuploid cells can have on the biopsies that are obtained, and show how a single biopsy may not accurately reflect the embryo from which it came."),
    p(
      "Tessera was created by",
      a("Dr Ben Skinner", href = "https://www.essex.ac.uk/people/skinn19306/benjamin-skinner"),
      "at the",
      a("University of Essex", href = "https://www.essex.ac.uk")
    ),
    p(
      "Download the source code and run tessera locally from",
      a("https://github.com/bmskinner/tessera", href = "https://github.com/bmskinner/tessera")
    )
  )

  ui <- fluidPage(
    theme = shinythemes::shinytheme("lumen"),
    titlePanel("Tessera: Effect of embryo mosaicism on trophectoderm biopsies for PGT-A"),
    tabsetPanel(
      type = "tabs",
      tabPanel("Model", biopsyUI("biop1")),
      tabPanel("Ranks", rankingUI("rank1")),
      tabPanel("About", aboutUi)
    )
  )

  server <- function(input, output, session) {
    rankingServer("rank1")
    biopsyServer("biop1")
  }

  shiny::shinyApp(ui = ui, server = server)
}
