# Server and UI for ranking module

#' Create the ranking server
#'
#' @import shiny
#' @import shinythemes
#' @rawNamespace import(plotly, except = "last_plot")
#' @import ggplot2
#' @return the server
rankingServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Create and store seeds for embryos in the ranking pool only when button is clicked
    pool.embryo.seeds <- shiny::eventReactive(input$rank.embryos, {
      sample.int(.Machine$integer.max, size = input$n.pool)
    })


    # Create a pool of embryos for ranking
    calculateRanks <- shiny::reactive({
      aneuploidies <- runif(input$n.pool, min = 0, max = 1)

      # Create an embryo with the given aneuploidy and seed, and take one random biopsy
      embryos <- mapply(function(a, s) {
        tessera::Embryo(n.cells = 200, n.chrs = 1, prop.aneuploid = a, dispersal = input$dispersal, rng.seed = s)
      }, a = aneuploidies, s = pool.embryo.seeds())

      biopsies <- sapply(embryos, tessera::takeBiopsy, biopsy.size = input$biopsy.size)

      # Rank the embryos by aneuploidy
      best.from.biopsy <- which(biopsies %in% sort(biopsies)[1:input$n.transfer])
      best.from.reality <- which(aneuploidies %in% sort(aneuploidies)[1:input$n.transfer])


      is.real.best <- aneuploidies %in% sort(aneuploidies)[1:input$n.transfer]
      is.biop.best <- biopsies %in% sort(biopsies)[1:input$n.transfer]

      # There may be ties in ranks; in this case, increase the number of transferred embryos
      is.transferrable <- sum(is.biop.best)
      pct.correct <- sum(is.real.best & is.biop.best) / is.transferrable * 100

      # Calculate the percent of embryos correctly selected for transfer by biopsy
      # pct.corrrect <- sum(best.from.biopsy %in% best.from.reality) / input$n.transfer * 100

      list(
        "embryos" = embryos,
        "is.real.best" = is.real.best,
        "is.biop.best" = is.biop.best,
        "best.from.reality" = best.from.reality,
        "best.from.biopsy" = best.from.biopsy,
        "aneuploidies" = aneuploidies,
        "biopsies" = biopsies,
        "pct.correct" = pct.correct
      )
    })

    # Create text for rank results
    output$pool.data <- shiny::renderText(paste(calculateRanks()[["pct.correct"]]))

    # Render the ranking output plot
    output$pool.aneuplodies <- shiny::renderPlot({
      result <- calculateRanks()

      pct.correct <- round(result[["pct.correct"]], digits = 2)

      best.embryo.cutoff <- max(result[["aneuploidies"]][result[["best.from.reality"]]]) + 0.005

      best.biopsy.cutoff <- max(result[["biopsies"]][result[["best.from.biopsy"]]]) + 0.05

      colours <- dplyr::case_when(
        result[["is.biop.best"]] & result[["is.real.best"]] ~ "Correctly transferred",
        result[["is.biop.best"]] ~ "Wrongly transferred",
        result[["is.real.best"]] ~ "Wrongly rejected",
        T ~ "Correctly rejected"
      )

      ggplot2::ggplot() +
        annotate("rect", xmin = -Inf, xmax = best.embryo.cutoff, ymin = -Inf, ymax = best.biopsy.cutoff, fill = "darkgreen", alpha = 0.3) +
        annotate("rect", xmin = -Inf, xmax = best.embryo.cutoff, ymin = best.biopsy.cutoff, ymax = Inf, fill = "orange", alpha = 0.3) +
        annotate("rect", xmin = best.embryo.cutoff, xmax = Inf, ymin = -Inf, ymax = best.biopsy.cutoff, fill = "red", alpha = 0.3) +
        geom_point(aes(
          x = result[["aneuploidies"]],
          y = result[["biopsies"]],
          col = colours
        ),
        size = 3
        ) +
        geom_vline(xintercept = best.embryo.cutoff, col = "black") +
        geom_hline(yintercept = best.biopsy.cutoff, col = "black") +
        geom_text(aes(x = 0, y = input$biopsy.size + 1, label = "Wrongly rejected", hjust = 0)) +
        geom_text(aes(x = 1, y = -1, label = "Wrongly transferred", hjust = 1)) +
        geom_text(aes(x = 1, y = input$biopsy.size + 1, label = "Correctly rejected", hjust = 1)) +
        geom_text(aes(x = 0, y = -1, label = "Correctly transferred", hjust = 0)) +
        labs(
          x = "True embryo aneuploidy", y = "Aneuploid cells in biopsy",
          col = "Biopsy selected correcly",
          title = paste0(pct.correct, "% of the best embryos would be selected for transfer")
        ) +
        scale_y_continuous(breaks = seq(0, input$biopsy.size, 1)) +
        scale_colour_manual(values = c(
          "Correctly rejected" = "darkgrey",
          "Correctly transferred" = "darkgreen",
          "Wrongly rejected" = "orange",
          "Wrongly transferred" = "red"
        )) +
        theme_bw() +
        theme(legend.position = "none")
    })
  })
}

#' Create the ranking UI
#'
#' @return the ui
#'
rankingUI <- function(id) {
  sidebarLayout(
    sidebarPanel(
      # Number of embryos in the pool
      numericInput(
        inputId = NS(id, "n.pool"),
        label = strong("Embryos in pool"),
        value = 6,
        min = 2,
        max = 10,
        step = 1
      ),
      numericInput(
        inputId = NS(id, "n.transfer"),
        label = strong("Embryos to select"),
        value = 3,
        min = 1,
        max = 10,
        step = 1
      ),
      numericInput(
        inputId = NS(id, "biopsy.size"),
        label = strong("Cells per biopsy"),
        value = 5,
        min = 1,
        max = 20,
        step = 1
      ),
      numericInput(
        inputId = NS(id, "dispersal"),
        label = strong("Dispersal of aneuploid cells (0-1)"),
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.05
      ),
      actionButton(inputId = NS(id, "rank.embryos"), "Rank by biopsy")
    ),
    mainPanel(
      p("Generate a pool of embryos with random aneuploidies, and take a single biopsy from each. The biopsies to rank the embryos from best to worst, and select embryos for transfer. The chart shows the relationship between the embryos and the biopsies. Biopsies that were able to select the correct embryos are shown in green."),
      plotOutput(outputId = NS(id, "pool.aneuplodies"))
    )
  )
}
