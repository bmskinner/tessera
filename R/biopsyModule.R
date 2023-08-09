# Server and UI for biopsy module

#' Create the biopsy server
#'
#' @param id the namespace id
#' @import shiny
#' @import shinythemes
#' @rawNamespace import(plotly, except = "last_plot")
#' @import ggplot2
#' @return the server
biopsyServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Fixed definitions
    pgdis.classes <- factor(c("Euploid", "Low level", "High level", "Aneuploid"),
      levels = c("Euploid", "Low level", "High level", "Aneuploid")
    )

    to.pgdis.class <- function(prop.aneuploidy) {
      if (prop.aneuploidy < 0.20) {
        return(pgdis.classes[1])
      }
      if (prop.aneuploidy < 0.40) {
        return(pgdis.classes[2])
      }
      if (prop.aneuploidy <= 0.80) {
        return(pgdis.classes[3])
      }
      return(pgdis.classes[4])
    }


    # Create and store seed value for model embryo only when button is clicked
    seedVals <- shiny::eventReactive(input$new.embryo, {
      sample.int(.Machine$integer.max, size = 1)
    })

    # Generate an embryo with the desired parameters and take all possible biopsies
    calculateData <- reactive({
      all.chr <- T

      props <- if (all.chr) rep(input$proportion, times = 23) else input$proportion
      disps <- if (all.chr) rep(input$dispersal, times = 23) else input$dispersal

      embryo <- Embryo(
        n.cells = input$n.cells,
        n.chrs = 1,
        prop.aneuploid = props,
        dispersal = disps,
        concordance = 1,
        rng.seed = seedVals()
      )

      result <- takeAllBiopsies(embryo,
        biopsy.size = input$n.samples,
        chromosome = 1
      ) # input$chr.to.view

      biopsy.classes <- sapply(result, function(x) to.pgdis.class(x / input$n.samples))

      return(list(
        "embryo" = embryo,
        "biopsies" = result,
        "classes" = biopsy.classes,
        "true.class" = to.pgdis.class(input$proportion)
      ))
    })



    # Plot the embryo
    output$embryo.model <- plotly::renderPlotly({
      embryo <- calculateData()[["embryo"]]
      plot(embryo)
    })

    # Create plot of accuracies
    output$biopsy.accuracy <- plotly::renderPlotly({
      true.class <- calculateData()[["true.class"]]

      # Ensure zero values are still on the chart
      zero.data <- data.frame("Class" = pgdis.classes, "Count" = 0)

      classes <- data.frame("Class" = calculateData()[["classes"]]) %>%
        dplyr::group_by(.data$Class) %>%
        dplyr::summarise(Count = dplyr::n()) %>%
        dplyr::bind_rows(.data$., zero.data) %>%
        dplyr::group_by(.data$Class) %>%
        dplyr::summarise(Count = sum(.data$Count)) %>%
        dplyr::mutate(Pct = (.data$Count / input$n.cells) * 100, Colour = .data$Class == true.class)

      # print(classes)
      p <- plot_ly(
        x = classes$Class,
        y = classes$Pct,
        type = "bar",
        color = classes$Colour,
        colors = c("#555555", "#009806")
      ) %>%
        layout(
          showlegend = F,
          xaxis = list(title = "Biopsy class (embryo's class highlighted)"),
          yaxis = list(
            title = "Percentage of biopsies",
            range = c(0, 100)
          )
        )
      return(p)
    })

    # Make the histogram
    output$biopsy.histogram <- plotly::renderPlotly({
      result <- calculateData()[["biopsies"]]


      n.euploids <- length(result[result == 0])
      n.aneuploids <- length(result[result == input$n.samples])
      euploid.ratio <- (n.euploids / length(result)) * 100
      aneuploid.ratio <- (n.aneuploids / length(result)) * 100
      mosaic.ratio <- (length(result) - n.euploids - n.aneuploids) / length(result) * 100

      result <- data.frame(values = result)
      result$colour <- factor(ifelse(result$values == 0, "#009806",
        ifelse(result$values == input$n.samples, "red", "orange")
      ),
      levels = c("#009806", "orange", "red")
      )

      p <- plot_ly(
        x = result$values,
        type = "histogram",
        color = result$colour,
        colors = c("#009806", "orange", "red")
      ) %>%
        layout(
          showlegend = F,
          xaxis = list(
            title = "Number of aneuploid cells in biopsy",
            range = c(-0.5, input$n.samples + 0.5)
          ),
          yaxis = list(
            title = "Number of biopsies",
            range = c(0, input$n.cells)
          )
        ) %>%
        add_annotations(
          text = paste0(format(euploid.ratio, nsmall = 1, digits = 3), "%\neuploid"),
          x = 0,
          y = input$n.cells * 0.9,
          align = "center",
          showarrow = F
        ) %>%
        add_annotations(
          text = paste0(format(aneuploid.ratio, nsmall = 1, digits = 3), "%\naneuploid"),
          x = input$n.samples,
          y = input$n.cells * 0.9,
          align = "center",
          showarrow = F
        )

      # Only draw the line and display mosaic value if there are mosaic cells
      if (mosaic.ratio > 0) {
        p <- p %>%
          add_segments(
            x = 1,
            xend = input$n.samples - 1,
            y = input$n.cells * 0.8,
            yend = input$n.cells * 0.8,
            hoverinfo = "none",
            line = list(color = "grey", width = 0.5)
          ) %>%
          add_annotations(
            text = paste0(format(mosaic.ratio, nsmall = 1, digits = 3), "%\nmosaic"),
            x = input$n.samples / 2,
            y = input$n.cells * 0.9,
            align = "center",
            showarrow = F
          )
      }




      return(p)
    })
  })
}

#' Create the biopsy UI
#'
#' @param id the namespace id
#' @return the ui
#'
biopsyUI <- function(id) {
  # UI for the biopsy plots
  sidebarLayout(
    sidebarPanel(
      numericInput(
        inputId = NS(id, "n.cells"),
        label = strong("Total cells"),
        value = 200,
        min = 10,
        max = 300,
        step = 1
      ),
      numericInput(
        inputId = NS(id, "proportion"),
        label = strong("Proportion of aneuploid cells (0-1)"),
        value = 0.2,
        min = 0,
        max = 1,
        step = 0.05
      ),
      numericInput(
        inputId = NS(id, "dispersal"),
        label = strong("Dispersal of aneuploid cells (0-1)"),
        value = 0.1,
        min = 0,
        max = 1,
        step = 0.05
      ),
      numericInput(
        inputId = NS(id, "n.samples"),
        label = strong("Cells per biopsy"),
        value = 5,
        min = 1,
        max = 20,
        step = 1
      ),
      actionButton(inputId = NS(id, "new.embryo"), "New")
    ),
    mainPanel(
      p("Click ", strong("'New'")," to create a new embryo. The embryo will contain euploid (green) and aneuploid (grey) cells."),
      p("Every time you click", strong("'New'")," a different randomly generated embryo will be shown."),
      p("Use the settings on the left to adjust the proportion of aneuploid
                           cells in the embryo, and their dispersal (low dispersal means they
                           are found mostly in clusters, high dispersal means individual cells are more
                           likely). The blastocyst trophectoderm generated will be shown below."),
      p(""),

      # Embryo
      plotlyOutput(outputId = NS(id, "embryo.model"), width = 500, height = 400),
      p(""),
      p("The plot below shows the distribution of all possible biopsies for the embryo:"),
      p(""),

      # Biopsy results
      plotlyOutput(outputId = NS(id, "biopsy.histogram"), height = 300),
      p(""),
      p("Now we group the biopsies as euploid (0-20%), low level (20-40%), high-level (40-80%) or aneuploid (>80-100%). 
        The bar in green shows the actual class of the embryo:"),
      p(""),

      # Summary of biopsy accuracy
      plotlyOutput(outputId = NS(id, "biopsy.accuracy"), height = 300)
    )
  )
}
