# Control the UI
library(shiny)
library(shinythemes)
library(plotly)

# Define the UI
fluidPage(theme = shinytheme("lumen"),
            titlePanel("Tessera: Effect of embryo mosaicism on trophectoderm biopsies"),
            sidebarLayout(
              sidebarPanel(

                numericInput(inputId = "n.cells",
                             label = strong("Total cells"),
                             value = 200,
                             min=10,
                             max = 300,
                             step = 1),

                numericInput(inputId = "proportion",
                             label = strong("Proportion of aneuploid cells (0-1)"),
                             value = 0.2,
                             min=0,
                             max = 1,
                             step = 0.05),

                numericInput(inputId = "dispersal",
                             label = strong("Dispersal of aneuploid cells (0-1)"),
                             value = 0.1,
                             min=0,
                             max = 1,
                             step = 0.05),

                numericInput(inputId = "n.samples",
                             label = strong("Cells per biopsy"),
                             value = 5,
                             min= 1,
                             max = 20,
                             step = 1),

                # radioButtons(inputId = "aneu.type",
                #              label = strong("Model all chromosomes?"),
                #              choices = c("One chr", "All chrs"),
                #              selected = "One chr"),
                #
                # conditionalPanel(
                #   condition = "input['aneu.type'] == 'All chrs'",
                #
                #   numericInput(inputId = "concordance",
                #                label = strong("Concordance between chromosomes (0-1)"),
                #                value = 1,
                #                min= 0,
                #                max = 1,
                #                step = 0.01),
                #
                #   numericInput(inputId = "chr.to.view",
                #                label = strong("Chromosome to view (0 to see all)"),
                #                value = 1,
                #                min= 0,
                #                max = 23,
                #                step = 1)
                #   ),

                actionButton("new.embryo", "New")
              ),

              mainPanel(
                p("Click 'New' to create a new embryo. The embryo will contain euploid (green) and aneuploid (grey) cells."),
                p("Every time you click 'New' a different randomly generated embryo will be shown."),
                p("Use the settings on the left to adjust the proportion of aneuploid
                           cells in the embryo, and their dispersal (low dispersal means they
                           are found mostly in clusters, high dispersal means individual cells are more
                           likely). The blastocyst trophectoderm generated will be shown below with a histogram
                  showing all possible biopsies for the embryo."),
                p(""),

                plotlyOutput("biopsyPlot",  width = 500, height = 400),
                plotlyOutput(outputId = "iterationSummary", height = 300)
              ))
)
