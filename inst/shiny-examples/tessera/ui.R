# Control the UI
library(shiny)
library(shinythemes)
library(plotly)

# Define the UI
fluidPage(theme = shinytheme("lumen"),
            titlePanel("Tessera: Aneuploidies seen in embryo biopsies"),
            sidebarLayout(
              sidebarPanel(

                numericInput(inputId = "n.cells",
                             label = strong("Total cells"),
                             value = 100,
                             min=10,
                             max = 300,
                             step = 1),

                radioButtons(inputId = "aneu.type",
                             label = "Type:",
                             choices = c("All chrs", "Per chr"),
                             selected = "All chrs"),

                numericInput(inputId = "proportion",
                             label = strong("Proportion of aneuploid cells (0-1)"),
                             value = 0.1,
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

                numericInput(inputId = "chr.to.view",
                             label = strong("Chromosome to view"),
                             value = 1,
                             min= 0,
                             max = 31,
                             step = 1)


              ),
              mainPanel(
                p("The embryo below contains euploid (green) and aneuploid (red) cells.
                           Use the settings on the left to adjust the proportion of aneuploid
                           cells in the data, as well as their dispersal (low dispersal means they
                           are found mostly in clumps, high dispersal means individual cells are more
                           likely). The blastocyst generated will be shown below, and a histogram
                  showing all possible biopsies of a given size for that embryo"),
                p(""),

                plotlyOutput("biopsyPlot",  width = 500, height = 400),
                plotOutput(outputId = "iterationSummary", height = "300px")
              ))
)
