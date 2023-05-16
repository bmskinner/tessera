# Control the UI
library(shiny)
library(shinythemes)
library(plotly)


# UI for the biopsy plots
modelUi <-  sidebarLayout(
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
    
    # Embryo
    plotlyOutput(outputId = "embryo.model",  width = 500, height = 400),
    
    # Biopsy results
    plotlyOutput(outputId = "biopsy.histogram", height = 300),
    
    # Summary of biopsy accuracy
    plotlyOutput(outputId = "biopsy.accuracy", height = 300)
    
    
  ))  

# UI for the ranking plots
rankingUi <-  sidebarLayout(
  sidebarPanel(
    # Number of embryos in the pool
               numericInput(inputId = "n.pool",
                            label = strong("Embryos in pool"),
                            value = 6,
                            min=2,
                            max = 10,
                            step = 1),
               
               numericInput(inputId = "n.transfer",
                            label = strong("Embryos to select"),
                            value = 3,
                            min=1,
                            max = 10,
                            step = 1),
               
               numericInput(inputId = "biopsy.size",
                            label = strong("Cells per biopsy"),
                            value = 5,
                            min= 1,
                            max = 20,
                            step = 1),      
               
               numericInput(inputId = "dispersal",
                            label = strong("Dispersal of aneuploid cells (0-1)"),
                            value = 0.1,
                            min=0,
                            max = 1,
                            step = 0.05),
               
               actionButton("rank.embryos", "Rank by biopsy")
               
               ),
  
  mainPanel(
    p("Generate a pool of embryos with random aneuploidies, and take a single biopsy from each. The biopsies to rank the embryos from best to worst, and select embryos for transfer. The chart shows the relationship between the embryos and the biopsies. Biopsies that were able to select the correct embryos are shown in green."),
    plotOutput("pool.aneuplodies")
    
  )
)

aboutUi <- mainPanel(
  p("This model visualises the possible tropectoderm biopsies obtained from a blastocyst in Preimplantation Genetic Testing for Aneuploidy (PGT-A) during IVF."),
  p("It is intended to help demonstrate the impact that the proportion and distribution of aneuploid cells can have on the biopsies that are obtained, and show how a single biopsy may not accurately reflect the embryo from which it came."),
  p("Tessera was created by", 
    a("Dr Ben Skinner", href="https://www.essex.ac.uk/people/skinn19306/benjamin-skinner"), 
    "at the",
    a("University of Essex", href="https://www.essex.ac.uk")
  ),
  
  p("Download the source code and run tessera locally from",
    a("https://github.com/bmskinner/tessera", href="https://github.com/bmskinner/tessera")
  )
)


# Define the overall UI
fluidPage(theme = shinytheme("lumen"),
          titlePanel("Tessera: Effect of embryo mosaicism on trophectoderm biopsies for PGT-A"),
          
          tabsetPanel(type = "tabs",
                      tabPanel("Model", modelUi),
                      tabPanel("Ranking",rankingUi),
                      tabPanel("About",aboutUi)
          )
)    

