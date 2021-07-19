# Control the server
library(shiny)
library(shinythemes)
library(plotly)
# source("tessera.R")

# Define the server
function(input, output, session){


  calculateData = reactive({
    create.embryo(n.cells        = input$n.cells,
                  prop.aneuploid = input$proportion,
                  dispersion     = input$dispersal)
  })

  output$biopsyPlot = renderPlotly({
    d = calculateData()
    plot_ly(x=d$x, y=d$y, z=d$z,
            type="scatter3d",
            mode="markers",
            color=d$isSeed,
            colors = c("#00FF00", "#FF0000")) %>%
      layout(showlegend = FALSE) %>%
      layout(title = "Click and drag to rotate the chart")
  })

  output$iterationSummary = renderPlot({

    d = calculateData()
    result = make.samples(d, input$n.samples)

    n.euploids = length(result[result==0])
    n.aneuploids  = length(result[result==input$n.samples])
    euploid.ratio = (n.euploids / length(result))*100
    aneuploid.ratio = (n.aneuploids / length(result))*100
    mosaic.ratio    = (length(result)-n.euploids-n.aneuploids)/length(result)*100

    # Plot the results
    hist(result, xlim=c(-0.5,input$n.samples+0.5),
         breaks=seq(-0.5, input$n.samples+0.5, 1),
         xlab = "Number of aneuploid cells in biopsy",
         ylab = "Fraction of biopsies",
         col = c("green", rep("orange", input$n.samples-1), "red"),
         xaxt="n",
         ylim = c(0,1),
         freq=F,
         main = paste("Biopsying",input$n.samples,
                      "cells from this blastocyst"))
    axis(1, at = seq(0, input$n.samples, 1))
    text( paste0(format(euploid.ratio, nsmall=1, digits = 3),"%\nall euploid"), x=0, y=0.9)
    text( paste0(format(mosaic.ratio, nsmall=1, digits = 3),"%\nmosaic"), x=input$n.samples/2, y=0.9)
    text( paste0(format(aneuploid.ratio, nsmall=1, digits = 3),"%\nall aneuploid"), x=input$n.samples, y=0.9)
    if(mosaic.ratio>0) segments(1, 0.7, input$n.samples-1, 0.7)
  })


}
