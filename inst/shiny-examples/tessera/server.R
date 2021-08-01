# Control the server
library(shiny)
library(shinythemes)
library(plotly)

# Define the server
function(input, output, session){

  calculateData = reactive({

    all.chr = input$aneu.type=="All chrs"

    props = if(all.chr) rep(input$proportion, times=22) else input$proportion
    disps = if(all.chr) rep(input$dispersal, times=22)  else input$dispersal

    create.embryo(n.cells        = input$n.cells,
                  prop.aneuploid = props,
                  dispersion     = disps)
  })

  output$biopsyPlot = renderPlotly({
    embryo = calculateData()

    show.legend = F
    if(input$chr.to.view==0){

      # Show number of aneuploid chromosomes in cell
      cell.list = c()
      for(cell in 1:nrow(embryo)){
        total = 0
        for(i in 1:31){
          if(bitwAnd(embryo$isAneuploid[cell], 2^(i-1))==2^(i-1)){
            total = total + 1
          }
        }
        cell.list = c(cell.list, total)
      }

      colours = factor(cell.list)
      show.legend = T

    } else {
      # just show the state of the chromosome of interest
      colours = factor(bitwAnd(embryo$isAneuploid, 2^(input$chr.to.view-1)),
                       levels = c(0, 2^(input$chr.to.view-1)))

    }

    plot_ly(x=embryo$x, y=embryo$y, z=embryo$z,
            type="scatter3d",
            mode="markers",
            color=colours,
            colors = c("#00FF00", "#FF0000")) %>%
      layout(showlegend = show.legend) %>%
      layout(title = "Click and drag to rotate the chart")
  })

  output$iterationSummary = renderPlot({

    embryo = calculateData()
    result = take.all.biopsies(embryo, input$n.samples, input$chr.to.view)

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
         main = paste("Chr", input$chr.to.view, ": Biopsying",input$n.samples,
                      "cells from this embryo"))
    axis(1, at = seq(0, input$n.samples, 1))
    text( paste0(format(euploid.ratio, nsmall=1, digits = 3),"%\nall euploid"), x=0, y=0.9)
    text( paste0(format(mosaic.ratio, nsmall=1, digits = 3),"%\nmosaic"), x=input$n.samples/2, y=0.9)
    text( paste0(format(aneuploid.ratio, nsmall=1, digits = 3),"%\nall aneuploid"), x=input$n.samples, y=0.9)
    if(mosaic.ratio>0) segments(1, 0.7, input$n.samples-1, 0.7)
  })


}
