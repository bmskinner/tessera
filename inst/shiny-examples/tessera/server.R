# Control the server
library(shiny)
library(shinythemes)
library(plotly)

# Define the server
function(input, output, session){


  seedVals = eventReactive(input$new.embryo, {
    sample.int(.Machine$integer.max, size=1)
  })

  calculateData = reactive({

    all.chr = input$aneu.type=="All chrs"

    print(paste("Seed", seedVals()))

    props = if(all.chr) rep(input$proportion, times=22) else input$proportion
    disps = if(all.chr) rep(input$dispersal, times=22)  else input$dispersal

    create.embryo(n.cells        = input$n.cells,
                  prop.aneuploid = props,
                  dispersion     = disps,
                  seed           = seedVals())
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

  output$iterationSummary = renderPlotly({

    embryo = calculateData()
    result = take.all.biopsies(embryo, input$n.samples, input$chr.to.view)

    n.euploids      = length(result[result==0])
    n.aneuploids    = length(result[result==input$n.samples])
    euploid.ratio   = (n.euploids / length(result))*100
    aneuploid.ratio = (n.aneuploids / length(result))*100
    mosaic.ratio    = (length(result)-n.euploids-n.aneuploids)/length(result)*100

    result = data.frame(values=result)
    result$colour = factor(ifelse(result$values==0, "green",
                           ifelse(result$values==input$n.samples, "red", "orange")),
                           levels = c("green", "orange", "red"))

    p = plot_ly(x = result$values,
            type="histogram",
            color = result$colour,
            colors = c("green", "orange", "red")) %>%
      layout(showlegend = F,
             xaxis = list( title = "Number of aneuploid cells in biopsy",
                           range=c(-0.5,input$n.samples+0.5)),
             yaxis = list( title = "Number of biopsies",
                           range=c(0, input$n.cells))) %>%
      add_annotations(text=paste0(format(euploid.ratio, nsmall=1, digits = 3),"%\nall euploid"),
                      x=0,
                      y=input$n.cells*0.9,
                      align="center",
                      showarrow=F) %>%
      add_annotations(text=paste0(format(aneuploid.ratio, nsmall=1, digits = 3),"%\nall aneuploid"),
                      x=input$n.samples,
                      y=input$n.cells*0.9,
                      align="center",
                      showarrow=F)

    # Only draw the line and display mosaic value if there are mosaic cells
    if(mosaic.ratio>0){

      p = p %>% add_segments(x=1,
                             xend=input$n.samples-1,
                             y=input$n.cells*0.8,
                             yend=input$n.cells*0.8,
                             hoverinfo="none",
                             line = list(color="grey", width=0.5)) %>%
        add_annotations(text=paste0(format(mosaic.ratio, nsmall=1, digits = 3),"%\nmosaic"),
                        x=input$n.samples/2,
                        y=input$n.cells*0.9,
                        align="center",
                        showarrow=F)
    }



    return(p)
  })
}
