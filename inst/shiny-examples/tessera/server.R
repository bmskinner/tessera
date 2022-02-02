# Control the server
library(shiny)
library(shinythemes)
library(plotly)

# Define the server
function(input, output, session){


  seedVals = eventReactive(input$new.embryo, {
    sample.int(.Machine$integer.max, size=1)
  })

  observeEvent(input$aneu.type, {
    if(input$aneu.type=="One chr"){
      updateNumericInput(session, "chr.to.view", value = 1)
    }
  })

  calculateData = reactive({

    all.chr = input$aneu.type == "All chrs"

    print(paste("Seed", seedVals()))

    props = if(all.chr) rep(input$proportion, times=23) else input$proportion
    disps = if(all.chr) rep(input$dispersal, times=23)  else input$dispersal

    Embryo(nCells = input$n.cells,
           nChrs   = 23,
           prop.aneuploid = props,
           dispersal = disps,
           concordance = input$concordance,
           rng.seed = seedVals())
  })

  output$biopsyPlot = renderPlotly({
    embryo = calculateData()

    plot(embryo)

  })

  output$iterationSummary = renderPlotly({

    embryo = calculateData()
    result = takeAllBiopsies(embryo, biopsy.size = input$n.samples,
                             chromosome = input$chr.to.view)

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
