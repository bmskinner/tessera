# Control the server
library(shiny)
library(shinythemes)
library(plotly)
library(tessera)

# Define the server
function(input, output, session){

  pgdis.classes = factor(c("Euploid", "Low level", "High level", "Aneuploid"),
                         levels = c("Euploid", "Low level", "High level", "Aneuploid"))

  to.pgdis.class = function(prop.aneuploidy){
    if(prop.aneuploidy<0.20) return(pgdis.classes[1])
    if(prop.aneuploidy<0.40) return(pgdis.classes[2])
    if(prop.aneuploidy<=0.80) return(pgdis.classes[3])
    return(pgdis.classes[4])
  }


  seedVals = eventReactive(input$new.embryo, {
    sample.int(.Machine$integer.max, size=1)
  })

  # observeEvent(input$aneu.type, {
  #   if(input$aneu.type=="One chr"){
  #     updateNumericInput(session, "chr.to.view", value = 1)
  #   }
  # })

  # Generate an embryo with the desired parameters and take all possible biopsies
  calculateData = reactive({

    # all.chr = input$aneu.type == "All chrs"
    all.chr = T

    # print(paste("Seed", seedVals()))

    props = if(all.chr) rep(input$proportion, times=23) else input$proportion
    disps = if(all.chr) rep(input$dispersal, times=23)  else input$dispersal

    embryo = Embryo(n.cells = input$n.cells,
           n.chrs   = 1,
           prop.aneuploid = props,
           dispersal = disps,
           concordance = 1,
           rng.seed = seedVals())

    result = takeAllBiopsies(embryo, biopsy.size = input$n.samples,
                             chromosome = 1) #input$chr.to.view

    biopsy.classes = sapply(result, function(x) to.pgdis.class(x/input$n.samples) )


    # embryo.class = to.pgdis.class(input$proportion)


    return(list("embryo"   = embryo,
                "biopsies" = result,
                "classes"  = biopsy.classes,
                "true.class"= to.pgdis.class(input$proportion)))
  })

  # Plot the embryo
  output$embryo.model = renderPlotly({
    embryo = calculateData()[['embryo']]
    plot(embryo)

  })

  output$biopsy.accuracy = renderPlotly({

    true.class = calculateData()[['true.class']]

    # Ensure zero values are still on the chart
    zero.data = data.frame("Class" = pgdis.classes, "Count" = 0)

    classes =  data.frame("Class" = calculateData()[['classes']]) %>%
      dplyr::group_by(Class) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::bind_rows(., zero.data) %>%
      dplyr::group_by(Class) %>%
      dplyr::summarise(Count = sum(Count)) %>%
      dplyr::mutate(Pct = (Count / input$n.cells)*100, Colour = Class==true.class)

    # print(classes)
    p = plot_ly(x = classes$Class,
                y = classes$Pct,
                type="bar",
                color = classes$Colour,
                colors = c("#555555", "#009806")) %>%
      layout(showlegend = F,
             xaxis = list( title = "Biopsy PGDIS class (embryo's class highlighted)"),
             yaxis = list( title = "Percentage of biopsies",
                           range=c(0, 100)))
    return(p)
  })

  # Make the histogram
  output$biopsy.histogram = renderPlotly({

    result = calculateData()[['biopsies']]


    n.euploids      = length(result[result==0])
    n.aneuploids    = length(result[result==input$n.samples])
    euploid.ratio   = (n.euploids / length(result))*100
    aneuploid.ratio = (n.aneuploids / length(result))*100
    mosaic.ratio    = (length(result)-n.euploids-n.aneuploids)/length(result)*100

    result = data.frame(values=result)
    result$colour = factor(ifelse(result$values==0, "#009806",
                           ifelse(result$values==input$n.samples, "red", "orange")),
                           levels = c("#009806", "orange", "red"))

    p = plot_ly(x = result$values,
            type="histogram",
            color = result$colour,
            colors = c("#009806", "orange", "red")) %>%
      layout(showlegend = F,
             xaxis = list( title = "Number of aneuploid cells in biopsy",
                           range=c(-0.5,input$n.samples+0.5)),
             yaxis = list( title = "Number of biopsies",
                           range=c(0, input$n.cells))) %>%
      add_annotations(text=paste0(format(euploid.ratio, nsmall=1, digits = 3),"%\neuploid"),
                      x=0,
                      y=input$n.cells*0.9,
                      align="center",
                      showarrow=F) %>%
      add_annotations(text=paste0(format(aneuploid.ratio, nsmall=1, digits = 3),"%\naneuploid"),
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
