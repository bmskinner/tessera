#' Create the tessera server
#'
#' @import shiny
#' @import shinythemes
#' @import plotly
#' @import ggplot2
#' @return the server
tesseraServer <- function(input, output, session){

  # Fixed definitions
  pgdis.classes = factor(c("Euploid", "Low level", "High level", "Aneuploid"),
                         levels = c("Euploid", "Low level", "High level", "Aneuploid"))

  to.pgdis.class = function(prop.aneuploidy){
    if(prop.aneuploidy<0.20) return(pgdis.classes[1])
    if(prop.aneuploidy<0.40) return(pgdis.classes[2])
    if(prop.aneuploidy<=0.80) return(pgdis.classes[3])
    return(pgdis.classes[4])
  }


  # Create and store seed value for model embryo only when button is clicked
  seedVals = shiny::eventReactive(input$new.embryo, {
    sample.int(.Machine$integer.max, size=1)
  })
  
  
  # Create and store seeds for embryos in the ranking pool only when button is clicked
  pool.embryo.seeds <- shiny::eventReactive(input$rank.embryos, {
    sample.int(.Machine$integer.max, size=input$n.pool)
  })
  
  
  # Create a pool of embryos for ranking
  calculateRanks = shiny::reactive({
    aneuploidies <- runif(input$n.pool, min = 0, max = 1)
    
    # Create an embryo with the given aneuploidy and seed, and take one random biopsy
    embryos <- mapply(function(a, s){
      tessera::Embryo(n.cells = 200, n.chrs = 1, prop.aneuploid = a, dispersal = input$dispersal, rng.seed = s)
    }, a= aneuploidies,  s = pool.embryo.seeds())
    
    biopsies <- sapply(embryos, tessera::takeBiopsy, biopsy.size=input$biopsy.size)

    # Rank the embryos by aneuploidy
    best.from.biopsy <- which(biopsies %in% sort(biopsies)[1:input$n.transfer])
    best.from.reality <- which(aneuploidies %in% sort(aneuploidies)[1:input$n.transfer])
    

    is.real.best = aneuploidies %in% sort(aneuploidies)[1:input$n.transfer]
    is.biop.best = biopsies %in% sort(biopsies)[1:input$n.transfer]
    
    # There may be ties in ranks; in this case, increase the number of transferred embryos
    is.transferrable = sum(is.biop.best)
    pct.correct = sum(is.real.best&is.biop.best) / is.transferrable * 100
    
    # Calculate the percent of embryos correctly selected for transfer by biopsy
    # pct.corrrect <- sum(best.from.biopsy %in% best.from.reality) / input$n.transfer * 100
    
    list("embryos"=embryos,
         "is.real.best" = is.real.best,
         "is.biop.best" = is.biop.best,
         "best.from.reality" = best.from.reality,
         "best.from.biopsy" = best.from.biopsy,
         "aneuploidies"=aneuploidies,
        "biopsies"=biopsies,
         "pct.correct"=pct.correct)
  })
  



  # Generate an embryo with the desired parameters and take all possible biopsies
  calculateData = reactive({

    all.chr = T

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

    return(list("embryo"   = embryo,
                "biopsies" = result,
                "classes"  = biopsy.classes,
                "true.class"= to.pgdis.class(input$proportion)))
  })
  
  

  # Plot the embryo
  output$embryo.model = plotly::renderPlotly({
    embryo = calculateData()[['embryo']]
    plot(embryo)

  })

  # Create plot of accuracies
  output$biopsy.accuracy = plotly::renderPlotly({

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
  output$biopsy.histogram = plotly::renderPlotly({

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

  
  # Create text for rank results
  output$pool.data <- shiny::renderText(paste(calculateRanks()[['pct.correct']]))
  
  # Render the ranking output plot
  output$pool.aneuplodies <- shiny::renderPlot({
    
    result = calculateRanks()
    
    pct.correct = round(result[['pct.correct']], digits = 2)
    
    best.embryo.cutoff = max(result[['aneuploidies']][result[['best.from.reality']]])+0.005
    
    best.biopsy.cutoff =  max(result[['biopsies']][result[['best.from.biopsy']]])+0.05
    
    colours = dplyr::case_when( result[['is.biop.best']]&result[['is.real.best']] ~ "Correctly transferred",
                                result[['is.biop.best']] ~ "Wrongly transferred",
                                result[['is.real.best']] ~ "Wrongly rejected",
                                T ~ "Correctly rejected")
    
    ggplot2::ggplot()+
      annotate("rect", xmin = -Inf, xmax = best.embryo.cutoff, ymin = -Inf, ymax = best.biopsy.cutoff, fill = 'darkgreen', alpha=0.3) +
      annotate("rect", xmin = -Inf, xmax = best.embryo.cutoff, ymin = best.biopsy.cutoff, ymax = Inf, fill = 'orange', alpha=0.3) +
      annotate("rect", xmin = best.embryo.cutoff, xmax = Inf, ymin = -Inf, ymax= best.biopsy.cutoff, fill = 'red', alpha=0.3) +
      geom_point(aes( x=result[['aneuploidies']],  
                    y =result[['biopsies']], 
                    col = colours),
                 size=3)+
      geom_vline(xintercept = best.embryo.cutoff, col='black')+
      geom_hline(yintercept = best.biopsy.cutoff, col="black")+
      geom_text(aes(x=0, y=input$biopsy.size+1, label="Wrongly rejected", hjust=0))+
      geom_text(aes(x=1, y=-1, label="Wrongly transferred", hjust=1))+
      geom_text(aes(x=1, y=input$biopsy.size+1, label="Correctly rejected", hjust=1))+
      geom_text(aes(x=0, y=-1, label="Correctly transferred", hjust=0))+

      labs(x = "True embryo aneuploidy", y = "Aneuploid cells in biopsy",
           col="Biopsy selected correcly",
           title = paste0(pct.correct, "% of the best embryos would be selected for transfer"))+
      scale_y_continuous(breaks = seq(0, input$biopsy.size, 1))+
      scale_colour_manual(values = c("Correctly rejected"="darkgrey", 
                                     "Correctly transferred"="darkgreen", 
                                     "Wrongly rejected"="orange", 
                                     "Wrongly transferred"="red"))+
      theme_bw()+
      theme(legend.position = "none")
    
  })
}
