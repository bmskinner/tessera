# Demo of how to create an S4 class
setClass("Embryo",

         # Define fields
         representation(
           x = "numeric",
           y = "numeric",
           z = "numeric",
           dists = "data.frame",
           ploidy = "data.frame"
         ),

         # Define defaults
         prototype(
           x = NA_real_,
           y = NA_real_,
           z = NA_real_,
           dists = data.frame(),
           ploidy = data.frame()

         )
)

# Constructor is not part of the class definition
Embryo <- function(nCells, nChrs){
  .N_NEIGHBOURS = 6
  # Make a sphere of evenly spaced points using the Fibonacci lattice
  indices = seq(0, nCells - 1, 1) + 0.5
  phi = acos(pmin(pmax( 1 - 2*indices/nCells,-1.0),1.0)) # constrain to avoid rounding errors
  theta = pi * (1 + sqrt(5)) * indices

  x = cos(theta)*sin(phi)
  y = sin(theta)*sin(phi)
  z = cos(phi)
  d = as.data.frame(cbind(x, y, z))

  # Create distance matrix for each point
  # Set the .N_NEIGHBOURS closest points to be neighbours
  for(i in 1:nrow(d)){
    dist = sqrt( (d$x - x[i])**2 + (d$y - y[i])**2 + (d$z - z[i])**2) # distance between points
    d[[paste0("d", i)]] = dist # create a column to store the distances
    # A point is a neighbour if it is not this point, and it is in the list of closest points
    d[[paste0("n", i)]] = dist > 0 & dist <= max(head(sort(dist), n = .N_NEIGHBOURS + 1))
  }

  ploidy = data.frame(matrix(data = 2, nrow = nCells, ncol = nChrs))
  colnames(ploidy) = paste0("chr", 1:nChrs)
  new("Embryo",
      x = d[,1], y = d[,2], z = d[,3],
      dists = d[,-(1:3)],
      ploidy = ploidy)
}

# Override plot function for an Embryo object
setMethod("plot", "Embryo", function(x){

  colours = factor(rowSums(x@ploidy)==ncol(x@ploidy*2), levels = c(F, T))

  plotly::plot_ly( type="scatter3d",
           mode="markers",
           colors = c("#22FF22", "#222222"),
           color=colours,
           marker = list(
             size = 25,
             line = list(
               color = '#111111',
               width = 1
             )
           ),
           hoverinfo="none") %>%
    plotly::add_trace(
      x = x@x,
      y = x@y,
      z = x@z,
      showlegend = F,
      hoverinfo = 'skip'
    ) %>%
    plotly::add_annotations( text="Click and drag to rotate",
                     xref="paper", yref="paper",
                     x=0.0, xanchor="left",
                     y=1.0, yanchor="top",
                     legendtitle=TRUE, showarrow=FALSE ) %>%
    plotly::config(displayModeBar = FALSE, scrollZoom = F) %>%
    plotly::layout(scene = list(
      xaxis = list(autorange = F,
                   fixedrange = TRUE,
                   showgrid = F,
                   showline = F,
                   showticklabels = F,
                   showaxeslabels = F,
                   title = "",
                   zeroline = F,
                   range = list(-1, 1)),
      yaxis = list(fixedrange = TRUE,
                   autorange = F,
                   showgrid = F,
                   showline = F,
                   showaxeslabels = F,
                   showticklabels = F,
                   title = "",
                   zeroline = F,
                   range = list(-1, 1)),
      zaxis = list(autorange = F,
                   fixedrange = TRUE,
                   showgrid = F,
                   showline = F,
                   showaxeslabels = F,
                   showticklabels = F,
                   zeroline = F,
                   title = "",
                   range = list(-1, 1))))
})

# Add a method override for existing generic function length to get the number of cells
setMethod("length", "Embryo", function(x) { length(x@x) } )

# Make a generic function
setGeneric(name="getCell",
           def = function(x, cell) { standardGeneric("getCell")})



# Provide implementation of the function for an Embryo
setMethod("getCell", signature = "Embryo", function(x, cell){ paste(x@x[cell], x@y[cell], x@z[cell])} )

em = Embryo(300, 1)
