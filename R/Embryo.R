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
  create.blank.sphere = function(n.points){
    # Make a sphere of evenly spaced points using the Fibonacci spiral
    indices = seq(0, n.points - 1, 1) + 0.5
    phi = acos(pmin(pmax( 1 - 2*indices/n.points,-1.0),1.0)) # constrain to avoid rounding errors
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
    return(d)
  }
  d = create.blank.sphere(nCells)
  ploidy = data.frame(matrix(data = 2, nrow = nCells, ncol = nChrs))
  colnames(ploidy) = paste0("chr", 1:nChrs)
  new("Embryo",
      x = d[,1], y = d[,2], z = d[,3],
      dists = (d[,-(1:3)]),
      ploidy = ploidy)
}

# Add a generic method for length to get the number of cells
setMethod("length", "Embryo", function(x) { length(x@x) } )

em = Embryo(10, 5)
