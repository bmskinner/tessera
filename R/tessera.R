# Create an embryo using the Fibonacci lattice model
# library(assertthat)

# The number of cells to be considered neighbours to a given cell
.N_NEIGHBOURS=6

#' Create a sphere of evenly spaced cells
#'
#' @param n.points the number of cells in the embryo
#'
#' @return a data frame with coordinates
#'
.create.blank.sphere = function(n.points){
  # Make a sphere of evenly spaced points using the Fibonacci spiral
  indices = seq(0, n.points-1, 1)+0.5
  phi = acos(pmin(pmax( 1-2*indices/n.points,-1.0),1.0)) # constrain to avoid rounding errors
  theta = pi * (1 + sqrt(5)) * indices

  x = cos(theta)*sin(phi)
  y = sin(theta)*sin(phi)
  z = cos(phi)
  d = as.data.frame(cbind(x, y, z))

  # Create distance matrix for each point
  # Set the .N_NEIGHBOURS closest points to be neighbours
  for(i in 1:nrow(d)){
    dist = sqrt( (d$x-x[i])**2 + (d$y-y[i])**2 + (d$z-z[i])**2) # distance between points
    d[[paste0("d", i)]] = dist # create a column to store the distances
    # A point is a neighbour if it is not this point, and it is in the list of closest points
    d[[paste0("isNeighbour", i)]] = dist > 0 & dist <= max(head(sort(dist), n=.N_NEIGHBOURS+1))
  }
  d$isAneuploid = 0

  # assertthat::assert_that(nrow(d)==n.points, msg=paste("Expected", n.points, "cells, found", nrow(d)))
  return(d)
}

# Test if any of the neighbouring cells have an aneuploidy
# d - the blastocyst
# index - the cell to test
# Returns true if any of the closest cells are aneuploid
.has.adjacent.aneuploid = function(d, index){
  adj.list = d[[paste0("isNeighbour", index)]]
  return(any(adj.list & d$isAneuploid>0))
}

#' Count the number of potential seed locations
#'
#' Cells that are not adjacent to a seed
#'
#' @param d the embryo
#' @param n.sampled.cells the number of cells to biopsy
#' @param  seed.sample the index of the cell to begin biopsying
#'
#' @return the number of potential seeds
.count.empty.blocks = function(d){
  n = 0
  for(i in 1:nrow(d)){
    if(d$isAneuploid[i]>0 & !.has.adjacent.aneuploid(d, i)) n = n+1
  }
  return(n)
}

#' Create an embryo
#'
#' A sphere of cells is created with the given proportion of aneuploidies.
#' Aneuploid cells are either adjacent or dispersed
#'
#' @param n.cells the number of cells in the embryo
#' @param prop.aneuploid the proportion of aneuploid cells (0-1)
#' @param dispersion the dispersion of the aneuploid cells (0-1)
#'
#' @return an embryo data frame
#'
#' @examples
#' embryo <- create.embryo(20, 0.1, 0.9)
create.embryo = function(n.cells, prop.aneuploid, dispersion){

  d = .create.blank.sphere(n.cells)

  # Shortcut the easy cases
  if(prop.aneuploid==0) return(d)

  if(prop.aneuploid==1){
    d$isAneuploid=1
    return(d)
  }

  # We must have an integer value of at least one aneuploid cell
  n.aneuploid = ceiling(max(1, n.cells * prop.aneuploid))

  # The approach for dispersal is to set seed cells which will
  # grow into separate aneuploid patches. The more dispersion, the more
  # initial seeds.

  # Choose number of seeds for aneuploid regions
  n.seeds = ceiling(max(1, n.aneuploid * dispersion))
  n.to.make = n.seeds

  # We can disperse up to a certain number of initial blocks with
  # no aneuploid neighbours. After this, every cell will have at least
  # one aneuploid neighbour. We stop a bit before this to make the maths simpler.
  initial.blocks = max(1,floor(n.cells/.N_NEIGHBOURS))

  # Disperse seeds as much as possible
  while(initial.blocks>0 & n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isAneuploid[seed]>0) next
    if(.has.adjacent.aneuploid(d, seed)) next # spread seeds out
    d$isAneuploid[seed] = 1
    n.to.make = n.to.make-1L
    initial.blocks = initial.blocks-1L
  }

  # When all dispersed seeds have been added, add the remaining seeds randomly
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isAneuploid[seed]>0) next
    d$isAneuploid[seed] = 1
    n.to.make = n.to.make-1L
  }
  # assertthat::assert_that(sum(d$isAneuploid)==n.seeds,
                          # msg = paste("Expected", n.seeds, "seeds, found", sum(d$isAneuploid)))

  # Grow the seeds into neighbouring cells for remaining aneuploid cells
  n.to.make = n.aneuploid - n.seeds
  while(n.to.make>0){
    seed = sample.int(n.cells, 1)
    if(d$isAneuploid[seed]>0) next # skip cells already aneuploid
    if(!.has.adjacent.aneuploid(d, seed)) next # only grow next to existing aneuploid
    d$isAneuploid[seed] = 1
    n.to.make = n.to.make-1
  }
  # assertthat::assert_that(sum(d$isAneuploid)==n.aneuploid,
  #                         msg = paste("Expected", n.aneuploid, "aneuploids, found", sum(d$isAneuploid)))
  return(d)
}

#' Take a sample from an embryo
#'
#' The cell at the given index is taken,
#' plus the closest n neighbouring cells where n = n.sampled.cells-1.
#'
#' @param embryo an embryo as created by \code{create.embryo}
#' @param n.sampled.cells the number of cells to biopsy
#' @param index.cell the index of the cell to begin biopsying. Must be a value
#'  between 1 and \code{nrow(embryo)}
#'
#' @return the number of aneuploid cells in the biopsy
#'
#' @examples
#' e <- create.embryo(100, 0.1, 0.1)
#' take.one.biopsy(e, 5, 1)
take.one.biopsy = function(embryo, n.sampled.cells, index.cell){
  if(index.cell < 1 | index.cell > nrow(embryo)){
    warning(paste("index.cell (", index.cell ,") must be between 1 and", nrow(embryo)))
    return()
  }
  sample.list = embryo[[paste0("d", index.cell)]]

  embryo$isSampled = embryo[[paste0("d", index.cell)]] <= max(head(sort(sample.list), n=n.sampled.cells))
  return(sum(embryo[embryo$isSampled,]$isAneuploid>0))
}


#' Find the number of aneuploid cells
#'
#' Take all possible biopsies of the given size from the given
#' embryo
#'
#' @param embryo an embryo as created by \code{create.embryo}
#' @param n.cells.per.sample the number of cells to take in each biopsy
#'
#' @return an integer vector of the number of aneuploid cells in each biopsy
#' @export
#'
#' @examples
#' e <- create.embryo(100, 0.1, 0.1)
#' take.all.biopsies(e, 5)
take.all.biopsies = function(embryo, n.cells.per.sample){
  result = c()
  for(i in 1:nrow(embryo)){ # sample each cell in turn, so we get every cell
    f = take.one.biopsy(embryo, n.cells.per.sample, i)
    result = c(result, f)
  }
  return(result)
}


