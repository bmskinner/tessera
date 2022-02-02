#' # Create an embryo using the Fibonacci lattice model
#' # library(assertthat)
#'
#' # The number of cells to be considered neighbours to a given cell
#' .N_NEIGHBOURS = 6
#'
#' #' Create a sphere of evenly spaced cells
#' #'
#' #' @param n.points the number of cells in the embryo
#' #'
#' #' @return a data frame with coordinates
#' #'
#' .create.blank.sphere = function(n.points){
#'   # Make a sphere of evenly spaced points using the Fibonacci spiral
#'   indices = seq(0, n.points-1, 1)+0.5
#'   phi = acos(pmin(pmax( 1-2*indices/n.points,-1.0),1.0)) # constrain to avoid rounding errors
#'   theta = pi * (1 + sqrt(5)) * indices
#'
#'   x = cos(theta)*sin(phi)
#'   y = sin(theta)*sin(phi)
#'   z = cos(phi)
#'   d = as.data.frame(cbind(x, y, z))
#'
#'   # Create distance matrix for each point
#'   # Set the .N_NEIGHBOURS closest points to be neighbours
#'   for(i in 1:nrow(d)){
#'     dist = sqrt( (d$x-x[i])**2 + (d$y-y[i])**2 + (d$z-z[i])**2) # distance between points
#'     d[[paste0("d", i)]] = dist # create a column to store the distances
#'     # A point is a neighbour if it is not this point, and it is in the list of closest points
#'     d[[paste0("isNeighbour", i)]] = dist > 0 & dist <= max(head(sort(dist), n=.N_NEIGHBOURS+1))
#'   }
#'   d$isAneuploid = 0
#'
#'   # assertthat::assert_that(nrow(d)==n.points, msg=paste("Expected", n.points, "cells, found", nrow(d)))
#'   return(d)
#' }
#'
#' # Test if any of the neighbouring cells have an aneuploidy
#' # d - the blastocyst
#' # index - the cell to test
#' # chromosome - the chromosome to test
#' # Returns true if any of the closest cells are aneuploid
#' .has.adjacent.aneuploid = function(d, index, chromosome){
#'   adj.list = d[[paste0("isNeighbour", index)]]
#'   isAenu = bitwAnd(d$isAneuploid, .bit.value(chromosome))==.bit.value(chromosome)
#'
#'   return(any(adj.list & isAenu))
#' }
#'
#' #' Count the number of potential seed locations
#' #'
#' #' Cells that are not adjacent to a seed
#' #'
#' #' @param d the embryo
#' #' @param chromosome the chromosome to test
#' #'
#' #' @return the number of potential seeds
#' .count.empty.blocks = function(d, chromosome){
#'   n = 0
#'   for(i in 1:nrow(d)){
#'     if(d$isAneuploid[i]>0 & !.has.adjacent.aneuploid(d, i, chromosome)) n = n+1
#'   }
#'   return(n)
#' }
#'
#' #' Set a cell to contain an aneuploid chromosome
#' #'
#' #' @param embryo the embryo
#' #' @param cell.index the cell to affect
#' #' @param chromosome the chromosome to make aneuploid (1-25)
#' #'
#' #' @return the modified embryo
#' #' @export
#' #'
#' #' @examples
#' set.aneuploid = function(embryo, cell.index, chromosome){
#'   if(chromosome < 1 | chromosome>31) {
#'     warning("Chromosome must be in range 1-31")
#'     return(embryo)
#'   }
#'   embryo$isAneuploid[cell.index] = bitwOr(embryo$isAneuploid[cell.index],
#'                                           .bit.value(chromosome))
#'   return(embryo)
#' }
#'
#' #' Test if the given chromosome in the given cell is aneuploid
#' #'
#' #' @param embryo the embryo
#' #' @param cell.index the cell to test (0 for all cells)
#' #' @param chromosome the chromosome to test
#' #'
#' #' @return if cell.index is >0, return true if the chromosome is aneuploid,
#' #' false otherwise. If cell.index is 0, return a boolean vector of all cells.
#' #' @export
#' #'
#' #' @examples
#' is.aneuploid = function(embryo, cell.index, chromosome){
#'   if(chromosome < 1 | chromosome>31) {
#'     warning("Chromosome must be in range 1-31")
#'     return(F)
#'   }
#'
#'   if(cell.index > nrow(embryo)){
#'     warning("Cell index must be in range 0-nrow(embryo)")
#'     return(F)
#'   }
#'
#'   # Return a vector if all cells requested
#'   if(cell.index==0){
#'     return(bitwAnd(embryo$isAneuploid,
#'                    .bit.value(chromosome)) == .bit.value(chromosome))
#'   }
#'
#'   # Otherwise return just the single cell value
#'   return(bitwAnd(embryo$isAneuploid[cell.index],
#'                  .bit.value(chromosome)) == .bit.value(chromosome))
#' }
#'
#' #' Get the bit value for the given chromosome
#' #'
#' #' @param chromosome the chromosome to test
#' #'
#' #' @return the bit value
#' #'
#' #' @examples
#' .bit.value = function(chromosome){
#'   return(2^(chromosome-1))
#' }
#'
#'
#' #' count the number of aneuploid cells for the given chromosome
#' #'
#' #' @param embryo the embryo
#' #' @param chromosome the chromosome to test
#' #'
#' #' @return the number of aneuploid cells
#' #' @export
#' #'
#' #' @examples
#' count.aneuploid = function(embryo, chromosome){
#'   return(sum(bitwAnd(embryo$isAneuploid,
#'                      .bit.value(chromosome)) == .bit.value(chromosome)))
#' }
#'
#' #' Create an embryo
#' #'
#' #' A sphere of cells is created with the given proportion of aneuploidies.
#' #' Aneuploid cells are either adjacent or dispersed. The concordance of aneuploid
#' #' cells for each chromosome can be specificed; if fully concordant, a cell aneuploid
#' #' for chr1 will also be aneuploid for chr2 etc.
#' #'
#' #' @param n.cells the number of cells in the embryo
#' #' @param prop.aneuploid the proportion vector of aneuploid cells (0-1) per chromosome
#' #' @param dispersion the dispersion vector of the aneuploid cells (0-1)
#' #' @param concordance the concordance between aneuploid cells for each chromosome (0-1).
#' #' @param seed the seed for the RNG. Defaults to NULL
#' #'
#' #' @return an embryo data frame
#' #'
#' #' @examples
#' #' embryo <- create.embryo(20, c(0.1, 0, 0, 0.4), 0.9)
#' create.embryo = function(n.cells, prop.aneuploids, dispersions, concordance=1, seed=NULL){
#'
#'   if(length(prop.aneuploids)>31){
#'     warning("Trying to set aneuploidies for more than 31 chromosomes")
#'     return(NULL)
#'   }
#'
#'   if(length(dispersions)>31){
#'     warning("Trying to set dispersions for more than 31 chromosomes")
#'     return(NULL)
#'   }
#'   if(length(dispersions)!=length(prop.aneuploids)){
#'     warning("Must have the same dimensions for input vectors")
#'     return(NULL)
#'   }
#'
#'   if(concordance<0 | concordance>1){
#'     warning("Concordance must be in the range 0-1")
#'     return(NULL)
#'   }
#'
#'   set.seed(seed)
#'
#'   embryo = .create.blank.sphere(n.cells)
#'
#'   for(chr in 1:length(prop.aneuploids)){
#'     embryo = set.aneuploidies(embryo,
#'                          chr,
#'                          prop.aneuploids[chr],
#'                          dispersions[chr],
#'                          concordance)
#'   }
#'
#'   return(embryo)
#' }
#'
#' #' Set aneuploidies in an embryo
#' #'
#' #' A sphere of cells is created with the given proportion of aneuploidies.
#' #' Aneuploid cells are either adjacent or dispersed
#' #'
#' #' @param embryo an embryo as created by \code{create.embryo}
#' #' @param chromosome the chromosome to set aneuploidies for
#' #' @param prop.aneuploid the proportion of aneuploid cells (0-1)
#' #' @param dispersion the dispersion of the aneuploid cells (0-1)
#' #' @param concordance the concordance between aneuploid cells for each chromosome (0-1).
#' #'
#' #' @return the embryo with aneuploidies
#' #'
#' #' @examples
#' #' embryo <- set.aneuploidies(embryo, 1, 0.1, 0.9)
#' set.aneuploidies = function(embryo, chromosome, prop.aneuploid, dispersion, concordance){
#'
#'   # Shortcut the easy cases
#'   if(prop.aneuploid==0) return(embryo)
#'
#'   if(prop.aneuploid==1){
#'     for(i in 1:nrow(embryo)){
#'       embryo = set.aneuploid(embryo, i, chromosome)
#'     }
#'     return(embryo)
#'   }
#'
#'   n.cells = nrow(embryo)
#'   # cat("Embryo has", n.cells, "cells\n")
#'
#'   # We must have an integer value of at least one aneuploid cell
#'   n.aneuploid = ceiling(max(1, n.cells * prop.aneuploid))
#'   # cat("Creating", n.aneuploid, "aneuploid cells for chr", chromosome,"\n")
#'
#'   # Decide how many cells need to be concordant with the previous
#'   # chromosome (if we are above chromosome 1)
#'   n.concordant = 0
#'   concordant.cells = rep(F, nrow(embryo))
#'   if(chromosome>1){
#'     prev.chr = chromosome - 1 # to be used when chr>1 only
#'     concordant.cells = is.aneuploid(embryo, 0, prev.chr)
#'     # cat("Prev chr", prev.chr, "has", length(concordant.cells[concordant.cells==T]), "aneuploid cells\n")
#'     n.concordant = length(concordant.cells[concordant.cells==T])*concordance
#'     # if there are more aneuploids in the prev chromosome, we can't match all
#'     n.concordant = min(n.aneuploid, n.concordant)
#'     # cat("Expecting", n.concordant, "concordant cells with chr", prev.chr, "\n")
#'   }
#'
#'   # The approach for dispersal is to set seed cells which will
#'   # grow into separate aneuploid patches. The more dispersion, the more
#'   # initial seeds.
#'
#'   # Choose number of seeds for aneuploid regions
#'   n.seeds = ceiling(max(1, n.aneuploid * dispersion))
#'   n.to.make = n.seeds
#'
#'   # We can disperse up to a certain number of initial blocks with
#'   # no aneuploid neighbours. After this, every cell will have at least
#'   # one aneuploid neighbour. We stop a bit before this to make the maths simpler.
#'   initial.blocks = max(1,floor(n.cells/.N_NEIGHBOURS))
#'
#'   # cat("Creating up to", initial.blocks, "initial seed positions\n")
#'   # cat("Creating", n.to.make, "initial seeds\n")
#'
#'   # Disperse seeds as much as possible
#'   while(initial.blocks>0 & n.to.make>0){
#'     seed = sample.int(n.cells, 1)
#'     if(is.aneuploid(embryo, seed, chromosome)) next
#'     if(.has.adjacent.aneuploid(embryo, seed, chromosome)) next # spread seeds out
#'     if(n.concordant>0 & !concordant.cells[seed]) next # skip non concordant cells
#'     embryo = set.aneuploid(embryo, seed, chromosome)
#'     n.to.make = n.to.make-1L
#'     initial.blocks = initial.blocks-1L
#'     n.concordant = max(0, n.concordant - 1L)
#'   }
#'
#'   # cat(n.concordant, "concordant cells remaining to place\n")
#'   # cat("Creating", n.to.make, "supplementary seeds\n")
#'
#'   # When all dispersed seeds have been added, add the remaining seeds randomly
#'   while(n.to.make>0){
#'     seed = sample.int(n.cells, 1)
#'     if(is.aneuploid(embryo, seed, chromosome)) next
#'     if(n.concordant>0 & !concordant.cells[seed]) next # skip non concordant cells
#'     embryo = set.aneuploid(embryo, seed, chromosome)
#'     n.to.make = n.to.make-1L
#'     n.concordant = max(0, n.concordant - 1L)
#'   }
#'   # assertthat::assert_that(sum(d$isAneuploid)==n.seeds,
#'   # msg = paste("Expected", n.seeds, "seeds, found", sum(d$isAneuploid)))
#'
#'   # Grow the seeds into neighbouring cells for remaining aneuploid cells
#'   # cat(n.concordant, "concordant cells remaining to place\n")
#'   n.to.make = n.aneuploid - n.seeds
#'   # cat("Placing", n.to.make, "final aneuploid cells\n")
#'
#'   # Handle any remaining concordant cells first - the placement rules don't apply
#'   while(n.to.make>0){
#'     seed = sample.int(n.cells, 1)
#'     if(n.concordant>0){
#'       # cat("Placing ", n.concordant, "concordant aneuploid cells\n")
#'       if(!concordant.cells[seed]) next
#'       if(is.aneuploid(embryo, seed, chromosome)) next # skip cells already aneuploid
#'       embryo = set.aneuploid(embryo, seed, chromosome)
#'       n.concordant = max(0, n.concordant - 1L)
#'
#'     } else {
#'       # cat("Placing ", n.to.make, "non-concordant aneuploid cells\n")
#'       if(is.aneuploid(embryo, seed, chromosome)) next # skip cells already aneuploid
#'       if(!.has.adjacent.aneuploid(embryo, seed, chromosome)) next # only grow next to existing aneuploid
#'       embryo = set.aneuploid(embryo, seed, chromosome)
#'     }
#'
#'     n.to.make = n.to.make-1
#'
#'   }
#'
#'   # cat("Finished placing chr", chromosome, "\n")
#'   # assertthat::assert_that(sum(d$isAneuploid)==n.aneuploid,
#'   #                         msg = paste("Expected", n.aneuploid, "aneuploids, found", sum(d$isAneuploid)))
#'   return(embryo)
#' }
#'
#' #' Take a sample from an embryo
#' #'
#' #' The cell at the given index is taken,
#' #' plus the closest n neighbouring cells where n = n.sampled.cells-1.
#' #'
#' #' @param embryo an embryo as created by \code{create.embryo}
#' #' @param n.sampled.cells the number of cells to biopsy
#' #' @param index.cell the index of the cell to begin biopsying. Must be a value
#' #'  between 1 and \code{nrow(embryo)}
#' #'  @param chromosome the chromosome to test
#' #'
#' #' @return the number of aneuploid cells in the biopsy
#' #'
#' #' @examples
#' #' e <- create.embryo(100, 0.1, 0.1)
#' #' take.one.biopsy(e, 5, 1)
#' take.one.biopsy = function(embryo, n.sampled.cells, index.cell, chromosome){
#'   if(index.cell < 1 | index.cell > nrow(embryo)){
#'     warning(paste("index.cell (", index.cell ,") must be between 1 and", nrow(embryo)))
#'     return(NULL)
#'   }
#'
#'   if(chromosome < 0 | chromosome>31){
#'     warning(paste("Chromosome (", chromosome ,") must be between 0 and 31"))
#'     return(NULL)
#'   }
#'
#'   sample.list = embryo[[paste0("d", index.cell)]]
#'
#'   isSampled = embryo[[paste0("d", index.cell)]] <= max(head(sort(sample.list), n=n.sampled.cells))
#'
#'   # count all chromsomes; don't care which chromosome is aneuploid
#'   # just is aneuploid or is not aneuploid
#'   if(chromosome==0){
#'     return(sum(embryo[isSampled,]$isAneuploid>0))
#'   }
#'
#'   return(count.aneuploid(embryo[isSampled,], chromosome))
#' }
#'
#'
#' #' Find the number of aneuploid cells
#' #'
#' #' Take all possible biopsies of the given size from the given
#' #' embryo
#' #'
#' #' @param embryo an embryo as created by \code{create.embryo}
#' #' @param n.cells.per.sample the ideal number of cells to take in each biopsy
#' #' @param chromosome the chromosome to test, or 0 for all chromosomes
#' #' @param n.cells.fixed true to take the same number of cells in each biopsy, false to
#' #' use a distribution model
#' #' @param n.cells.sd the standard deviation of the normal distribution used to model
#' #' the cell biopsy size
#' #'
#' #' @return an integer vector of the number of aneuploid cells in each biopsy
#' #' @export
#' #'
#' #' @examples
#' #' e <- create.embryo(100, 0.1, 0.1)
#' #' take.all.biopsies(e, 5, 1)
#' take.all.biopsies = function(embryo, n.cells.per.sample, chromosome, n.cells.fixed=T, n.cells.sd = 1) {
#'
#'   if(chromosome < 0 | chromosome > 31) {
#'     warning(paste("Chromosome (", chromosome ,") must be between 0 and 31"))
#'     return()
#'   }
#'
#'
#'   #' Model the number of biopsied cells in a sample.
#'   #'
#'   #' When biopsying cells, we may not get exactly the target number; there
#'   #' may be one too many or too few. We model the number of cells to take in
#'   #' a biopsy as a normal distribution with a mean around the desired number of
#'   #' cells and a standard deviation provided.
#'   create.n.cells.function = function(){
#'
#'     if (n.cells.fixed) {
#'       # If we are keeping a fixed number of cells in each biopsy, we don't need a
#'       # model
#'       return( function(){ n.cells.per.sample })
#'     } else {
#'       # model the number of biopsied cells as a distribution
#'       # Ensure sd is at least 1
#'       return( function() { max(1, ceiling(rnorm(1,
#'                                          mean = n.cells.per.sample,
#'                                          sd = max(1, n.cells.sd))))  })
#'     }
#'   }
#'
#'   fn = create.n.cells.function()
#'
#'
#'   # If just one chromosome sampled
#'   result = c()
#'   for(i in 1:nrow(embryo)) { # sample each cell in turn, so we get every cell
#'     f = take.one.biopsy(embryo, fn(), i, chromosome)
#'     result = c(result, f)
#'   }
#'   return(result)
#' }
