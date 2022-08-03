#' An S4 class to model an embryo
#'
#' @slot x numeric. x coordinates of cells
#' @slot y numeric.  y coordinates of cells
#' @slot z numeric.  z coordinates of cells
#' @slot aneu numeric. fraction of aneuploid cells
#' @slot disp numeric. fraction of dispersal of cells
#' @slot dists data.frame. pairwise distances between cells
#' @slot euploidy numeric. number of chromosomes considered euploid
#' @slot ploidy data.frame. number of chromosomes per cell
#'
#' @return
#' @export
setClass("Embryo",
         
         # Define fields
         representation(
           x = "numeric",
           y = "numeric",
           z = "numeric",
           aneu = "numeric",
           disp = "numeric",
           dists = "data.frame",
           neighbours = "data.frame",
           euploidy = "numeric",
           ploidy = "data.frame"
         ),
         
         # Define defaults
         prototype(
           x = NA_real_,
           y = NA_real_,
           z = NA_real_,
           aneu = NA_real_,
           disp = NA_real_,
           dists = data.frame(),
           neighbours = data.frame(),
           euploidy = NA_real_,
           ploidy = data.frame()
           
         )
)

#' Create an embryo
#'
#' A sphere of cells is created with the given proportion of aneuploidies.
#' Aneuploid cells are either adjacent or dispersed. The concordance of aneuploid
#' cells for each chromosome can be specificed; if fully concordant, a cell aneuploid
#' for chr1 will also be aneuploid for chr2 etc.
#'
#' @param n.cells the number of cells in the embryo
#' @param n.chrs the number of chromosome pairs per cell
#' @param prop.aneuploid the proportion vector of aneuploid cells (0-1) per chromosome
#' @param dispersal the dispersion vector of the aneuploid cells (0-1)
#' @param concordance the concordance between aneuploid cells for each chromosome (0-1).
#' @param embryo.size.fixed if true, the embryo is exactly the size in \code{n.cells}. If false, the embryo
#' size can vary according to \code{embryo.size.sd}.
#' @param embryo.size.sd the standard deviation of cell number if \code{embryo.size.fixed} is true. 
#' The actual embryo size will be sampled from a normal distribution with mean of \code{n.cells} and 
#' standard deviation \code{embryo.size.sd}.
#' @param euploidy the number of copies of a chromosome to consider euploid. For a diploid embryo this should be 2.
#' @param rng.seed the seed for the RNG. Defaults to NULL. Use this to get the same embryo each time
#'
#' @return an Embryo object
#'
#' @examples
#' Create an embryo with 200 cells, 20% aneuploid and a single pair of chromosomes
#' per cell. Aneuploid cells are highly dispersed
#' embryo <- Embryo(n.cells = 200, n.chrs = 1,  prop.aneuploid = 0.2,
#'                  dispersal =  0.9)
#'
#' Create the embryo above, but using a fixed seed for the random number generator
#' so the resulting embryo is reproducible.
#' embryo <- Embryo(n.cells = 200, n.chrs = 1,  prop.aneuploid = 0.2,
#'                  dispersal =  0.9, rng.seed = 42)
#'
#' Create an embryo with 3 pairs of chromosomes per cell, with all chromosome pairs
#' aneuploid in the same cells.
#' embryo <- Embryo(n.cells = 200, n.chrs = 3,  prop.aneuploid = 0.2,
#'                  dispersal =  0.9, concordance = 1)
#'
#' As above, but specifying a different aneuploidy level for each chromosome pair.
#' embryo <- Embryo(n.cells = 200, n.chrs = 3,  prop.aneuploid = c(0.2, 0.1, 0.4),
#'                  dispersal =  0.9)
Embryo <- function(n.cells = 200, n.chrs = 1, prop.aneuploid = 0.2, dispersal = 0.1,
                   concordance = 1, embryo.size.fixed = T, embryo.size.sd = 5, euploidy = 2,
                   rng.seed = NULL) {
  set.seed(rng.seed)

  if (embryo.size.sd <= 0) {
    stop(paste0("Number of cells sd (", embryo.size.sd, ") must be greater than 0"))
  }

  if (n.cells <= 1) {
    stop(paste0("Number of cells (", n.cells, ") must be greater than 1"))
  }

  if (euploidy <= 0) {
    stop(paste0("Number of chromosomes considered euploid (", euploidy, ") must be greater than 0"))
  }

  # Calculate the number of cells to use if not fixed
  if (!embryo.size.fixed) {
    n.cells <- max(1, ceiling(rnorm(1,
      mean = n.cells,
      sd = embryo.size.sd
    )))
  }

  if (n.chrs < 1) {
    stop(paste0("Number of chromosome pairs (", n.chrs, ") must be greater than 0"))
  }

  if (any(prop.aneuploid < 0) | any(prop.aneuploid > 1)) {
    stop(paste0(
      "prop.aneuploid (", paste(prop.aneuploid, collapse = ", "),
      ") must be between 0 and 1 inclusive"
    ))
  }

  if (any(dispersal < 0) | any(dispersal > 1)) {
    stop(paste0(
      "dispersal (", paste(dispersal, collapse = ", "),
      ") must be between 0 and 1 inclusive"
    ))
  }

  if (concordance < 0 | concordance > 1) {
    stop(paste0("Concordance (", concordance, ") must be between 0 and 1 inclusive"))
  }

  if (n.chrs > 1 & length(prop.aneuploid) == 1) prop.aneuploid <- rep(prop.aneuploid, n.chrs)
  if (n.chrs > 1 & length(dispersal) == 1) dispersal <- rep(dispersal, n.chrs)

  if (n.chrs > 1 & length(prop.aneuploid) != n.chrs) {
    stop(paste0(
      "Length of prop.aneuploid (", length(prop.aneuploid),
      ") must match the number of chromosomes in n.chrs (", n.chrs, ")"
    ))
  }


  .N_NEIGHBOURS <- 6
  # Make a sphere of evenly spaced points using the Fibonacci lattice
  indices <- seq(0, n.cells - 1, 1) + 0.5
  phi <- acos(pmin(pmax(1 - 2 * indices / n.cells, -1.0), 1.0)) # constrain to avoid rounding errors
  theta <- pi * (1 + sqrt(5)) * indices

  x <- cos(theta) * sin(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(phi)
  d <- as.data.frame(cbind(x, y, z))
  n <- as.data.frame(cbind(x, y, z))

  # Create distance matrix for each point
  # Set the .N_NEIGHBOURS closest points to be neighbours
  for (i in 1:nrow(d)) {
    dist <- sqrt((d$x - x[i])**2 + (d$y - y[i])**2 + (d$z - z[i])**2) # distance between points
    d[[paste0("d", i)]] <- dist # create a column to store the distances
    # A point is a neighbour if it is not this point, and it is in the list of closest points
    n[[paste0("n", i)]] <- dist > 0 & dist <= max(head(sort(dist), n = .N_NEIGHBOURS + 1))
  }

  # Make a column for each chromosome
  ploidy <- data.frame(matrix(data = euploidy, nrow = n.cells, ncol = n.chrs))
  colnames(ploidy) <- paste0("chr", 1:n.chrs)

  # Set a cell to contain an aneuploid chromosome
  #
  # @param ploidy the embryo
  # @param cell.index the cell to affect
  # @param chromosome the chromosome to make aneuploid
  #
  # @return the modified ploidy table
  set.aneuploid <- function(ploidy, cell.index, chromosome) {
    if (chromosome < 1 | chromosome > n.chrs) {
      stop(paste0("Chromosome must be in range 1-", n.chrs))
    }
    ploidy[cell.index, chromosome] <- 0 # For now, we just model all aneuploidy as nullisomy
    return(ploidy)
  }

  # Test if any of the neighbouring cells have an aneuploidy
  # d - the distance matrix for cells
  # ploidy - the ploidy matrix for cells and chromosomes
  # index - the cell to test
  # chromosome - the chromosome to test
  # Returns true if any of the closest cells are aneuploid
  .has.adjacent.aneuploid <- function(n, ploidy, index, chromosome, euploidy) {
    adj.list <- n[[paste0("n", index)]]
    isAenu <- ploidy[, chromosome] != euploidy

    return(any(adj.list & isAenu))
  }

  # Test if the given chromosome in the given cell is aneuploid
  #
  # @param embryo the embryo
  # @param cell.index the cell to test (0 for all cells)
  # @param chromosome the chromosome to test
  #
  # @return if cell.index is >0, return true if the chromosome is aneuploid,
  # false otherwise. If cell.index is 0, return a boolean vector of all cells.
  is.aneuploid <- function(ploidy, cell.index, chromosome, euploidy) {
    if (chromosome < 1 | chromosome > n.chrs) {
      stop(paste0("Chromosome must be in range 1-", n.chrs))
    }

    if (cell.index > n.cells) {
      stop(paste0("Cell index must be in range 0-", n.cells))
    }

    # cat("Testing aneuploidy of cell", cell.index, "chr", chromosome, "\n")

    # Return a vector if all cells requested
    if (cell.index == 0) {
      return(ploidy[, chromosome] != euploidy)
    }

    # Otherwise return just the single cell value
    return(ploidy[cell.index, chromosome] != euploidy)
  }

  # Set aneuploidies in an embryo
  #
  # Aneuploid cells are either adjacent or dispersed
  #
  # @param ploidy the ploidy matrix from the embryo
  # @param chromosome the chromosome to set aneuploidies for
  # @param prop.aneuploid the proportion of aneuploid cells (0-1)
  # @param dispersion the dispersion of the aneuploid cells (0-1)
  # @param concordance the concordance between aneuploid cells for each chromosome (0-1).
  #
  # @return the ploidy matrix with aneuploidies
  set.aneuploidies <- function(ploidy, chromosome, prop.aneuploid, dispersion, concordance) {

    # Shortcut the easy cases
    if (prop.aneuploid == 0) {
      return(ploidy)
    }

    if (prop.aneuploid == 1) {
      for (i in 1:nrow(ploidy)) {
        ploidy <- set.aneuploid(ploidy, i, chromosome)
      }
      return(ploidy)
    }

    # We must have an integer value of at least one aneuploid cell
    n.aneuploid <- ceiling(max(1, n.cells * prop.aneuploid))
    # cat("Creating", n.aneuploid, "aneuploid cells for chr", chromosome,"\n")

    # Decide how many cells need to be concordant with the previous
    # chromosome (if we are above chromosome 1)
    n.concordant <- 0
    concordant.cells <- rep(F, nrow(ploidy))
    if (chromosome > 1) {
      prev.chr <- chromosome - 1 # to be used when chr>1 only
      concordant.cells <- is.aneuploid(ploidy, 0, prev.chr, euploidy)
      # cat("Prev chr", prev.chr, "has", length(concordant.cells[concordant.cells==T]), "aneuploid cells\n")
      n.concordant <- length(concordant.cells[concordant.cells == T]) * concordance
      # if there are more aneuploids in the prev chromosome, we can't match all
      n.concordant <- min(n.aneuploid, n.concordant)
      # cat("Expecting", n.concordant, "concordant cells with chr", prev.chr, "\n")
    }

    # The approach for dispersal is to set seed cells which will
    # grow into separate aneuploid patches. The more dispersion, the more
    # initial seeds.

    # Choose number of seeds for aneuploid regions
    n.seeds <- ceiling(max(1, n.aneuploid * dispersion))
    n.to.make <- n.seeds

    # Disperse seeds as much as possible initially
    # We create a list of possible seed locations, and as each seed is assigned
    # we remove it and its neighbours
    # If there are enough seeds required, we will exhaust the possible locations
    # and finish seeding in the next step
    # cat("Creating", n.to.make, "seeds\n")
    possible.initial.seeds <- 1:n.cells # vector of indexes we can sample
    while (n.to.make > 0 & length(possible.initial.seeds) > 0) {

      # Take the next seed from those available. NOTE: if the vector is length 1,
      # sample will sample from 1:n, so we need to correct for this
      seed <- ifelse(length(possible.initial.seeds) == 1, possible.initial.seeds[1], sample(possible.initial.seeds, 1))
      if (n.concordant > 0 & !concordant.cells[seed]) {
        possible.initial.seeds <- possible.initial.seeds[!possible.initial.seeds %in% c(seed)]
        next # skip non concordant cells
      }

      if (.has.adjacent.aneuploid(n, ploidy, seed, chromosome, euploidy)) {
        possible.initial.seeds <- possible.initial.seeds[!possible.initial.seeds %in% c(seed)]
        next # spread seeds out
      }

      ploidy <- set.aneuploid(ploidy, seed, chromosome)

      # Remove seed and its neighbours from the possible seed list
      possible.initial.seeds <- possible.initial.seeds[!possible.initial.seeds %in% c(seed, which(n[[paste0("n", seed)]]))]
      n.to.make <- n.to.make - 1L
      n.concordant <- max(0, n.concordant - 1L)
    }

    # When all dispersed seeds have been added, add the remaining seeds randomly
    while (n.to.make > 0) {
      seed <- sample.int(n.cells, 1)
      if (is.aneuploid(ploidy, seed, chromosome)) next
      if (n.concordant > 0 & !concordant.cells[seed]) next # skip non concordant cells
      ploidy <- set.aneuploid(ploidy, seed, chromosome)
      n.to.make <- n.to.make - 1L
      n.concordant <- max(0, n.concordant - 1L)
    }

    # Grow the seeds into neighbouring cells for remaining aneuploid cells
    n.to.make <- n.aneuploid - n.seeds

    # Handle any remaining concordant cells first - the placement rules don't apply
    while (n.to.make > 0) {
      seed <- sample.int(n.cells, 1)
      if (n.concordant > 0) {
        # cat("Placing ", n.concordant, "concordant aneuploid cells\n")
        if (!concordant.cells[seed]) next
        if (is.aneuploid(ploidy, seed, chromosome, euploidy)) next # skip cells already aneuploid
        ploidy <- set.aneuploid(ploidy, seed, chromosome)
        n.concordant <- max(0, n.concordant - 1L)
      } else {

        # Grow from an existing seed.
        if (is.aneuploid(ploidy, seed, chromosome, euploidy)) next # skip cells already aneuploid
        if (!.has.adjacent.aneuploid(n, ploidy, seed, chromosome, euploidy)) next # only grow next to existing aneuploid
        ploidy <- set.aneuploid(ploidy, seed, chromosome)
      }

      n.to.make <- n.to.make - 1
    }
    return(ploidy)
  }


  # Set the aneuploid cells in the ploidy matrix
  for (chr in 1:n.chrs) {
    ploidy <- set.aneuploidies(
      ploidy,
      chr,
      prop.aneuploid[chr],
      dispersal[chr],
      concordance
    )
  }


  new("Embryo",
    x = d[, 1], y = d[, 2], z = d[, 3],
    aneu = prop.aneuploid,
    disp = dispersal,
    dists = d[, -(1:3)],
    neighbours = n[, -(1:3)],
    euploidy = euploidy,
    ploidy = ploidy
  )
}


# Override show function for an Embryo object
setMethod("show", "Embryo", function(object) {
  cat("Embryo with ", length(object@x), " cells\n",
      ncol(object@ploidy), " chromosome pairs per cell\n",
      object@euploidy, " copies of each euploid chromosome per cell\n",
      "Aneuploidy: ", object@aneu, " dispersal ", object@disp, "\n",
      sep = ""
  )
})

# Override plot function for an Embryo object
setMethod("plot", "Embryo", function(x) {
  if (!require(magrittr)) stop("Package magrittr is not installed. Install using 'install.packages'")
  library(magrittr)

  colours <- factor(sapply(1:length(x), function(i) all(x@ploidy[i, ] == x@euploidy)),
    levels = c(T, F)
  )

  plotly::plot_ly(
    type = "scatter3d",
    mode = "markers",
    colors = c("#22FF22", "#222222"),
    color = colours,
    marker = list(
      size = 25,
      line = list(
        color = "#111111",
        width = 1
      )
    ),
    hoverinfo = "none"
  ) %>%
    plotly::add_trace(
      x = x@x,
      y = x@y,
      z = x@z,
      showlegend = F,
      hoverinfo = "skip"
    ) %>%
    plotly::add_annotations(
      text = "Click and drag to rotate",
      xref = "paper", yref = "paper",
      x = 0.0, xanchor = "left",
      y = 1.0, yanchor = "top",
      legendtitle = TRUE, showarrow = FALSE
    ) %>%
    plotly::config(displayModeBar = FALSE, scrollZoom = F) %>%
    plotly::layout(scene = list(
      xaxis = list(
        autorange = F,
        fixedrange = TRUE,
        showgrid = F,
        showline = F,
        showticklabels = F,
        showaxeslabels = F,
        title = "",
        zeroline = F,
        range = list(-1, 1)
      ),
      yaxis = list(
        fixedrange = TRUE,
        autorange = F,
        showgrid = F,
        showline = F,
        showaxeslabels = F,
        showticklabels = F,
        title = "",
        zeroline = F,
        range = list(-1, 1)
      ),
      zaxis = list(
        autorange = F,
        fixedrange = TRUE,
        showgrid = F,
        showline = F,
        showaxeslabels = F,
        showticklabels = F,
        zeroline = F,
        title = "",
        range = list(-1, 1)
      )
    ))
})

# Add a method override for existing generic function length to get the number of cells
setMethod("length", "Embryo", function(x) { length(x@x) } )


#' Take a sample from an embryo
#'
#' The cell at the given index is taken,
#' plus the closest n neighbouring cells where n = n.sampled.cells-1.
#'
#' @param embryo an embryo
#' @param biopsy.size the number of cells to biopsy
#' @param index.cell the index of the cell to begin biopsying. Must be a value
#'  between 1 and \code{nrow(embryo)}
#'  @param chromosome the chromosome to test
#'
#' @return the number of aneuploid cells in the biopsy
#'
#' @examples
#' e <- Embryo()
#' takeBiopsy(e, 5, 1)
setGeneric(name="takeBiopsy",
           def = function(embryo, ...) { standardGeneric("takeBiopsy")})


setMethod("takeBiopsy", signature = "Embryo", function(embryo, biopsy.size = 5,
                                                       index.cell = 1, chromosome = 0) {
  if (index.cell < 1 | index.cell > length(embryo@x)) {
    stop(paste("index.cell (", index.cell, ") must be between 1 and", length(embryo@x)))
  }

  if (chromosome < 0 | chromosome > ncol(embryo@ploidy)) {
    stop(paste("Chromosome (", chromosome, ") must be between 0 and", ncol(embryo@ploidy)))
  }

  # Get the distance list for the index cell
  sample.list <- embryo@dists[[paste0("d", index.cell)]]

  # Choose the cells to join the biopsy based on distance
  isSampled <- embryo@dists[[paste0("d", index.cell)]] <= max(head(sort(sample.list), n = biopsy.size))

  # count all chromsomes; don't care which chromosome is aneuploid
  # just is aneuploid or is not aneuploid
  if (chromosome == 0) {
    return(sum(embryo@ploidy[isSampled, ] != embryo@euploidy))
  }
  return(sum(embryo@ploidy[isSampled, chromosome] != embryo@euploidy))
})

#' Take all possible biopsies from an embryo
#'
#' Take a biopsy starting from each cell in turn of the given size from the given
#' embryo. The biopsy size is fixed by default; use the \code{n.cells.fixed} to choose biopsy
#' size with a normal distribution with mean biopsy.size and standard deviation specified
#' by \code{n.cells.sd}.
#'
#' @param embryo an embryo as created by \code{tessera::Embryo()}
#' @param biopsy.size the ideal number of cells to take in each biopsy
#' @param chromosome the chromosome to test, or 0 for all chromosomes
#' @param biopsy.size.fixed true to take the same number of cells in each biopsy, false to
#' use a distribution model
#' @param biopsy.size.sd the standard deviation of the normal distribution used to model
#' the cell biopsy size if \code{n.cells.fixed = F}.
#' @param calc.percent return the percentage of aneuploid cells in the biopsy, instead
#' of the absolute number. Note, this will use \code{biopsy.size} for the calculation
#' even if \code{biopsy.size.fixed=F}
#' @param summarise if true, summarise the biopsy data as a tibble of counts rather than a vector.
#' This parameter also works with \code{calc.percent}.
#'
#' @return an integer vector containing the number of aneuploid cells in each biopsy, or
#' if \code{summarise = T}, a tibble containing the number and percentage of biopsies
#' grouped by the number of aneuploid cells.
#'
#' Note that if you are simulating lots of embryos, it is more efficent to use
#'  \code{summarise = F} and perform aggregation later than to aggregate at this
#'  stage.
#' @export
#'
#' @examples
#' Create an embryo with default parameters
#' e <- Embryo()
#'
#' Biopsy size fixed at 5 cells
#' takeAllBiopsies(e, biopsy.size = 5, chromosome = 1)
#'
#' Biopsy size varies with a mean of 6 and sd of 1.5
#' takeAllBiopsies(e, biopsy.size = 6, n.cells.fixed = F, n.cells.sd = 1.5)
#'
#' Calculate percentage aneuploidy in each biopsy instead of absolute number of cells
#' takeAllBiopsies(e, biopsy.size = 5, calc.percent = T)
#'
#' Calculate a summary tibble instead of absolute counts
#' takeAllBiopsies(e, biopsy.size = 5, summarise = T)
setGeneric(name="takeAllBiopsies",
           def = function(embryo, ...) { standardGeneric("takeAllBiopsies")})


setMethod("takeAllBiopsies",
  signature = "Embryo",
  function(embryo, biopsy.size = 5,
           chromosome = 0, biopsy.size.fixed = T, biopsy.size.sd = 1, calc.percent = F, summarise = F) {
    if (chromosome < 0 | chromosome > ncol(embryo@ploidy)) {
      stop(paste("Chromosome (", chromosome, ") must be between 0 and", ncol(embryo@ploidy)))
    }


    #' Model the number of biopsied cells in a sample.
    #'
    #' When biopsying cells, we may not get exactly the target number; there
    #' may be one too many or too few. We model the number of cells to take in
    #' a biopsy as a normal distribution with a mean around the desired number of
    #' cells and a standard deviation provided.
    create.n.cells.function <- function() {
      if (biopsy.size.fixed) {
        # If we are keeping a fixed number of cells in each biopsy, we don't need a
        # model
        return(function() {
          biopsy.size
        })
      } else {
        # model the number of biopsied cells as a distribution
        # Ensure sd is at least 1
        return(function() {
          max(1, ceiling(rnorm(1,
            mean = biopsy.size,
            sd = max(1, biopsy.size.sd)
          )))
        })
      }
    }

    fn <- create.n.cells.function()


    # If just one chromosome sampled
    result <- c()
    for (i in 1:nrow(embryo@ploidy)) { # sample each cell in turn, so we get every cell
      f <- takeBiopsy(embryo,
        biopsy.size = fn(),
        index.cell = i, chromosome = chromosome
      )
      result <- c(result, f)
    }


    if (calc.percent) {
      result <- (result / biopsy.size) * 100
    }

    if (summarise) {
      cnames <- c("AneuploidCells" = "result")

      if (calc.percent) {
        cnames <- c("PctAneuploid" = "result")
      }

      result <- as.data.frame(result) %>%
        dplyr::group_by(result) %>%
        dplyr::rename(!!!cnames) %>%
        dplyr::summarise(Biopsies = dplyr::n()) %>%
        dplyr::mutate(PctBiopsies = (Biopsies / sum(Biopsies)) * 100)
    }

    return(result)
  }
)



#' Find neighbouring cell indexes
#'
#' From the given embryo, find the cell indexes that are neighbours of the given
#' cell. This will return an integer vector of neighbouring cells, excluding the
#' initially requested cell index.
#'
#' @param embryo an embryo as created by \code{Embryo()}
#' @param cell.index the cell whose neighbours to find; must be an integer between 1 and embryo size
#'
#' @return an integer vector of the cell indexes of the neighbouring cells
#' @export
#'
#' @examples
#' e <- Embryo(100, 1, 0.1, 0.1)
#' tessera::getNeighbouringCellIndexes(e, cell.index = 1)
setGeneric(name="getNeighbouringCellIndexes",
           def = function(embryo, cell.index, ...) { standardGeneric("getNeighbouringCellIndexes")})


setMethod("getNeighbouringCellIndexes",
  signature = "Embryo",
  function(embryo, cell.index) {
    if (cell.index < 1 | cell.index > nrow(embryo@ploidy)) {
      stop(paste("Cell index (", cell.index, ") must be between 1 and", nrow(embryo@ploidy)))
    }

    return(which(embryo@dists[[paste0("n", cell.index)]]))
  }
)


#' count aneuploid cells in an embryo
#'
#' From the given embryo, count the number of aneuploid cells
#'
#' @param embryo an embryo as created by \code{Embryo()}
#'
#' @return an integer number of aneuploid cells
#' @export
#'
#' @examples
#' e <- Embryo(100, 1, 0.1, 0.1)
#' tessera::countAneuploid(e) # 10
setGeneric(name="countAneuploidCells", def = function(embryo, ...) {standardGeneric("countAneuploidCells")})

setMethod("countAneuploidCells",
  signature = "Embryo",
  function(embryo) {
    # Check all chromosomes for each cell match euploid value
    length(embryo) - sum(sapply(1:length(embryo), function(i) all(embryo@ploidy[i, ] == embryo@euploidy)))
  }
)
