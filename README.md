# Tessera

This is a tool to visualise the impact of aneuploidy on embryo biopsies in pre-implantation genetic testing.

Tessera is a Shiny app that runs in RStudio. You will need both R and RStudio installed.

## Install
```
install.packages("devtools")
devtools::install_github("bmskinner/tessera")
```

## Run the web app

```
library(tessera)
runTessera()
```

## The parameters

- number of cells: the size of the embryo
- proportion: the proportion of cells in the embryo that are aneuploid. Values between 0 - 1.
- dispersal: low dispersal means the aneuploid cells are clumped, high means they are scattered. Values between 0 - 1.
- number of samples: the number of cells in a single biopsy

- model all chromosomes: if selected, 23 chromosome pairs will be modelled separately, rather than a single
assessment of whether the cell is aneuploid. This also reveals two new controls in the ui, detailed below.
- concordance: if concordance is 1, all chromosomes will be aneuploid in the same cells. The lower the concordance,
the more chance cells will differ in which chromosomes are aneuploid. Values between 0 - 1.
- chromosome to view: choose which chromosome to display in the chart, or set as 0 to see the overall number
of aneuploid chromosomes in each cell. Values between 0 - 23

## Compute values yourself

If you want to run simulations computationally, you can create embryos using the `create.embryo` function with desired parameters, and count the number of aneuploid cells in all possible biopsies with the `take.all.biopsies` function:

```
e <- create.embryo(n.cells = 100, 
                   prop.aneuploid = rep(0.1, 23), # vector with proportion for each chromosome
                   dispersal = rep(0.2, 23),      # vector with dispersal for each chromosome
                   concordance = 1)
take.all.biopsies(e, chromosome = 5)
# Output is a vector of the number of aneuploid cells from all possible biopsies
```

Note that since proportion and dispersal are vectors, you can give each chromosome a different value if you wish.
