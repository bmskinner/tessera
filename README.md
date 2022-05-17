# Tessera

[![R version](https://img.shields.io/github/r-package/v/bmskinner/tessera)]()

Tessera is a tool to visualise the impact of aneuploidy on embryo biopsies in pre-implantation genetic testing.

It includes a Shiny app that can be run from base R or RStudio.

## Install
```
install.packages("devtools")
devtools::install_github("bmskinner/tessera")
```

## Run the web app

The Shiny app is available at [reproduction.essex.ac.uk](reproduction.essex.ac.uk), and is also included in the package for you to run locally if you prefer:

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
of aneuploid chromosomes in each cell.

## Compute values yourself

If you want to run simulations computationally, you can create embryos using the `Embryo` function with desired parameters, and count the number of aneuploid cells in all possible biopsies with the `takeAllBiopsies` function:

```

e <- Embryo(n.cells = 200, 
            n.chr   = 23,
            prop.aneuploid = 0.2,
            dispersal = 0.1,
            concordance = 1)
            
takeAllBiopsies(e, chromosome = 1, biopsy.size = 5)
# Output is a vector of the number of aneuploid cells from all possible biopsies
```

Note that since `prop.aneuploid` and `dispersal` are vectors, you can give each chromosome a different value if you wish.
