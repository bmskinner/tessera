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

## Compute values yourself

If you want to run simulations computationally, you can create embryos using the `create.embryo` function with desired parameters, and count the number of aneuploid cells in all possible biopsies with the `take.all.biopsies` function:

```
e <- create.embryo(100, 0.1, 0.2)
take.all.biopsies(e, 5)
# Output is a vector of the number of aneuploid cells from all possible biopsies
```
