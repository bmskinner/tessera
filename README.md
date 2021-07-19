# Tessera

This is a tool to visualise the impact of aneuploidy on embryo biopsies in pre-implantation genetic testing.

This is a Shiny app that runs in RStudio.

## Install
```
install.packages("devtools")
devtools::install_github("bmskinner/tessera")
```

## Run

```
library(tessera)
runTessera()
```

If you want to run simulations computationally, you can create embryos using the `create.embryo` function with desired parameters, and count the number of aneuploid cells in all possible biopsies with the `take.all.biopsies` function