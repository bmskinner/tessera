# Demo of how to create an S4 class
setClass("Embryo",
         slots = c(
           cell = "numeric",
           pos = "matrix",
           isAneuploid = "logical"
         ),

         prototype = list(
           cell = NA_real_,
           pos = NA,
           isAneuploid = F

         )
)

# Constructor is not part of the class definition
Embryo <- function(cellCount, x, y, z){
  new("Embryo", cell=cellCount, pos = matrix(cbind(x, y, z)), isAneuploid = F)
}
