#' genomicSimulation: Stochastic simulation of breeding programs
#'
#' Functions to simulate breeding programs on (user-provided) sets of genetic markers.
#' 
#'
#' @docType package
#' @name genomicSimulation
#' @useDynLib genomicSimulation, .registration = TRUE
NULL
#> NULL


sim.data <- new.env(parent = emptyenv())
sim.data$p <- NULL

.onUnload <- function (libpath) {
  library.dynam.unload("genomicSimulation", libpath)
}

