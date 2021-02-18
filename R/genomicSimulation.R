#' genomicSimulation: Crossing and genomic selection simulation tool
#'
#' Functions to simulate breeding programs on SNP marker genomes, and 
#' functions to calculate or perform selection based on GEBVs from provided effect values.
#' 
#'
#' @docType package
#' @name genomicSimulation
#' @useDynLib genomicSimulation
NULL
#> NULL


sim.data <- new.env(parent = emptyenv())
sim.data$p <- NULL

.onUnload <- function (libpath) {
  library.dynam.unload("genomicSimulation", libpath)
}

