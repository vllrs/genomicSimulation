#' genomicSimulation: Stochastic simulation of breeding programs
#'
#' Functions to simulate breeding programs on (user-provided) sets of genetic markers.
#' 
#' Latest package vignette: https://vllrs.github.io/genomicSimulation/doc/gSvignette.html
#' 
#' Older overview of package (some details may be out of date): https://www.biorxiv.org/content/10.1101/2021.12.12.472291v2
#'
#' Package sources: https://github.com/vllrs/genomicSimulation
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

