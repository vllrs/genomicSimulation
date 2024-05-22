#' Create a new SimData object from data loaded from files
#'
#' \code{load.data} sets up the simulation object with initial data
#' from the provided files. Either a genotype file or a map file, or 
#' both, should be provided, as a set of marker effects cannot be loaded
#' without a reference set of markers from either of the former two
#' files. 
#' 
#' The list of markers tracked by the simulation is immutable after this
#' call. Therefore, the map file (or genotype file, if no map file 
#' is provided) should contain all markers you will wish to simulate. 
#' 
#' Additional founder genotypes, recombination maps, and marker effect
#' sets can be loaded after this call using @seealso \link{load.genotypes}, 
#' @seealso \link{load.map},
#' @seealso \link{load.effects}
#' 
#' See package vignette for file formats that the package can load.
#'
#' @param genotype.file A string containing a filename. The file should
#' contain a matrix of markers and alleles
#' @param map.file A string containing a filename. The file should contain
#' a linkage map for the markers loaded from allele.file
#' @param effect.file A string containing a filename. The
#' file should contain effect values for calculating GEBVs of a trait
#' @return A list with entries $groupNum, $mapID, and $effectID, representing
#' the group number of the genotypes loaded from the genotype file, the 
#' recombination map identifier of the map loaded from the map.file, and 
#' the marker effect set identifier of the marker effects loaded from 
#' the effect.file respectively.
#'
#' @family loader functions
#' @export
load.data <- function(genotype.file=NULL, map.file=NULL, effect.file=NULL) {
 	sim.data$p <- .Call(SXP_load_data, genomicSimulation:::expand.path(genotype.file), 
 	                    genomicSimulation:::expand.path(map.file), 
 	                    genomicSimulation:::expand.path(effect.file))
	gn <- ifelse(is.null(genotype.file), NA, 1L)
	mi <- ifelse(is.null(map.file) || is.null(gn), NA, 1L)
	ei <- ifelse(is.null(effect.file) || is.null(mi), NA, 1L)
	return(list(groupNum=gn,mapID=mi,effectID=ei))
}

#' Load genotypes to the existing SimData object from a file
#'
#' \code{load.genotypes} returns the group number of the group
#' that the new genotypes were loaded into. Only the alleles at markers
#' already tracked by the SimData are saved.
#'
#' @inheritParams load.data
#' @return The group number of the genotypes loaded from allele.file.
#'
#' @family loader functions
#' @export
load.genotypes <- function(allele.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_load_genotypes, sim.data$p, genomicSimulation:::expand.path(allele.file))) 
}

#' OLD NAME | Load more genotypes to the existing SimData object from a file
#' 
#' ! This is the old name for \code{load.genotypes}. From genomicSimulation v0.2.5,
#' \code{load.genotypes} is the recommended name over \code{load.more.genotypes}. 
#' \code{load.more.genotypes} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{load.genotypes}
#' 
#' @keywords internal 
#' @export
load.more.genotypes <- function(allele.file) {
  return(load.genotypes(allele.file))
}

#' Load an additional recombination map from a file
#'
#' \code{load.map} returns the identifier of the new recombination
#' map. Only the markers already tracked by the SimData are included
#' in this map.
#'
#' @inheritParams load.data
#' @return The identifier of the recombination map loaded from map.file
#'
#' @family loader functions
#' @export
load.map <- function(map.file) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_load_map, sim.data$p, genomicSimulation:::expand.path(map.file))) 
}

#' Add a new set of marker effect values
#'
#' Only the alleles at SNPs
#' already tracked by the SimData are saved. 
#'
#' @param effect.file A string containing a filename. The file should
#' contain effect values for calculating GEBV for a trait
#' @return The effect set ID. 1 for the first set loaded, 2 for the next, and so on.
#' This should be passed to breeding-value-calculating functions to tell them to 
#' use this set of marker effects.
#'
#' @family loader functions
#' @export
load.effects <- function(effect.file) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_load_effects, sim.data$p, genomicSimulation:::expand.path(effect.file))) 
}

#' OLD NAME | Add a new set of marker effect values
#' 
#' ! This is the old name for \code{load.effects}. From genomicSimulation v0.2.5,
#' \code{load.effects} is the recommended name over \code{load.different.effects}. 
#' \code{load.different.effects} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{load.effects}
#' 
#' @keywords internal 
#' @export
load.different.effects <- function(effect.file) {
  return(load.effects(effect.file))
}

#' Create a custom label
#'
#' \code{make.label} creates a new custom label with the given default value.
#'
#' @param default the default value of the new custom label. All current genotypes
#' and all genotypes generated in future will have this value for this label.
#' @return The label number of the new label created.
#'
#' @family label functions
#' @export
create.new.label <- function(default) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_create_new_label, sim.data$p, default))
}

#' OLD NAME | Create a custom label
#' 
#' ! This is the old name for \code{create.new.label}. From genomicSimulation v0.2.5,
#' \code{create.new.label} is the recommended name over \code{make.label}. 
#' \code{make.label} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{create.new.label}
#' 
#' @keywords internal 
#' @export
make.label <- function(default) {
  return(create.new.label(default))
}
