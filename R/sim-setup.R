#' Create a new SimData object from data loaded from files
#'
#' \code{load.data} returns the group number of the group the starting
#' data set was loaded into. It also sets up the markers, linkage map,
#' and GEBV calculators.
#'
#' @param allele.file A string containing a filename. The file should
#' contain a matrix of markers and alleles
#' @param map.file A string containing a filename. The file should contain
#' a linkage map for the markers loaded from allele.file
#' @param effect.file (optional) A string containing a filename. The
#' file should contain effect values for calculating GEBVs of a trait
#' @return The group number of the genotypes loaded from allele.file. This is
#' always 1 in the current implementation. If an effect file was also loaded, the 
#' return value will be a list with two entries: groupNum, for the group number of the
#' first group loaded, and EffectID, for the identifier for the set of marker effects
#' loaded from the effect file.
#'
#' @family loader functions
#' @export
load.data <- function(allele.file, map.file, effect.file=NULL) {
	#the group number of the first group is always 1
	if (is.null(effect.file)) {
	  sim.data$p <- .Call(SXP_load_data, fs::path_expand(allele.file), 
	                      fs::path_expand(map.file), 
	                      NULL)
		return(c(groupNum=1L)) 
	} else {
	  sim.data$p <- .Call(SXP_load_data, fs::path_expand(allele.file), 
	                      fs::path_expand(map.file), 
	                      fs::path_expand(effect.file))
		return(list(groupNum=1L,effectID=1L))
	}
}

#' Load more genotypes to the existing SimData object from a file
#'
#' \code{load.more.genotypes} returns the group number of the group
#' that the new genotypes were loaded into. Only the alleles at SNPs
#' already tracked by the SimData are saved.
#'
#' @inheritParams load.data
#' @return The group number of the genotypes loaded from allele.file.
#'
#' @family loader functions
#' @export
load.more.genotypes <- function(allele.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_load_more_genotypes, sim.data$p, fs::path_expand(allele.file))) 
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
load.more.effects <- function(effect.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_load_more_effects, sim.data$p, fs::path_expand(effect.file))) 
}

#' OLD NAME | Add a new set of marker effect values
#' 
#' ! This is the old name for \code{load.more.effects}. From genomicSimulation v0.2.5,
#' \code{load.more.effects} is the recommended name over \code{load.different.effects}. 
#' \code{load.different.effects} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{load.more.effects}
#' 
#' @keywords internal 
#' @export
load.different.effects <- function(effect.file) {
  return(load.more.effects(effect.file))
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
