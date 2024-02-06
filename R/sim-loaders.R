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
	if (is.null(effect.file)) {
		sim.data$p <- .Call(SXP_load_data, allele.file, map.file)
		#the group number of the first group is always 1
		return(c(groupNum=1L)) 
	} else {
		sim.data$p <- .Call(SXP_load_data_weff, allele.file, map.file, effect.file)
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
	return(.Call(SXP_load_more_genotypes, sim.data$p, allele.file)) 
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
load.different.effects <- function(effect.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_load_new_effects, sim.data$p, effect.file)) 
}

#' Clear the internal storage of all data
#'
#' \code{clear.simdata} frees and deletes all data stored in the
#' package's internal SimData struct (sim.data$p). Returns 0 on success.
#'
#' @family loader functions
#' @export
clear.simdata <- function() {
	if (is.null(sim.data$p)) { return(0L) }
	sim.data$p <- NULL
	gc()
	return(0L)
}



#' Set the names of genotypes in simulation memory
#'
#' \code{change.names.to.values} changes the names of members of group(s) or
#'genotypes in general to a sequence of values.
#' 
#' The function orders the genotypes selected by the `group` parameter by their
#' index, finds the `startIndex`th genotype (counting up from 0), then pastes
#' the vector `values` in order into the names of 
#' the genotypes from there on.
#' 
#' If the vector of values provided is longer than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' values are ignored. 
#' If the vector of values provided is shorter than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' genotype names are left unchanged. 
#'
#' @param values A vector of strings to copy into the genotypes' name slots.
#' @param group NA or 0 to apply names to all genotypes in simulation memory,
#' or a group number to apply the names to only members of that group.
#' @param startIndex The position within the list of genotypes identified by `group`
#' and sorted by their index that the first entry of `values` should be copied into.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.names.to.values <- function(values, group=NA, startIndex=0) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_name_values, sim.data$p, values, group, startIndex))    
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
make.label <- function(default) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_create_label, sim.data$p, default))
}


#' Change the default value of a custom label
#'
#' \code{change.label.default} changes the default (birth) value of a custom
#' label(s) to the given integer value(s).
#' 
#' If `defaults` is shorter than `labels` 
#'
#' @param labels the label(s) of which to change the default values.
#' @param defaults the new default value of the new custom label. All genotypes
#' generated in future will have this value for this label. If there are multiple
#' labels, then `defaults` should be a vector of the same length, with the new
#' defaults for each label in `labels` at corresponding positions.
#' @return 0 on success
#'
#' @family grouping functions
#' @export
change.label.default <- function(labels, defaults) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_default, sim.data$p, labels, defaults))
}


#' Set values of a custom label
#'
#' \code{change.label.to.values} changes the values of a custom label for all 
#' genotypes or for members of group(s) to a sequence of values.
#' 
#' The function orders the genotypes selected by the `group` parameter by their
#' index, finds the `startIndex`th genotype (counting up from 0), then copies
#' the vector `values` entry-by-entry into the appropriate label of 
#' the genotypes from there on.
#' 
#' If the vector of values provided is longer than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' values are ignored. 
#' If the vector of values provided is shorter than the number of genotypes after
#' the `startIndex`th that are selected by the `group` parameter, then trailing 
#' genotypes are left unchanged. 
#'
#' @param label the label number of the label to change.
#' @param values A vector of integer values to copy into the label. 
#' @param group NA or 0 to apply the values to all genotypes in simulation memory,
#' or a group number to apply the values to only members of that group.
#' @param startIndex The position within the list of genotypes identified by `group`
#' and sorted by their index that the first entry of `values` should be copied into.
#' The first position in the list is considered index 1.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.to.values <- function(label, values, group=NA, startIndex=1) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_values, sim.data$p, label, 
               values, group, startIndex))    
}


#' Set a flat value for a custom label
#'
#' \code{change.label.to} changes the value of a custom label for all 
#' genotypes or for members of group(s) to a certain value.
#' 
#' Multiple group input supported.
#'
#' @param label the label number of the label to set.
#' @param value an integer to which the label value of all genotypes selected 
#' by the `group` parameter will be changed.
#' @param group NA or 0 to change all genotypes in simulation memory,
#' or a group number or vector of group numbers to change only
#' members of those groups.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.to.this <- function(label, value, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_const, sim.data$p, label, value, group))  
}


#' Add/subtract from a custom label
#'
#' \code{change.label.by.amount} changes the values of a custom
#' label for all genotypes or for members of group(s) by a given increment amount, 
#' positive or negative.
#' 
#' Multiple group input supported.
#'
#' @param label the label number of the label to change.
#' @param amount the integer change to modify the label value for the chosen
#' genotypes. eg. 1 for incrementing label values by 1, -3 to subtract 3 from
#' all label values.
#' @param group NA or 0 to apply this increment to all genotypes in simulation memory,
#' or a group number or vector of group numbers to apply the incrememnt to only
#' members of those groups.
#' @return 0 on success
#'
#' @family label functions
#' @export
change.label.by.amount <- function(label, amount, group=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_change_label_amount, sim.data$p, label, amount, group))
}


#' Delete custom label(s)
#'
#' \code{delete.label} removes the given custom label(s) from all genotypes
#' and deletes them from simulation memory.
#'
#' @param labels an vector containing the label numbers of the labels to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family label functions
#' @family grouping functions
#' @export
delete.label <- function(labels) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_delete_label, sim.data$p, labels))
}


#' Delete marker effect set(s)
#'
#' \code{delete.effect.set} removes the given marker effect set(s) from simulation memory.
#'
#' @param effect.sets an vector containing the label numbers of the labels to be deleted
#' @return 0 on success. An error is raised on failure.
#'
#' @family grouping functions
#' @export
delete.effect.set <- function(effect.sets) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_delete_eff_set, sim.data$p, effect.sets))
}



#' Identify historical crossover events
#'
#' \code{find.crossovers} identifies alleles in genotypes that must
#' have come from one parent or another, and so estimates the sections
#' of the genotype that must have come from each parent.
#'
#' @param parentage.file A string containing a filename. The file should
#' contain tab-separated, the name of an line and the name of its parents.
#' Parents do not have to match the parents tracked by its pedigree, but instead
#' are the ones from which the recombination frequencies will be calculated.
#' @param out.file A string containing a filename. The results of the analysis 
#' will be saved here as tab-separated matrix with child lines as the columns,
#' SNPs as the rows, and either a 0 (for undefined ancestry) or the id of the 
#' parent from which its alleles were inherited as the value in each cell. 
#' @param window.size An odd integer representing the size of the window centered at
#' each SNP in the matrix over which to check which parent the alleles come from. 
#' @param certainty A boolean. If TRUE, only SNPs where we can tell which parent the
#' alleles came from are filled in. If FALSE, uncertain SNPs are filled with the id of
#' the parent most recently identified as having given an allele.
#' @return 0 on success. The result will be printed to the file of name out.file
#'
#' @seealso \code{\link{find.plot.crossovers}}, which performs the same calculations
#' but also produces a plot. 
#' @export
find.crossovers <- function(parentage.file, out.file, window.size=1, certainty=TRUE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_find_crossovers, sim.data$p, parentage.file, out.file, window.size, certainty))
}

#' Identify and plot historical crossover events
#'
#' \code{plot.crossovers} identifies alleles in genotypes that must
#' have come from one parent or another, and so estimates the sections
#' of the genotype that must have come from each parent. It plots the result using ggplot
#'
#' @inheritParams find.crossovers
#' @return 0 on success. The result will be printed to the file of name out.file
#'
#' @seealso \code{\link{find.crossovers}}, which does not produce the plot and so will
#' be faster. 
#' @export
find.plot.crossovers <- function(parentage.file, out.file, window.size=1, certainty=TRUE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	.Call(SXP_find_crossovers, sim.data$p, parentage.file, out.file, window.size, certainty)
	
	my.palette = c("#cccccc", "#0077AA", "#66ccee", "#228833", "#DDCC77", "#EE6677", "#aa3377", "yellow",
           'green3', 'orange')
	
	parent <- genome <- SNP <- NULL #because CHECK is worried about the use of these in ggplot
	
	#load the result matrix with a parent combination string associated with each line
	crossovers <- utils::read.table(out.file, sep='\t', header=T)
	families <- utils::read.table(parentage.file, sep='\t', header=F, col.names=c("X","P1","P2"))
	families$P <- sprintf("%s.%s", families$P1, families$P2)
	families <- families[,c("X","P")]
	crossovers <- merge(crossovers, families)
	
	#load the map
	map <- send.map()
	
	# Melt the crossovers and match cells in the matrix to genome linkage map positions
	crossovers <- crossovers[ , ! colnames(crossovers) == "P"]
	crossovers <- reshape2::melt(crossovers, id.vars="X")
	colnames(crossovers) <- c("genome", "SNP", "parent")
	crossovers$n.snp <- match(crossovers$SNP, map$SNP) 
	#assumes the map is ordered which it is since it was loaded from SimData
	crossovers <- merge(crossovers, map)
	
	#plot the results
	ggplot2::ggplot(crossovers, ggplot2::aes(fill=parent, y=genome, x=rep(1, times=length(SNP)))) +
	  ggplot2::facet_grid(~crossovers$chr, scales = 'free_x', space = 'free_x', switch = 'x') + #divide into chromosomes
	  ggplot2::geom_bar(stat='identity') +
	  ggplot2::theme_classic() + ggplot2::scale_x_continuous(expand = c(0, 0)) + #no bg and nicer chr divides
	  ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank()) +
	  ggplot2::labs(x ="Chromosome", y = "", fill="Origin") + #Titles
	  ggplot2::scale_fill_gradientn(colors=my.palette, guide = ggplot2::guide_legend(),
						breaks=sort(unique(crossovers$parent)))
}

send.map <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	m <- data.frame(.Call(SXP_send_map, sim.data$p))
	colnames(m) <- c("SNP","chr","pos")
	return(m)
}