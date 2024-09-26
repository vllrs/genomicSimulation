#' Get a string containing the ultimate/highest-scoring set of alleles.
#'
#' Allows you to extract the optimal haplotype 
#' according to the current loaded effect values as an R string. An error is 
#' raised if no effect values have been loaded. The optimal genotype (assuming
#' only additive allele effects) is just the doubled version of this haplotype.
#'
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' haplotypes and haplotype scores.
#' @return A string containing the allele out of available alleles at that
#' marker that has the highest effect value. The string will be ordered in
#' genome order (lowest chromosome and lowest position to highest) according
#' to the map that was included on initialisation.
#'
#' @family data access functions
#' @export
see.optimal.haplotype <- function(eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_see_optimal_haplotype, sim.data$p, eff.set))
}

#' Get a string containing the ultimate/highest-scoring set of alleles available
#' in the current group.
#' 
#' That is, consider the pool of alleles that exist in the group for each locus,
#' and take the highest-scoring allele at each locus.
#'
#' An error is raised if no effect values have been loaded.
#'
#' @param group an integer: the group number of the group to have a look at. 
#' Can be a vector of groups, in which case the calculation will be run independently
#' for each group.
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' optimal haplotype and haplotype breeding values
#' @return A string containing the allele out of all alleles in the group
#' that has the highest effect value at that locus. The string will be ordered in
#' genome order (lowest chromosome and lowest position to highest) according
#' to the map that was included on initialisation.
#'
#' @family data access functions
#' @export
see.optimal.possible.haplotype <- function(group, eff.set=1L) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (!is.integer(group)) {
    tmp <- group
    group <- as.integer(group)
    if (!isTRUE(all(tmp==group))) { stop("Group identifiers must be integers.") }
  }
  return(.Call(SXP_get_optimal_possible_haplotype, sim.data$p, group, eff.set))  
}

#' Get the ultimate/highest-possible GEBV given the current loaded effect values.
#'
#' Allows you to extract the optimal GEBV
#' according to the current loaded effect values. An error is 
#' raised if no effect values have been loaded.
#' 
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' maximal GEBV
#'
#' @return The highest possible GEBV
#'
#' @family data access functions
#' @export
#' @aliases see.optimal.gebv
see.optimal.GEBV <- function(eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_get_optimal_GEBV, sim.data$p, eff.set))
}

#' Get the ultimate/highest-possible GEBV given the pool of alleles available in
#' the current group.
#' 
#' Effectively, given the additive model of trait effects, this is 2 times the 
#' breeding value score of \code{see.optimal.possible.haplotype}. 
#'
#' An error is raised if no effect values have been loaded.
#'
#' @param group an integer: the group number of the group to have a look at.
#' Can be a vector of groups, in which case the calculation will be run independently
#' for each group.
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' optimal possible GEBV
#' @return The highest possible GEBV that can be created from alleles available
#' in the group.
#'
#' @family data access functions
#' @export
#' @aliases see.optimal.possible.gebv
see.optimal.possible.GEBV <- function(group, eff.set=1L) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (!is.integer(group)) {
    tmp <- group
    group <- as.integer(group)
    if (!isTRUE(all(tmp==group))) { stop("Group identifiers must be integers.") }
  }
  return(.Call(SXP_get_optimal_possible_GEBV, sim.data$p, group, eff.set))    
}

#' Get the lowest-possible GEBV given the current loaded effect values.
#'
#' Allows you to extract the lowest GEBV score
#' possible according to the current loaded effect values. An error is 
#' raised if no effect values have been loaded.
#' 
#' @param eff.set identifier for the set of marker effects to use to calculate the
#' minimum GEBV.
#'
#' @return The lowest possible GEBV
#'
#' @family data access functions
#' @export
#' @aliases see.minimal.gebv
see.minimal.GEBV <- function(eff.set=1L) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_get_minimal_GEBV, sim.data$p, eff.set))
}

#' OLD NAME | Get the lowest-possible GEBV given the current loaded effect values.
#' 
#' ! This is the old name for \code{see.minimal.GEBV}. From genomicSimulation v0.2.5,
#' \code{see.minimal.GEBV} is the recommended name over \code{see.minimum.GEBV}. 
#' \code{see.minimum.GEBV} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{see.minimal.GEBV}
#' 
#' @keywords internal 
#' @export
#' @aliases see.minimum.gebv
see.minimum.GEBV <- function(eff.set=1L) {
  return(see.minimal.GEBV(eff.set))
}


#' Identify historical crossover events
#'
#' \code{find.crossovers} identifies alleles in genotypes that must
#' have come from one parent or another, and so estimates the sections
#' of the genotype that must have come from each parent.
#' 
#' **Experimental function: may be changed or removed at any time, and validity of results is not guaranteed.**
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
#' **Experimental function: may be changed or removed at any time, and validity of results is not guaranteed.**
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

send.map <- function(mapID=0L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	m <- data.frame(.Call(SXP_send_map, sim.data$p,mapID))
	colnames(m) <- c("SNP","chr","pos")
	return(m)
}