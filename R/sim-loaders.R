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
#' always 1 in the current implementation.
#'
#' @family loader functions
#' @useDynLib genomicSimulation load_data
#' @useDynLib genomicSimulation load_data_weff
#' @export
load.data <- function(allele.file, map.file, effect.file=NULL) {
	if (is.null(effect.file)) {
		sim.data$p <- .Call(load_data, allele.file, map.file)
	} else {
		sim.data$p <- .Call(load_data_weff, allele.file, map.file, effect.file)
	}
	#the group number of the first group is always 1
	return(1L) 
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
#' @useDynLib genomicSimulation load_more_genotypes
#' @export
load.more.genotypes <- function(allele.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(load_more_genotypes, sim.data$p, allele.file)) 
}

#' Replace effect values
#'
#' \code{load.different.effects} returns 0. Only the alleles at SNPs
#' already tracked by the SimData are saved. The old effect values,
#' if they exist, will be replaced with these new ones.
#'
#' @param effect.file A string containing a filename. The file should
#' contain effect values for calculating GEBV for a trait
#' @return 0 on success.
#'
#' @family loader functions
#' @useDynLib genomicSimulation load_new_effects
#' @export
load.different.effects <- function(effect.file) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(load_new_effects, sim.data$p, effect.file)) 
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
#' @useDynLib genomicSimulation find_crossovers
#' @export
find.crossovers <- function(parentage.file, out.file, window.size=1, certainty=TRUE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(find_crossovers, sim.data$p, parentage.file, out.file, window.size, certainty))
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
	.Call(find_crossovers, sim.data$p, parentage.file, out.file, window.size, certainty)
	
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

#' @useDynLib genomicSimulation send_map
send.map <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	m <- data.frame(.Call(send_map, sim.data$p))
	colnames(m) <- c("SNP","chr","pos")
	return(m)
}