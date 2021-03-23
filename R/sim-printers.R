#' Save the details of the genome that the SimData uses.
#'
#' \code{save.genome.model} saves the SNP names, their linkage map positions, 
#' and the effect values for calculating GEBVs (if applicable) to a file.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_simdata
#' @export
save.genome.model <- function(filename) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_simdata, sim.data$p, filename))
}

#' Save genotypes/lines currently tracked by the simulation.
#'
#' \code{save.genotypes} saves the alleles at each SNP of the
#' lines stored in the SimData to a file.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group. Non-integers and negatives raise an error. Nonexistent 
#' groups result in empty files.
#' @param type The printing format. Use a string starting with 'R' or 'r' 
#' to save in regular format (SNPs as columns, lines as rows). Use a string
#' starting with 'T' or 't' to save in transposed format (SNPs as rows, lines
#' as columns). Both formats produce tab-separated matrices. The regular
#' format may be slightly faster because reads from the SimData are more
#' contiguous.
#' @return 0 on success. On failure an error will be raised.
#' 
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_genotypes
#' @export
save.genotypes <- function(filename, group=NULL, type="R") {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_genotypes, sim.data$p, filename, group, type))
}

#' Generate and save count matrices for alleles in the SimData
#'
#' \code{save.allele.counts} counts the number of occurences of a given
#' allele at each SNP of each genotype in the selected set, and prints
#' the output to a file.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will count and print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group.
#' @param allele The first character of this string will be used as the 
#' allele to count. 
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_counts
#' @export
save.allele.counts <- function(filename, group=NULL, allele) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_counts, sim.data$p, filename, group, allele))
}

#' Save the pedigrees of genotypes currently saved to the SimData
#'
#' \code{save.pedigrees} saves a record of the ancestors of genotypes
#' for which this data was tracked. 
#'
#' Two formats are available for printing pedigrees. 
#'
#' The recursive format (type="R") recursively traces back and prints
#' all known parents in the ancestry of the genotype. It first prints the
#' ID of the genotype, followed by a tab, then recursively prints the pedigree
#' with format [name of genotype, or ID if no name exists]=([parent1],[parent2]).
#' (each of the [parent1] and [parent2] will follow a similar format if they
#' are not founder genotypes or other genotypes with unknown/untracked pedigrees.
#' One row in the produced file corresponds to the pedigree of one line. 
#' Genotypes produced by selfing or doubling haploids only show one parent 
#' (eg. F2=(F1), not F2=(F1,F1)).
#' 
#' Parents-only format (type="P") prints the line's name (or ID if it has none), 
#' then the two parents' names (or IDs if they have none) to a tab-separated file.
#' Even if both parents are the same this format prints two columns of parent names. 
#' Each row of the produced  file is the pedigree of a different genotype. 
#'
#' Note that the parents-only format prints only the most immediate generation of 
#' parents of a line. Also note that these 'immediate parents'
#' are the genotypes from which the ones to be printed were generated (i.e. 
#' the immediate parents of the result of a call to \code{\link{self.n.times}} are 
#' the group passed to the \code{\link{self.n.times}} function, not the (n-1)th 
#' generation of selfing).
#'
#' The file produced by this function when type="P" can be used as the 
#' parental.file parameter to \code{\link{find.crossovers}}.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group.
#' @param type The printing format. Use a string starting with 'R' or 'r' 
#' to save in recursive format. Use a string
#' starting with 'P' or 'p' to save in parents-only format.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_pedigrees
#' @export
save.pedigrees <- function(filename, group=NULL, type="R") {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_pedigrees, sim.data$p, filename, group, type))
}

#' Save the GEBVs of genotypes currently saved to the SimData.
#'
#' \code{save.GEBVs} calculates GEBVs for a set of genotypes using 
#' the current effect values saved to the SimData, then saves those GEBVs 
#' to a file.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_GEBVs
#' @export
save.GEBVs <- function(filename, group=NULL) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_GEBVs, sim.data$p, filename, group))
}

#' Save the local GEBVs of each block in each selected line's haplotypes to a file. 
#'
#' \code{save.local.GEBVs} calculates GEBVs for each block of markers listed in 
#' `block.file` for the maternal and paternal halves of its genotype. The results
#' are saved to a file. The output file is formatted with the blocks as unlabelled columns
#' and two rows for each genotype, named [line name]_1 and [line name]_2. The entries
#' in this matrix are the calculated local GEBVs/block effects.
#'
#' Note that this function was written to run on an HPC cluster and so aims for speed rather
#' than considering potential memory constraints. If the function crashes, it is likely
#' that the system on which it is being run is not giving it enough memory to make a matrix
#' of all lines and blocks. Try running the function on smaller groups of lines and concatenating
#' the output files, or contact the package maintainer to ask them to have another look at
#' this function and fix it up.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param block.file A string containing a filename from which the blocks should be read.
#' It should have five columns and a header row. The third column (block name) and 
#' fifth column (semicolon-separated marker names for the markers in the block) are 
#' the only columns that will be used. Designed for the output of the function 
#' `st.def.hblocks` from the package SelectionTools
#' @param group Save only lines that belong to this group.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @useDynLib genomicSimulation SXP_save_block_effects
#' @export
save.local.GEBVs <- function(filename, block.file, group=NULL) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_block_effects, sim.data$p, filename, block.file, group))
}