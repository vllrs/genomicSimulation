#' Save the details of the genome that the SimData uses.
#'
#' **DEFUNCT**
#' > This function was never correctly implemented. 
#' An implementation may return at some future date.
#' Original documentation below.
#'
#' Saves the SNP names, their linkage map positions, 
#' and the effect values for calculating GEBVs (if applicable) to a file.
#'
#' The output is printed with first a tab-separated header row, containing
#' the column names "name", "chr", "pos", then all alleles for which the saved
#' SimData has effect values, if it has any. The subsequent lines, also tab-separated,
#' contain the name of a SNP, the chromosome on which it is found, the positions
#' at which it is found, and then the contribution to the GEBV of every allele
#' at that marker for which effects are loaded.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @return 0 on success. On failure an error will be raised.
#'
#' @rdname genomicSimulation-defunct
#' @family saving functions
#' @keywords internal
#' @export
save.genome.model <- function(filename) {
  .Defunct("save.genotypes",msg="Defunct. This function was never correctly implemented. An implementation may return at some future date")
  #if (is.null(sim.data$p)) { stop("Please load.data first.") }
  #return(.Call(SXP_save_simdata, sim.data$p, filename))
}

#' Save genotypes/lines currently tracked by the simulation.
#'
#' Saves the alleles at each SNP of the
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
#' @export
save.genotypes <- function(filename, group=NULL, type="R") {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_genotypes, sim.data$p, genomicSimulation:::expand.path(filename), group, type))
}

#' Generate and save count matrices for alleles in the SimData
#'
#' Counts the number of occurences of a given
#' allele at each SNP of each genotype in the selected set, and prints
#' the output to a file.
#'
#' The function prints its results to a tab-separated matrix file. The first line/
#' column header gives the name (or id if the genotype does not have a name allocated)
#' of each genotype that is having its counts printed. After that header row, all rows 
#' begin with the name of a SNP, followed by the number of the given allele at that
#' marker in the corresponding genotype from the header row.
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
#' @export
save.allele.counts <- function(filename, group=NULL, allele) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_allele_counts, sim.data$p, genomicSimulation:::expand.path(filename), group, allele))
}

#' Save the pedigrees of genotypes currently saved to the SimData
#'
#' Saves a record of the ancestors of genotypes
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
#' @export
save.pedigrees <- function(filename, group=NULL, type="R") {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_pedigrees, sim.data$p, genomicSimulation:::expand.path(filename), group, type))
}

#' Save the GEBVs of genotypes currently saved to the SimData.
#'
#' Calculates GEBVs for a set of genotypes using 
#' the current effect values saved to the SimData, then saves those GEBVs 
#' to a file.
#'
#' The function prints its results to a tab-separated file, where each
#' row contains, in order, the genotype's ID, name, and calculated GEBV
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group.
#' @param eff.set identifier for the set of marker effects with which to calculate
#' GEBVs.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @export
#' @aliases save.gebvs
save.GEBVs <- function(filename, group=NULL, eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_GEBVs, sim.data$p, genomicSimulation:::expand.path(filename), group, eff.set))
}

#' Save the local GEBVs of each block in each selected line's haplotypes to a file, using
#' blocks taken from a file.
#'
#' Calculates GEBVs for each block of markers listed in 
#' `block.file` for the maternal and paternal halves of its genotype. The results
#' are saved to a file. The output file is formatted with the blocks as unlabelled columns
#' and two rows for each genotype, named [line name]_1 and [line name]_2. The entries
#' in this matrix are the calculated local GEBVs/block effects.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param block.file A string containing a filename from which the blocks should be read.
#' It should have five columns and a header row. The third column (block name) and 
#' fifth column (semicolon-separated marker names for the markers in the block) are 
#' the only columns that will be used. Designed for the output of the function 
#' `st.def.hblocks` from the package SelectionTools
#' @param group Save only lines that belong to this group.
#' @param eff.set identifier for the set of marker effects with which to calculate 
#' local GEBVs.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @export
#' @aliases save.local.gebvs.blocks.from.file
save.local.GEBVs.blocks.from.file <- function(filename, block.file, group=NULL, eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_local_GEBVs_blocks_from_file, sim.data$p, genomicSimulation:::expand.path(filename), 
	             genomicSimulation:::expand.path(block.file), group, eff.set))
}

#' OLD NAME | Save the local GEBVs of each block in each selected line's 
#' haplotypes to a file, using blocks taken from a file.
#' 
#' ! This is the old name for \code{save.local.GEBVs.blocks.from.file}. From genomicSimulation v0.2.5,
#' \code{save.local.GEBVs.blocks.from.file} is the recommended name over \code{save.local.GEBVs.by.file}. 
#' \code{save.local.GEBVs.by.file} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{save.local.GEBVs.blocks.from.file}
#' 
#' @keywords internal 
#' @export
#' @aliases save.local.gebvs.by.file
save.local.GEBVs.by.file <- function(filename, block.file, group=NULL, eff.set=1L) {
  return(save.local.GEBVs.blocks.from.file(filename,block.file,group,eff.set))
}

#' Save the local GEBVs of each block in each selected line's haplotypes to a file,
#' using blocks created by slicing chromosomes into segments. 
#'
#' \code{save.local.GEBVs.blocks.from.chrsplit} calculates GEBVs for each block of markers created by
#' splitting each chromosome into `n.blocks.per.chr` same-length blocks for the maternal
#' and paternal halves of its genotype. The results are saved to a file. The output file 
#' is formatted with the blocks as unlabelled columns and two rows for each genotype, named 
#' [line name]_1 and [line name]_2. The entries in this matrix are the calculated local 
#' GEBVs/block effects.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param n.blocks.per.chr An integer containing the number of blocks each chromosome
#' will be divided into.
#' @param group Save only lines that belong to this group.
#' @param map identifier for the recombination map to use to calculate genetic 
#' distances with which to split the map. By default uses the oldest loaded map
#' currently active in simulation.
#' @param eff.set identifier for the set of marker effects with which to calculate 
#' local GEBVs. By default uses the first loaded set of effect values.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @export
#' @aliases save.local.gebvs.blocks.from.chrsplit
save.local.GEBVs.blocks.from.chrsplit <- function(filename, n.blocks.per.chr, group=NULL, map=0L, eff.set=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_local_GEBVs_blocks_from_chrsplit, sim.data$p, 
	             genomicSimulation:::expand.path(filename), n.blocks.per.chr, group, map, eff.set))
}

#' OLD NAME | Save the local GEBVs of each block in each selected line's haplotypes to a file,
#' using blocks created by slicing chromosomes into segments. 
#' 
#' ! This is the old name for \code{save.local.GEBVs.blocks.from.chrsplit}. From genomicSimulation v0.2.5,
#' \code{save.local.GEBVs.blocks.from.chrsplit} is the recommended name over \code{save.local.GEBVs.by.chr}. 
#' \code{save.local.GEBVs.by.chr} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{save.local.GEBVs.blocks.from.chrsplit}
#' 
#' @keywords internal 
#' @export
#' @aliases save.local.gebvs.by.chr
save.local.GEBVs.by.chr <- function(filename, n.blocks.per.chr, group=NULL, eff.set=1L) {
  return(save.local.GEBVs.blocks.from.chrsplit(filename,n.blocks.per.chr,group,eff.set))
}