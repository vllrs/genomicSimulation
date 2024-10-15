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

#' Save genotypes of lines in the simulation to a file
#'
#' Save the genotypes of current candidates to a file, as a 
#' tab-separated matrix of candidates x markers.
#' 
#' From version 0.2.6 of the package, the use of parameter `markers.as.rows` 
#' for defining the orientation of the matrix is preferred over use of the old 
#' parameter `type`.
#' 
#' Candidate lines are represented in the output file by their names, if they
#' have names, or by their pedigree IDs if they are unnamed.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group. Non-integers and negatives raise an error. Nonexistent 
#' groups result in empty files.
#' @param type The printing format. It is now preferred that you use `markers.as.rows` 
#' instead of `type`. For `type`, use a string starting with 'R' or 'r' 
#' to save in regular format (SNPs as columns, lines as rows). Use a string
#' starting with 'T' or 't' to save in transposed format (SNPs as rows, lines
#' as columns). Both formats produce tab-separated matrices. The regular
#' format may be slightly faster because reads from the SimData are more
#' contiguous.
#' @param markers.as.rows A logical representing the orientation of the matrix
#' in the output file. If TRUE, SNP markers will be rows of the matrix, and candidate 
#' lines will be columns. If FALSE, SNP markers will be columns of the matrix, and 
#' candidate lines will be rows.
#' @return 0 on success. On failure an error will be raised.
#' 
#' @family saving functions
#' @export
save.genotypes <- function(filename, group=NULL, type=NULL, markers.as.rows=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	
  if (!is.null(type)) {
	  if (is.logical(type)) { # in this case the third parameter is intended to be markers.as.rows
	    markers.as.rows <- type
	  } else {
	    .Deprecated(msg="Logical parameter `markers.as.rows` is now preferred over string parameter `type`.\ntype=\"R\" can be replaced with markers.as.rows=FALSE, and type=\"T\" with markers.as.rows=TRUE", 
	                old="save.genotypes(filename,group,type)", new="save.genotypes(filename,group,markers.as.rows)")
	    typeupper <- toupper(type)
	    if (startsWith(typeupper,"T")) {
	      markers.as.rows <- TRUE
	    } else if (startsWith(typeupper,"R")) {
	      markers.as.rows <- FALSE
	    } else {
	      stop("`type` parameter not recognised")
	    }
	  }
	}
  
  return(.Call(SXP_save_genotypes, sim.data$p, genomicSimulation:::expand.path(filename), 
	             group, markers.as.rows))
}

#' Generate and save count matrices of alleles in the simulation to a file
#'
#' Counts the number of occurrences of a given
#' allele at each SNP for each candidate in the simulation or for all 
#' candidates in a group. Prints
#' the output to a file in a tab-separated matrix.
#'
#' Candidate lines are represented in the output file by their names, if they
#' have names, or by their pedigree IDs if they are unnamed.
#'
#' @param filename A string containing a filename to which the output will
#' be written
#' @param group If not set/set to NULL, will count and print all genotypes.
#' Otherwise, if a group of that number exists, save only lines that belong
#' to that group.
#' @param allele The first character of this string will be used as the 
#' allele to count. 
#' @param markers.as.rows A logical representing the orientation of the matrix
#' in the output file. If TRUE, SNP markers will be rows of the matrix, and candidate 
#' lines will be columns. If FALSE, SNP markers will be columns of the matrix, and 
#' candidate lines will be rows.
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @export
save.allele.counts <- function(filename, group=NULL, allele, markers.as.rows=TRUE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_save_allele_counts, sim.data$p, genomicSimulation:::expand.path(filename), 
	             group, allele, markers.as.rows))
}

#' Save pedigrees from the simulation to a file
#'
#' Saves a record of the ancestors of genotypes
#' for which this data was tracked. 
#'
#' Two formats are available for printing pedigrees. 
#'
#' The recursive format (recursive.format=TRUE) recursively traces back and prints
#' all known parents in the ancestry of the genotype. It first prints the
#' ID of the genotype, followed by a tab, then recursively prints the pedigree
#' with format [name of genotype, or ID if no name exists]=([parent1],[parent2]).
#' (each of the [parent1] and [parent2] will follow a similar format if they
#' are not founder genotypes or other genotypes with unknown/untracked pedigrees.
#' One row in the produced file corresponds to the pedigree of one line. 
#' Genotypes produced by selfing or doubling haploids only show one parent 
#' (eg. F2=(F1), not F2=(F1,F1)). This is the default printing format.
#' 
#' Parents-only format (recursive.format=FALSE) prints the line's name (or ID if it has none), 
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
#' @param type The printing format. It is now preferred you use `recursive.format`
#' instead of `type`. For `type`, use a string starting with 'R' or 'r' 
#' to save in recursive format. Use a string
#' starting with 'P' or 'p' to save in parents-only format.
#' @param recursive.format A logical representing the format of the output pedigree.
#' If TRUE, print the recursive pedigree that contains all ancestors that genomicSimulation 
#' can link to a candidate in a nested-bracket formatted file. If FALSE, 
#' print the immediate parents of the candidate in a three-column file. 
#' @return 0 on success. On failure an error will be raised.
#'
#' @family saving functions
#' @export
save.pedigrees <- function(filename, group=NULL, type=NULL, recursive.format=TRUE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	
  if (!is.null(type)) {
    if (is.logical(type)) { # in this case the third parameter is intended to be recursive.format
      recursive.format <- type
    } else {
      .Deprecated(msg="Logical parameter `recursive.format` is now preferred over string parameter `type`.\ntype=\"R\" can be replaced with recursive.format=TRUE, and type=\"P\" with recursive.format=FALSE", 
                  old="save.pedigrees(filename,group,type)", new="save.pedigrees(filename,group,recursive.format)")
      typeupper <- toupper(type)
      if (startsWith(typeupper,"R")) {
        recursive.format <- TRUE
      } else if (startsWith(typeupper,"P")) {
        recursive.format <- FALSE
      } else {
        stop("`type` parameter not recognised")
      }
    }
  }
  
  return(.Call(SXP_save_pedigrees, sim.data$p, genomicSimulation:::expand.path(filename), group, recursive.format))
}

#' Calculate and save breeding values from the simulation to a file
#'
#' The function prints its results to a tab-separated file with three columns 
#' (the genotype's ID, its name, and its calculated breeding value/GEBV).
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