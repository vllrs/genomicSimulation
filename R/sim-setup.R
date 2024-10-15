#' Create a new SimData object from data loaded from files
#'
#' Sets up the simulation object with initial data
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
#' @param allele.file A string containing a filename. The file should
#' contain a matrix of markers and alleles
#' @param genotype.file An alternate parameter name for `allele.file`. `allele.file`
#' is preferred, for consistency with \link{load.genotypes}
#' @param map.file A string containing a filename. The file should contain
#' a linkage map for the markers loaded from allele.file
#' @param effect.file A string containing a filename. The
#' file should contain effect values for calculating GEBVs of a trait
#' @param format Output of \link{define.matrix.format.details} to manually specify
#' file format, or an empty list by default to detect file formats automatically
#' @return A list with entries $groupNum, $mapID, and $effectID, representing
#' the group number of the genotypes loaded from the genotype file, the 
#' recombination map identifier of the map loaded from the map.file, and 
#' the marker effect set identifier of the marker effects loaded from 
#' the effect.file respectively.
#'
#' @family loader functions
#' @export
load.data <- function(allele.file=NULL, map.file=NULL, effect.file=NULL, format=list(), genotype.file=NULL) {
  if (!is.null(genotype.file)) {
    if (!is.null(allele.file)) {
      stop("Cannot use both `allele.file` and `genotype.file` parameters in the same call.")
    } else {
      .Deprecated(msg="For consistency with load.genotypes, prefer parameter name `allele.file` over `genotype.file` in load.data.",
                  old="load.data(genotype.file, map.file, effect.file)", new="load.data(allele.file, map.file, effect.file)")
      allele.file <- genotype.file
    }
  }
  
 	sim.data$p <- .Call(SXP_load_data, genomicSimulation:::expand.path(allele.file), 
 	                    genomicSimulation:::expand.path(map.file), 
 	                    genomicSimulation:::expand.path(effect.file),
 	                    format)
	gn <- ifelse(is.null(allele.file), NA, 1L)
	mi <- ifelse(is.null(map.file) || is.null(gn), NA, 1L)
	ei <- ifelse(is.null(effect.file) || is.null(mi), NA, 1L)
	return(list(groupNum=gn,mapID=mi,effectID=ei))
}

#' Load genotypes to the existing SimData object from a file
#'
#' Returns the group number of the group
#' that the new genotypes were loaded into. Only the alleles at markers
#' already tracked by the SimData are saved.
#' 
#' @details
#' The simplest genotype matrix file is formatted as follows:
#'   ```
#' name	G01	G02	G03	G04	G05	G06
#' m1	TT	TT	TT	TA	TT	AT
#' m3	TT	TT	TA	TA	TT	TT
#' m2	AA	AA	AA	AA	CC	AA
#' ```
#' where G01, G02, ..., are names of the founder genotypes; 
#' m1, m2, ..., are the genetic markers; and entries in the matrix are the 
#' alleles that the founder genotypes have at those markers. 
#' 
#' Other valid genotype matrix files might include:
#'   ```
#' m100, m101, m102
#' cand1,0,0,1
#' cand2,1,2,2
#' cand4,2,1,2
#' ```
#' or
#' ```
#' marker1	T/T	T/T	T/T	T/A	T/T	A/T
#' marker3	T/T	T/T	T/A	T/A	T/T	T/T
#' marker2	A/A	A/A	A/A	A/A	T/T	A/A
#' ```
#' 
#' The genotype matrix can be row-major or column-major (that is, the genetic markers 
#' may be rows, or columns). The program will determine the orientation by attempting
#' to match row and column headers with names of markers tracked by the simulation 
#' (extracted from a genetic map file). If the simulation does not yet have
#' a list of tracked markers (that is, no genetic map has been simultaneously or
#' previously provided), then it defaults to assuming that rows represent genetic markers.  
#' 
#' The order in which genetic markers are presented in the file does not matter. 
#' Candidate genotypes, however, will be saved internally in the simulation in 
#' the order that they appear in the file. Candidate genotype names do not need
#' to be unique. Candidate genotype names are optional, if candidates are stored 
#' as columns. genomicSimulation cannot read a gentoype matrix with no row headers,
#' so if candidates are stored as rows, then they are required to have names.
#' 
#' Genetic markers' names are required. 
#' 
#' When there are both row and column headers, the first/corner cell (the one 
#' that contained the value "name" in the first example above) can be deleted 
#' or can contain any text. Its value will be ignored. 
#' 
#' The table cells in the genotype matrix file may be separated by spaces, tabs,
#' commas, or any combination thereof. Cell spacers do not need to be consistent across the file. 
#' 
#' There are several options for how the alleles at each marker can be presented. 
#' All allele pair cells in the same genotype matrix must be in the same format. 
#' The format of the allele pairs in the genotype matrix is automatically detected. 
#' There are four acceptable formats for allele pairs in the genotype matrix:
#' 
#' 1. Any pair of characters (eg. "AA", "TA", "nW"). Each character is an allele. 
#'    This is a format that specifies allele phase (ie. "AT" and "TA" are different genotypes).
#' 2. Any pair of characters, separated by a forwards slash "/" character (eg. "A/A", "T/A", "n/W"). 
#'    The two characters either side of the slash are the alleles. 
#'    This is a format that specifies allele phase (ie. "AT" and "TA" are different genotypes).
#' 3. Alternate allele counts ("0", "1", "2"). This is a format that does not 
#'    specify allele phase, so phase of heterozygotes will be randomised when loaded. 
#'    The counts represent the number of copies of the alternate allele (stored inside
#'    genomicSimulation as "A", while the reference allele is stored as "T"). "0" 
#'    then corresponds to "TT", "2" to "AA", and "1" will be randomised as either 
#'    "TA" or "AT". Corresponding marker effect files for calculating breeding 
#'    values must use allele "A" consistently to represent the alternate allele, 
#'    and "T" to represent the reference. (alternatively, `change.allele.symbol` 
#'    can be used after loading to change T and/or A to other symbols.)
#' 4. (subset of) IUPAC encodings of DNA bases ("G", "A", "T", "C", "R", "Y", "M", "K", "S", "W", "N").
#'    This is a format that does not specify allele phase, so phase of heterozygotes 
#'    will be randomised when loaded. See package vignette for more details.
#'
#' **Note you might have a genotype matrix of only "AA", "AT", and "TT". This 
#' uses "alternate allele counts"-style encoding (like format 3) but presents it 
#' in a format that looks like pairs of alleles (format 1).** genomicSimulation 
#' expects allele pair encodings to include haplotype phase, (that is, to have 
#' four possible values for genotypes of two alleles, not three: eg. "AA", "AT", 
#' "TA", and "TT" instead of just "AA", "AT", "TT"). 
#'
#' Two options for loading a dataset with non-phased "AA"/"AT"/"TT" allele pairs are:
#' * Use “haplotyping”/”haplotype phasing”/”haplotype inference” software to infer 
#' whether heterozygotes are "AT" or "TA", before loading into genomicSimulation.
#' * Find-and-replace "TT" with "0", "AT" with "1", and "AA" with "2" before loading 
#' into genomicSimulation. genomicSimulation will then randomise the phase of each haplotype.
#'
#' @inheritParams load.data
#' @return The group number of the genotypes loaded from allele.file.
#'
#' @family loader functions
#' @export
#' 
#' @md
load.genotypes <- function(allele.file, format=list()) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_load_genotypes, sim.data$p, genomicSimulation:::expand.path(allele.file), format)) 
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
#' Returns the identifier of the new recombination
#' map. Only the markers already tracked by the SimData are included
#' in this map.
#' 
#' The simplest genetic map file is formatted as follows:
#' ```
#' marker chr pos
#' m3 3 15
#' m2 1 8.3
#' m1 1 5.2
#' ```
#' Other valid genetic map files might include:
#'   ```
#' chr marker pos
#' 1A 1243509 173.2
#' 1A 2350898 462.2
#' 1B 4360992 32.009
#' 2A 1243556 243.5
#' ```
#' or
#' ```
#' gene 10 3.24
#' othergene 10 8.3e-1
#' etc 15 1.203e2
#' ```
#' 
#' The header line is optional. If there is a header line, the three columns 
#' "marker" "chr" and "pos" may be rearranged into any ordering. If the header is 
#' not provided, the order is assumed to be marker name in the first column, 
#' followed by chromosome in the second column, followed by position in the third.
#' 
#' Regarding marker names: For all genomicSimulation features to work, genetic 
#' markers must have names. There is no issue with purely numeric names.
#' 
#' Regarding chromosomes/linkage group: any alphanumeric combination may be used 
#' to denote a chromosome/linkage group, for example '9' or '1A'. There is no 
#' limit on the number of unique chromosomes/linkage groups in the map.
#' 
#' Regarding marker position: this represents the position in centimorgans of 
#' each marker along its linkage group. (Distance of 1 cM = expected probability
#' of 0.01 chromosomal crossovers in that range.)
#' 
#' The cells in the map file may be separated by spaces, tabs, commas, or any 
#' combination thereof. Cell spacers do not need to be consistent across the file, 
#' but therefore marker names/linkage group names cannot contain spaces, tabs, or commas. 
#' 
#' The order in which genetic markers are presented in the file does not matter. 
#' If a marker is duplicated in the file, it will be recorded twice as two 
#' separate markers in the same position and with the same name.  However, 
#' because loading genotypes depends on searching for markers by name, if two 
#' markers have the same name, genotypes may not be loaded correctly. Therefore 
#' it is suggested that all markers have unique names.
#'
#' @inheritParams load.data
#' @return The identifier of the recombination map loaded from map.file
#'
#' @family loader functions
#' @export
#' 
#' @md
load.map <- function(map.file) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  return(.Call(SXP_load_map, sim.data$p, genomicSimulation:::expand.path(map.file))) 
}

#' Add a new set of marker effect values
#'
#' Only the alleles at SNPs
#' already tracked by the SimData are saved. 
#' 
#' Loading marker effect file(s) is optional for running the simulation. The simplest marker effect file is formatted as follows:
#'   ```
#' marker allele eff
#' m1 A -0.8
#' m2 A -0.1
#' m3 A 0.1
#' m1 T 0.9
#' m3 T -0.1
#' ```
#' Other valid marker effect files might include:
#' ```
#' marker eff allele
#' m1243509 0.1 A
#' m1243509 -0.1 T
#' m2350898 0.15 T
#' m2350898 -0.1 A
#' ```
#' or
#' ```
#' specialgene G 1.0
#' ```
#' 
#' A single marker effect file represents the allele effects for one trait. 
#' Multiple files like this can be loaded in order to calculate breeding values 
#' for multiple traits separately.
#' 
#' The header line is optional. If there is a header line, the three columns 
#' "marker" "allele" and "eff" may be rearranged into any ordering. If the header
#'  is not provided, the order is assumed to be marker name in the first column, 
#'  followed by allele, followed by additive effect value.
#' 
#' The order of the rows in the file does not matter. If multiple rows exist in 
#' the file for the same marker name/allele combination, only the last row's 
#' additive effect value will be saved. 
#' 
#' Regarding marker names: Only rows in this file whose marker name matches a 
#' marker name from a previously- or simultaneously-loaded map file or genotype 
#' matrix file will be loaded. 
#' 
#' Regarding allele: This column should be the allele (non-space character, eg "A")
#' that this effect value corresponds to. 
#' 
#' Regarding allele effect: This column should be a decimal representing the 
#' additive effect value of the same-row allele for the same-row marker. It can 
#' be positive or negative or zero, and can be represented by an integer, a 
#' decimal, or scientific notation. 
#' 
#' Cells may be separated by spaces, tabs, commas, or any combination thereof. 
#' Cell spacers do not need to be consistent across the file. 
#'
#' @param effect.file A string containing a filename. The file should
#' contain effect values for calculating GEBV for a trait
#' @return The effect set ID. 1 for the first set loaded, 2 for the next, and so on.
#' This should be passed to breeding-value-calculating functions to tell them to 
#' use this set of marker effects.
#'
#' @family loader functions
#' @export
#' 
#' @md
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
#' Creates a new custom label with the given default value.
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

#' Manually specify the format for a genotype matrix file that will be loaded
#' 
#' For use in generating the \code{format} parameter of
#' \link{load.genotypes} or \link{load.data}, when the genotype file 
#' being loaded is a genotype matrix.
#' 
#' It is not required that you manually specify all parameters, or none.
#' Whatever details of the file format are not manually specified
#' (that are not passed to this function or that are defined as NA) 
#' will be automatically detected by genomicSimulation when a genotype file is loaded
#' with this file format specification.
#'
#' @param has.header TRUE if the first line of the target file is a header row, 
#' FALSE if the first line of the target file is not a header row.
#' @param markers.as.rows TRUE if the rows of the genotype matrix represent genetic markers,
#' FALSE if the columns of the gentoype matrix represent genetic markers.
#' @param cell.style a string that will be used to identify the encoding of alleles in 
#' the body of the genotype matrix. String parsing is case-insensitive.
#' \itemize{
#' \item If the first character is "P", the target file stores alleles as ordered pairs of characters
#' (eg. "AA", "AT", "TA", "TT", "xy", "++")
#' \item If the first character is "/", the target file stores alleles as ordered pairs of 
#' characters separated by a forward slash character (eg. "A/A", "A/T", "T/A", "T/T", "x/y", "+/+")
#' \item If the first character is "C", the target file stores alleles as alternate allele counts
#' (ie. "0", "1", "2")
#' \item If the first character is "I", the target file stores alleles as IUPAC-encoded unordered pairs
#' (eg. "A"(=AA),"T"(=TT),"Y"(=TC or CT))
#' }
#' @returns a list suitable to be passed to the \code{format} parameter of \code{load.genotypes}
#' or \code{load.data}.
#' 
#' @export
define.matrix.format.details <- function(has.header=NULL,markers.as.rows=NULL,cell.style=NULL) {
  format <- list()
  
  if (is.logical(has.header)) {       format[["has.header"]] <- has.header 
  } else if (!is.null(has.header)) {
    warning("has.header should be logical") 
  }
  
  if (is.logical(markers.as.rows)) {  format[["markers.as.rows"]] <- markers.as.rows
  } else if (!is.null(markers.as.rows)) { 
    warning("markers.as.rows should be logical") 
  }
  
  if (is.character(cell.style)) {     format[["cell.style"]] <- toupper(cell.style) 
  } else if (!is.null(cell.style)) {
    warning("cell.style should be a string")
  }
  return(format)
}
