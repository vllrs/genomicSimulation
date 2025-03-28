#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names, IDs, indexes, 
#' genotypes (as strings), breeding values, or parentage of group members as an R vector.
#' As long as the group is not modified in between calls, accessing any datatype 
#' will return a vector using the same ordering of group members, 
#' so this function can be used to construct tables of data, eg.
#' \code{data.frame("i"=see.group.data(3,"X"),"GEBV"=see.group.data(3,"BV"))}.
#' 
#' Multiple groups can be passed to this function. The function will still 
#' return a single vector of data. Data from each group will be concatenated in 
#' the order that the groups were provided to the function. That is, a call to
#' \code{see.group.data(c(4,8,1),"N")} will return a single vector of names, 
#' containing the names of members of group 4, followed by the names of members 
#' of group 8, followed by the names of members of group 1.
#'
#' @section Warning about indexes:
#' Deleting a group (\code{\link{delete.group}}) may cause the indexes of group members to 
#' change (as the positions/indexes that they are stored in the internal matrix may be
#' rearranged). Therefore, you must ensure that there is no call to \code{\link{delete.group}}
#' between loading indexes with this function and passing them to a function that takes indexes
#' as inputs (eg \code{\link{make.group}}).
#'
#' @param group a vector of integers: the group number(s) of the group(s) to
#' access data from
#' @param data.type a string that will be used to identify the type of data
#' to be extracted. String parsing is case-insensitive.
#' \itemize{
#' \item If the first character is 'N', the 'N'ames of group members will be extracted. 
#' \item If the first character is 'D', the pedigree I'D's of group members will be extracted. 
#' \item If the first character is 'X', the Inde'X'es of group members will be extracted. 
#' \item If the first character is 'G', the 'G'enotypes of the group members will be extracted.
#' \item If the first character is 'B', the 'b'reeding values/GEBVs of group members will be returned.
#' \item If it starts with "L", the values of a custom 'l'abel will be returned. The values 
#' will be shown for the custom label whose identifier is provided by the label parameter
#' (by default, the first custom label created).
#' The breeding values will be calculated according to the effect set whose identifier
#' is provided by the effect.set parameter (by default, the first set of marker effects loaded).
#' \item If it starts with "P1" (but has no third letter or the third letter is not D), 
#' the name of the first parent of each group member is returned.
#' \item If it starts with "P1D", the pedigree ID of the first parent of each 
#' group member is returned.
#' \item If it starts with "P2" (but has no third letter or the third letter is not D), 
#' the name of the second parent of each group member is returned.
#' \item If it starts with "P2D", the pedigree ID of the second parent of each 
#' group member is returned.
#' \item If it starts with "ped", the full pedigree log of each group member is returned.
#' }
#' @param effect.set identifier for the set of marker effects to be used to calculate breeding values.
#' This parameter is only used if the data.type corresponds to a request for breeding values/GEBVs.
#' @param label identifier for the custom label to be viewed. This parameter is only used 
#' if the data.type corresponds to a request for viewing custom label values.
#' @return A vector containing the desired data output for each member of the group.
#'
#' @family grouping functions
#' @family data access functions
#' @export
see.group.data <- function(group, data.type, effect.set=1L, label=1L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (!is.integer(group)) {
    tmp <- group
    group <- as.integer(group)
    if (!isTRUE(all(tmp==group))) { stop("Group identifiers must be integers.") }
  }
	return(.Call(SXP_see_group_data, sim.data$p, group, toupper(data.type), effect.set, label))
}

#' Get genotypes or allele counts of a group as a matrix.
#'
#' \code{see.group.gene.data} allows you to extract a (markers)x(group members)
#' matrix of genotypes in simulation.
#' 
#' Multiple group numbers can be passed to this function. The function will still 
#' return a single matrix. Column representing all group members from the first
#' group will be first in the matrix, followed by columns representing all group
#' members in the second group, and so on. Groups are concatenated in the same order
#' as they were ordered in the input parameter passed to this function.
#' 
#' @param group a vector of integers: the group number(s) of the group(s) to
#' access data from
#' @param count.allele a character that will be used to identify the type of data
#' to be extracted.
#' If it is NA, as it is by default, then return the pair of alleles as a string eg "AT"
#' If it is a character, then return the count of that allele eg. 'A' -> 0 if 
#' no copies of A at that marker, 1 if heterozygous for A, 2 if "AA" at that marker.
#' @param unknown.allele The single character that will be used to represent alleles with 
#' the value '\0' in the table, if @a count.allele is NA. It is unused if @a count.allele
#' is not NA. '\0' alleles represent genetic marker/candidate combinations that were 
#' missing data/failed to load data in the original input file. They propagate 
#' through the generations the same as non-missing alleles. 
#' 
#' @family grouping functions
#' @family data access functions
#' @export
see.group.gene.data <- function(group, count.allele=NA_character_, unknown.allele="-") {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (!is.integer(group)) {
    tmp <- group
    group <- as.integer(group)
    if (!isTRUE(all(tmp==group))) { stop("Group identifiers must be integers.") }
  }
  m <- .Call(SXP_see_group_gene_data, sim.data$p, group, count.allele, unknown.allele)
  colnames(m) <- suppressWarnings(see.group.data(group,"Names"))
  rownames(m) <- .Call(SXP_see_marker_names, sim.data$p)
  return(m)
} 


#' Get a list of the groups currently existing in the SimData
#'
#' Scans the saved data for groups that currently
#' have members and returns their group numbers and number of members.
#'
#' @return A dataframe containing two columns, the first, named "Group", 
#' being the group numbers of every group in the SimData that currently 
#' has members, the second, named "GroupSize", being 
#' the number of genotypes currently allocated to that group
#'
#' @family grouping functions
#' @export
see.existing.groups <- function() {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	d <- data.frame(.Call(SXP_see_existing_groups, sim.data$p))
	colnames(d) <- c("Group", "GroupSize")
	return(d)
}


#' View one of the genetic maps being used by the simulation
#'
#' Access the list of markers of one of the genetic maps loaded in 
#' the simulation. The row order of the returned dataframe will 
#' match the order of the markers as stored by the simulation (eg.
#' will match the marker order in genotypes from see.group.data(data.type="G")).
#' 
#' Note that genomicSimulation does not store chromosome names after initial
#' import. The "chr" column of the output contains the index of a chromosome
#' as stored in the simulation. 
#' 
#' Similarly, genomicSimulation does not store the marker positions exactly as 
#' they were in the input file. The first marker in a chromosome is designated 
#' to be at position 0, and the other markers are stored as their cM offset 
#' from the first marker's position. Chromosomes containing only one marker
#' do not track that marker's position, so that marker's position will be NaN 
#' in the output dataframe.
#'
#' @param mapID the identifier of the specific genetic map to view, or zero to 
#' view the first-loaded/default genetic map.
#' @return A dataframe containing three columns: the first, named "marker", 
#' containing marker names; the second, named "chr", containing chromosome numbers;
#' and the third, named "pos", containing positions in cM.
#'
#' @family data access functions
#' @export
see.genetic.map <- function(mapID=0L) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  m <- data.frame(.Call(SXP_see_map, sim.data$p,mapID))
  colnames(m) <- c("marker","chr","pos")
  return(m)
}


#' View one of the sets of marker effects being used by the simulation
#' 
#' The function 
#' will return a dataframe with three columns: 
#' the first, named "marker", containing marker names; 
#' the second, named "allele", containing an allele symbol;
#' and the third, named "eff", containing the contribution of that allele at
#' that marker to the estimated breeding value.
#' 
#' If the set of marker effects has centering values, then rows containing the
#' special value "(centre)" in the allele column will contain that marker's 
#' centring value.
#'
#' @param effectID the identifier of the specific marker effect set to view, or 
#' zero to view the first-loaded/default marker effect set.
#' @return A dataframe with columns "marker", "allele" and "eff"
#'
#' @family data access functions
#' @export
see.marker.effects <- function(effectID=0L) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  m <- data.frame(.Call(SXP_see_effects, sim.data$p, effectID))
  if (is.data.frame(m)) {
    colnames(m) <- c("marker","allele","eff")
  }
  return(m)
}


#' Set the names of genotypes in simulation memory
#'
#' Changes the names of members of group(s) or
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
#'
#' @family label functions
#' @export
change.names.to.values <- function(values, group=NA, startIndex=0) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  .Call(SXP_change_name_to_values, sim.data$p, values, group, startIndex)   
}


#' Change the visual representation of an allele
#'
#' genomicSimulation represents alleles as single-character symbols. 
#' This function swaps occurrences of the allele represented
#' by one character with a different character, in a particular marker 
#' or across all markers.
#' 
#' All genotypes present in simulation will be affected.
#' 
#' If @a to is an allele that already exists, all distinction between 
#' alleles @a from and @a to will be irreversibly lost.
#' 
#' @param from Character symbol that will be replaced. eg. "A"
#' @param to Character symbol that @a from will be replaced with. eg. "b"
#' @param marker Name of the genetic marker across which to perform this
#' replacement operation, or NA if the allele should be replaced
#' across all genetic markers. 
#' 
#' @family saving functions 
#' @export
change.allele.symbol <- function(from, to, marker=NA) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  .Call(SXP_change_allele_symbol, sim.data$p, marker, from, to) 
} 