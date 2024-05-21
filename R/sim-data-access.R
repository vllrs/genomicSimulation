#' Get data of all individuals that belong to a group.
#'
#' \code{see.group.data} allows you to extract the names, IDs, indexes, 
#' genotypes (as strings), breeding values, or parentage of group members as an R vector.
#' As long as the group is not modified in between calls, the same ordering of group
#' members will be used, so this function can be used to construct tables of data, eg.
#' \code{data.frame("i"=see.group.data(3,"X"),"GEBV"=see.group.data(3,"BV"))}.
#'
#' @section Warning about indexes:
#' Deleting a group (\code{\link{delete.group}}) may cause the indexes of group members to 
#' change (as the positions/indexes that they are stored in the internal matrix may be
#' rearranged). Therefore, you must ensure that there is no call to \code{\link{delete.group}}
#' between loading indexes with this function and passing them to a function that takes indexes
#' as inputs (eg \code{\link{make.group}}).
#'
#' @param group an integer: the group number of the group to have a look at
#' @param data.type a string that will be used to identify the type of data
#' to be extracted. String parsing is case-insensitive.
#' \itemize{
#' \item If the first character is 'N', the 'N'ames of group members will be extracted. 
#' \item If the first character is 'D', the I'D's of group members will be extracted. 
#' \item If the first character is 'X', the Inde'X'es of group members will be extracted. 
#' \item If the first character is 'G', the 'G'enotypes of the group members will be extracted.
#' \item If the first character is 'C', the allele 'c'ounts of the group members for each marker
#' will be extracted (that is, 0 for no copies of `count.allele` at that marker,
#'1 for 1 copy, 2 for 2 copies, etc). A (group members)x(markers) matrix.
#' \item If the first character is 'B', the 'b'reeding values/GEBVs of group members will be returned.
#' The breeding values will be calculated according to the effect set whose identifier
#' is provided by the effect.set parameter (by default, the first set of marker effects loaded).
#' \item If it starts with "P1", the first parent of each group member is returned. As usual in 
#' genomicSimulation, if the parent has a name, its name will be what is returned, otherwise,
#' it will be its ID.
#' \item If it starts with "P2", the second parent of each group member is returned. As usual in 
#' genomicSimulation, if the parent has a name, its name will be what is returned, otherwise,
#' it will be its ID.
#' \item If it starts with "ped", the full pedigree log of each group member is returned.
#' \item If it starts with "L", the values of a custom 'l'abel will be returned. The values 
#' will be shown for the custom label whose identifier is provided by the label parameter
#' (by default, the first custom label created).
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
	return(.Call(SXP_see_group_data, sim.data$p, group, toupper(data.type), effect.set, label))
}

#' Get genotypes or allele counts of a group as a matrix.
#'
#' \code{see.group.gene.data} allows you to extract (markers)x(group members)
#' matrices of the genetics at each marker of some of the genotypes in simulation
#' For the moment only implemented in R, not C
#' 
#' @param group an integer: the group number of the group to have a look at
#' @param count.allele a character that will be used to identify the type of data
#' to be extracted.
#' If it is NA, then return the pair of alleles as a string eg "AT"
#' If it is a character, then return the count of that allele eg. 'A' -> 0 if 
#' no copies of A at that marker, 1 if heterozygous for A, 2 if "AA" at that marker.
#' 
#' @family grouping functions
#' @family data access functions
#' @export
see.group.gene.data <- function(group, count.allele=NA_character_) {
  if (is.null(sim.data$p)) { stop("Please load.data first.") }
  
  return(.Call(SXP_see_group_gene_data, sim.data$p, group, count.allele))
} 


#' Get a list of the groups currently existing in the SimData
#'
#' \code{see.existing.groups} scans the saved data for groups that currently
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
  return(.Call(SXP_change_name_to_values, sim.data$p, values, group, startIndex))    
}


#' Change the visual representation of an allele
#'
#' genomicSimulation represents alleles as single-character symbols. 
#' \code{change.allele.symbol} swaps occurrences of the allele represented
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
  return(.Call(SXP_change_allele_symbol, sim.data$p, marker, from, to))  
} 