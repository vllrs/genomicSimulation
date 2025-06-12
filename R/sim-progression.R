#' Performs random crosses between members of a group.
#'
#' Performs random crosses between members of a group. Selfing is not a possible operation.
#' Returns the group number of the group
#' that the new genotypes were loaded into. 
#'
#' For random crossing, the n.crosses parameter represents the number of 
#' random crosses to perform. Each cross will be repeated a number of times
#' determined by the offspring parameter. Each member of a group becomes 
#' a parent with equal probability.
#'
#' Random crosses can only be performed within a group, not within the
#' whole set of genotypes tracked by the simulation. 
#' If multiple groups are passed to this function, random crossing occurs 
#' independently within each of the input groups (i.e. there are no random crosses
#' where one parent is from one group, and the other parent is from a different group).
#'
#' @param group The group number from which to draw the parents. If a vector, 
#' the parameters will be used for each group number from the vector, independently,
#' therefore the output of the function will be a vector of the same length 
#' where each group number in the output is the offspring produced from the 
#' corresponding group in the input vector.
#' @param cap If nonzero, this is the maximum number of times each member of 
#' the group can be used as the parent of a cross.
#' Set to 0 for no restriction on the number of offspring per parent.
#' @param map The identifier for the recombination map with which gametes 
#' from members of \code{group} will be generated. By default uses the oldest 
#' loaded map currently active in simulation. If a vector of groups is provided,
#' this parameter can be a single map to be used for crossing operations on all groups, 
#' or this parameter can be a vector of the same length as `group`, in which case 
#' each map will be used for simulating the crossing operations on the group which
#' is at the corresponding position in the vector.
#' @param n.crosses The function will pick this many random pairs of parents
#' to cross.
#' @param offspring The number of times to cross each randomly-chosen pair. 
#' @param retain A logical value, repesenting whether to save the generated
#' genotypes to memory or discard them. You may wish to discard them
#' but save to file if you are generating too many crosses to save into
#' memory.
#' @param give.names A logical value representing whether or not to produce names
#' for the new genotypes generated. The names produced would have format [name.prefix][id]
#' @param name.prefix A string. If give.names is TRUE, the id is concatenated to 
#' this to produce the name of each new genotype.
#' @param track.pedigree A logical value representing whether or not to save the ids
#' of the parents of each new genotype to the new genotype's pedigree. If this is 
#' false, the new genotype's pedigree is unknown.
#' @param give.ids A logical value representing whether or not to allocate each new 
#' genotype an id. If this is FALSE, the new genotype is 'invisible' to pedigree trackers and
#' even if the pedigree of its offspring is supposedly tracked, the pedigree trackers
#' will not be able to identify the progenitors of its offspring. Furthermore, if it is
#' false and names are generated using give.names, all names generated in the same group will
#' be the same. Probably you'd only have this FALSE if you were discarding the results or worried
#' about id overflow.
#' @param file.prefix A string representing the prefix of files produced if save.pedigree=TRUE, 
#' save.gebv>0, or save.genotype=TRUE. If no file.prefix is provided, but at least one of 
#' those settings is TRUE, the default file prefix is "out".
#' save.pedigree, save.gebv, and save.genotype are known
#' as the "save-as-you-go" settings, which allow you to output data on the genotypes at the same 
#' time as they are being generated. The save-as-you-go settings can be used in conjunction with
#' retain=FALSE, if you want to generate results from a greater number of offspring than your 
#' system memory can hold. 
#' @param save.pedigree A logical value. If TRUE, saves the pedigrees of each
#' generated genotype to the file with filename "[file.prefix]-pedigree.txt". 
#' Pedigrees are saved progressively (1000 at a time). They are saved in recursive 
#' predigree format - see function \code{save.pedigrees} or the package vignette for more details.
#' @param save.gebv An integer representing a marker effect set identifier, or 0.
#' If FALSE, 0, or negative, does not save GEBVs. If positive and corresponding to the identifier
#' of a loaded set of marker effects, saves the GEBVs of each generated genotype to the
#' file with filename "[file.prefix]-bv.txt". This is a tab-separated text file. 
#' Breeding values for generated genotypes are calculated and saved progressively (1000 at a 
#' time). See function \code{save.GEBVs} for more details on the output format. A reminder that 
#' if only one set of marker effects have been loaded, the identifier for that set of 
#' marker effects is 1.
#' @param save.genotype A logical value. If TRUE, saves the SNP matrix of the generated 
#' genotypes to the file with filename "[file.prefix]-genotype.txt". 
#' Generated genotypes are saved progressively (1000 at a time). The output file 
#' uses genetic markers as columns and genotypes as rows. See \code{save.genotypes}
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.random.crosses <- function(group, n.crosses=5, cap=0, map=0L, offspring=1, retain=TRUE, 
		give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  group <- convert.to.integer(group,"Group identifiers")
  map <- convert.to.integer(map,"Map identifiers")
	return(.Call(SXP_make_random_crosses, sim.data$p, group, n.crosses, cap, map, give.names, 
	             name.prefix, offspring, track.pedigree, give.ids, my.expand.path(file.prefix), 
	             save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs random crosses between members of a group.
#' 
#' ! This is the old name for \code{make.random.crosses}. From genomicSimulation v0.2.5,
#' \code{make.random.crosses} is the recommended name over \code{cross.randomly}. 
#' \code{cross.randomly} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.random.crosses}
#' 
#' @keywords internal 
#' @export
cross.randomly <- function(group, n.crosses=5, cap=0, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.random.crosses")
  return(make.random.crosses(group=group,n.crosses=n.crosses,cap=cap,offspring=offspring,
                             retain=retain,give.names=give.names,name.prefix=name.prefix,
                             track.pedigree=track.pedigree,give.ids=give.ids,
                             file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs random crosses between two groups.
#'
#' Performs crosses where the first parent comes 
#' from one group and the second from another. It returns the group number of the group
#' that the new genotypes were loaded into.
#'
#' For random crossing, the n.crosses parameter represents the number of 
#' random crosses to perform. Each cross will be repeated a number of times
#' determined by the offspring parameter. Each member of a group becomes 
#' a parent with equal probability.
#' 
#' If the user only wants one parent to be randomly chosen, the flag `set.parent1`
#' can be set to TRUE, in which case `group1` will instead be interpreted as the 
#' index of an individual genotype, to which members of the other group will be crossed.
#' (and similarly the second parent can be set using `set.parent2`. A reminder 
#' that parent order affects nothing but pedigree printing, so pick whichever.)
#' 
#' Parameters set.parent1 and set.parent2 are deprecated and removed!
#' Use make.group and combine.groups to temporarily move an individual to their own group
#' if you wish to cross randomly from a group to an individual.
#'
#' @param group1 The group number from which to draw the first parent of the random 
#' crosses.
#' @param group2 The group number from which to draw the second parent of the random 
#' crosses.
#' @param cap1 If nonzero, this is the maximum number of times each member of 
#' group1 can be used as the (first) parent of a cross.
#' Set to 0 for no restriction on the number of offspring per parent.
#' @param cap2 If nonzero, this is the maximum number of times each member of 
#' group2 can be used as the (second) parent of a cross.
#' Set to 0 for no restriction on the number of offspring per parent.
#' @param map1 The identifier for the recombination map with which gametes 
#' from members of \code{group1} will be generated. By default uses the oldest 
#' loaded map currently active in simulation.
#' @param map2 The identifier for the recombination map with which gametes 
#' from members of \code{group2} will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @param map If \code{map1} and \code{map2} are not specified 
#' but \code{map} is, then the map specified by parameter \code{map} will be used for 
#' both map1 and map2.
#' @param cap If \code{cap1} and \code{cap2} are not specified 
#' but \code{cap} is, then the map specified by parameter \code{cap} will be used for 
#' both cap1 and cap2.
#' @inheritParams make.random.crosses
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.random.crosses.between <- function(group1, group2, n.crosses=5, cap1=0, cap2=0, 
	  map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE, map=0L, cap=0L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (map1[1] == 0L && map2[1] == 0L && map[1] != 0L) {
    map1 = map; map2 = map;
  }
  if (cap1[1] == 0L && cap2[1] == 0L && cap[1] != 0L) {
    cap1 = cap; cap2 = cap;
  }
  
	return(.Call(SXP_make_random_crosses_between, sim.data$p, as.integer(group1), 
	             as.integer(group2), cap1, cap2, map1, map2,
				 n.crosses, give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs random crosses between two groups.
#' 
#' ! This is the old name for \code{make.random.crosses.between}. From genomicSimulation v0.2.5,
#' \code{make.random.crosses.between} is the recommended name over \code{cross.randomly.between}. 
#' \code{cross.randomly.between} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.random.crosses.between}
#' 
#' @keywords internal 
#' @export
cross.randomly.between <- function(group1, group2, cap1=0, cap2=0, 
		n.crosses=5, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.random.crosses.between")
  return(make.random.crosses.between(group1=group1,group2=group2,cap1=cap1,cap2=cap2,n.crosses=n.crosses,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}


#' Performs defined crosses as passed in as R vectors.
#'
#' Performs specific crosses between defined pairs of parents. Returns the group number of the group
#' that the offspring genotypes were loaded into. 
#'
#' The offspring parameter represents the number of times each cross in the file
#' is carried out.
#'
#' Parents of each cross can be identified by name or by index. The first 
#' entry in the first.parents vector and the first entry in the second.parents 
#' vector are crossed to make the first offspring, then the next offspring comes 
#' from crossing the parent identified at first.parents[2] with the
#' parent identified at second.parents[2], etc.
#'
#' @param first.parents a vector identifying the first parent in each cross to
#' be carried out. Can be a vector of names or of indexes.
#' @param second.parents a vector identifying the second parent in each cross to 
#' be carried out. Can be a vector of names or of indexes.
#' @param map1 The identifier for the recombination map with which gametes 
#' from members of \code{group1} will be generated. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the second parent of all the targeted crosses.
#' @param map2 The identifier for the recombination map with which gametes 
#' from members of \code{group2} will be generated. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the second parent of all the targeted crosses.
#' @param map If \code{map1} and \code{map2} are not specified 
#' but \code{map} is, then the map specified by parameter \code{map} will be used for 
#' both map1 and map2.
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.targeted.crosses <- function(first.parents, second.parents, map1=0L, map2=0L,
		offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE, map=0L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (is.numeric(first.parents)) {
    first.parents <- convert.to.integer(first.parents,"Parent 1 indexes")
  } 
  if (is.numeric(second.parents)) {
    second.parents <- convert.to.integer(second.parents,"Parent 2 indexes")
  }
  if (map1[1] == 0L && map2[1] == 0L && map[1] != 0L) {
    map1 = map; map2 = map;
  }
	return(.Call(SXP_make_targeted_crosses, sim.data$p, first.parents, second.parents,
				 map1, map2, give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs defined crosses as passed in as R vectors.
#' 
#' ! This is the old name for \code{make.targeted.crosses}. From genomicSimulation v0.2.5,
#' \code{make.targeted.crosses} is the recommended name over \code{cross.combinations}. 
#' \code{cross.combinations} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.targeted.crosses}
#' 
#' @keywords internal 
#' @export
cross.combinations <- function(first.parents, second.parents,
		offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.targeted.crosses")
  return(make.targeted.crosses(first.parents=first.parents,second.parents=second.parents,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs defined crosses as laid out in a file.
#'
#' Performs specific crosses between defined pairs of parents taken from the input file. 
#' Returns the group number of the group
#' that the offspring genotypes were loaded into. 
#'
#' The offspring parameter represents the number of times each cross in the file
#' is carried out.
#'
#' The function searches for the parents of each cross by name, so parents must
#' have loaded names for this function to work. 
#'
#' @param cross.file a string containing a filename. The file should be available
#' to read and contain a tab-separated pair of names on each line. Each line
#' represents a cross to make.
#' @param map1 The identifier for the recombination map with which gametes 
#' will be generated from the first parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the first parent of all the targeted crosses.
#' @param map2 The identifier for the recombination map with which gametes 
#' will be generated from the second parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the second parent of all the targeted crosses.
#' @param map If \code{map1} and \code{map2} are not specified 
#' but \code{map} is, then the map specified by parameter \code{map} will be used for 
#' both map1 and map2.
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.crosses.from.file <- function(cross.file, map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE, map=0L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (map1[1] == 0L && map2[1] == 0L && map[1] != 0L) {
    map1 = map; map2 = map;
  }
	return(.Call(SXP_make_crosses_from_file, sim.data$p, cross.file, map1, map2, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs defined crosses as laid out in a file.
#' 
#' ! This is the old name for \code{make.crosses.from.file}. From genomicSimulation v0.2.5,
#' \code{make.crosses.from.file} is the recommended name over \code{cross.combinations.file}. 
#' \code{cross.combinations.file} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.crosses.from.file}
#' 
#' @keywords internal 
#' @export
cross.combinations.file <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.crosses.from.file")
  return(make.crosses.from.file(cross.file=cross.file,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=my.expand.path(file.prefix),save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs defined crosses between children of known parents as 
#' laid out in a file.
#'
#' Performs specific crosses between grandparents and their offspring 
#' to create genotypes with certain grandparents as defined in the input file. 
#' Returns the group number of the group
#' that the offspring genotypes were loaded into. 
#'
#' The function is designed for use in a situation where you wish to cross
#' unnamed lines or make crosses between simulated lines without having to check
#' what they were named. Each line of the file has four names, and the function
#' searches for an offspring of the first two that already exists in the simulation, 
#' and an offspring of the last two that already exists in the simulation,
#' then crosses these offspring. The grandparent generation, therefore, must have
#' names, so that they can be located in the pedigree.
#'
#' @param cross.file a string containing a filename. The file should be available
#' to read and contain a tab-separated quartet of names on each line. Each line
#' represents a cross to make.
#' @param map1 The identifier for the recombination map with which gametes 
#' will be generated from the first parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the first parent of all the targeted crosses.
#' @param map2 The identifier for the recombination map with which gametes 
#' will be generated from the second parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation. The same map will be used for simulating
#' the gametes of the second parent of all the targeted crosses.
#' @param map If \code{map1} and \code{map2} are not specified 
#' but \code{map} is, then the map specified by parameter \code{map} will be used for 
#' both map1 and map2.
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.double.crosses.from.file <- function(cross.file, map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE, map=0L) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (map1[1] == 0L && map2[1] == 0L && map[1] != 0L) {
    map1 = map; map2 = map;
  }
	return(.Call(SXP_make_double_crosses_from_file, sim.data$p, cross.file, map1, map2, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs defined crosses between children of known parents as 
#' laid out in a file.
#' 
#' ! This is the old name for \code{make.double.crosses.from.file}. From genomicSimulation v0.2.5,
#' \code{make.double.crosses.from.file} is the recommended name over \code{cross.dc.combinations.file}. 
#' \code{cross.dc.combinations.file} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.double.crosses.from.file}
#' 
#' @keywords internal 
#' @export
cross.dc.combinations.file <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.double.crosses.from.file")
  return(make.double.crosses.from.file(cross.file=cross.file,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs crosses between every line and every other line in a group
#' in one direction
#'
#' Simulates a cross for each pair in the half-diallel of parents in a group.
#' Returns the group number of the group
#' that the offspring genotypes were loaded into. 
#' 
#' If multiple groups are passed to this function, the half-diallele of each 
#' group is crossed and the different half-diallel outputs are returned as separate groups 
#' (i.e. there are no crosses made with one parent from one group and the other 
#' parent from another).
#'
#' The function performs one cross between each pair of genotypes in the
#' group. Order of parents is not considered because in the simulation method
#' this has no effect.
#'
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination of group members is crossed.
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.all.unidirectional.crosses <- function(group, map=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  group <- convert.to.integer(group,"Group identifiers")
  map <- convert.to.integer(map,"Map identifiers")
	return(.Call(SXP_make_all_unidirectional_crosses, sim.data$p, group, map, give.names, name.prefix,
	             offspring, track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' OLD NAME | Performs crosses between every line and every other line in a group
#' in one direction
#' 
#' ! This is the old name for \code{make.all.unidirectional.crosses}. From genomicSimulation v0.2.5,
#' \code{make.all.unidirectional.crosses} is the recommended name over \code{cross.all.pairs}. 
#' \code{cross.all.pairs} may become deprecated in future, when the package reaches 
#' stability.
#'
#' @seealso \link{make.all.unidirectional.crosses}
#' 
#' @keywords internal 
#' @export
cross.all.pairs <- function(group, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
  .Deprecated("make.all.unidirectional.crosses")
  return(make.all.unidirectional.crosses(group=group,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs n selfing steps on the lines in a group
#'
#' Crosses each genotype with itself repeatedly for a given number of generations.
#' Returns the group number of the group
#' that the selfed progeny genotypes were loaded into. 
#'
#' The function performs a selfing step (crosses each genotype
#' with itself), then performs a selfing step on those offspring, and repeats for a total of 
#' n generations of selfing. 
#'
#' The offspring parameter represents the number of genotypes produced from each 
#' genotype in the original group. These genotypes from the same family are 
#' independently descended from the original genotype.
#'
#' @param n An integer representing the number of steps of selfing to perform.
#' @inheritParams make.random.crosses
#' @param offspring This many offspring of each group member will be produced at the 
#' first selfing step. After that step, exactly one selfed offspring from each
#' parent will be progressed, no matter this value.
#' @return The group number of the new genotypes produced, or 0 if none could be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
self.n.times <- function(group, n, map=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  group <- convert.to.integer(group,"Group identifiers")
  map <- convert.to.integer(map,"Map identifiers")
	return(.Call(SXP_self_n_times, sim.data$p, group, n, map, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}			
			
#' Creates doubled haploids from each genotype in a group
#'
#' \code{make.doubled.haploids} returns the group number of the group
#' that the new genotypes were loaded into.
#'
#' The function generates a gamete from each genome in the original group. 
#' Output genomes are therefore perfectly homozygous.
#'
#' The offspring parameter represents the number of genotypes produced from each 
#' genotype in the original group. These genotypes from the same family are 
#' independently produced from the original genotype (i.e. a different gamete
#' to duplicate is generated for each of them)
#'
#' @inheritParams make.random.crosses
#' @param offspring The number of doubled haploids to make from each group member.
#' @return The group number of the new genotypes produced, or 0 if none could be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.doubled.haploids <- function(group, map=0L, offspring=1, retain=TRUE, 
		give.names=FALSE, name.prefix=NULL, track.pedigree=TRUE, give.ids=TRUE, 
		file.prefix="", save.pedigree=FALSE, save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  group <- convert.to.integer(group,"Group identifiers")
  map <- convert.to.integer(map,"Map identifiers")
	return(.Call(SXP_make_doubled_haploids, sim.data$p, group, map, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
}

#' Creates a genetically identical copy of each member of a group
#'
#' Copy the genotypes in a particular group.
#' Returns the group number of the group
#' that the clone genotypes were loaded into. 
#'
#' The function copies over the genotypes of the original group members exactly,
#' and applies metadata according to the optional parameters chosen.
#'
#' If track.pedigree is true, each clone is labelled as having one parent, the individual
#' of which it is a copy.
#'
#' The offspring parameter represents the number of clones produced from each 
#' genotype in the original group.
#'
#' If inherit.names is true, each clone shares the name of its progenitor, regardless
#' of the setting of give.names. If inherit.names is false, then clones are either named 
#' or left nameless according to the usual give.names/name.prefix settings.
#'
#' @inheritParams make.random.crosses
#' @param offspring The number of clones to make from each group member.
#' @param inherit.names A logical value, whether to give each clone the same name as the group
#' member from which it was created. Overrides give.names.
#' @return The group number of the new genotypes produced, or 0 if none could be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.clones <- function(group, offspring=1, retain=TRUE, inherit.names=TRUE, 
		give.names=FALSE, name.prefix=NULL, track.pedigree=TRUE, give.ids=TRUE, 
		file.prefix="", save.pedigree=FALSE, save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  group <- convert.to.integer(group,"Group identifiers")
	return(.Call(SXP_make_clones, sim.data$p, group, inherit.names, give.names, 
				 name.prefix, offspring, track.pedigree, give.ids, my.expand.path(file.prefix), save.pedigree,
				 save.gebv, save.genotype, retain))
}
