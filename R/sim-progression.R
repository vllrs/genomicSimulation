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
#'
#' @param group The group number from which to draw the parents of these
#' random crosses. If a vector, the parameters and breeding method 
#' (eg. random crossing) will be applied to each group number in the vector.
#' @param cap If nonzero, this is the maximum number of times each member of 
#' the group can be used as the parent of a cross.
#' Set to 0 for no restriction on the number of offspring per parent.
#' @param map The identifier for the recombination map with which gametes 
#' from members of @a group will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @param n.crosses The function will pick this many random pairs of parents
#' to cross.
#' @param offspring The number of times to cross each randomly-chosen pair. 
#' @param retain A boolean, repesenting whether to save the generated
#' genotypes to memory or discard them. You may wish to discard them
#' but save to file if you are generating too many crosses to save into
#' memory.
#' @param give.names A boolean representing whether or not to produce names
#' for the new genotypes generated. The names produced would have format [name.prefix][id]
#' @param name.prefix A string. If give.names is TRUE, the id is concatenated to 
#' this to produce the name of each new genotype.
#' @param track.pedigree A boolean representing whether or not to save the ids
#' of the parents of each new genotype to the new genotype's pedigree. If this is 
#' false, the new genotype's pedigree is unknown.
#' @param give.ids A boolean representing whether or not to allocate each new 
#' genotype an id. If this is FALSE, the new genotype is 'invisible' to pedigree trackers and
#' even if the pedigree of its offspring is supposedly tracked, the pedigree trackers
#' will not be able to identify the progenitors of its offspring. Furthermore, if it is
#' false and names are generated using give.names, all names generated in the same group will
#' be the same. Probably you'd only have this FALSE if you were discarding the results or worried
#' about id overflow.
#' @param file.prefix A string representing the prefix of files produced if save.pedigree=TRUE, 
#' save.gebv=TRUE, or save.genotype=TRUE.
#' @param save.pedigree A boolean. If TRUE, saves the pedigree in recursive format of each
#' generated genotype to the file with filename "[file.prefix]-pedigree". This is a text file. 
#' Generated genotypes are saved progressively (up to 1000 at a time), so if the 
#' full result of a crosser function call will not fit in memory, this setting can allow
#' you to still get results.
#' @param save.gebv A boolean. If TRUE, saves the GEBVs of each generated genotype to the
#' file with filename "[file.prefix]-eff". This is a tab-separated text file. 
#' Valuse for generated genotypes are calculated and saved progressively (up to 1000 at a 
#' time), so if the full result of a crosser function call will not fit in memory, this 
#' setting can allow you to still get results.
#' @param save.genotype A boolean. If TRUE, saves the SNP/line matrix in regular format
#' (generated genotypes as rows, SNPs as columns) to the file with filename 
#' "[file.prefix]-genome". This is a tab-separated text file. 
#' Generated genotypes are saved progressively (up to 1000 at a time), so if the 
#' full result of a crosser function call will not fit in memory, this setting can allow
#' you to still get results.
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
	return(.Call(SXP_make_random_crosses, sim.data$p, group, n.crosses, cap, map, give.names, 
	             name.prefix, offspring, track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), 
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
#' from members of @a group1 will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @param map2 The identifier for the recombination map with which gametes 
#' from members of @a group2 will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @inheritParams make.random.crosses
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.random.crosses.between <- function(group1, group2, n.crosses=5, cap1=0, cap2=0, 
	  map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  
	return(.Call(SXP_make_random_crosses_between, sim.data$p, as.integer(group1), 
	             as.integer(group2), cap1, cap2, map1, map2,
				 n.crosses, give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
#' from members of @a first.parents will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @param map2 The identifier for the recombination map with which gametes 
#' from members of @a second.parents will be generated. By default uses the oldest 
#' loaded map currently active in simulation
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.targeted.crosses <- function(first.parents, second.parents, map1=0L, map2=0L,
		offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  if (is.numeric(first.parents) && !is.integer(first.parents)) {
    numeric.fp <- first.parents
    first.parents <- as.integer(first.parents)
    if (!all.equal(numeric.fp,first.parents)) { stop("Parent 1 indexes must be integers.") }
  }
  if (is.numeric(second.parents) && !is.integer(second.parents)) {
    numeric.sp <- second.parents
    second.parents <- as.integer(second.parents)
    if (!all.equal(numeric.sp,second.parents)) { stop("Parent 2 indexes must be integers.") }
  }
	return(.Call(SXP_make_targeted_crosses, sim.data$p, first.parents, second.parents,
				 map1, map2, give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
#' loaded map currently active in simulation
#' @param map2 The identifier for the recombination map with which gametes 
#' will be generated from the second parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.crosses.from.file <- function(cross.file, map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_make_crosses_from_file, sim.data$p, cross.file, map1, map2, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
  return(make.crosses.from.file(cross.file=cross.file,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=genomicSimulation:::expand.path(file.prefix),save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
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
#' searches for an offspring of the first two, and an offspring of the last two,
#' then crosses these offspring. The grandparent generation, therefore, must have
#' names.
#'
#' @param cross.file a string containing a filename. The file should be available
#' to read and contain a tab-separated quartet of names on each line. Each line
#' represents a cross to make.
#' @param map1 The identifier for the recombination map with which gametes 
#' will be generated from the first parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation
#' @param map2 The identifier for the recombination map with which gametes 
#' will be generated from the second parent in each cross. By default uses the oldest 
#' loaded map currently active in simulation
#' @inheritParams make.random.crosses
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
make.double.crosses.from.file <- function(cross.file, map1=0L, map2=0L, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix="", save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_make_double_crosses_from_file, sim.data$p, cross.file, map1, map2, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
  return(make.double.crosses.from.file(cross.file=cross.file,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs crosses between every line and every other line in a group
#' in one direction
#'
#' Simulates a cross for each pair in the half-diallel of parents in a group.
#' Returns the group number of the group
#' that the offspring genotypes were loaded into. 
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
	return(.Call(SXP_make_all_unidirectional_crosses, sim.data$p, group, map, give.names, name.prefix,
	             offspring, track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
  return(make.all.unidirectional.crosses(group=group,offspring=offspring,retain=retain,give.names=give.names,name.prefix=name.prefix,track.pedigree=track.pedigree,give.ids=give.ids,file.prefix=file.prefix,save.pedigree=save.pedigree,save.gebv=save.gebv,save.genotype=save.genotype))		
}

#' Performs n selfing steps on the lines in a group
#'
#' Crosses each genotype with itself repeatedly for a given number of generations.
#' Returns the group number of the group
#' that the selfed progeny genotypes were loaded into. 
#'
#' The function performs a selfing step (crosses each genotype
#' with itself), then another on the output of that setp, etc., for a total of 
#' n steps of selfing. 
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
	return(.Call(SXP_self_n_times, sim.data$p, group, n, map, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
	return(.Call(SXP_make_doubled_haploids, sim.data$p, group, map, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree, save.gebv, save.genotype, retain))
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
#' @param inherit.names Boolean, whether to give each clone the same name as the group
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
	return(.Call(SXP_make_clones, sim.data$p, group, inherit.names, give.names, 
				 name.prefix, offspring, track.pedigree, give.ids, genomicSimulation:::expand.path(file.prefix), save.pedigree,
				 save.gebv, save.genotype, retain))
}
