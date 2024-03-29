#' Performs random crosses between members of a group.
#'
#' \code{cross.randomly} returns the group number of the group
#' that the new genotypes were loaded into. Selfing is not permitted.
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
cross.randomly <- function(group, n.crosses=5, cap=0, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_cross_randomly, sim.data$p, group, n.crosses, cap, give.names, 
	             name.prefix, offspring, track.pedigree, give.ids, file.prefix, 
	             save.pedigree, save.gebv, save.genotype, retain))
}


#' Performs random crosses between two groups.
#'
#' \code{cross.randomly.between} performs crosses where the first parent comes 
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
#' @inheritParams cross.randomly
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
cross.randomly.between <- function(group1, group2, cap1=0, cap2=0, 
		n.crosses=5, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
  
	return(.Call(SXP_cross_randomly_btwn, sim.data$p, as.integer(group1), 
	             as.integer(group2), cap1, cap2,
				 n.crosses, give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}


#' Performs defined crosses as passed in as R vectors.
#'
#' \code{cross.combinations} returns the group number of the group
#' that the new genotypes were loaded into. 
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
#' @inheritParams cross.randomly
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
cross.combinations <- function(first.parents, second.parents,
		offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_cross_Rcombinations, sim.data$p, first.parents, second.parents,
				 give.names, name.prefix, offspring, track.pedigree, give.ids, 
				 file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs defined crosses as laid out in a file.
#'
#' \code{cross.combinations.file} returns the group number of the group
#' that the new genotypes were loaded into. 
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
#' @inheritParams cross.randomly
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
cross.combinations.file <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_cross_combinations, sim.data$p, cross.file, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs defined crosses between children of known parents as 
#' laid out in a file.
#'
#' \code{cross.dc.combinations.file} returns the group number of the group
#' that the new genotypes were loaded into. 
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
#' @inheritParams cross.randomly
#' @param offspring The number of times each combination in the file is crossed.
#' @return The group number of the new crosses produced
#'
#' @family crossing functions
#' @export
cross.dc.combinations.file <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_dcross_combinations, sim.data$p, cross.file, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs crosses between every line and every other line in a group
#' in one direction
#'
#' \code{cross.all.pairs} returns the group number of the group
#' that the new genotypes were loaded into.
#'
#' The function performs one cross between each pair of genotypes in the
#' group. Order of parents is not considered because in the simulation method
#' this has no effect.
#'
#' @inheritParams cross.randomly
#' @param offspring The number of times each combination of group members is crossed.
#' @return The group number of the new crosses produced, or 0 if they could not be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
cross.all.pairs <- function(group, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_cross_unidirectional, sim.data$p, group, give.names, name.prefix,
	             offspring, track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs n selfing steps on the lines in a group
#'
#' \code{self.n.times} returns the group number of the group
#' that the new genotypes were loaded into.
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
#' @inheritParams cross.randomly
#' @param offspring This many offspring of each group member will be produced at the 
#' first selfing step. After that step, exactly one selfed offspring from each
#' parent will be progressed, no matter this value.
#' @return The group number of the new genotypes produced, or 0 if none could be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
self.n.times <- function(group, n, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_selfing, sim.data$p, group, n, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
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
#' @inheritParams cross.randomly
#' @param offspring The number of doubled haploids to make from each group member.
#' @return The group number of the new genotypes produced, or 0 if none could be
#' produced due to an invalid parent group number being provided.
#'
#' @family crossing functions
#' @export
make.doubled.haploids <- function(group, offspring=1, retain=TRUE, 
		give.names=FALSE, name.prefix=NULL, track.pedigree=TRUE, give.ids=TRUE, 
		file.prefix=NULL, save.pedigree=FALSE, save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_doubled, sim.data$p, group, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Creates a genetically identical copy of each member of a group
#'
#' \code{make.clones} returns the group number of the group
#' that the new clone genotypes were loaded into.
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
#' @inheritParams cross.randomly
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
		file.prefix=NULL, save.pedigree=FALSE, save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_clone, sim.data$p, group, inherit.names, give.names, 
				 name.prefix, offspring, track.pedigree, give.ids, file.prefix, save.pedigree,
				 save.gebv, save.genotype, retain))
}
