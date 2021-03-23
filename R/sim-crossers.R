#' Performs random crosses between members of a group.
#'
#' \code{cross.randomly} returns the group number of the group
#' that the new genotypes were loaded into. Selfing is not permitted.
#'
#' For random crossing, the offspring parameter represents the number of 
#' random crosses to perform.
#'
#' Random crosses can only be performed within a group, not within the
#' whole set of genotypes tracked by the simulation.
#'
#' @param group The group number from which to draw the parents of these
#' random crosses. If a vector, the parameters and breeding method 
#' (eg. random crossing) will be applied to each group number in the vector.
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
#' @useDynLib genomicSimulation cross_randomly
#' @export
cross.randomly <- function(group, n.crosses=5, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(cross_randomly, sim.data$p, length(group), group, n.crosses, give.names, name.prefix, 
	             offspring, track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs defined crosses as laid out in a file.
#'
#' \code{cross.combinations} returns the group number of the group
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
#' @useDynLib genomicSimulation cross_combinations
#' @export
cross.combinations <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(cross_combinations, sim.data$p, cross.file, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Performs defined crosses between children of known parents as 
#' laid out in a file.
#'
#' \code{cross.dc.combinations} returns the group number of the group
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
#' @useDynLib genomicSimulation dcross_combinations
#' @export
cross.dc.combinations <- function(cross.file, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(dcross_combinations, sim.data$p, cross.file, give.names, name.prefix, offspring, 
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
#' @useDynLib genomicSimulation cross_unidirectional
#' @export
cross.all.pairs <- function(group, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(cross_unidirectional, sim.data$p, length(group), group, give.names, name.prefix,
	             offspring, track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

# Performs random crosses between a high-scoring subset of the group of genotypes
#
# \code{cross.from.top.pc} returns the group number of the group
# that the new genotypes were loaded into. Only the alleles at SNPs
# already tracked by the SimData are saved.
#
# The function calculates the GEBVs of the genotypes in the given group,
# then takes the genotypes that are in the top [threshold]% of the group
# by GEBV and performs random crosses betwee them.
#
#
# The offspring parameter represents the number of random crosses performed.
#
# @param threshold The percent threshold with the highest GEBVs to subset and do 
# random crosses of. Should not be the decimal version of a percentage: eg for top 5%, 
# threshold should be 5, not 0.05
# @inheritParams cross.randomly
# @return The group number of the new genotypes produced
#
# @family crossing functions
# @useDynLib genomicSimulation cross_top
# @export
#cross.from.top.pc <- function(group, threshold, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
#		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
#		save.gebv=FALSE, save.genotype=FALSE) {
#	if (is.null(sim.data$p)) { stop("Please load.data first.") }
#	return(.Call(cross_top, sim.data$p, group, threshold, give.names, name.prefix, offspring, track.pedigree, 
#				 give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
#}

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
#' @useDynLib genomicSimulation SXP_selfing
#' @export
self.n.times <- function(group, n, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_selfing, sim.data$p, length(group), group, n, give.names, name.prefix, offspring, 
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
#' @useDynLib genomicSimulation SXP_doubled
#' @export
make.doubled.haploids <- function(group, offspring=1, retain=TRUE, give.names=FALSE, name.prefix=NULL, 
		track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_doubled, sim.data$p, length(group), group, give.names, name.prefix, offspring, 
				 track.pedigree, give.ids, file.prefix, save.pedigree, save.gebv, save.genotype, retain))
}

#' Perform a cross between two specific lines.
#'
#' \code{cross} allows you to perform crosses between specific lines. It 
#' can be called in a for loop, looping over the rows of a dataframe, to 
#' achieve similar functionality to \code{\link{cross.combinations}} without
#' needing to save the combinations to be performed to a file.
#'
#' The offspring parameter represents the number of offspring to make from
#' this pair of parents.
#' 
#' @inheritParams cross.randomly
#' @param offspring The number of offspring to produce from this pair of parents
#' @return The group number of the group that the generated offspring were loaded into.
#'
#' @family crossing functions
#' @useDynLib genomicSimulation SXP_one_cross
#' @export
cross <- function(parent1.index, parent2.index, offspring=1, retain=TRUE, give.names=FALSE, 
		name.prefix=NULL, track.pedigree=TRUE, give.ids=TRUE, file.prefix=NULL, save.pedigree=FALSE, 
		save.gebv=FALSE, save.genotype=FALSE) {
	if (is.null(sim.data$p)) { stop("Please load.data first.") }
	return(.Call(SXP_one_cross, sim.data$p, parent1.index, parent2.index, give.names, 
				 name.prefix, offspring, track.pedigree, give.ids, file.prefix, save.pedigree, 
				 save.gebv, save.genotype, retain))
}

