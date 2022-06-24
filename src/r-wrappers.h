#include "sim-operations.h"

/*-------------------------- Loaders -------------------------*/

SEXP SXP_load_data(SEXP s_alleleFile, SEXP s_mapFile);
SEXP SXP_load_data_weff(SEXP s_alleleFile, SEXP s_mapFile, SEXP s_effectFile);
SEXP SXP_load_more_genotypes(SEXP exd, SEXP s_alleleFile);
SEXP SXP_load_new_effects(SEXP exd, SEXP s_effectFile);

/*----------------Crossing-----------------*/
GenOptions create_genoptions(SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_cross_randomly(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_crosses, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_cross_randomly_btwn(SEXP exd, SEXP s_group1, SEXP s_group2, SEXP s_cap1, SEXP s_cap2,
		SEXP s_crosses, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_cross_Rcombinations(SEXP exd, SEXP s_firstparents, SEXP s_secondparents,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_cross_combinations(SEXP exd, SEXP s_filename, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_dcross_combinations(SEXP exd, SEXP s_filename, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_cross_unidirectional(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
/*SEXP cross_top(SEXP exd, SEXP group, SEXP percent, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);*/
SEXP SXP_selfing(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_ngen, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_doubled(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_one_cross(SEXP exd, SEXP s_parent1_index, SEXP s_parent2_index, SEXP s_name, 
		SEXP s_namePrefix, SEXP s_familySize, SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, 
		SEXP s_savePedigree, SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);

/*-----------------------------------Groups----------------------------------*/
SEXP SXP_combine_groups(SEXP exd, SEXP s_len, SEXP s_groups);
SEXP SXP_split_individuals(SEXP exd, SEXP s_group);
SEXP SXP_split_familywise(SEXP exd, SEXP s_group);
SEXP SXP_split_halfsibwise(SEXP exd, SEXP s_group, SEXP s_parent);
SEXP SXP_split_randomly(SEXP exd, SEXP s_group, SEXP s_n);
SEXP SXP_split_evenly(SEXP exd, SEXP s_group, SEXP s_n);
SEXP SXP_split_buckets(SEXP exd, SEXP s_group, SEXP s_buckets);
SEXP SXP_split_probabilities(SEXP exd, SEXP s_group, SEXP s_probs);

SEXP SXP_split_out(SEXP exd, SEXP s_len, SEXP s_indexes);

/*-----------------Fitness-------------------*/
SEXP SXP_group_eval(SEXP exd, SEXP s_group);
SEXP SXP_simple_selection(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_number, SEXP s_bestIsLow);
SEXP SXP_simple_selection_bypercent(SEXP exd, SEXP s_glen, SEXP s_groups, SEXP s_percent, SEXP s_bestIsLow);


/*-----------------Data access---------------*/
SEXP SXP_get_best_haplotype(SEXP exd);
SEXP SXP_get_best_GEBV(SEXP exd);
SEXP SXP_get_worst_GEBV(SEXP exd);

SEXP SXP_find_crossovers(SEXP exd, SEXP s_parentFile, SEXP s_outFile, SEXP s_windowSize, SEXP s_certainty);
SEXP SXP_send_map(SEXP exd);

SEXP SXP_get_groups(SEXP exd);
SEXP SXP_get_group_data(SEXP exd, SEXP s_group, SEXP s_whatData);

/*--------------------------------Printing-----------------------------------*/
SEXP SXP_save_simdata(SEXP exd, SEXP s_filename);
SEXP SXP_save_genotypes(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_type);
SEXP SXP_save_counts(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_allele);
SEXP SXP_save_pedigrees(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_type);
SEXP SXP_save_GEBVs(SEXP exd, SEXP s_filename, SEXP s_group);
SEXP SXP_save_file_block_effects(SEXP exd, SEXP s_filename, SEXP block_file, SEXP s_group);
SEXP SXP_save_chrsplit_block_effects(SEXP exd, SEXP s_filename, SEXP s_nslices, SEXP s_group);

/*--------------------------------Deletors------------------------------------*/
SEXP clear_simdata(SEXP exd);
void SXP_delete_simdata(SEXP sd);
SEXP SXP_delete_group(SEXP exd, SEXP s_len, SEXP s_groups);