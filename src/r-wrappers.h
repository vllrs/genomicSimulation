#define R_Typeless_Calloc(size) R_Calloc(size,char)
#define GSC_MALLOC R_Typeless_Calloc
#define GSC_FREE R_Free
#include "sim-operations.h"

/*-------------------------- Setup -------------------------*/

FileFormatSpec SXP_parse_filespec_list(SEXP s_fileSpec);
SEXP SXP_load_data(SEXP s_alleleFile, SEXP s_mapFile, SEXP s_effectFile, SEXP s_fileSpec);
SEXP SXP_load_genotypes(SEXP exd, SEXP s_alleleFile, SEXP s_fileSpec);
SEXP SXP_load_effects(SEXP exd, SEXP s_effectFile);
SEXP SXP_load_map(SEXP exd, SEXP s_mapFile);
SEXP SXP_create_new_label(SEXP exd, SEXP s_default);

/*----------------------- Calculators ---------------------*/
SEXP SXP_see_optimal_haplotype(SEXP exd, SEXP s_eff_set, SEXP s_unknown_char);
SEXP SXP_get_optimal_possible_haplotype(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_unknown_char);
SEXP SXP_get_optimal_GEBV(SEXP exd, SEXP s_eff_set);
SEXP SXP_get_optimal_possible_GEBV(SEXP exd, SEXP s_groups, SEXP s_eff_set);
SEXP SXP_get_minimal_GEBV(SEXP exd, SEXP s_eff_set);

SEXP SXP_find_crossovers(SEXP exd, SEXP s_parentFile, SEXP s_outFile, SEXP s_windowSize, SEXP s_certainty);
SEXP SXP_send_map(SEXP exd, SEXP s_map);

/*----------------------- Data Access ---------------------*/
SEXP SXP_see_group_data(SEXP exd, SEXP s_group, SEXP s_whatData, SEXP s_eff_set_id, SEXP s_label_id);
SEXP SXP_see_marker_names(SEXP exd);
SEXP SXP_see_group_gene_data(SEXP exd, SEXP s_group, SEXP s_countAllele, SEXP s_unknownAllele);
SEXP SXP_see_existing_groups(SEXP exd);
SEXP SXP_change_name_to_values(SEXP exd, SEXP s_values, SEXP s_group, SEXP s_start);
SEXP SXP_change_allele_symbol(SEXP exd, SEXP s_markername, SEXP s_from, SEXP s_to);

/*------------------------ Deletors ------------------------*/
SEXP SXP_clear_simdata(SEXP exd);
void SXP_delete_simdata(SEXP sd);
SEXP SXP_delete_group(SEXP exd, SEXP s_groups);
SEXP SXP_delete_label(SEXP exd, SEXP s_labels);
SEXP SXP_delete_eff_set(SEXP exd, SEXP s_eff_sets);
SEXP SXP_delete_recombination_map(SEXP exd, SEXP s_maps);

/*------------------------- Groups -------------------------*/
SEXP SXP_combine_groups(SEXP exd, SEXP s_groups);
SEXP SXP_make_group_from(SEXP exd, SEXP s_indexes);
SEXP SXP_break_group_into_individuals(SEXP exd, SEXP s_group);
SEXP SXP_break_group_into_families(SEXP exd, SEXP s_group);
SEXP SXP_break_group_into_halfsib_families(SEXP exd, SEXP s_group, SEXP s_parent);
SEXP SXP_break_group_randomly(SEXP exd, SEXP s_group, SEXP s_n);
SEXP SXP_break_group_evenly(SEXP exd, SEXP s_group, SEXP s_n);
SEXP SXP_break_group_into_buckets(SEXP exd, SEXP s_group, SEXP s_buckets);
SEXP SXP_break_group_into_probabilities(SEXP exd, SEXP s_group, SEXP s_probs);
SEXP SXP_break_group_by_label_value(SEXP exd, SEXP s_label, SEXP s_value, SEXP s_groups);
SEXP SXP_break_group_by_label_range(SEXP exd, SEXP s_label, SEXP s_lowbound, 
									SEXP s_highbound, SEXP s_groups);

SEXP SXP_change_label_default(SEXP exd, SEXP s_labels, SEXP s_defaults);
SEXP SXP_change_label_to_values(SEXP exd, SEXP s_label, SEXP s_values, SEXP s_group, SEXP s_skip);
SEXP SXP_change_label_to_this(SEXP exd, SEXP s_label, SEXP s_const, SEXP s_groups);
SEXP SXP_change_label_by_amount(SEXP exd, SEXP s_label, SEXP s_incr, SEXP s_groups);

SEXP SXP_break_group_by_GEBV_num(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_number, SEXP s_bestIsLow);
SEXP SXP_break_group_by_GEBV_percent(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_percent, SEXP s_bestIsLow);

/*---------------------- Progression ----------------------*/
GenOptions SXP_create_genoptions(SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_random_crosses(SEXP exd, SEXP s_groups, SEXP s_crosses, SEXP s_cap, SEXP s_map,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_random_crosses_between(SEXP exd, SEXP s_group1, SEXP s_group2, SEXP s_cap1, SEXP s_cap2, SEXP s_map1, SEXP s_map2,
		SEXP s_crosses, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_targeted_crosses(SEXP exd, SEXP s_firstparents, SEXP s_secondparents, SEXP s_map1, SEXP s_map2,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_crosses_from_file(SEXP exd, SEXP s_filename, SEXP s_map1, SEXP s_map2, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_double_crosses_from_file(SEXP exd, SEXP s_filename, SEXP s_map1, SEXP s_map2, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_all_unidirectional_crosses(SEXP exd, SEXP s_groups, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_self_n_times(SEXP exd, SEXP s_groups, SEXP s_ngen, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_doubled_haploids(SEXP exd, SEXP s_groups, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);
SEXP SXP_make_clones(SEXP exd, SEXP s_groups, SEXP s_inherit_name, SEXP s_name,
        SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain);

/*----------------------- Printers ----------------------*/
SEXP SXP_save_genotypes(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_markersasrows);
SEXP SXP_save_allele_counts(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_allele, SEXP s_markersasrows);
SEXP SXP_save_pedigrees(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_fullpedigree);
SEXP SXP_save_GEBVs(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_eff_set);
SEXP SXP_save_local_GEBVs_blocks_from_file(SEXP exd, SEXP s_filename, SEXP block_file, SEXP s_group, SEXP s_eff_set);
SEXP SXP_save_local_GEBVs_blocks_from_chrsplit(SEXP exd, SEXP s_filename, SEXP s_nslices, SEXP s_group, SEXP s_map, SEXP s_eff_set);

