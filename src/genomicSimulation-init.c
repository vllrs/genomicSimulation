#include <R_ext/Rdynload.h>
#include "r-wrappers.h"

/*------------------------ R setup functions---------------------------------*/

R_CallMethodDef calledMethods[] = {
	{"SXP_change_label_by_amount", (DL_FUNC) &SXP_change_label_by_amount, 4},
	{"SXP_change_label_to_this", (DL_FUNC) &SXP_change_label_to_this, 4},
	{"SXP_change_label_default", (DL_FUNC) &SXP_change_label_default, 3},
	{"SXP_change_label_to_values", (DL_FUNC) &SXP_change_label_to_values, 5},
	{"SXP_change_name_to_values", (DL_FUNC) &SXP_change_name_to_values, 4},
	{"SXP_change_allele_symbol", (DL_FUNC) &SXP_change_allele_symbol, 4},
	{"SXP_combine_groups", (DL_FUNC) &SXP_combine_groups, 2},
	{"SXP_create_new_label", (DL_FUNC) &SXP_create_new_label, 2},
	{"SXP_delete_group", (DL_FUNC) &SXP_delete_group, 2},
	{"SXP_delete_label", (DL_FUNC) &SXP_delete_label, 2},
	{"SXP_delete_eff_set", (DL_FUNC) &SXP_delete_eff_set, 2},
	{"SXP_delete_recombination_map", (DL_FUNC) &SXP_delete_recombination_map, 2},
	{"SXP_get_optimal_possible_haplotype", (DL_FUNC) &SXP_get_optimal_possible_haplotype, 3},
	{"SXP_get_optimal_possible_GEBV", (DL_FUNC) &SXP_get_optimal_possible_GEBV, 3},
	{"SXP_see_optimal_haplotype", (DL_FUNC) &SXP_see_optimal_haplotype, 2},
	{"SXP_get_optimal_GEBV", (DL_FUNC) &SXP_get_optimal_GEBV, 2},
	{"SXP_get_minimal_GEBV", (DL_FUNC) &SXP_get_minimal_GEBV, 2},
	{"SXP_see_group_data", (DL_FUNC) &SXP_see_group_data, 5},
	{"SXP_see_marker_names", (DL_FUNC) &SXP_see_marker_names, 1},
	{"SXP_see_group_gene_data", (DL_FUNC) &SXP_see_group_gene_data, 3},
	{"SXP_see_existing_groups", (DL_FUNC) &SXP_see_existing_groups, 1},
	{"SXP_save_GEBVs", (DL_FUNC) &SXP_save_GEBVs, 4},
	{"SXP_save_local_GEBVs_blocks_from_chrsplit", (DL_FUNC) &SXP_save_local_GEBVs_blocks_from_chrsplit, 6},
	{"SXP_save_allele_counts", (DL_FUNC) &SXP_save_allele_counts, 4},
	{"SXP_save_local_GEBVs_blocks_from_file", (DL_FUNC) &SXP_save_local_GEBVs_blocks_from_file, 5},
	{"SXP_save_genotypes", (DL_FUNC) &SXP_save_genotypes, 4},
	{"SXP_save_pedigrees", (DL_FUNC) &SXP_save_pedigrees, 4},
	{"SXP_self_n_times", (DL_FUNC) &SXP_self_n_times, 14},
	{"SXP_break_group_by_GEBV_num", (DL_FUNC) &SXP_break_group_by_GEBV_num, 5},
	{"SXP_break_group_by_GEBV_percent", (DL_FUNC) &SXP_break_group_by_GEBV_percent, 5},
	{"SXP_break_group_into_buckets", (DL_FUNC) &SXP_break_group_into_buckets, 3},
	{"SXP_break_group_by_label_range", (DL_FUNC) &SXP_break_group_by_label_range, 5},
	{"SXP_break_group_by_label_value", (DL_FUNC) &SXP_break_group_by_label_value, 4},
	{"SXP_break_group_evenly", (DL_FUNC) &SXP_break_group_evenly, 3},
	{"SXP_break_group_into_families", (DL_FUNC) &SXP_break_group_into_families, 2},
	{"SXP_break_group_into_halfsib_families", (DL_FUNC) &SXP_break_group_into_halfsib_families, 3},
	{"SXP_break_group_into_individuals", (DL_FUNC) &SXP_break_group_into_individuals, 2},
	{"SXP_make_group_from", (DL_FUNC) &SXP_make_group_from, 2},
	{"SXP_break_group_into_probabilities", (DL_FUNC) &SXP_break_group_into_probabilities, 3},
	{"SXP_break_group_randomly", (DL_FUNC) &SXP_break_group_randomly, 3},
	{"SXP_make_clones", (DL_FUNC) &SXP_make_clones, 13},
	{"SXP_make_crosses_from_file", (DL_FUNC) &SXP_make_crosses_from_file, 14},
	{"SXP_make_random_crosses", (DL_FUNC) &SXP_make_random_crosses, 15},
	{"SXP_make_random_crosses_between", (DL_FUNC) &SXP_make_random_crosses_between, 18},
	{"SXP_make_targeted_crosses", (DL_FUNC) &SXP_make_targeted_crosses, 15},
	{"SXP_make_all_unidirectional_crosses", (DL_FUNC) &SXP_make_all_unidirectional_crosses, 13},
	{"SXP_make_double_crosses_from_file", (DL_FUNC) &SXP_make_double_crosses_from_file, 14},
	{"SXP_make_doubled_haploids", (DL_FUNC) &SXP_make_doubled_haploids, 13},
	{"SXP_find_crossovers", (DL_FUNC) &SXP_find_crossovers, 5},
	{"SXP_load_data", (DL_FUNC) &SXP_load_data, 3},
	{"SXP_load_genotypes", (DL_FUNC) &SXP_load_genotypes, 2},
	{"SXP_load_effects", (DL_FUNC) &SXP_load_effects, 2},
	{"SXP_load_map", (DL_FUNC) &SXP_load_map, 2},
	{"SXP_send_map", (DL_FUNC) &SXP_send_map, 2},
	{NULL}
};

void R_init_genomicSimulation(DllInfo *info) {
	// R_RegisterCCallable doesn't seem to have any effect, so currently commented out.
	/*R_RegisterCCallable("genomicSimulation", "SXP_combine_groups", (DL_FUNC) &SXP_combine_groups);
	R_RegisterCCallable("genomicSimulation", "SXP_delete_group", (DL_FUNC) &SXP_delete_group);
	R_RegisterCCallable("genomicSimulation", "SXP_doubled", (DL_FUNC) &SXP_doubled);
	R_RegisterCCallable("genomicSimulation", "SXP_get_best_genotype", (DL_FUNC) &SXP_get_best_genotype);
	R_RegisterCCallable("genomicSimulation", "SXP_get_group_data", (DL_FUNC) &SXP_get_group_data);
	R_RegisterCCallable("genomicSimulation", "SXP_get_groups", (DL_FUNC) &SXP_get_groups);
	R_RegisterCCallable("genomicSimulation", "SXP_group_eval", (DL_FUNC) &SXP_group_eval);
	R_RegisterCCallable("genomicSimulation", "SXP_one_cross", (DL_FUNC) &SXP_one_cross);
	R_RegisterCCallable("genomicSimulation", "SXP_save_GEBVs", (DL_FUNC) &SXP_save_GEBVs);
	R_RegisterCCallable("genomicSimulation", "SXP_save_block_effects", (DL_FUNC) &SXP_save_block_effects);
	R_RegisterCCallable("genomicSimulation", "SXP_save_counts", (DL_FUNC) &SXP_save_counts);
	R_RegisterCCallable("genomicSimulation", "SXP_save_genotypes", (DL_FUNC) &SXP_save_genotypes);
	R_RegisterCCallable("genomicSimulation", "SXP_save_pedigrees", (DL_FUNC) &SXP_save_pedigrees);
	R_RegisterCCallable("genomicSimulation", "SXP_save_simdata", (DL_FUNC) &SXP_save_simdata);
	R_RegisterCCallable("genomicSimulation", "SXP_selfing", (DL_FUNC) &SXP_selfing);
	R_RegisterCCallable("genomicSimulation", "SXP_simple_selection", (DL_FUNC) &SXP_simple_selection);
	R_RegisterCCallable("genomicSimulation", "SXP_simple_selection_bypercent", (DL_FUNC) &SXP_simple_selection_bypercent);
	R_RegisterCCallable("genomicSimulation", "SXP_split_familywise", (DL_FUNC) &SXP_split_familywise);
	R_RegisterCCallable("genomicSimulation", "SXP_split_individuals", (DL_FUNC) &SXP_split_individuals);
	R_RegisterCCallable("genomicSimulation", "SXP_split_out", (DL_FUNC) &SXP_split_out);
	R_RegisterCCallable("genomicSimulation", "SXP_cross_combinations", (DL_FUNC) &SXP_cross_combinations);
	R_RegisterCCallable("genomicSimulation", "SXP_cross_randomly", (DL_FUNC) &SXP_cross_randomly);
	R_RegisterCCallable("genomicSimulation", "SXP_cross_unidirectional", (DL_FUNC) &SXP_cross_unidirectional);
	R_RegisterCCallable("genomicSimulation", "SXP_dcross_combinations", (DL_FUNC) &SXP_dcross_combinations);
	R_RegisterCCallable("genomicSimulation", "SXP_find_crossovers", (DL_FUNC) &SXP_find_crossovers);
	R_RegisterCCallable("genomicSimulation", "SXP_load_data", (DL_FUNC) &SXP_load_data);
	R_RegisterCCallable("genomicSimulation", "SXP_load_data_weff", (DL_FUNC) &SXP_load_data_weff);
	R_RegisterCCallable("genomicSimulation", "SXP_load_more_genotypes", (DL_FUNC) &SXP_load_more_genotypes);
	R_RegisterCCallable("genomicSimulation", "SXP_load_new_effects", (DL_FUNC) &SXP_load_new_effects);
	R_RegisterCCallable("genomicSimulation", "SXP_send_map", (DL_FUNC) &SXP_send_map);*/

	R_registerRoutines(info, NULL, calledMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
