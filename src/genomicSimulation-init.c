#include <R_ext/Rdynload.h>
#include "r-wrappers.h"

/*------------------------ R setup functions---------------------------------*/

R_CallMethodDef calledMethods[] = {
	{"SXP_change_label_amount", (DL_FUNC) &SXP_change_label_amount, 4},
	{"SXP_change_label_const", (DL_FUNC) &SXP_change_label_const, 4},
	{"SXP_change_label_default", (DL_FUNC) &SXP_change_label_default, 3},
	{"SXP_change_label_values", (DL_FUNC) &SXP_change_label_values, 5},
	{"SXP_change_name_values", (DL_FUNC) &SXP_change_name_values, 4},
	{"SXP_combine_groups", (DL_FUNC) &SXP_combine_groups, 2},
	{"SXP_create_label", (DL_FUNC) &SXP_create_label, 2},
	{"SXP_delete_group", (DL_FUNC) &SXP_delete_group, 2},
	{"SXP_delete_label", (DL_FUNC) &SXP_delete_label, 2},
	{"SXP_delete_eff_set", (DL_FUNC) &SXP_delete_eff_set, 2},
	{"SXP_get_best_available_haplotype", (DL_FUNC) &SXP_get_best_available_haplotype, 3},
	{"SXP_get_best_available_GEBV", (DL_FUNC) &SXP_get_best_available_GEBV, 3},
	{"SXP_get_best_haplotype", (DL_FUNC) &SXP_get_best_haplotype, 2},
	{"SXP_get_best_GEBV", (DL_FUNC) &SXP_get_best_GEBV, 2},
	{"SXP_get_worst_GEBV", (DL_FUNC) &SXP_get_worst_GEBV, 2},
	{"SXP_get_group_data", (DL_FUNC) &SXP_get_group_data, 4},
	{"SXP_get_group_gene_data", (DL_FUNC) &SXP_get_group_gene_data, 3},
	{"SXP_get_groups", (DL_FUNC) &SXP_get_groups, 1},
	{"SXP_group_eval", (DL_FUNC) &SXP_group_eval, 3},
	{"SXP_save_GEBVs", (DL_FUNC) &SXP_save_GEBVs, 4},
	{"SXP_save_chrsplit_block_effects", (DL_FUNC) &SXP_save_chrsplit_block_effects, 5},
	{"SXP_save_counts", (DL_FUNC) &SXP_save_counts, 4},
	{"SXP_save_file_block_effects", (DL_FUNC) &SXP_save_file_block_effects, 5},
	{"SXP_save_genotypes", (DL_FUNC) &SXP_save_genotypes, 4},
	{"SXP_save_pedigrees", (DL_FUNC) &SXP_save_pedigrees, 4},
	{"SXP_save_simdata", (DL_FUNC) &SXP_save_simdata, 2},
	{"SXP_selfing", (DL_FUNC) &SXP_selfing, 13},
	{"SXP_simple_selection", (DL_FUNC) &SXP_simple_selection, 5},
	{"SXP_simple_selection_bypercent", (DL_FUNC) &SXP_simple_selection_bypercent, 5},
	{"SXP_split_buckets", (DL_FUNC) &SXP_split_buckets, 3},
	{"SXP_split_by_label_range", (DL_FUNC) &SXP_split_by_label_range, 5},
	{"SXP_split_by_label_value", (DL_FUNC) &SXP_split_by_label_value, 4},
	{"SXP_split_evenly", (DL_FUNC) &SXP_split_evenly, 3},
	{"SXP_split_familywise", (DL_FUNC) &SXP_split_familywise, 2},
	{"SXP_split_halfsibwise", (DL_FUNC) &SXP_split_halfsibwise, 3},
	{"SXP_split_individuals", (DL_FUNC) &SXP_split_individuals, 2},
	{"SXP_split_out", (DL_FUNC) &SXP_split_out, 2},
	{"SXP_split_probabilities", (DL_FUNC) &SXP_split_probabilities, 3},
	{"SXP_split_randomly", (DL_FUNC) &SXP_split_randomly, 3},
	{"SXP_clone", (DL_FUNC) &SXP_clone, 13},
	{"SXP_cross_combinations", (DL_FUNC) &SXP_cross_combinations, 12},
	{"SXP_cross_randomly", (DL_FUNC) &SXP_cross_randomly, 14},
	{"SXP_cross_randomly_btwn", (DL_FUNC) &SXP_cross_randomly_btwn, 16},
	{"SXP_cross_Rcombinations", (DL_FUNC) &SXP_cross_Rcombinations, 13},
	{"SXP_cross_unidirectional", (DL_FUNC) &SXP_cross_unidirectional, 12},
	{"SXP_dcross_combinations", (DL_FUNC) &SXP_dcross_combinations, 12},
	{"SXP_doubled", (DL_FUNC) &SXP_doubled, 12},
	{"SXP_find_crossovers", (DL_FUNC) &SXP_find_crossovers, 5},
	{"SXP_load_data", (DL_FUNC) &SXP_load_data, 2},
	{"SXP_load_data_weff", (DL_FUNC) &SXP_load_data_weff, 3},
	{"SXP_load_more_genotypes", (DL_FUNC) &SXP_load_more_genotypes, 2},
	{"SXP_load_new_effects", (DL_FUNC) &SXP_load_new_effects, 2},
	{"SXP_send_map", (DL_FUNC) &SXP_send_map, 1},
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
