#include <R_ext/Rdynload.h>
#include "r-wrappers.h"

/*------------------------ R setup functions---------------------------------*/

R_CallMethodDef calledMethods[] = {
	{"SXP_combine_groups", (DL_FUNC) &SXP_combine_groups, 3},
	{"SXP_delete_group", (DL_FUNC) &SXP_delete_group, 2},
	{"SXP_get_best_genotype", (DL_FUNC) &SXP_get_best_genotype, 1},
	{"SXP_get_best_GEBV", (DL_FUNC) &SXP_get_best_GEBV, 1},
	{"SXP_get_group_data", (DL_FUNC) &SXP_get_group_data, 3},	
	{"SXP_get_groups", (DL_FUNC) &SXP_get_groups, 1},	
	{"SXP_group_eval", (DL_FUNC) &SXP_group_eval, 2},	
	{"SXP_one_cross", (DL_FUNC) &SXP_one_cross, 13},	
	{"SXP_save_GEBVs", (DL_FUNC) &SXP_save_GEBVs, 3},	
	{"SXP_save_chrsplit_block_effects", (DL_FUNC) &SXP_save_chrsplit_block_effects, 4},
	{"SXP_save_counts", (DL_FUNC) &SXP_save_counts, 4},	
	{"SXP_save_file_block_effects", (DL_FUNC) &SXP_save_file_block_effects, 4},
	{"SXP_save_genotypes", (DL_FUNC) &SXP_save_genotypes, 4},	
	{"SXP_save_pedigrees", (DL_FUNC) &SXP_save_pedigrees, 4},	
	{"SXP_save_simdata", (DL_FUNC) &SXP_save_simdata, 2},	
	{"SXP_selfing", (DL_FUNC) &SXP_selfing, 14},	
	{"SXP_simple_selection", (DL_FUNC) &SXP_simple_selection, 5},	
	{"SXP_simple_selection_bypercent", (DL_FUNC) &SXP_simple_selection_bypercent, 5},	
	{"SXP_split_familywise", (DL_FUNC) &SXP_split_familywise, 2},	
	{"SXP_split_individuals", (DL_FUNC) &SXP_split_individuals, 2},	
	{"SXP_split_out", (DL_FUNC) &SXP_split_out, 3},	
	{"SXP_cross_combinations", (DL_FUNC) &SXP_cross_combinations, 12},
	{"SXP_cross_randomly", (DL_FUNC) &SXP_cross_randomly, 14},
	{"SXP_cross_Rcombinations", (DL_FUNC) &SXP_cross_Rcombinations, 13},
	{"SXP_cross_unidirectional", (DL_FUNC) &SXP_cross_unidirectional, 13},
	{"SXP_dcross_combinations", (DL_FUNC) &SXP_dcross_combinations, 12},
	{"SXP_doubled", (DL_FUNC) &SXP_doubled, 13},
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
