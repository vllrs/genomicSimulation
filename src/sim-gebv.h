#ifndef SIM_GEBV_H
#define SIM_GEBV_H

#include "utils.h"

/* Fitness calculators */
int split_group_by_fitness(SimData* d, int group, int top_n, int lowIsBest);
DecimalMatrix calculate_fitness_metric_of_group(SimData* d, int group);
DecimalMatrix calculate_fitness_metric( AlleleMatrix* m, EffectMatrix* e);
DecimalMatrix calculate_count_matrix_of_allele_for_ids( AlleleMatrix* m, unsigned int* for_ids, unsigned int n_ids, char allele);
DecimalMatrix calculate_full_count_matrix_of_allele( AlleleMatrix* m, char allele);

/** A struct used to store a set of blocks of markers. 
 *
 * @param num_blocks the number of blocks whose details are stored here.
 * @param num_markers_in_block pointer to a heap array of length num_blocks
 * containing the number of markers that make up each block
 * @param markers_in_block pointer to a heap array of length num_blocks, each
 * entry in which is a pointer to a heap array with length corresponding to 
 * the value of the corresponding entry in num_markers_in_block whose values
 * are the indexes in the SimData of the markers that make up that block.
 */
struct markerBlocks {
	int num_blocks;
	int* num_markers_in_block;
	int** markers_in_block;
};

void calculate_group_block_effects(SimData* d, const char* block_file, const char* output_file, int group);
struct markerBlocks read_block_file(SimData* d, const char* block_file);
void calculate_all_block_effects(SimData* d, const char* block_file, const char* output_file);

char* calculate_ideal_genotype(SimData* d);
SEXP SXP_get_best_genotype(SEXP exd);

SEXP SXP_group_eval(SEXP exd, SEXP group);
SEXP SXP_simple_selection(SEXP exd, SEXP glen, SEXP groups, SEXP number, SEXP bestIsLow);
SEXP SXP_simple_selection_bypercent(SEXP exd, SEXP glen, SEXP groups, SEXP percent, SEXP bestIsLow);

#endif