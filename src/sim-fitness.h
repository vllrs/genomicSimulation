#ifndef SIM_FITNESS_H
#define SIM_FITNESS_H

#include "sim-utils.h"
#include "sim-groups.h"

/* Fitness calculators */
int split_by_bv(SimData* d, int group, int top_n, int lowIsBest);
DecimalMatrix calculate_group_bvs(SimData* d, int group);
DecimalMatrix calculate_bvs( AlleleMatrix* m, EffectMatrix* e);
DecimalMatrix calculate_count_matrix_of_allele_for_ids( AlleleMatrix* m, unsigned int* for_ids, unsigned int n_ids, char allele);
DecimalMatrix calculate_full_count_matrix_of_allele( AlleleMatrix* m, char allele);

MarkerBlocks create_n_blocks_by_chr(SimData* d, int n);
MarkerBlocks read_block_file(SimData* d, const char* block_file);
void calculate_group_local_bvs(SimData* d, MarkerBlocks b, const char* output_file, int group);
void calculate_local_bvs(SimData* d, MarkerBlocks b, const char* output_file);

char* calculate_optimal_alleles(SimData* d);
double calculate_optimum_bv(SimData* d);
double calculate_minimum_bv(SimData* d);

/* Recombination calculators */
int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2, 
		unsigned int p2num, char* offspring, int certain); // forward filling, window size 1
int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2, 
		unsigned int p2num, char* offspring, int window_size, int certain); // forward filling, window size n
int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file, 
		int window_len, int certain);
		
// int has_same_alleles(char* p1, char* p2, int i);
// int has_same_alleles_window(char* g1, char* g2, int start, int w);
/** Simple operator to determine if at marker i, two genotypes share at least
 * one allele. Checks only 3 of four possible permutations because assumes
 * there cannot be more than two alleles at a given marker.
 *
 * @param p1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param p2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param i index of the marker at which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles(char* p1, char* p2, int i) {
	return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}
// w is window length, i is start value
/** Simple operator to determine if at markers with indexes i to i+w inclusive, two genotypes 
 * share at least one allele. Checks only 3 of four possible permutations at each marker 
 * because assumes there cannot be more than two alleles at a given marker. For the return value
 * to be true, there must be at least one match at every one of the markers in the window.
 *
 * @param g1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param g2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param start index of the first marker in the window over which to perform the check
 * @param w length of the window over which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles_window(char* g1, char* g2, int start, int w) {
	int same = TRUE;
	int i;
	for (int j = 0; j < w; ++j) {
		i = start + j;
		same = same && (g1[i<<1] == g2[i<<1] || g1[(i<<1) + 1] == g2[i] || g1[i] == g2[(i<<1) + 1]);
	}
	return same;
}



#endif