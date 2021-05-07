#ifndef SIM_LOADERS_H
#define SIM_LOADERS_H

#include "utils.h"
#include "sim-printers.h"

SEXP SXP_load_data(SEXP alleleFile, SEXP mapFile);
SEXP SXP_load_data_weff(SEXP alleleFile, SEXP mapFile, SEXP effectFile);
SEXP SXP_load_more_genotypes(SEXP exd, SEXP alleleFile);
SEXP SXP_load_new_effects(SEXP exd, SEXP effectFile);

/* Loaders */
int load_transposed_genes_to_simdata(SimData* d, const char* filename);
int load_more_transposed_genes_to_simdata(SimData* d, const char* filename);
//int load_genes_to_simdata(SimData* d, const char* filename); //@ add
//int load_more_genes_to_simdata(SimData* d, const char* filename); //@ add
int load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename);
void load_genmap_to_simdata(SimData* d, const char* filename);
void get_sorted_markers(SimData* d, int actual_n_markers);
void get_chromosome_locations(SimData *d);
void load_effects_to_simdata(SimData* d, const char* filename);
int load_all_simdata(SimData* d, const char* data_file, const char* map_file, const char* effect_file);

/* Recombination calculators */
int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2, 
		unsigned int p2num, char* offspring, int certain); // forward filling, window size 1
int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2, 
		unsigned int p2num, char* offspring, int window_size, int certain); // forward filling, window size n
		

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

int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file, 
		int window_len, int certain);
		
SEXP SXP_find_crossovers(SEXP exd, SEXP parentFile, SEXP outFile, SEXP windowSize, SEXP certainty);
SEXP SXP_send_map(SEXP exd);

#endif