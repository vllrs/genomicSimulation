#ifndef SIM_LOADERS_H
#define SIM_LOADERS_H

#include "utils.h"
#include "sim-printers.h"

SEXP load_data(SEXP alleleFile, SEXP mapFile, SEXP groupName);
SEXP load_data_weff(SEXP alleleFile, SEXP mapFile, SEXP effectFile, SEXP groupName);
SEXP load_more_genotypes(SEXP exd, SEXP alleleFile, SEXP groupName);
SEXP load_new_effects(SEXP exd, SEXP effectFile);

/* Loaders */
int load_transposed_genes_to_simdata(SimData* d, const char* filename, const char* groupName);
int load_more_transposed_genes_to_simdata(SimData* d, const char* filename, const char* groupName);
//int load_genes_to_simdata(SimData* d, const char* filename); //@ add
//int load_more_genes_to_simdata(SimData* d, const char* filename); //@ add
int load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename, const char* groupName);
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
static inline int has_same_alleles(char* p1, char* p2, int i) {
	return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}
// w is window length, i is start value
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
		
SEXP find_crossovers(SEXP exd, SEXP parentFile, SEXP outFile, SEXP windowSize, SEXP certainty);
SEXP send_map(SEXP exd);

#endif