#ifndef SIM_LOADERS_H
#define SIM_LOADERS_H

#include "sim-utils.h"
#include "sim-printers.h"

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

#endif