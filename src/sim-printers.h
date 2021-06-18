#ifndef SIM_PRINTERS_H
#define SIM_PRINTERS_H

#include "sim-utils.h"
#include "sim-fitness.h"

/*--------------------------------Printing-----------------------------------*/
void save_simdata(FILE* f, SimData* m);

void save_marker_blocks(FILE* f, SimData* d, MarkerBlocks b);

void save_allele_matrix(FILE* f, AlleleMatrix* m, char** markers);
void save_transposed_allele_matrix(FILE* f, AlleleMatrix* m, char** markers);
void save_group_alleles(FILE* f, SimData* d, int group_id);
void save_transposed_group_alleles(FILE* f, SimData* d, int group_id);

void save_group_one_step_pedigree(FILE* f, SimData* d, int group); 
void save_one_step_pedigree(FILE* f, SimData* d); 
void save_group_full_pedigree(FILE* f, SimData* d, int group);
void save_full_pedigree(FILE* f, SimData* d);
void save_AM_pedigree(FILE* f, AlleleMatrix* m, SimData* parents);
void save_parents_of(FILE* f, AlleleMatrix* m, unsigned int p1, unsigned int p2);

void save_group_fitness(FILE* f, SimData* d, int group);
void save_fitness(FILE* f, DecimalMatrix* e, unsigned int* ids, char** names);
void save_all_fitness(FILE* f, SimData* d);

void save_count_matrix(FILE* f, SimData* d, char allele);
void save_count_matrix_of_group(FILE* f, SimData* d, char allele, int group);

#endif