#ifndef SIM_CROSSERS_H
#define SIM_CROSSERS_H

#include "utils.h"
#include "sim-gebv.h"
#include "sim-printers.h"

GenOptions create_genoptions(SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP cross_randomly(SEXP exd, SEXP glen, SEXP groups, SEXP crosses, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP cross_combinations(SEXP exd, SEXP filename, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP dcross_combinations(SEXP exd, SEXP filename, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP cross_unidirectional(SEXP exd, SEXP glen, SEXP groups, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
/*SEXP cross_top(SEXP exd, SEXP group, SEXP percent, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);*/
SEXP SXP_selfing(SEXP exd, SEXP glen, SEXP groups, SEXP n, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP SXP_doubled(SEXP exd, SEXP glen, SEXP groups, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain);
SEXP SXP_one_cross(SEXP exd, SEXP parent1_index, SEXP parent2_index, SEXP name, 
		SEXP namePrefix, SEXP familySize, SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, 
		SEXP savePedigree, SEXP saveEffects, SEXP saveGenes, SEXP retain);


/* Crossers */
void generate_gamete(SimData* d, char* parent_genome, char* output);
void generate_cross(SimData* d, char* parent1_genome, char* parent2_genome, char* output);
void generate_doubled_haploid(SimData* d, char* parent_genome, char* output);

int cross_this_pair(SimData* d, int parent1_index, int parent2_index, GenOptions g);
int cross_random_individuals(SimData* d, int from_group, int n_crosses, GenOptions g);
int cross_these_combinations(SimData* d, int n_combinations, int combinations[2][n_combinations],  GenOptions g);
int self_n_times(SimData* d, int n, int group, GenOptions g);
int make_doubled_haploids(SimData* d, int group, GenOptions g); //@add

int make_all_unidirectional_crosses(SimData* d, int from_group, GenOptions g);
int make_n_crosses_from_top_m_percent(SimData* d, int n, int m, int group, GenOptions g);
int make_crosses_from_file(SimData* d, const char* input_file, GenOptions g);
int make_double_crosses_from_file(SimData* d, const char* input_file, GenOptions g);

#endif