#ifndef SIM_GROUPS_H
#define SIM_GROUPS_H

#include "sim-utils.h"
#include "sim-fitness.h" // for get_group_bvs only
#include "sim-printers.h" // for get_group_pedigrees only

/* Group Modification */
int combine_groups( SimData* d, int list_len, int group_ids[list_len]);
void split_into_individuals( SimData* d, int group_id);
void split_into_families(SimData* d, int group_id);
int split_from_group( SimData* d, int n, int indexes_to_split[n]);

/* Group data access */
int get_group_size( SimData* d, int group_id);
char** get_group_genes( SimData* d, int group_id, int group_size);
char** get_group_names( SimData* d, int group_id, int group_size);
unsigned int* get_group_ids( SimData* d, int group_id, int group_size);
unsigned int* get_group_indexes(SimData* d, int group_id, int group_size);
double* get_group_bvs( SimData* d, int group_id, int group_size);
unsigned int* get_group_parent_IDs( SimData* d, int group_id, int group_size, int parent);
char** get_group_parent_names( SimData* d, int group_id, int group_size, int parent);
char** get_group_pedigrees( SimData* d, int group_id, int group_size);

#endif