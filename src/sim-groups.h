#ifndef SIM_GROUPS_H
#define SIM_GROUPS_H

#include "sim-utils.h"

/* Group Modification */
int combine_groups( SimData* d, int list_len, int group_ids[list_len]);
void split_into_individuals( SimData* d, int group_id);
void split_into_families(SimData* d, int group_id);
int split_from_group( SimData* d, int n, int indexes_to_split[n]);

#endif