#include "sim-groups.h"

/*-----------------------------------Groups----------------------------------*/
/** Combine a set of groups into one group.
 *
 * The function does so by setting the group membership of every genotype 
 * belonging to one of the groups to the same group number. 
 *
 * @param d the SimData struct on which to perform the operation
 * @param list_len the number of groups to be combined
 * @param group_ids an array of group numbers containing the groups that
 * are to be combined.
 * @returns the group number of the new combined group.
 */
int combine_groups( SimData* d, int list_len, int group_ids[list_len]) {
	int outGroup = group_ids[0];
	if (list_len < 2) {
		return outGroup;
	} else if (list_len == 2) {
		AlleleMatrix* m = d->m;
		int i;
		while (1) {
			// for each subject, check all group numbers.
			for (i = 0; i < m->n_subjects; ++i) {
				if (m->groups[i] == group_ids[1]) {
					m->groups[i] = outGroup;
				}
			}		
			
			if (m->next == NULL) {
				return outGroup;
			} else {
				m = m->next;
			}
		}	
		
		
	} else {
		AlleleMatrix* m = d->m;
		int i, j;
		while (1) {
			// for each subject, check all group numbers.
			for (i = 0; i < m->n_subjects; ++i) {
				
				// check the group number against all of those in our list
				for (j = 1; j < list_len; ++j) {
					if (m->groups[i] == group_ids[j]) {
						m->groups[i] = outGroup;
					}
				}
			}		
			
			if (m->next == NULL) {
				return outGroup;
			} else {
				m = m->next;
			}
		}
	}
}

/** Give every individual in the group a new group number that does not
 * belong to any other existing group (thereby allocating each genotype 
 * in the group to a new group of 1).
 *
 * Note: this does not return the group numbers of all the newly created 
 * groups-of-one.
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 */
void split_into_individuals( SimData* d, int group_id) {
	// get pre-existing numbers
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	// have another variable for the next id we can't allocate so we can still free the original.
	int level = 0;
	
	// looping through all entries
	AlleleMatrix* m = d->m;
	int i;
	int next_id = 0;
	while (1) {
		// check all subjects to see if this one belongs to the group.
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				// change it to a new unique group
				// first, find the next unused group;
				/*while (1) {
					if (next_id > existing_groups[level]) { // need to deal with when level is not a valid index anymore.
						++level;
					} else {
						++next_id;
						if (level >= n_groups || next_id < existing_groups[level]) {
							break;
						}
					}
				}*/
				++next_id;
				while (level < n_groups) {
					if (next_id < existing_groups[level]) {
						break;
					}
					
					++level;
					++next_id;
				}
				
				// save this entry as a member of that new group
				m->groups[i] = next_id;
			}
		}		
		
		if (m->next == NULL) {
			free(existing_groups);
			return;
		} else {
			m = m->next;
		}
	}
}

/** Split a group into a set of smaller groups, each of which contains the 
 * genotypes in the original group that share a particular pair of parents.
 * The number of new groups produced depends on the number of parent-combinations
 * in the set of genotypes in the provided group.
 *
 * Individuals with both parents unknown will be grouped together.
 *
 * Note: this does not return the group numbers of all the newly created 
 * groups-of-one.
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 */
void split_into_families(SimData* d, int group_id) {
	// get pre-existing numbers
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	// have another variable for the next id we can't allocate so we can still free the original.
	int level = 0;
	
	int families_found = 0;
	unsigned int family_groups[1000];
	unsigned int family_identities[2][1000];
	
	// looping through all entries
	AlleleMatrix* m = d->m;
	int i;
	int next_id = 0;
	while (1) {
		// check all subjects to see if this one belongs to the group.
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				// First, see if it is a member of a family we already know.
				for (int j = 0; j < families_found; ++j) {
					if (m->pedigrees[0][i] == family_identities[0][j] && 
							m->pedigrees[1][i] == family_identities[1][j]) {
						m->groups[i] = family_groups[j];
						break;
					}
				}
				
				// if the group number has not been updated in the above loop
				// (so we don't know this family yet)
				if (m->groups[i] == group_id) {
					// find the next unused group;
					++next_id;
					while (level < n_groups) {
						if (next_id < existing_groups[level]) {
							break;
						}
						
						++level;
						++next_id;
					}
					
					// save this entry as a member of that new group
					m->groups[i] = next_id;
					family_groups[families_found] = next_id;
					family_identities[0][families_found] = m->pedigrees[0][i];
					family_identities[1][families_found] = m->pedigrees[1][i];
					++families_found;
				}
			}
		}		
		
		if (m->next == NULL) {
			free(existing_groups);
			return;
		} else {
			m = m->next;
		}
	}
}

/** Take a list of indexes and allocate the genotypes at those indexes
 * to a new group.
 *
 * Does not check if all the indexes are valid/if all indexes have successfully
 * had their groups changed.
 *
 * @param d the SimData struct on which to perform the operation
 * @param n the number of indexes provided
 * @param indexes_to_split an array containing the indexes (0-based, starting
 * at the first entry at `d->m`) of the genotypes to allocate to the new group.
 * @returns the group number of the new group to which the provided indexes
 * have been allocated. 
 */
int split_from_group( SimData* d, int n, int indexes_to_split[n]) {
	int new_group = get_new_group_num(d);
	
	// Order the indexes
	qsort(indexes_to_split, n, sizeof(int), _ascending_int_comparer);
	
	AlleleMatrix* m = d->m;
	int total_i = 0;
	
	for (int i = 0; i < n; ++i) {
		while (indexes_to_split[i] >= total_i + m->n_subjects) {
			if (m->next == NULL) {
				warning("Only found %d out of %d indexes\n", i, n);
				return new_group;
			}
			total_i += m->n_subjects;
			m = m->next;
		}
		
		m->groups[indexes_to_split[i] - total_i] = new_group;
	}
	return new_group;
}

