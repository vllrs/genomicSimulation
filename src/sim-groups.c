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




//-----------------------------------Data Access-----------------------------------

/** Function to count the number of genotypes that currently belong to
 * the specified group.
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * This goes through and checks every genotype in the SimData, because 
 * there is currently no centralised tracking of groups.
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group whose members need to
 * be counted.
 * @returns the number of genotypes currently belonging to this group.
 */
int get_group_size( SimData* d, int group_id) {
	AlleleMatrix* m = d->m;
	int size = 0;
	int i;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				++size;
			}
		}		
		
		if (m->next == NULL) {
			return size;
		} else {
			m = m->next;
		}
	}
}

/** Gets a shallow copy of the genes/alleles of each member of the group. Not 
 * guaranteed to be null-terminated. 
 * @see get_group_size()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group we want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing pointers to the genes of each member of the group.
 * The vector itself is on the heap and should be freed, but its contents are only
 * shallow copies that should not be freed.
 */
char** get_group_genes( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
	char** genes;
	if (group_size > 0) {
		genes = get_malloc(sizeof(char*) * group_size);
	} else {
		genes = get_malloc(sizeof(char*) * get_group_size( d, group_id ));
	}
	int i, genes_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				genes[genes_i] = m->alleles[i];
				++genes_i;
			}
		}		
		
		if (m->next == NULL) {
			return genes;
		} else {
			m = m->next;
		}
	}	
}

/** Gets a shallow copy of the names of each member of the group.
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing pointers to the names of each member of the group.
 * The vector itself is on the heap and should be freed, but its contents are only
 * shallow copies that should not be freed.
 */
char** get_group_names( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
	char** names;
	if (group_size > 0) {
		names = get_malloc(sizeof(char*) * group_size);
	} else {
		names = get_malloc(sizeof(char*) * get_group_size( d, group_id ));
	}
	int i, names_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				names[names_i] = m->subject_names[i];
				++names_i;
			}
		}		
		
		if (m->next == NULL) {
			return names;
		} else {
			m = m->next;
		}
	}	
}

/** Gets the ids of each member of the group.
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the ids of each member of the group.
 * The vector itself is on the heap and should be freed.
 */
unsigned int* get_group_ids( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
	unsigned int* gids;
	if (group_size > 0) {
		gids = get_malloc(sizeof(unsigned int) * group_size);
	} else {
		gids = get_malloc(sizeof(unsigned int) * get_group_size( d, group_id ));
	}
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				gids[ids_i] = m->ids[i];
				++ids_i;
			}
		}		
		
		if (m->next == NULL) {
			return gids;
		} else {
			m = m->next;
		}
	}	
}

/** Gets the indexes (0-based, from the start of the linked list in the SimData)
 * of each member of the group.
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the indexes of each member of the group.
 * The vector itself is on the heap and should be freed.
 */
unsigned int* get_group_indexes(SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
	unsigned int* gis;
	if (group_size > 0) {
		gis = get_malloc(sizeof(unsigned int) * group_size);
	} else {
		gis = get_malloc(sizeof(unsigned int) * get_group_size( d, group_id ));
	}
	int i, total_i = 0, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i, ++total_i) {
			if (m->groups[i] == group_id) {
				gis[ids_i] = total_i;
				++ids_i;
			}
		}		
		
		if (m->next == NULL) {
			return gis;
		} else {
			m = m->next;
		}
	}	
}

/** Gets the breeding values/GEBVs/fitnesses of each member of the group
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the breeding values of each member of the group.
 * The vector itself is on the heap and should be freed.
 */
double* get_group_bvs( SimData* d, int group_id, int group_size) {
	double* bvs;
	if (group_size > 0) {
		bvs = get_malloc(sizeof(double) * group_size);
	} else {
		bvs = get_malloc(sizeof(double) * get_group_size( d, group_id ));
	}
	
	DecimalMatrix dm_bvs = calculate_group_bvs(d, group_id);
	
	for (int i = 0; i < dm_bvs.cols; ++i) {
		bvs[i] = dm_bvs.matrix[0][i];
	}
	
	delete_dmatrix(&dm_bvs);
	
	return bvs;
}

/** Gets the ids of either the first or second parent of each member
 * of the group.
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @param parent 1 to get the first parent of each group member, 2 to get the second.
 * Raises an error if this parameter is not either of those values.
 * @returns a vector containing the id of either the first or second parent of
 * each member of the group.
 * The vector itself is on the heap and should be freed.
 */
unsigned int* get_group_parent_IDs( SimData* d, int group_id, int group_size, int parent) {
	if (!(parent == 1 || parent == 2)) {
		error("Value error: `parent` must be 1 or 2.");
	}
	--parent;
	
	AlleleMatrix* m = d->m;
	unsigned int* pids;
	if (group_size > 0) {
		pids = get_malloc(sizeof(unsigned int) * group_size);
	} else {
		pids = get_malloc(sizeof(unsigned int) * get_group_size( d, group_id ));
	}
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				pids[ids_i] = m->pedigrees[parent][i];
				++ids_i;
			}
		}		
		
		if (m->next == NULL) {
			return pids;
		} else {
			m = m->next;
		}
	}		
}
 
/** Gets the names of either the first or second parent of each member
 * of the group. Names may be NULL.
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @param parent 1 to get the first parent of each group member, 2 to get the second.
 * Raises an error if this parameter is not either of those values.
 * @returns  a vector containing pointers to the names either the first or second
 * parent of each member of the group. The vector itself is on the heap and should be 
 * freed, but its contents are only shallow copies that should not be freed.
 */
char** get_group_parent_names( SimData* d, int group_id, int group_size, int parent) {
	if (!(parent == 1 || parent == 2)) {
		error("Value error: `parent` must be 1 or 2.");
	}
	--parent;
	
	AlleleMatrix* m = d->m;
	char** pnames;
	if (group_size > 0) {
		pnames = get_malloc(sizeof(char*) * group_size);
	} else {
		pnames = get_malloc(sizeof(char*) * get_group_size( d, group_id ));
	}
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				pnames[ids_i] = get_name_of_id(d->m, m->pedigrees[parent][i]);
				++ids_i;
			}
		}		
		
		if (m->next == NULL) {
			return pnames;
		} else {
			m = m->next;
		}
	}	
}

/** Gets the full pedigree string (as per save_group_full_pedigree() )
 * of each member of the group.
 *
 * This function is not fast, or particularly memory-efficient. To avoid recursive
 * string-building/concatenation because that's pretty tricky in C, this just 
 * calls save_group_full_pedigree() to a temporary file, then reads that file 
 * back in to make a list of strings. Contact the package maintainers if you'd 
 * benefit from a faster version of this function, otherwise I probably won't bother 
 * improving this. 
 *
 * If you just want to have the full pedigrees for further analysis, consider
 * save_group_full_pedigree(). With the current implementation that one will 
 * be faster than this anyway.
 *
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_IDs()
 * @see get_group_parent_names()
 * @see save_group_full_pedigree()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns  a vector containing pointers to the full pedigree string of each 
 * member of the group. Both the vector and its contents are dynamically 
 * allocated on the heap, and so should both be freed. This is because, unlike
 * other data-getter functions, this data is not stored but must be generated 
 * specifically to answer this function call.
 */
char** get_group_pedigrees( SimData* d, int group_id, int group_size) {
	char* fname = "gS_gpptmp";
	FILE* fp = fopen(fname, "w");
	save_group_full_pedigree(fp, d, group_id);
	fclose(fp);
	
	FILE* fp2;
	if ((fp2 = fopen(fname, "r")) == NULL) {
		error( "Failed to use temporary file.\n");
	}
	
	// Create the list that we will return
	if (group_size <= 0) {
		group_size = get_group_size( d, group_id );
	}
	char** gp_ped = get_malloc(sizeof(char*) * group_size);
	
	// read one line at a time
	//size_t n;
	//int line_len;
	unsigned int size;
	unsigned int index;
	int nextc;
	for (int i = 0; i < group_size; ++i) {
		// getline is not available in MinGW it looks like (AUG 2021)
		/*gp_ped[i] = NULL;
		if ((line_len = getline(&(gp_ped[i]), &n, fp2)) == -1) {
			error("Failed to get %d th pedigree\n", i);
		} 
		// remove the newline character
		if (gp_ped[i][line_len - 1] == '\n') {
			gp_ped[i][line_len - 1] = '\0';
		}*/
		
		// a not-very-size-efficient, fgets-based line getter
		size = 50;
		index = 0;
		gp_ped[i] = get_malloc(sizeof(char) * size);
		while ((nextc = fgetc(fp)) != '\n' && nextc != EOF) {
			gp_ped[i][index] = nextc;
			++index;

			if (index >= size) {
				size *= 2;
				char* temp = realloc(gp_ped[i], size);
				if (temp == NULL) {
					free(gp_ped[i]);
					fprintf(stderr, "Memory allocation of size %u failed.\n", size);
				} else {
					gp_ped[i] = temp;
				}
			}
		}
		gp_ped[i][index] = '\0';
	}
	
	fclose(fp2);
	remove(fname);
	
	return gp_ped;
}