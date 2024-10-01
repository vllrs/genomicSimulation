#ifndef SIM_OPERATIONS
#define SIM_OPERATIONS
#include "sim-operations.h"
/* genomicSimulationC v0.2.5.04 - last edit 1 October 2024 */

/** Default parameter values for GenOptions, to help with quick scripts and prototypes.
 *
 * @shortnamed{BASIC_OPT}
 */
const gsc_GenOptions GSC_BASIC_OPT = {
	.will_name_offspring = GSC_FALSE,
	.offspring_name_prefix = NULL,
	.family_size = 1,
	.will_track_pedigree = GSC_FALSE,
	.will_allocate_ids = GSC_TRUE,
	.filename_prefix = NULL,
	.will_save_pedigree_to_file = GSC_FALSE,
    .will_save_bvs_to_file = GSC_NO_EFFECTSET,
	.will_save_alleles_to_file = GSC_FALSE,
	.will_save_to_simdata = GSC_TRUE
};

/** Replace calls to malloc direct with this function.
 *
 * @param size size of space to be allocated. Treat like the parameter of a
 * regular malloc call.
 * @param exitonfail if truthy, if memory allocation fails it will terminate execution.
 * if falsy, on memory allocation failure it will only print an error message and
 * return NULL (useful, for example, if the calling function wants to add its own error
 * message, or return partial results even if it fails here).
 * @returns pointer to the allocated space.
 */
static void* gsc_malloc_wrap(const unsigned int size, char exitonfail) {
    if (size == 0) {
        warning( "0 memory allocation requested.\n");
        return NULL;
    }
	void* v = GSC_MALLOC(size);
	if (v == NULL) {
        if (exitonfail) {
            error( "Memory allocation failed. Exiting.\n");
        } else {
            warning( "Memory allocation failed.\n");
        }
	}
	return v;
}

/** Creator for an empty gsc_AlleleMatrix object of a given size. Includes memory
 * allocation for `n_genotypes` worth of `.alleles`.
 *
 * @param n_markers number of rows/markers to create
 * @param n_labels number of custom labels to create
 * @param labelDefaults an array of length [n_labels] of the default value to pre-fill for each custom label.
 * Can be null if n_labels == 0.
 * @param n_genotypes number of individuals to create. This includes filling the first
 * n_genotypes entries of .alleles with heap char* of length n_markers, so that the
 * alleles for these can be added without further memory allocation.
 * @returns pointer to the empty created gsc_AlleleMatrix
 */
gsc_AlleleMatrix* gsc_create_empty_allelematrix(const int n_markers, const int n_labels, const int* labelDefaults, const int n_genotypes) {
    gsc_AlleleMatrix* m = gsc_malloc_wrap(sizeof(gsc_AlleleMatrix),GSC_TRUE);

	m->n_genotypes = n_genotypes;
	m->n_markers = n_markers;
    m->n_labels = n_labels;
	//m->alleles = gsc_malloc_wrap(sizeof(char*) * CONTIG_WIDTH);
	for (int i = 0; i < n_genotypes; ++i) {
        m->alleles[i] = gsc_malloc_wrap(sizeof(char) * (n_markers<<1),GSC_TRUE);
		memset(m->alleles[i], 0, sizeof(char) * (n_markers<<1));
		//m->ids[i] = 0;
	}
	memset(m->alleles + n_genotypes, 0, sizeof(char*) * (CONTIG_WIDTH - n_genotypes)); // setting the pointers to NULL

    if (n_labels > 0) {
        m->labels = gsc_malloc_wrap(sizeof(int*) * n_labels,GSC_TRUE);
        for (int i = 0; i < n_labels; ++i) {
            m->labels[i] = gsc_malloc_wrap(sizeof(int) * CONTIG_WIDTH,GSC_TRUE);
            for (int j = 0; j < CONTIG_WIDTH; ++j) {
                m->labels[i][j] = labelDefaults[i];
            }
        }
    } else if (n_labels == 0) {
        m->labels = NULL;
    } else {
        warning( "Invalid negative number of labels provided to gsc_create_empty_allelematrix");
        m->labels = NULL;
    }

    memset(m->ids, 0, sizeof(gsc_PedigreeID) * CONTIG_WIDTH);
    memset(m->pedigrees[0], 0, sizeof(gsc_PedigreeID) * CONTIG_WIDTH);
    memset(m->pedigrees[1], 0, sizeof(gsc_PedigreeID) * CONTIG_WIDTH);
    memset(m->groups, 0, sizeof(gsc_GroupNum) * CONTIG_WIDTH);
	memset(m->names, 0, sizeof(char*) * CONTIG_WIDTH); // setting the pointers to NULL

	m->next = NULL;

	return m;
}

/** Creator for an empty gsc_SimData object on the heap. This is the main struct
 * that will contain/manage simulation data.
 *
 * @shortnamed{create_empty_simdata}
 *
 * @returns pointer to the empty created gsc_SimData
 */
gsc_SimData* gsc_create_empty_simdata() {
    gsc_SimData* d = gsc_malloc_wrap(sizeof(gsc_SimData),GSC_TRUE);
    d->n_labels = 0;
    d->label_ids = NULL;
    d->label_defaults = NULL;
    d->genome.n_markers = 0;
    d->genome.marker_names = NULL;
    d->genome.names_alphabetical = NULL;
    d->genome.n_maps = 0;
    d->genome.map_ids = NULL;
    d->genome.maps = NULL;
	d->m = NULL;
    d->n_eff_sets = 0;
    d->e = NULL;
    ;
    d->current_id = GSC_NO_PEDIGREE;
    d->n_groups = 0;
	return d;
}

/** Clear a gsc_SimData object on the heap.
 *
 *  Has the effects of gsc_delete_simdata followed by gsc_create_empty_simdata,
 *  but guarantees the use of the same memory location for the
 *  new gsc_SimData.
 *
 * @shortnamed{clear_simdata}
 *
 *  @param d pointer to the gsc_SimData to be cleared.
 */
void gsc_clear_simdata(gsc_SimData* d) {
    // Free label defaults
    if (d->n_labels > 0) {
        if (d->label_ids != NULL) {
            GSC_FREE(d->label_ids);
        }
        if (d->label_defaults != NULL) {
            GSC_FREE(d->label_defaults);
        }
    }

    // Free other details
    gsc_delete_genome(&(d->genome));
    for (int i = 0; i < d->n_eff_sets; ++i) {
        gsc_delete_effect_matrix(&(d->e[i]));
    }
    if (d->n_eff_sets > 0) {
        GSC_FREE(d->eff_set_ids);
        GSC_FREE(d->e);
    }
    gsc_delete_allele_matrix(d->m);

    // Clear all values
    d->n_labels = 0;
    d->label_ids = NULL;
    d->label_defaults = NULL;
    d->genome.n_markers = 0;
    d->genome.marker_names = NULL;
    d->genome.n_maps = 0;
    d->genome.map_ids = NULL;
    d->genome.maps = NULL;
    d->m = NULL;
    d->n_eff_sets = 0;
    d->e = NULL;
    d->current_id = GSC_NO_PEDIGREE;
    d->n_groups = 0;
}
/*------------------------Supporter Functions--------------------------------*/

/** Opens a table file and reads the number of columns and rows
 * (including headers) separated by `sep` into a gsc_TableSize struct that is
 * returned.
 *
 * If the file fails to open, the simulation exits.
 *
 * Rows must be either empty or have same number of columns as the first.
 * Empty rows are not counted.
 *
 * If a row with an invalid number of columns is found, the number of
 * columns in the return value is set to 0. number of rows in the return
 * value in this case is arbitrary.
 *
 * @param filename the path/name to the table file whose dimensions we want
 * @param sep the character that separates columns in the file eg tab
 * @returns gsc_TableSize struct with .num_columns and .num_rows filled. These
 * counts include header rows/columns and exclude blank rows.
 */
struct gsc_TableSize gsc_get_file_dimensions(const char* filename, const char sep) {
	struct gsc_TableSize details;
	details.num_columns = 0;
	details.num_rows = 0;

	FILE* fp;
	int c; // this is used to store the output of fgetc i.e. the next character in the file
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	c = fgetc(fp);

	while (c != EOF && c != '\n') {
		R_CheckUserInterrupt();
		if (c == sep) {
			details.num_columns += 1; // add count for columns of form [colname]sep
		}
		c = fgetc(fp);
	}

	details.num_columns += 1; // add another column that was bounded by sep[colname][EOF or \n]
	details.num_rows = 1; // we successfully got the first row

	// now get all the rows. What we care about in the rows is the number of them
	c = fgetc(fp);
	int sep_count = 0; // for each row, count the columns to make sure they match and the file is valid
	int has_length = GSC_FALSE;
	while (c != EOF) {
		R_CheckUserInterrupt();
		if (c == '\n') {
			details.num_rows += 1; // add count for columns of form [colname]sep

			// check we have right number of columns and reset counter
			if (has_length && sep_count != details.num_columns-1) {
				// we have a bad number of columns
				details.num_columns = 0;
				fclose(fp);
				error( "Bad columns on row %d\n", details.num_rows + 1);
			}
			sep_count = 0;
			has_length = GSC_FALSE;

		} else if (c == sep) {
			sep_count += 1;
		} else if (has_length == GSC_FALSE) {
			has_length = GSC_TRUE;
		}
		c = fgetc(fp);
	}
	if (has_length) {
		details.num_rows += 1; // for the last row before EOF
	}

	fclose(fp);
	return details;
}

/** Binary search through list of unsigned integers
 *
 * Returns the located index in an array of integers where the integer
 * is `target`. Returns -1 if no match was found.
 * @see gsc_get_from_unordered_str_list()
 * @see gsc_get_from_ordered_pedigree_list()
 * @see gsc_get_from_ordered_str_list()
 *
 * The list is assumed to be sorted in ascending order. Only integers
 * >0 are considered valid; entries of 0 are considered empty and can
 * be located at any point in the list.
 *
 * It uses a binary search method, but has to widen its search
 * both directions if the desired midpoint has value 0.
 *
 * @param target the integer to be located
 * @param list an array of integers to search, with at least [list_len] entries
 * @param list_len length of the array of integers to search
 * @returns Index in `list` where we find the same integer as
 * `target`, or -1 if no match is found.
 */
int gsc_get_from_ordered_uint_list(const unsigned int target, const unsigned int listLen, const unsigned int* list) {
    unsigned int first = 0, last = listLen - 1;
    int index = (first + last) / 2;
    while (list[index] != target && first <= last) {
        if (list[index] == 0) {
            int lookahead = 1;
            while(1) {
                if (index+lookahead <= last && list[index+lookahead] != 0) {
                    if (list[index+lookahead] == target) {
                        return index+lookahead;
                    } else if (list[index+lookahead] < target) {
                        first = index+lookahead + 1;
                        break;
                    } else {
                        last = index - 1;
                        break;
                    }
                } else if (index-lookahead <= last && list[index-lookahead] != 0) {
                    if (list[index-lookahead] == target) {
                        return index-lookahead;
                    } else if (list[index-lookahead] < target) {
                        first = index + 1;
                        break;
                    } else {
                        last = index-lookahead - 1;
                        break;
                    }
                }
                ++lookahead;
                if (index+lookahead <= last || index-lookahead >= first) {
                    // failed to find any nonzeros between first and last
                    return -1;
                }
            }

        } else { // No need to dodge 0. Normal binary search.
            if (list[index] == target) {
                return index;
            } else if (list[index] < target) {
                first = index + 1;
            } else {
                last = index - 1;
            }

        }
        // index has been updated, no matter the branch.
        index = (first + last) / 2;
    }

    if (first > last) {
        return -1;
    }
    return index;
}

/** Binary search through a list of PedigreeIDs
 *
 * Returns the located index in an array of gsc_PedigreeIDs where the gsc_PedigreeID
 * is `target`. Returns -1 if no match was found.
 * @see gsc_get_from_unordered_str_list()
 * @see gsc_get_from_ordered_uint_list()
 * @see gsc_get_from_ordered_str_list()
 *
 * The list is assumed to be sorted in ascending order. Only IDs
 * >0 are considered valid; entries of 0 are considered empty and can
 * be located at any point in the list.
 *
 * It uses a binary search method, but has to widen its search
 * both directions if the desired midpoint has value 0.
 *
 * @param target the integer to be located
 * @param list the array of integers to search, with at least [list_len] entries
 * @param list_len length of the array of PedigreeIDs to search
 * @returns Index in `list` where we find the same integer as
 * `target`, or -1 if no match is found.
 */
int gsc_get_from_ordered_pedigree_list(const gsc_PedigreeID target, const unsigned int listLen, const gsc_PedigreeID* list) {
    unsigned int first = 0, last = listLen - 1;
    int index = (first + last) / 2;
    while (list[index].id != target.id && first <= last) {
        if (list[index].id == GSC_NO_PEDIGREE.id) {
            int lookahead = 1;
            while(1) {
                if (index+lookahead <= last && list[index+lookahead].id != GSC_NO_PEDIGREE.id) {
                    if (list[index+lookahead].id == target.id) {
                        return index+lookahead;
                    } else if (list[index+lookahead].id < target.id) {
                        first = index+lookahead + 1;
                        break;
                    } else {
                        last = index - 1;
                        break;
                    }
                } else if (index-lookahead <= last && list[index-lookahead].id != GSC_NO_PEDIGREE.id) {
                    if (list[index-lookahead].id == target.id) {
                        return index-lookahead;
                    } else if (list[index-lookahead].id < target.id) {
                        first = index + 1;
                        break;
                    } else {
                        last = index-lookahead - 1;
                        break;
                    }
                }
                ++lookahead;
                if (index+lookahead <= last || index-lookahead >= first) {
                    // failed to find any nonzeros between first and last
                    return -1;
                }
            }

        } else { // No need to dodge 0. Normal binary search.
            if (list[index].id == target.id) {
                return index;
            } else if (list[index].id < target.id) {
                first = index + 1;
            } else {
                last = index - 1;
            }

        }
        // index has been updated, no matter the branch.
        index = (first + last) / 2;
    }

    if (first > last) {
        return -1;
    }
    return index;
}

/** Linear search through a list of strings
 *
 * Returns the first located index in an array of strings where the string
 * is the same as the string `target`. Returns -1 if no match was found.
 * @see gsc_get_from_ordered_uint_list()
 * @see gsc_get_from_ordered_pedigree_list()
 * @see gsc_get_from_ordered_str_list()
 *
 * The list of strings is not assumed to be sorted.
 *
 * @param target the string to be located
 * @param list the array of strings to search, with at least [list_len] entries
 * @param list_len length of the array of strings to search
 * @returns Index in `list` where we find the same string as
 * `target`, or -1 if no match is found.
 */
int gsc_get_from_unordered_str_list(const char* target, const int listLen, const char** list) {
    for (int i = 0; i < listLen; ++i) {
		if (strcmp(list[i], target) == 0) {
			return i;
		}
	}
	return -1; // did not find a match.
}

/** Binary search through a list of strings
 *
 * Returns the first located index in an array of strings where the string
 * is the same as the string `target`. Returns -1 if no match was found.
 * @see gsc_get_from_ordered_uint_list()
 * @see gsc_get_from_ordered_pedigree_list()
 * @see gsc_get_from_unordered_str_list()
 *
 * The list of strings is assumed to be sorted in alphabetical order.
 *
 * @param target the string to be located
 * @param list the array of strings to search, with at least [list_len] entries
 * @param list_len length of the array of strings to search
 * @returns Index in `list` where we find the same string as
 * `target`, or -1 if no match is found.
 */
int gsc_get_from_ordered_str_list(const char* target, const int listLen, const char** list) {
    unsigned int first = 0, last = listLen - 1;
    int index = (first + last) / 2;
    int comparison = strcmp(target,list[index]);
    while (comparison != 0 && first <= last) {
        if (comparison == 0) {
            return index;
        } else if (comparison < 0) {
            first = index + 1;
        } else {
            last = index - 1;
        }

        // index has been updated, no matter the branch.
        index = (first + last) / 2;
        comparison = strcmp(target, list[index]);
    }

    if (first > last) {
        return -1;
    }
    return index;
}


/** Produce a random ordering of the first n elements in an array of integers
 * using a (partial) Fisher-Yates shuffle.
 *
 * Modified from https://benpfaff.org/writings/clc/shuffle.html
 *
 * @param d gsc_SimData, only used for pointer to random number generator
 * @param sequence the array
 * @param total_n length of the array
 * @param n_to_shuffle After calling this function, the first n_to_shuffle
 * integers in the array will be randomly ordered by a Fischer-Yates shuffle. 
 * Every entry in the array could end up at any position, but the post-shuffle
 * positions have only been calculated for the first 'n_to_shuffle' entries.
 */
void shuffle_up_to( unsigned int* sequence, const unsigned int total_n, const unsigned int n_to_shuffle) {
	if (n_to_shuffle > 1) {
        unsigned int maxi = total_n > n_to_shuffle ? n_to_shuffle - 1 : total_n - 1;
		unsigned int i;
        for (i = 0; i < maxi; ++i) {
			// items before i are already shuffled
            unsigned int j = i + round(unif_rand() * (total_n - i - 1));

			// add the next chosen value to the end of the shuffle
            unsigned int t = sequence[j];
			sequence[j] = sequence[i];
			sequence[i] = t;
		}
	}
}

/** Fills the designated section of the `.names` array in an
 * gsc_AlleleMatrix with the pattern "`prefix`index".
 *
 * In future this function could be expanded to allow for different naming formats.
 *
 * The index is padded with zeros depending on the size of
 * `a->n_genotypes`.
 *
 * @param a pointer to the gsc_AlleleMatrix whose `.names` to modify
 * @param prefix the prefix to add to the suffix to make the new genotype name
 * @param suffix suffixes start at this value and increment for each additional name
 * @param from_index the new names are added to this index and all those following it in this
 * gsc_AlleleMatrix.
*/
static void gsc_set_names(gsc_AlleleMatrix* a, const char* prefix, const int suffix, const int from_index) {
	char sname[NAME_LENGTH];
	char format[NAME_LENGTH];
	if (prefix == NULL) {
		// make it an empty string instead, so it is not displayed as (null)
		prefix = "";
	}
	// use sname to save the number of digits to pad by:
	sprintf(sname, "%%0%dd", gsc_get_integer_digits(a->n_genotypes - from_index));  // Creates: %0[n]d
	sprintf(format, "%s%s", prefix, sname);

    int livingsuffix = suffix;
    ++livingsuffix;
	for (int i = from_index; i < a->n_genotypes; ++i) {
		// clear name if it's pre-existing
		if (a->names[i] != NULL) {
			GSC_FREE(a->names[i]);
		}

		// save new name
        sprintf(sname, format, livingsuffix);
        a->names[i] = gsc_malloc_wrap(sizeof(char) * (strlen(sname) + 1),GSC_TRUE);
		strcpy(a->names[i], sname);

        ++livingsuffix;
	}
}

/** Initialises a new custom label.
 *
 * Creates a new custom label on every genotype currently and in future belonging
 * to the gsc_SimData. The value of the label is set as `setTo` for every genotype.
 *
 * @shortnamed{create_new_label}
 *
 * @param d pointer to the `gsc_SimData` whose child `gsc_AlleleMatrix`s will be given the new label.
 * @param setTo the value to which every genotype's label is initialised.
 * @returns the label id of the new label
*/
gsc_LabelID gsc_create_new_label(gsc_SimData* d, const int setTo) {
    // Add new label default
    if (d->n_labels == 0) {
        d->label_ids = gsc_malloc_wrap(sizeof(gsc_LabelID) * 1,GSC_TRUE);
        d->label_ids[0] = (gsc_LabelID){.id=1};

        d->label_defaults = gsc_malloc_wrap(sizeof(int) * 1,GSC_TRUE);
        d->label_defaults[0] = setTo;

    } else if (d->n_labels > 0) {

        gsc_LabelID* new_label_ids;
        if (d->label_ids != NULL) {
            new_label_ids = gsc_malloc_wrap(sizeof(gsc_LabelID) * (d->n_labels + 1),GSC_TRUE);
            memcpy(new_label_ids,d->label_ids,sizeof(gsc_LabelID)*d->n_labels);
            new_label_ids[d->n_labels] = gsc_get_new_label_id(d);
            GSC_FREE(d->label_ids);

        } else { // d->label_ids == NULL
            // If the other labels do not have identifiers, they're corrupted and
            // deserve to be destroyed.
            new_label_ids = gsc_malloc_wrap(sizeof(gsc_LabelID) * 1,GSC_TRUE);
            d->n_labels = 0;
            new_label_ids[d->n_labels] = gsc_get_new_label_id(d);
        }
        d->label_ids = new_label_ids;

        int* new_label_defaults = gsc_malloc_wrap(sizeof(int) * (d->n_labels + 1),GSC_TRUE);
        if (d->label_defaults != NULL) {
            for (int i = 0; i < d->n_labels; ++i) {
                new_label_defaults[i] = d->label_defaults[i];
            }
            GSC_FREE(d->label_defaults);
        } else if (d->n_labels > 0) {
            memset(new_label_defaults, 0, sizeof(int) * d->n_labels);
        }
        new_label_defaults[d->n_labels] = setTo;
        d->label_defaults = new_label_defaults;

    } else {
        warning( "Labels malformed; gsc_SimData may be corrupted.\n");
        return (gsc_LabelID){.id=GSC_UNINIT};
    }
    d->n_labels += 1;

    // Set all values of that label to the default
    gsc_AlleleMatrix* m = d->m;
    int warned = GSC_FALSE;
    do {
        // Do we need to destroy the extant label table? happens if label_ids were missing and we discarded them
        if (m->n_labels != d->n_labels - 1 && m->labels != NULL) {
            for (int i = 0; i < m->n_labels; ++i) {
                GSC_FREE(m->labels[i]);
            }
            GSC_FREE(m->labels);
            m->labels = NULL;
        }

        m->n_labels = d->n_labels;

        // Consider the case when we need to expand the label list
        if (m->n_labels > 1 && m->labels != NULL) {
            int newLabel = m->n_labels - 1;

            // Create label list
            int** oldLabelList = m->labels;
            m->labels = gsc_malloc_wrap(sizeof(int*) * m->n_labels,GSC_TRUE);
            for (int i = 0; i < m->n_labels - 1; ++i) {
                m->labels[i] = oldLabelList[i];
            }
            m->labels[newLabel] = gsc_malloc_wrap(sizeof(int) * CONTIG_WIDTH,GSC_TRUE);
            GSC_FREE(oldLabelList);

            // Set labels
            if (setTo == 0) {
                memset(m->labels[newLabel], 0, sizeof(int) * CONTIG_WIDTH);
            } else {
                for (int i = 0; i < CONTIG_WIDTH; ++i) {
                    m->labels[newLabel][i] = setTo;
                }
            }

        // Consider the case we need to initialise the label list
        } else if (m->n_labels == 1 && m->labels == NULL) {
            // Create the label list
            m->labels = gsc_malloc_wrap(sizeof(int*) * 1,GSC_TRUE);
            m->labels[0] = gsc_malloc_wrap(sizeof(int) * CONTIG_WIDTH,GSC_TRUE);

            // Set labels
            if (setTo == 0) {
                memset(m->labels[0], 0, sizeof(int) * CONTIG_WIDTH);
            } else {
                for (int i = 0; i < CONTIG_WIDTH; ++i) {
                    m->labels[0][i] = setTo;
                }
            }

        } else if (!warned) {
            warning( "Unable to create new label for all genotypes; gsc_SimData may be corrupted.\n");
            warned = GSC_TRUE;
        }

    } while ((m = m->next) != NULL);
    return d->label_ids[d->n_labels - 1];
}

/** Set the default value of a custom label
 *
 * Sets the default (birth) value of the custom label that has index `whichLabel` to the
 * value `newDefault`.
 *
 * @shortnamed{change_label_default}
 *
 * @param d pointer to the `gsc_SimData` containing the genotypes and labels to be relabelle
 * @param whichLabel the label id of the relevant label.
 * @param newDefault the value to which the appropriate label's default will be set.
*/
void gsc_change_label_default(gsc_SimData* d, const gsc_LabelID whichLabel, const int newDefault) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return;
    }
    d->label_defaults[labelIndex] = newDefault;
}

/** Set the values of a custom label
 *
 * Sets the values of the custom label that has index `whichLabel` to the
 * value `setTo`. Depending on the value of `whichGroup`, the function
 * will modify the label of every single genotype in the simulation, or
 * just modify the label of every member of a given group.
 *
 * @shortnamed{change_label_to}
 *
 * @param d pointer to the `gsc_SimData` containing the genotypes and labels to be relabelled
 * @param whichGroup 0 to modify the relevant labels of all extant genotypes,
 * or a positive integer to modify the relevant labels of all members of group `whichGroup`.
 * @param whichLabel the label id of the relevant label.
 * @param setTo the value to which the appropriate labels will be set.
*/
void gsc_change_label_to(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel, const int setTo) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return;
    }
    // Risks: if m->labels or m->labels[i] don't exist for labels where they should,
    // will get some out of bounds accesses.

    gsc_AlleleMatrix* m = d->m;
    if (whichGroup.num != GSC_NO_GROUP.num) { // set the labels of group members
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                if (m->groups[i].num == whichGroup.num) {
                    m->labels[labelIndex][i] = setTo;
                }
            }

        } while ((m = m->next) != NULL);

    } else { // whichGroup == 0 so set the labels of all genotypes
        do {

            if (setTo == 0) {
                memset(m->labels[labelIndex], 0, sizeof(int) * m->n_genotypes);
            } else {
                for (int i = 0; i < m->n_genotypes; ++i) {
                    m->labels[labelIndex][i] = setTo;
                }
            }

        } while ((m = m->next) != NULL);
    }
}

/** Increment the values of a custom label
 *
 * Increments the values of the custom label that has index `whichLabel` by the
 * value `byValue`. Depending on the value of `whichGroup`, the function
 * will modify the label of every single genotype in the simulation, or
 * just modify the label of every member of a given group.
 *
 * @shortnamed{change_label_by_amount}
 *
 * @param d pointer to the `gsc_SimData` containing the genotypes and labels to be relabelled
 * @param whichGroup 0 to modify the relevant labels of all extant genotypes,
 * or a positive integer to modify the relevant labels of all members of group `whichGroup.
 * @param whichLabel the label id of the relevant label.
 * @param byValue the value by which the appropriate labels will be incremented. For example,
 * a value of 1 would increase all relevant labels by 1, a value of -2 would subtract 2 from
 * each relevant label.
*/
void gsc_change_label_by_amount(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel, const int byValue) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return;
    }
    // Risks: if m->labels or m->labels[i] don't exist for labels where they should,
    // will get some out of bounds accesses.

    gsc_AlleleMatrix* m = d->m;
    if (whichGroup.num != GSC_NO_GROUP.num) { // set the labels of group members
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                if (m->groups[i].num == whichGroup.num) {
                    m->labels[labelIndex][i] += byValue;
                }
            }

        } while ((m = m->next) != NULL);

    } else { // whichGroup == 0 so set the labels of all genotypes
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                m->labels[labelIndex][i] += byValue;
            }

        } while ((m = m->next) != NULL);
    }

}

/** Copy a vector of integers into a custom label
 *
 * Sets values of the custom label that has index `whichLabel` to the
 * contents of the array `values`. The genotypes with the `n_values` contiguous
 * simulation indexes starting with `startIndex` in group `whichGroup` are the ones given
 * those labels. If `whichGroup` is 0, then it sets the labels of the genotypes with
 * global indexes starting at `startIndex`.
 *
 * If `n_values` is longer than the number of genotypes in the right group
 * after `startIndex`, then the extra values are ignored.
 *
 * @shortnamed{change_label_to_values}
 *
 * @param d pointer to the `gsc_SimData` containing the genotypes and labels to be relabelled
 * @param whichGroup 0 to set the label of the genotypes with global indexes between
 * `startIndex` and `startIndex + n_values`, or a positive integer to set the label
 * of the `startIndex`th to `startIndex + n_values`th members of group `whichGroup`.
 * @param startIndex the first index of the group to set to a value
 * @param whichLabel the label id of the relevant label.
 * @param n_values length (number of entries) of the array `values`
 * @param values vector of integers, of length at least [n_values],
 * to paste into the chosen custom label of the chosen genotypes.
*/
void gsc_change_label_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const int startIndex, const gsc_LabelID whichLabel,
                          const int n_values, const int* values) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return;
    }

    gsc_AlleleMatrix* m = d->m;
    int currentIndex = 0;
    if (whichGroup.num != GSC_NO_GROUP.num) { // set the labels of group members
        // First scan through to find firstIndex
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                if (m->groups[i].num == whichGroup.num) {
                    // Update label if it is between startIndex and startIndex + n_values
                    if (currentIndex >= startIndex) {
                        m->labels[labelIndex][i] = values[currentIndex - startIndex];
                    }
                    currentIndex++;
                    if (currentIndex - startIndex >= n_values) {
                        return;
                    }
                }
            }

        } while ((m = m->next) != NULL);

    } else { // whichGroup == 0 so set the labels of all genotypes
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                // Update label if it is between startIndex and startIndex + n_values
                if (currentIndex >= startIndex) {
                    m->labels[labelIndex][i] = values[currentIndex - startIndex];
                }
                currentIndex++;
                if (currentIndex - startIndex >= n_values) {
                    return;
                }
            }

        } while ((m = m->next) != NULL);
    }
}

/** Copy a vector of strings into the genotype name field
 *
 * Sets genotype names to the contents of the array `values`. The genotypes with
 * the `n_values` contiguous simulation indexes starting with `startIndex` in
 * group `whichGroup` are the ones given those names. If `whichGroup` is 0,
 * then it sets the names of the genotypes with
 * global indexes starting at `startIndex`.
 *
 * If `n_values` is longer than the number of genotypes in the right group
 * after `startIndex`, then the extra values are ignored.
 *
 * A deep copy of each name in `values` is made. This way, `values` can be a pointer to
 * names of existing genotypes (say, if you are creating clones and want them to have
 * the same names), and there will be no issues if only one of the genotypes sharing the
 * name is deleted.
 *
 * @shortnamed{change_names_to_values}
 *
 * @param d pointer to the `gsc_SimData` containing the genotypes to be renamed
 * @param whichGroup 0 to set the names of the genotypes with global indexes between
 * `startIndex` and `startIndex + n_values`, or a positive integer to set the names
 * of the `startIndex`th to `startIndex + n_values`th members of group `whichGroup`.
 * @param startIndex the first index of the group to set to a value
 * @param n_values length (number of entries) of `values`
 * @param values vector of strings containing at least [n_values] strings,
 *  to paste into the name field of the chosen genotypes.
*/
void gsc_change_names_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const int startIndex, const int n_values, const char** values) {
    // this will be much improved once we can hash our names.

    gsc_AlleleMatrix* m = d->m;
    int currentIndex = 0;
    if (whichGroup.num != GSC_NO_GROUP.num) { // set the names of group members
        // First scan through to find firstIndex
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                if (m->groups[i].num == whichGroup.num) {
                    // Update name if index is between startIndex and startIndex + n_values
                    if (currentIndex >= startIndex) {
                        // clear name if it's pre-existing
                        if (m->names[i] != NULL) {
                            GSC_FREE(m->names[i]);
                        }

                        // save new name
                        const int whichName = currentIndex - startIndex;
                        m->names[i] = gsc_malloc_wrap(sizeof(char) * (strlen(values[whichName]) + 1),GSC_TRUE);
                        strcpy(m->names[i], values[whichName]);
                    }
                    currentIndex++;
                    if (currentIndex > n_values) {
                        return;
                    }
                }
            }

        } while ((m = m->next) != NULL);

    } else { // whichGroup == 0 so set the names of all genotypes
        do {

            for (int i = 0; i < m->n_genotypes; ++i) {
                // Update name if it is between startIndex and startIndex + n_values
                if (currentIndex >= startIndex) {
                    // clear name if it's pre-existing
                    if (m->names[i] != NULL) {
                        GSC_FREE(m->names[i]);
                    }

                    // save new name
                    const int whichName = currentIndex - startIndex;
                    const int nameLen = strlen(values[whichName]);
                    m->names[i] = gsc_malloc_wrap(sizeof(char) * (nameLen + 1),GSC_TRUE);
                    strncpy(m->names[i], values[whichName], nameLen);
                }
                currentIndex++;
                if (currentIndex > n_values) {
                    return;
                }
            }

        } while ((m = m->next) != NULL);
    }
}

/** Replace all occurences of a given allele with a different symbol representation
 *
 *  Alleles in genomicSimulation are represented by any single character. This function
 *  allows you to replace every instance of some allele with a different character. It may
 *  be used to replace unprintable alleles before saving output (like the null character '\0'
 *  used for a loaded founder genotype's alleles at a marker where no data was provided.)
 *
 *  If the allele is changed to a character which already represents another allele at that marker,
 *  the distinction between those two alleles will be lost.
 *
 *  If no marker name is provided, changes the allele symbol in all markers tracked by the
 *  simulation.
 *
 * @param d SimData on which to perform the operation
 * @param which_marker if null, any occurences of that allele symbol in any marker tracked by the
 * simulation will be changed. If the name of a marker, replace all occurences of that allele symbol at
 * that marker with the new symbol, and do not replace that allele symbol anywhere else.
 * @param from character that currently represents the allele, whose representation is to be changed.
 * @param to character to be the new representation of that allele
 */
void gsc_change_allele_symbol(gsc_SimData* d, const char* which_marker, char from, char to) {
    unsigned int nmarkers = 0;
    unsigned int ngenos = 0;
    unsigned int nalleles = 0;

    unsigned int markeri;
    if (which_marker == NULL) {
        gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,NO_GROUP);
        gsc_GenoLocation loc = gsc_set_bidirectional_iter_to_start(&it);

        while (IS_VALID_LOCATION(loc)) {
            for (unsigned int m = 0; m < d->genome.n_markers; ++m) {
                if (from == loc.localAM->alleles[loc.localPos][m << 1]) {
                    loc.localAM->alleles[loc.localPos][m << 1] = to;
                    ++nalleles;
                    ++ngenos;
                }
                if (from == loc.localAM->alleles[loc.localPos][(m << 1) + 1]) {
                    loc.localAM->alleles[loc.localPos][(m << 1) + 1] = to;
                    ++nalleles;
                    if (loc.localAM->alleles[loc.localPos][m << 1] != loc.localAM->alleles[loc.localPos][(m << 1) + 1]) {
                        ++ngenos;
                    }
                }
            }

            loc = gsc_next_forwards(&it);
        }

        gsc_delete_bidirectional_iter(&it);


    } else if (gsc_get_index_of_genetic_marker(which_marker, d->genome, &markeri)) {
        nmarkers = 1;
        gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,NO_GROUP);
        gsc_GenoLocation loc = gsc_set_bidirectional_iter_to_start(&it);
        while (IS_VALID_LOCATION(loc)) {
            if (from == loc.localAM->alleles[loc.localPos][markeri << 1]) {
                loc.localAM->alleles[loc.localPos][markeri << 1] = to;
                ++nalleles;
                ++ngenos;
            }
            if (from == loc.localAM->alleles[loc.localPos][(markeri << 1) + 1]) {
                loc.localAM->alleles[loc.localPos][(markeri << 1) + 1] = to;
                ++nalleles;
                if (loc.localAM->alleles[loc.localPos][markeri << 1] != loc.localAM->alleles[loc.localPos][(markeri << 1) + 1]) {
                    ++ngenos;
                }
            }

            loc = gsc_next_forwards(&it);
        }

        gsc_delete_bidirectional_iter(&it);

    } else {
        nmarkers = 0;
        ngenos = 0;
    }

    Rprintf("Changed allele %c to %c %u times across %u markers and %u genotypes\n", from, to, nalleles, nmarkers, ngenos);
}

/** Count and return the number of digits in `i`.
 *
 * @param i the integer whose digits are to be counted.
 * @returns the number of digits to print `i`
 */
int gsc_get_integer_digits(const int i) {
    int digits = 0, ii = i;
    while (ii != 0) {
        ii = ii / 10;
		digits ++;
	}
	return digits;
}

/** Comparator function for qsort. Used to compare a pair of doubles* to sort
 * them in descending order of the doubles they point to.
 * @see gsc_split_by_bv()
 *
 * Sorts higher numbers before lower numbers. If they are equal, their
 * order after comparison is undefined.
 */
static int gsc_helper_descending_pdouble_comparer(const void* pp0, const void* pp1) {
    double d0 = **(double **)pp0;
    double d1 = **(double **)pp1;
	if (d0 > d1) {
		return -1;
	} else {
		return (d0 < d1); // 0 if equal, 1 if d0 is smaller
	}
}

/** Comparator function for qsort. Used to compare a pair of doubles to sort in ascending order
 * @see gsc_split_by_bv()
 *
 * Sorts lower numbers before higher numbers. If they are equal, their
 * order after comparison is undefined.
 */
static int gsc_helper_ascending_double_comparer(const void* pp0, const void* pp1) {
    double d0 = *(double *)pp0;
    double d1 = *(double *)pp1;
	if (d0 < d1) {
		return -1;
	} else {
		return (d0 > d1); // 0 if equal, 1 if d0 is smaller
	}
}

/** Comparator function for qsort. Used to compare a pair of doubles* to sort
 * them in ascending order of the doubles they point to.
 * @see gsc_split_by_bv()
 *
 * Sorts lower numbers before higher numbers. If they are equal, their
 * order after comparison is undefined.
 */
static int gsc_helper_ascending_pdouble_comparer(const void* pp0, const void* pp1) {
    double d0 = **(double **)pp0;
    double d1 = **(double **)pp1;
    if (d0 < d1) {
        return -1;
    } else {
        return (d0 > d1); // 0 if equal, 1 if d0 is smaller
    }
}

/** Comparator function for qsort. Used to compare a pair of pointers to strings
 * to sort them in alphabetical order of the strings.
 */
static int gsc_helper_indirect_alphabetical_str_comparer(const void* p0, const void* p1) {
    char* str1 = **(char***)p0;
    char* str2 = **(char***)p1;
    return strcmp(str1,str2);
}

/** Comparator function for qsort. Used to compare a pair of gsc_MapFileUnits
 * to sort them in ascending order of chromosome number.
 */
static int gsc_helper_mapfileunit_ascending_chr_comparer(const void* p0, const void* p1) {
    struct gsc_MapfileUnit s0 = *(struct gsc_MapfileUnit*)p0;
    struct gsc_MapfileUnit s1 = *(struct gsc_MapfileUnit*)p1;
    //return s0.ul - s1.ul;
    return (s0.chr < s1.chr) ? -1 : (s0.chr > s1.chr);
}

/** Comparator function for qsort. Used to compare a pair of gsc_MapFileUnits
 * to sort them in ascending order of position.
 */
static int gsc_helper_mapfileunit_ascending_d_comparer(const void* p0, const void* p1) {
    struct gsc_MapfileUnit s0 = *(struct gsc_MapfileUnit*)p0;
    struct gsc_MapfileUnit s1 = *(struct gsc_MapfileUnit*)p1;
    return (s0.pos < s1.pos) ? -1 : (s0.pos > s1.pos);
}

/** Move all details of the genotype at one gsc_GenoLocation to another gsc_GenoLocation.
 *
 * After the genotype at @a gsc_GenoLocation @a from has been copied to @a gsc_GenoLocation @a to,
 * the information is cleared from @a gsc_GenoLocation @a from.
 *
 * Copying a genotype that belongs to one @a gsc_SimData object to a @a gsc_GenoLocation in a
 * different @a gsc_SimData object is not recommended, as genetic maps and numbers of custom labels
 * may differ. This function is not suitable for that purpose.
 *
 * @param from copying source. The information at this GenoLocation will be cleared.
 * @param to copying destination. A prior genotype existing at this location will be overwritten with
 * a warning.
 * @param label_defaults array of default values for custom labels. It must have at least enough values
 * that @a to.localAM.n_labels values can be read from it.
 */
void gsc_move_genotype(gsc_GenoLocation from, gsc_GenoLocation to, int* label_defaults) {
    if (to.localAM == from.localAM && to.localPos == from.localPos) {
        return;
    }
    if (to.localAM->groups[to.localPos].num != GSC_NO_GROUP.num) {
        warning("In moving a genotype from %p:%d to %p:%d, the genotype at %p:%d will be overwritten\n", from.localAM, from.localPos, to.localAM, to.localPos, to.localAM, to.localPos);
        --to.localAM->n_genotypes;
    }
    to.localAM->alleles[to.localPos] = from.localAM->alleles[from.localPos];
    from.localAM->alleles[from.localPos] = NULL;

    to.localAM->names[to.localPos] = from.localAM->names[from.localPos];
    from.localAM->names[from.localPos] = NULL;

    to.localAM->ids[to.localPos] = from.localAM->ids[from.localPos];
    from.localAM->ids[from.localPos] = GSC_NO_PEDIGREE;

    to.localAM->pedigrees[0][to.localPos] = from.localAM->pedigrees[0][from.localPos];
    from.localAM->pedigrees[0][from.localPos] = GSC_NO_PEDIGREE;
    to.localAM->pedigrees[1][to.localPos] = from.localAM->pedigrees[1][from.localPos];
    from.localAM->pedigrees[1][from.localPos] = GSC_NO_PEDIGREE;

    to.localAM->groups[to.localPos] = from.localAM->groups[from.localPos];
    from.localAM->groups[from.localPos] = GSC_NO_GROUP;

    if (to.localAM->n_labels != from.localAM->n_labels) {
        warning("Origin and destination when copying genotype do not have the same number of custom labels (n_labels). The genotype now at %p:%d will have lost its label data.\n", to.localAM, to.localPos);
    } else {
        for (int i = 0; i < to.localAM->n_labels; ++i) {
            to.localAM->labels[i][to.localPos] = from.localAM->labels[i][from.localPos];
            from.localAM->labels[i][from.localPos] = label_defaults[i];
        }
    }

    if (from.localAM != to.localAM) {
        --from.localAM->n_genotypes;
        ++to.localAM->n_genotypes;
    }
}

/** Sets the current cursor position in a gsc_GappyIterator to the next valid position,
 *  if the cursor is not already a valid position.
 *
 * Valid positions in the linked list of AlleleMatrices may contain genotypes or not.
*/
static gsc_GenoLocation gsc_nextgappy_valid_pos(struct gsc_GappyIterator* it) {
    if (it->cursor.localAM == NULL) {
        it->cursor = GSC_INVALID_GENO_LOCATION;
    } else if (it->cursor.localPos >= CONTIG_WIDTH) {
        it->cursor.localPos = 0;
        it->cursor.localAM = it->cursor.localAM->next;
        ++it->cursorAMIndex;
        if (it->cursor.localAM == NULL) {
            it->cursor = GSC_INVALID_GENO_LOCATION;
        }
    }
    return it->cursor;
}

/** Sets the current cursor position in a gsc_GappyIterator to the next empty
 *  position, if the cursor is not already an empty position.
 *
 * Empty positions do not contain genotypes. This function does not change
 * the cursor position if the curstor is already on a gap.
*/
static gsc_GenoLocation gsc_nextgappy_get_gap(struct gsc_GappyIterator* it) {
    if (!GSC_IS_VALID_LOCATION(gsc_nextgappy_valid_pos(it))) {
        return GSC_INVALID_GENO_LOCATION;
    }

    while (it->cursor.localAM->groups[it->cursor.localPos].num != GSC_NO_GROUP.num) {

        // Trusts that n_genotypes is correct.
        if (it->cursor.localAM->n_genotypes == CONTIG_WIDTH) { // work-saver: skip this gsc_AlleleMatrix if it is already known to be full.
            it->cursor.localAM = it->cursor.localAM->next;
        } else {
            ++it->cursor.localPos;
        }

        if (!GSC_IS_VALID_LOCATION(gsc_nextgappy_valid_pos(it))) {
            return GSC_INVALID_GENO_LOCATION;
        }
    }

    return it->cursor;
}

/** Sets the current cursor position in a gsc_GappyIterator to the next filled
 *  position, if the cursor is not already a filled position.
 *
 * Non-gap positions do contain genotypes. This function does not change the
 *  cursor position if the curstor is already on a non-gap.
*/
static gsc_GenoLocation gsc_nextgappy_get_nongap(struct gsc_GappyIterator* it) {
    if (!GSC_IS_VALID_LOCATION(gsc_nextgappy_valid_pos(it))) {
        return GSC_INVALID_GENO_LOCATION;
    }

    while (it->cursor.localAM->groups[it->cursor.localPos].num == GSC_NO_GROUP.num) {
        ++it->cursor.localPos;
        if (!GSC_IS_VALID_LOCATION(gsc_nextgappy_valid_pos(it))) {
            return GSC_INVALID_GENO_LOCATION;
        }
    }

    return it->cursor;
}


/** A function to tidy the internal storage of genotypes after addition
 * or deletion of genotypes in the gsc_SimData. Not intended to be called by an
 * end user - functions which require it should be calling it already.
 *
 * Ideally, we want all gsc_AlleleMatrix structs in the gsc_SimData's linked list
 * to have no gaps. That is, if there are more than CONTIG_WIDTH genotypes, the
 * all gsc_AlleleMatrix structs except the last should be full (contain CONTIG_WIDTH genotypes),
 * and the last should have the remaining n genotypes at local indexes 0 to n-1
 * (so all at the start of the gsc_AlleleMatrix, with no gaps between them).
 *
 * We trust that each gsc_AlleleMatrix's n_genotypes is correct.
 *
 * This function achieves the cleanup by using two pointers: a checker out
 * the front that identifies a genotype that needs to be shifted back/that
 * occurs after a gap, and a filler that identifies each gap and copies
 * the genotype at the checker back into it.
 *
 * @param d The gsc_SimData struct on which to operate.
 */
void gsc_condense_allele_matrix( gsc_SimData* d) {
    // Find the first gap
    struct gsc_GappyIterator filler = {.cursor=(gsc_GenoLocation){.localAM=d->m, .localPos=0}, .cursorAMIndex=0};
    gsc_nextgappy_get_gap(&filler);

    if (!GSC_IS_VALID_LOCATION(filler.cursor)) {
        return; // no gaps found
    }

    struct gsc_GappyIterator checker = filler; // copy filler
    ++checker.cursor.localPos;
    gsc_nextgappy_get_nongap(&checker);

    while (GSC_IS_VALID_LOCATION(filler.cursor) && GSC_IS_VALID_LOCATION(checker.cursor)) {
        gsc_move_genotype(checker.cursor, filler.cursor, d->label_defaults);

        ++filler.cursor.localPos;
        gsc_nextgappy_get_gap(&filler);

        ++checker.cursor.localPos;
        gsc_nextgappy_get_nongap(&checker);
    }

    if (filler.cursor.localAM->next != NULL && filler.cursor.localAM->next->n_genotypes == 0) {
        gsc_delete_allele_matrix(filler.cursor.localAM->next);
        filler.cursor.localAM->next = NULL;
    }
}



/*----------------------------------Locators---------------------------------*/


/** Create a bidirectional iterator
 *
 *  A bidirectional iterator can be used to loop through members of a particular group
 *  or all genotypes in the simulation, forwards or backwards.
 *  The iterator is not initialised to any location at this point. The first call to
 *  a next* function will initialise it, or you can manually initialise using
 *  gsc_set_bidirectional_iter_to_start or gsc_set_bidirectional_iter_to_end. gsc_next_forwards will
 *  initialise it to the first in the group, or call gsc_next_backwards to initialise it to the
 *  last in the group.
 *
 * @shortnamed{create_bidirectional_iter}
 *
 *  @warning An initialised iterator is only valid if no genotypes have been added
 *  to the group and no genotypes in the simulation as a whole have been removed
 *  since its creation. Discard the old iterator and create a new one when the state
 *  of the simulation changes.
 *
 *  @param d pointer to the gsc_SimData containing the genotypes to iterate through
 *  @param group the group number of the group to iterate through, or 0 to iterate
 *  through all genotypes.
 * @returns uninitialised gsc_BidirectionalIterator for the provided group.
 */
gsc_BidirectionalIterator gsc_create_bidirectional_iter( gsc_SimData* d, const gsc_GroupNum group) {
    return (gsc_BidirectionalIterator) {
        .d = d,
        .group = group,
        .localPos = GSC_UNINIT,

        .cachedAM = d->m,
        .cachedAMIndex = 0,

        .atStart = GSC_FALSE,
        .atEnd = GSC_FALSE
    };
}

/** Create a Random Access Iterator
 *
 *  A random access iterator and the function gsc_next_get_nth() can be used to
 *  get any genotype in the simulation by index or any genotype in a group
 *  by within-group index.
 *
 *  The random access iterator gives you access to the (i)th group member in
 *  O(n) time (n being the number of genotypes in the simulation), but contains
 *  a persistent cache that allows subsequent O(1) access to group members.
 *
 *  The random access iterator created by this function is initialised. If the
 *  `groupSize` attribute of the returned iterator is set to 0, no members of
 *  that group could be found in the simulation, so the iterator will never
 *  return a valid position.
 *
 * @shortnamed{create_randomaccess_iter}
 *
 *  @warning An initialised iterator is only valid if no genotypes have been added
 *  to the group and no genotypes in the simulation as a whole have been removed
 *  since its creation. Discard the old iterator and create a new one when the state
 *  of the simulation changes.
 *
 *  @param d pointer to the gsc_SimData containing the genotypes to iterate through
 *  @param group the group number of the group to iterate through, or 0 to iterate
 *  through all genotypes.
 * @returns initialised Random Access iterator for the provided group.
*/
gsc_RandomAccessIterator gsc_create_randomaccess_iter( gsc_SimData* d, const gsc_GroupNum group) {
    unsigned long first = 0;
    gsc_AlleleMatrix* firstAM = d->m;
    char anyExist = GSC_TRUE;

    // Want to know:
    // - is this group empty? (randomAccess should know if group size is 0)
    // - what is the first genotype index in this group?

    if (group.num == GSC_NO_GROUP.num) { // scanning all genotypes
        while (firstAM->n_genotypes == 0) {
            if (firstAM->next == NULL) {
                // gsc_SimData is empty. Nowhere to go.
                anyExist = GSC_FALSE;

            } else { // Keep moving forwards through the list. Not polite enough to clean up the blank AM.

                firstAM = firstAM->next;
            }
        }

    } else { // scanning a specific group

        char exitNow = GSC_FALSE;
        while (exitNow == GSC_FALSE) {

            // Set first, firstAM, firstAMIndex if appropriate
            for (int i = 0; i < firstAM->n_genotypes; ++i) {
                if (firstAM->groups[i].num == group.num) {
                    first = i;
                    exitNow = GSC_TRUE;
                    break;
                }
            }

            // Move along and set anyExist if appropriate
            if (exitNow == GSC_FALSE) {
                firstAM = firstAM->next;
                if (firstAM == NULL) {
                    anyExist = GSC_FALSE;
                    exitNow = GSC_TRUE;
                }
            }
        }

    }

    gsc_GenoLocation* cache = NULL;
    unsigned int cacheSize = 0;
    if (anyExist) {
        cacheSize = 50;
        cache = gsc_malloc_wrap((sizeof(gsc_GenoLocation)*cacheSize),GSC_TRUE);
        cache[0] = (gsc_GenoLocation) {
                .localAM= firstAM,
                .localPos = first,
        };
        for (int i = 1; i < cacheSize; ++i) {
            cache[i] = GSC_INVALID_GENO_LOCATION;
        }

    }

    return (gsc_RandomAccessIterator) {
        .d = d,
        .group = group,

        .largestCached = anyExist ? 0 : GSC_UNINIT,
        .groupSize = anyExist ? GSC_UNINIT : 0,
        .cacheSize = cacheSize,
        .cache = cache
    };
}

/** Get an gsc_AlleleMatrix by index in the linked list
 *
 *  listStart is considered the 0th gsc_AlleleMatrix (0-based indexing)
 *
 *  @param listStart the 0th gsc_AlleleMatrix in the linked list
 *  @param n the index of the desired gsc_AlleleMatrix in the linked list
 *  @returns pointer to the nth gsc_AlleleMatrix, if it exists in the list, or
 *  NULL if the list is shorter than n
 */
gsc_AlleleMatrix* gsc_get_nth_AlleleMatrix(gsc_AlleleMatrix* listStart, const unsigned int n) {
    unsigned int currentIndex = 0;
    gsc_AlleleMatrix* am = listStart;
    while (currentIndex < n) {
        if (am->next == NULL) {
            return NULL;
        } else {
            am = am->next;
            currentIndex++;
        }
    }
    return am;
}

/** Initialise a Bidirectional iterator to the start of its sequence.
 *
 *  Can be used to reset a gsc_BidirectionalIterator so that it is pointing
 *  at the very first member of the group it is looping through.
 *
 * @shortnamed{set_bidirectional_iter_to_start}
 *
 * @param it BidirectioanlIterator to initialise
 * @return location of the first member of the sequence.
 */
gsc_GenoLocation gsc_set_bidirectional_iter_to_start(gsc_BidirectionalIterator* it) {
    unsigned int first = 0;
    gsc_AlleleMatrix* firstAM = it->d->m;
    unsigned int firstAMIndex = 0;
    char anyExist = GSC_TRUE;

    // Want to know:
    // - is this group empty? (iterator should know if it is at the end as well as at the start)
    // - what is the first genotype index in this group?

    if (it->group.num == GSC_NO_GROUP.num) {
        while (firstAM->n_genotypes == 0) {
            if (firstAM->next == NULL) {
                anyExist = GSC_FALSE; // gsc_SimData is empty.

            } else { // (Not polite enough to clean up the blank AM.)
                firstAM = firstAM->next;
                firstAMIndex++;
                // first += 0;
            }
        }

        // After this runs we have set firstAM, first, firstAMIndex, anyExist appropriately

    } else { // scanning a specific group

        char exitNow = GSC_FALSE;
        while (exitNow == GSC_FALSE) {

            // Set first, firstAM, firstAMIndex if appropriate
            for (int i = 0; i < firstAM->n_genotypes; ++i) {
                if (firstAM->groups[i].num == it->group.num) {
                    first = i;
                    exitNow = GSC_TRUE;
                    break;
                }
            }

            // Move along and set anyExist if appropriate
            if (exitNow == GSC_FALSE) {
                firstAM = firstAM->next;
                firstAMIndex++;
                if (firstAM == NULL) {
                    first = GSC_UNINIT;
                    anyExist = GSC_FALSE;
                    exitNow = GSC_TRUE;
                }
            }
        }
    }

    it->localPos = first;
    it->atStart = GSC_TRUE;
    it->atEnd = !anyExist;
    it->cachedAM = firstAM;
    it->cachedAMIndex = firstAMIndex;

    return (gsc_GenoLocation) {
        .localAM = firstAM,
        .localPos = first
    };
}

/** Initialise a Bidirectional iterator to the end of its sequence.
 *
 *  Can be used to reset a gsc_BidirectionalIterator so that it is pointing
 *  at the very last member of the group it is looping through.
 *
 * @shortnamed{set_bidirectional_iter_to_end}
 *
 * @param it BidirectioanlIterator to initialise
 * @return location of the last member of the sequence.
 */
gsc_GenoLocation gsc_set_bidirectional_iter_to_end(gsc_BidirectionalIterator* it) {
    unsigned long last = 0;
    gsc_AlleleMatrix* lastAM = it->d->m;
    unsigned int lastAMIndex = 0;
    char anyExist = GSC_TRUE;

    // Want to know:
    // - is this group empty? (iterator should know if it is at the end as well as at the start)
    // - what is the first genotype index in this group?

    if (it->group.num == GSC_NO_GROUP.num) {
        while (lastAM->next != NULL && lastAM->next->n_genotypes != 0) {
            lastAM = lastAM->next;
            lastAMIndex++;
        }
        if (lastAMIndex > 0 || lastAM->n_genotypes > 0) {
            last = lastAM->n_genotypes - 1;
        } else {
            anyExist = GSC_FALSE;
        }

    } else { // scanning a specific group

        // Find last AM
        while (lastAM->next != NULL && lastAM->next->n_genotypes != 0) {
            lastAM = lastAM->next;
            lastAMIndex++;
        }

        char exitNow = GSC_FALSE;
        while (exitNow == GSC_FALSE) {

            // Set first, firstAM, firstAMIndex if appropriate
            for (int i = lastAM->n_genotypes - 1; i >= 0; --i) {
                if (lastAM->groups[i].num == it->group.num) {
                    last = i;
                    exitNow = GSC_TRUE;
                    break;
                }
            }

            // Move along and set anyExist if appropriate
            if (exitNow == GSC_FALSE) {
                --lastAMIndex;
                lastAM = gsc_get_nth_AlleleMatrix(it->d->m, lastAMIndex);
                if (lastAM->n_genotypes == 0) {
                    last = GSC_UNINIT;
                    anyExist = GSC_FALSE;
                    exitNow = GSC_TRUE;
                }
            }
        }
    }

    it->localPos = last;
    it->atStart = !anyExist;
    it->atEnd = GSC_TRUE;
    it->cachedAM = lastAM;
    it->cachedAMIndex = lastAMIndex;

    return (gsc_GenoLocation) {
        .localAM = lastAM,
        .localPos = last
    };
}


/** Get the next location from a bidirectional iterator
 *
 * Moves the pointer of a gsc_BidirectionalIterator forwards by
 * one step, and returns the new location it points to. If
 * the gsc_BidirectionalIterator is not initialised, then initialises
 * it to the very first element.
 *
 * Returns GSC_INVALID_GENO_LOCATION
 * if the gsc_BidirectionalIterator is corrupted or if it is at the
 * end of the sequence. Test the return value of this function
 * with @ref GSC_IS_VALID_LOCATION().
 *
 * @shortnamed{next_forwards}
 *
 * @param it the gsc_BidirectionalIterator to iterate forwards
 * @returns the location of the next genotype in the sequence,
 * or GSC_INVALID_GENO_LOCATION if the iterator is corrupted or
 * the iterator's pointer is already at the last element.
 */
gsc_GenoLocation gsc_next_forwards(gsc_BidirectionalIterator* it) {
    if (it->localPos == GSC_UNINIT) {
        return gsc_set_bidirectional_iter_to_start(it);
    }

    if (it->atEnd) { // || validate_bidirectional_cache(it) == GSC_FALSE) { // can't use this because what if our iterator user is modifying group allocations?
        return GSC_INVALID_GENO_LOCATION;
    }

    if (it->group.num == GSC_NO_GROUP.num) {

        // Search for the next value.
        if (it->localPos + 1 < it->cachedAM->n_genotypes) {
            // The next value is in the same gsc_AlleleMatrix
            it->localPos++;
            it->atStart = GSC_FALSE;
            return (gsc_GenoLocation) {
                .localAM = it->cachedAM,
                .localPos = it->localPos
            };

        } else {
            // The next value is in the next gsc_AlleleMatrix
            gsc_AlleleMatrix* nextAM = it->cachedAM;
            int nextAMIndex = it->cachedAMIndex;
            do {
                nextAM = nextAM->next;
                nextAMIndex++;
            } while (nextAM != NULL && nextAM->n_genotypes == 0);

            if (nextAM == NULL) {
                // There is no further gsc_AlleleMatrix; we are at the end of the iterator.
                it->atEnd = GSC_TRUE;
                return GSC_INVALID_GENO_LOCATION;
            } else {
                it->cachedAM = nextAM;
                it->cachedAMIndex = nextAMIndex;
                it->localPos = 0;
                it->atStart = GSC_FALSE;
                return (gsc_GenoLocation) {
                    .localAM = it->cachedAM,
                    .localPos = 0
                };
            }
        }

    } else { // We are iterating through a specific group

        // Search for the next value
        while(1) {
            if (it->localPos + 1 < it->cachedAM->n_genotypes) {
                for (++it->localPos; it->localPos < it->cachedAM->n_genotypes; ++it->localPos) {
                    if (it->cachedAM->groups[it->localPos].num == it->group.num) {
                        it->atStart = GSC_FALSE;
                        return (gsc_GenoLocation) {
                            .localAM = it->cachedAM,
                            .localPos = it->localPos
                        };
                    }
                }
            }

            gsc_AlleleMatrix* nextAM = it->cachedAM;
            int nextAMIndex = it->cachedAMIndex;
            do {
                nextAM = nextAM->next;
                nextAMIndex++;
            } while (nextAM != NULL && nextAM->n_genotypes == 0);

            if (nextAM == NULL) {
                // There is no further gsc_AlleleMatrix; we are at the end of the iterator.
                it->atEnd = GSC_TRUE;
                return GSC_INVALID_GENO_LOCATION;
            } else {
                it->cachedAM = nextAM;
                it->cachedAMIndex = nextAMIndex;
                it->localPos = 0;
                if (it->cachedAM->groups[it->localPos].num == it->group.num) {
                    it->atStart = GSC_FALSE;
                    return (gsc_GenoLocation) {
                        .localAM = it->cachedAM,
                        .localPos = it->localPos
                    };
                }
            }
        }

    }
}


/** Get the previous location from a bidirectional iterator
 *
 * Moves the pointer of a gsc_BidirectionalIterator backwards by
 * one step, and returns the new location it points to. If
 * the gsc_BidirectionalIterator is not initialised, then initialises
 * it to the very last element.
 *
 * Slightly slower than gsc_next_forwards, because the gsc_AlleleMatrix linked list
 * is not bidirectional. To find the preceding gsc_AlleleMatrix, it needs to
 * count forwards from the beginning of the list to find the n-1th gsc_AlleleMatrix.
 *
 * Returns GSC_INVALID_GENO_LOCATION
 * if the gsc_BidirectionalIterator is corrupted or if it is at the
 * beginning of the sequence. Test the return value of this function
 * with @ref GSC_IS_VALID_LOCATION().
 *
 * @shortnamed{next_backwards}
 *
 * @param it the gsc_BidirectionalIterator to iterate backwards
 * @returns the location of the previous genotype in the sequence,
 * or GSC_INVALID_GENO_LOCATION if the iterator is corrupted or
 * the iterator's pointer is already at the first element.
 */
gsc_GenoLocation gsc_next_backwards(gsc_BidirectionalIterator* it) {
    if (it->localPos == GSC_UNINIT) {
        return gsc_set_bidirectional_iter_to_end(it);
    }

    if (it->atStart) { //|| validate_bidirectional_cache(it) == GSC_FALSE) {
        return GSC_INVALID_GENO_LOCATION;
    }

    if (it->group.num == GSC_NO_GROUP.num) {

        // Search for the previous value.
        if (it->localPos > 0) {
            // The previous value is in the same gsc_AlleleMatrix
            it->localPos--;
            it->atEnd = GSC_FALSE;
            return (gsc_GenoLocation) {
                .localAM = it->cachedAM,
                .localPos = it->localPos
            };

        } else {
            // The previous value is in the previous gsc_AlleleMatrix
            if (it->cachedAMIndex == 0) {
                it->atStart = GSC_TRUE;
                return GSC_INVALID_GENO_LOCATION;
            } else {
                gsc_AlleleMatrix* nextAM = it->cachedAM;
                int nextAMIndex = it->cachedAMIndex;
                do {
                    nextAMIndex--;
                    nextAM = gsc_get_nth_AlleleMatrix(it->d->m, nextAMIndex);
                } while (nextAM != NULL && nextAM->n_genotypes == 0);

                if (nextAM == NULL) {
                    it->atStart = GSC_TRUE;
                    return GSC_INVALID_GENO_LOCATION;
                } else {
                    it->cachedAM = nextAM;
                    it->cachedAMIndex = nextAMIndex;
                    it->localPos = it->cachedAM->n_genotypes - 1;
                    it->atEnd = GSC_FALSE;
                    return (gsc_GenoLocation) {
                        .localAM = it->cachedAM,
                        .localPos = it->localPos
                    };
                }
            }
        }

    } else { // We are iterating through a specific group

        // Search for the next value
        while(1) {
            if (it->localPos > 0) {
                for (--it->localPos; it->localPos >= 0; --it->localPos) {
                    if (it->cachedAM->groups[it->localPos].num == it->group.num) {
                        it->atEnd = GSC_FALSE;
                        return (gsc_GenoLocation) {
                            .localAM = it->cachedAM,
                            .localPos = it->localPos
                        };
                    }
                }
            }

            if (it->cachedAMIndex == 0) {
                it->atStart = GSC_TRUE;
                it->localPos = 0;
                return GSC_INVALID_GENO_LOCATION;
            } else {
                gsc_AlleleMatrix* nextAM = it->cachedAM;
                int nextAMIndex = it->cachedAMIndex;
                do {
                    nextAMIndex--;
                    nextAM = gsc_get_nth_AlleleMatrix(it->d->m, nextAMIndex);
                } while (nextAM != NULL && nextAM->n_genotypes == 0);

                if (nextAM == NULL) {
                    it->atStart = GSC_TRUE;
                    return GSC_INVALID_GENO_LOCATION;
                } else {
                    it->cachedAM = nextAM;
                    it->cachedAMIndex = nextAMIndex;
                    it->localPos = it->cachedAM->n_genotypes - 1;
                    if (it->cachedAM->groups[it->localPos].num == it->group.num) {
                        it->atEnd = GSC_FALSE;
                        return (gsc_GenoLocation) {
                            .localAM = it->cachedAM,
                            .localPos = it->localPos
                        };
                    }
                }
            }
        }
    }
}


/** Get a location by index using a gsc_RandomAccessIterator
 *
 * Gives the location of the provided global index (if `it->group == 0`)
 * or the location of the provided group index (if `it->group` is not 0),
 * by first searching the gsc_RandomAccessIterator's cache for it and if
 * not, searching the gsc_SimData for it and adding it and its predecessors
 * to the cache.
 *
 * Returns GSC_INVALID_GENO_LOCATION
 * if the iterator is corrupted or the index is invalid. Check the
 * return value with @ref GSC_IS_VALID_LOCATION().
 *
 * @shortnamed{next_get_nth}
 *
 * @param it the gsc_RandomAccessIterator to read and update the cache of.
 * @param n If the iterator is iterating through all genotypes, the global 
 * index in the simulation of the genotype you want to access. If the iterator
 * is iterating through the group, the within-group index of the genotype 
 * you want to access. In either case the first genotype is at index 0.
 * @returns the location of the nth genotype/nth group member, or
 * GSC_INVALID_GENO_LOCATION if the index is invalid.
 */
gsc_GenoLocation gsc_next_get_nth(gsc_RandomAccessIterator* it, const unsigned int n) {
    // Step 0: First check n is in our group index range as far as we know it
    // (If it->groupSize is negative then we don't know the range)
    if (it->groupSize > 0 && it->groupSize <= n) {
        return GSC_INVALID_GENO_LOCATION;
    }

    // Step 1: Check if we have it in the cache.
    if (n < it->cacheSize) {
        // 'n' is less than or equal to our current furthest cached group member.

        if (GSC_IS_VALID_LOCATION(it->cache[n])) {
            return it->cache[n];
        }

        // Otherwise we do not have it cached, but we will enter it into the cache in the next section

    } else {
        // We need to expand the cache to fit it.
        unsigned int newCacheSize = it->cacheSize;
        if (it->cacheSize == 0) {
            newCacheSize = 50;
        }
        while (newCacheSize < n+1) {
            newCacheSize = newCacheSize << 1;
        }
        gsc_GenoLocation* newCache = gsc_malloc_wrap(sizeof(gsc_GenoLocation)*newCacheSize,GSC_TRUE);
        int i = 0;
        for (; i < it->cacheSize; ++i) {
            newCache[i] = it->cache[i];
        }
        for (; i < newCacheSize; ++i) {
            newCache[i] = GSC_INVALID_GENO_LOCATION;
        }
        GSC_FREE(it->cache);
        it->cacheSize = newCacheSize;
        it->cache = newCache;
    }

    // Validity checks for a random access iterator: largestCached must exist,
    // is indeed cached and belongs to the same group
    if (it->largestCached < 0 || (!GSC_IS_VALID_LOCATION(it->cache[it->largestCached]) &&
            (it->group.num == GSC_NO_GROUP.num || it->group.num != gsc_get_group(it->cache[it->largestCached]).num))) {
        return GSC_INVALID_GENO_LOCATION;
    }

    // Step 2: Actually finding the nth group member.
    if (it->group.num == GSC_NO_GROUP.num) {
        // Assuming all non-end gsc_AlleleMatrix are filled to CONTIG_WIDTH
        gsc_GenoLocation expectedLocation = {
            .localAM = gsc_get_nth_AlleleMatrix(it->d->m, n / CONTIG_WIDTH),
            .localPos = n % CONTIG_WIDTH
        };
        // Check n was not too large
        if (expectedLocation.localAM == NULL ||
                expectedLocation.localAM->n_genotypes <= expectedLocation.localPos) {
            return GSC_INVALID_GENO_LOCATION;
        }
        return expectedLocation;

    } else { // searching for a particular group

        gsc_AlleleMatrix* currentAM;
        int groupN;
        int localPos;
        if (it->largestCached < n) {
            // Search forwards from largestCached
            currentAM = it->cache[it->largestCached].localAM;
            groupN = it->largestCached;
            localPos = it->cache[it->largestCached].localPos + 1;
            while (1) {
                for (; localPos < currentAM->n_genotypes; ++localPos) {
                    // If we found a group member, cache it and count upwards towards n
                    if (currentAM->groups[localPos].num == it->group.num) {
                        it->largestCached = ++groupN;
                        it->cache[groupN] = (gsc_GenoLocation) {
                            .localAM = currentAM,
                            .localPos = localPos
                        };
                        if (groupN == n) {
                            return it->cache[n];
                        }
                    }
                }

                if (currentAM->next == NULL || currentAM->next->n_genotypes == 0) {
                    // We are at the end of the iterator and have not found n
                    it->groupSize = groupN + 1;
                    return GSC_INVALID_GENO_LOCATION;
                } else {
                    currentAM = currentAM->next;
                    localPos = 0;
                }

           }

        } else {
            // With current method of filling cache, this branch will never be taken

            // Search backwards from largestCached for one gsc_AlleleMatrix, to take advantage of that pointer
            if (it->cache[it->largestCached].localAM != it->d->m) {
                currentAM = it->cache[it->largestCached].localAM;
                groupN = it->largestCached;
                localPos = it->cache[it->largestCached].localPos;
                for (; localPos >= 0; --localPos) {
                    // If we found a group member, cache it and count down towards n
                    if (currentAM->groups[localPos].num == it->group.num) {
                        --groupN;
                        it->cache[groupN] = (gsc_GenoLocation) {
                            .localAM = currentAM,
                            .localPos = localPos
                        };
                        if (groupN == n) {
                            return it->cache[n];
                        }
                    }
                }
            }

            // Then search forwards from start
            currentAM = it->d->m;
            localPos = 0;
            groupN = -1; // overwrite everything
            while (currentAM != it->cache[it->largestCached].localAM) {
                for (; localPos < currentAM->n_genotypes; ++localPos) {
                    // If we found a group member, cache it and count upwards towards n
                    if (currentAM->groups[localPos].num == it->group.num) {
                        ++groupN;
                        it->cache[groupN] = (gsc_GenoLocation) {
                            .localAM = currentAM,
                            .localPos = localPos
                        };
                        if (groupN == n) {
                            return it->cache[n];
                        }
                    }
                }

                if (currentAM->next == NULL || currentAM->next->n_genotypes == 0) {
                    // We are at the end of the iterator and have not found n. Also we didn't reach
                    // 'largestCached' yet. Something is wrong with the iterator
                    it->groupSize = groupN + 1;
                    return GSC_INVALID_GENO_LOCATION;
                } else {
                    currentAM = currentAM->next;
                    localPos = 0;
                }
            }

            // We were somehow unable to find the nth group member. Something is probably wrong
            // with the iterator
            return GSC_INVALID_GENO_LOCATION;

        }
    }
}

/** Returns the name of the genotype with a given id.
 *
 * The function uses a bisection search on the gsc_AlleleMatrix where it should be
 * found. Searching by id is therefore fast.
 * @see gsc_get_from_ordered_pedigree_list()
 *
 * This function assumes that ids are never reshuffled in the gsc_SimData. This is true
 * as long as gsc_condense_allele_matrix's process of moving genotypes while retaining
 * their order is the only genotype-rearranging function in use.
 * This function's algorithm will need to
 * be revisited if different genotype-rearranging processes are implemented.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose name is sought
 * @returns the name of the genotype that has id `id`, as a copy of the pointer
 * to the heap memory where the name is saved (so *don't* free the pointer returned
 * from this function). Returns NULL if the ID does not exist. Note the return value
 * might also be NULL if the genotype of this ID has no name.
 */
char* gsc_get_name_of_id( const gsc_AlleleMatrix* start, const gsc_PedigreeID id) {
    if (id.id == GSC_NO_PEDIGREE.id) {
        warning( "Invalid ID %d\n", id.id);
		return NULL;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
    const gsc_AlleleMatrix* m = start;

	while (1) {
        // try to find our id. Does this AM potentially have the right range for it?
        // If we're not sure, because either of the endpoints does not have its ID tracked,
        // check anyway
        if (m->n_genotypes != 0 && (id.id >= m->ids[0].id || m->ids[0].id == GSC_NO_PEDIGREE.id) &&
                (id.id <= m->ids[m->n_genotypes - 1].id || m->ids[m->n_genotypes - 1].id == GSC_NO_PEDIGREE.id)) {

            int index = gsc_get_from_ordered_pedigree_list(id, m->n_genotypes, m->ids);

            if (index < 0) {
                // search failed
                if (m->next == NULL) {
                    warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id.id);
                    return NULL;
                } else {
                    m = m->next;
                    continue;
                }
            }

            return m->names[index];

        }

        if (m->next == NULL) {
            warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id.id);
            return NULL;
        } else {
            m = m->next;
        }
    }
}

/** Saves the ids of the parents of a genotype with a particular id to
 * the output array `output`.
 *
 * The function uses a bisection search on the gsc_AlleleMatrix where it should be
 * found. Searching by id is therefore fast.
 * @see gsc_get_from_ordered_pedigree_list()
 *
 * This function assumes that ids are never reshuffled in the gsc_SimData. This is true
 * as long as gsc_condense_allele_matrix's process of moving genotypes while retaining
 * their order is the only genotype-rearranging function in use.
 * This function's algorithm will need to
 * be revisited if different genotype-rearranging processes are implemented.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose parents are sought
 * @param output An array which the calling function can access where this function
 * will put its results.
 * @returns 0 when the id is successfully identified and at least one parent's
 * id is known, 1 if neither parent is known, and 2 if the ID passed in does
 * not exist. The ids of both parents if at least one parent is
 * known/nonzero are saved to the array `output`.
 */
int gsc_get_parents_of_id( const gsc_AlleleMatrix* start, const gsc_PedigreeID id, gsc_PedigreeID output[static 2]) {
    if (id.id == GSC_NO_PEDIGREE.id) {
		return 1;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
    const gsc_AlleleMatrix* m = start;
	while (1) {
		// try to find our id. Does this AM have the right range for it?
        if (m->n_genotypes != 0 && id.id >= m->ids[0].id && id.id <= m->ids[m->n_genotypes - 1].id) {
			// perform binary search to find the exact index.
            int index = gsc_get_from_ordered_pedigree_list(id, m->n_genotypes, m->ids);

			if (index < 0) {
				// search failed
                /*if (m->next == NULL) {
                    warning( "Unable to locate ID %d in simulation memory (genotype has likely been deleted): pedigree past this point cannot be determined.\n", id.id);
					return 2;
				} else {
					m = m->next;
                }*/
                continue;
            } else {

                if (m->pedigrees[0][index].id != GSC_NO_PEDIGREE.id || m->pedigrees[1][index].id != GSC_NO_PEDIGREE.id) {
                    output[0] = m->pedigrees[0][index];
                    output[1] = m->pedigrees[1][index];
                    return 0;
                }
                return 1; // if neither parent's id is known
            }

		}

		if (m->next == NULL) {
            warning( "Unable to locate ID %d in simulation memory (genotype has likely been deleted): pedigree past this point cannot be determined.\n", id.id);
			return 2;
		} else {
			m = m->next;
		}
	}
}

/** Search for genotypes with certain names in a linked list of gsc_AlleleMatrix and
 * save the ids of those names. Exits if any name cannot be found.
 *
 * This function must check every name in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param n_names the length of the array of names which are being sought.
 * @param names an array of strings/names whose ids we want to find. It must have  
 * at least [n_names] entries.
 * @param output pointer to an array with at least enough space to store [n_names]
 * PedigreeIDs. The ID corresponding to each name in `names` is saved at the corresponding
 * index in `output` when it is found by the function.
 */
void gsc_get_ids_of_names( const gsc_AlleleMatrix* start, const int n_names, const char** names, gsc_PedigreeID* output) {
    if (start == NULL || (start->n_genotypes <= 0 && start->next == NULL)) {
        warning("Invalid start parameter: gsc_AlleleMatrix* `start` must exist\n");
        return;
    }
    if (n_names < 1) {
        warning("Invalid n_names parameter: Search list length must be positive\n");
        return;
    }

    int found;
	//int ids = malloc(sizeof(int) * n_names);
    const gsc_AlleleMatrix* m;
	int i, j;

	for (i = 0; i < n_names; ++i) {
		found = GSC_FALSE;
        output[i] = GSC_NO_PEDIGREE;
		m = start;
		while (1) {
			// try to identify the name in this AM
			for (j = 0; j < m->n_genotypes; ++j) {
				if (strcmp(m->names[j], names[i]) == 0) {
					found = GSC_TRUE;
					output[i] = m->ids[j];
					break;
				}
			}

			if (found) {
				break;
			}
			if ((m = m->next) == NULL) {
				warning( "Didn't find the name %s\n", names[i]);
			}
		}
	}
}

/** Search for a genotype with parentage matching two given parent ids in a linked
 * list of gsc_AlleleMatrix, and return its index. Exits if such a child cannot be found.
 *
 * This function must check every genotype in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * the child is assumed to be found.
 * @param parent1id one of the parents of the genotype must have this id.
 * @param parent2id the other parent of the genotype must have this id.
 * @returns the index (0-based, starting at the start of `start`) of the first sequentially
 * located genotype whose parents match the two parent ids provided
 */
int gsc_get_index_of_child( const gsc_AlleleMatrix* start, const gsc_PedigreeID parent1id, const gsc_PedigreeID parent2id) {
    if (start == NULL || (start->n_genotypes <= 0 && start->next == NULL)) {
        warning("Invalid start parameter: gsc_AlleleMatrix* `start` must exist\n");
        return GSC_UNINIT;
    }
    const gsc_AlleleMatrix* m = start;
	int j, total_j = 0;

	while (1) {
		// try to identify the child in this AM
		for (j = 0; j < m->n_genotypes; ++j, ++total_j) {
            if ((parent1id.id == m->pedigrees[0][j].id && parent2id.id == m->pedigrees[1][j].id) ||
                    (parent1id.id == m->pedigrees[1][j].id && parent2id.id == m->pedigrees[0][j].id)) {
				return total_j;
			}
		}

		if ((m = m->next) == NULL) {
            error( "Didn't find the child of %d & %d\n", parent1id.id, parent2id.id);
		}
	}
}

/** Search for a genotype with a particular name in a linked
 * list of gsc_AlleleMatrix, and return its index. Exits if such a child cannot be found.
 *
 * This function must check every genotype in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * the genotype is assumed to be found.
 * @param name a string to match to the name of the target
 * @returns the index (0-based, starting at the start of `start`) of the first sequentially
 * located genotype whose name is the same as the provided name.
 */
int gsc_get_index_of_name( const gsc_AlleleMatrix* start, const char* name) {
    if (start == NULL || (start->n_genotypes <= 0 && start->next == NULL)) {
        warning("Invalid start parameter: gsc_AlleleMatrix* `start` must exist\n");
        return GSC_UNINIT;
    }
    const gsc_AlleleMatrix* m = start;
	int j, total_j = 0;

	while (1) {
		// try to identify the child in this AM
		for (j = 0; j < m->n_genotypes; ++j, ++total_j) {
			if (strcmp(m->names[j], name) == 0) {
				return total_j;
			}
		}

		if ((m = m->next) == NULL) {
			error( "Didn't find the name %s\n", name);
		}
	}
}

/** Get the id of a genotype by its index. The index is assumed to be 0-based,
 * starting at the first entry of `start` and continuing through the linked list
 * to which `start` belongs. Exits if the linked list does not contain that index.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * to locate the id at a particular index.
 * @param index the index of the target. It is assumed to start at 0 at the start
 * of the matrix in `start` and be incremented up until the last entry in the matrix
 * of the last entry in the linked list.
 * @returns the lifetime-unique id of the genotype found at that index.
 */
gsc_PedigreeID gsc_get_id_of_index( const gsc_AlleleMatrix* start, const int index) {
    if (start == NULL) {
        warning("Invalid start parameter: gsc_AlleleMatrix* `start` must exist\n");
        return GSC_NO_PEDIGREE;
    }
    const gsc_AlleleMatrix* m = start;
	int total_j = 0;

	while (1) {
		if (total_j == index) {
			return m->ids[0];
		} else if (total_j < index && total_j + m->n_genotypes > index) {
			return m->ids[index - total_j];
		}
		total_j += m->n_genotypes;

		if ((m = m->next) == NULL) {
			warning( "Didn't find the index %d\n", index);
            return GSC_NO_PEDIGREE;
		}
	}
}

/** Get the alleles of a genotype by its index. The index is assumed to be 0-based,
 * starting at the first entry of `start` and continuing through the linked list
 * to which `start` belongs. Exits if the linked list does not contain that index.
 *
 * @param start Pointer to the first of a linked list of gsc_AlleleMatrixes in which
 * to locate the id at a particular index.
 * @param index the index of the target. It is assumed to start at 0 at the start
 * of the matrix in `start` and be incremented up until the last entry in the matrix
 * of the last entry in the linked list.
 * @returns the alleles of the genotype that has idnex `index`, as a copy of the pointer
 * to the heap memory where the genotype is saved (so *don't* free the pointer returned
 * from this function). It points to a sequence of characters, ordered according to
 * the markers in the gsc_SimData to which the gsc_AlleleMatrix belongs.
 */
char* gsc_get_genes_of_index( const gsc_AlleleMatrix* start, const int index) {
	if (index < 0) {
		warning( "Invalid negative index %d\n", index);
		return NULL;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
    const gsc_AlleleMatrix* m = start;
	int total_j = 0;

	while (1) {
		if (total_j == index) {
			return m->alleles[0];
		} else if (total_j < index && total_j + m->n_genotypes > index) {
			return m->alleles[index - total_j];
		}
		total_j += m->n_genotypes;

		if ((m = m->next) == NULL) {
			warning( "Didn't find the index %d\n", index);
			return NULL;
		}
	}
}



/*-----------------------------------Groups----------------------------------*/

/** Combine a set of groups into one group.
 *
 * The function does so by setting the group membership of every genotype
 * belonging to one of the groups to the same group number.
 *
 * @shortnamed{combine_groups}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param list_len the number of groups to be combined
 * @param grouplist an array of at least [list_len] group numbers,
 * representing the groups that are to be combined.
 * @returns the group number of the new combined group.
 */
gsc_GroupNum gsc_combine_groups( gsc_SimData* d, const int list_len, const gsc_GroupNum* grouplist) {

    // Find the first group in the list that exists. In most use cases this will be the
    // first group in the list, so not too much of a performance penalty.
    gsc_GroupNum outGroup = GSC_NO_GROUP;
    int i = 0;
    for (; i < list_len; ++i) {
        gsc_GroupNum candidate = grouplist[i];
        gsc_BidirectionalIterator testit = gsc_create_bidirectional_iter(d,candidate);
        gsc_GenoLocation testloc = gsc_next_forwards(&testit);
        if (GSC_IS_VALID_LOCATION(testloc)) {
            outGroup = candidate;
            break;
        }
    }

    int remaininglistlen = list_len - i;
    if (remaininglistlen < 2) {
        return outGroup;
    } else if (remaininglistlen == 2) {
        if (grouplist[i].num == grouplist[i+1].num) {
            return outGroup;
        }
        gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,grouplist[i+1]);
        gsc_GenoLocation loc = gsc_next_forwards(&it);
        int anyFound = GSC_IS_VALID_LOCATION(loc);

        while (GSC_IS_VALID_LOCATION(loc)) {
            gsc_set_group(loc,outGroup);
            loc = gsc_next_forwards(&it);
        }

        if (anyFound) {
            d->n_groups--;
        }
        return outGroup;

	} else {
        GSC_CREATE_BUFFER(isDuplicate,int,remaininglistlen);
        memset(isDuplicate, GSC_FALSE, sizeof(int)*remaininglistlen);
        for (int ii = i; ii < list_len; ++ii) {
            for (int jj = ii+1; jj < list_len; ++jj) {
                if (grouplist[ii].num == grouplist[jj].num) {
                    isDuplicate[jj-i] = GSC_TRUE;
                }
            }
        }

        GSC_CREATE_BUFFER(anyFound,int,remaininglistlen);
        memset(anyFound, GSC_FALSE, sizeof(int)*remaininglistlen);

        gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,GSC_NO_GROUP);
        gsc_GroupNum cachedgroup = GSC_NO_GROUP; // just for speedier lookups. Groups tend to be stored contiguous in most simulations.
        gsc_GenoLocation loc = gsc_next_forwards(&it);

        while (GSC_IS_VALID_LOCATION(loc)) {
            if (gsc_get_group(loc).num == cachedgroup.num) {
                gsc_set_group(loc,outGroup);
            } else {
                for (int k = i+1; k < list_len; ++k) {
                    if (gsc_get_group(loc).num == grouplist[k].num) {
                        gsc_set_group(loc,outGroup);
                        cachedgroup = grouplist[k];
                        anyFound[k-i] = GSC_TRUE;
                        break;
                    }
                }
            }

            loc = gsc_next_forwards(&it);
        }

        int groupsgone = 0;
        for (int j = 0; j < remaininglistlen; ++j) {
            if (!isDuplicate[j] && anyFound[j]) {
                groupsgone++;
            }
        }
        d->n_groups -= groupsgone;
        GSC_DELETE_BUFFER(anyFound);
        GSC_DELETE_BUFFER(isDuplicate);
        return outGroup;
    }
}

/** Take a list of indexes and allocate the genotypes at those indexes
 * to a new group.
 *
 * Does not check if all the indexes are valid/if all indexes have successfully
 * had their groups changed.
 *
 * @shortnamed{make_group_from}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param index_list_len the number of indexes provided
 * @param genotype_indexes an array containing the global indexes (0-based, starting
 * at the first entry at `d->m`) of the genotypes to allocate to the new group.
 * This array must have at least [index_list_len] entries.
 * @returns the group number of the new group to which the provided indexes
 * have been allocated.
 */
gsc_GroupNum gsc_make_group_from( gsc_SimData* d, const int index_list_len, const unsigned int* genotype_indexes) {
    if (index_list_len < 1) {
        warning("Invalid index_list_len value: length of allocation list must be positive.\n");
        return GSC_NO_GROUP;
    }
	
    gsc_GroupNum newGroup = gsc_get_new_group_num(d);
    gsc_RandomAccessIterator it = gsc_create_randomaccess_iter(d, GSC_NO_GROUP);
    unsigned int invalidLocations = 0;
    for (int i = 0; i < index_list_len; ++i) {
		gsc_GenoLocation loc = gsc_next_get_nth(&it, genotype_indexes[i]);
        if (GSC_IS_VALID_LOCATION(loc)) {
            gsc_set_group(loc,newGroup);
        } else {
            invalidLocations++;
        }
	}

    if (invalidLocations > 0) {
        warning("%d indexes were invalid.\n",invalidLocations);
    }
    if (invalidLocations < index_list_len) {
        d->n_groups++;
    }

	gsc_delete_randomaccess_iter(&it);
	return newGroup;
}

/** Allocates the genotypes with a particular value of a label to a new group.
 *
 * Searches through every genotype (or every member of the group, if a specific group is
 * given) to find all genotypes with value `valueToSplit` in the `whichLabel`th label,
 * and puts those genotypes into a new group.
 *
 * Returns 0 for invalid parameters, or if no genotypes were found that fit the
 * criteria to be moved to the new group.
 *
 * @shortnamed{split_by_label_value}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group If `group` > 0, then only genotypes with group number `group`
 * AND `whichLabel` value `valueToSplit` will be moved to the new group. Otherwise,
 * all genotypes with `whichLabel` value `valueToSplit` will be moved to the next group.
 * @param whichLabel the label id of the relevant label.
 * @param valueToSplit the value of the label that defines the genotypes that will be
 * moved to the new group.
 * @returns the group number of the new group to which the genotypes with that value
 * for that label were allocated, or 0 if no genotypes that fit the criteria were found.
 */
gsc_GroupNum gsc_split_by_label_value( gsc_SimData* d, const gsc_GroupNum group, const gsc_LabelID whichLabel, const int valueToSplit) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return GSC_NO_GROUP;
    }
	
	gsc_GroupNum newGroup = gsc_get_new_group_num(d);
    int anyFound = GSC_FALSE;
	
    gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d, group);
	gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
		if (gsc_get_label_value(loc,labelIndex) == valueToSplit) {
			gsc_set_group(loc,newGroup);
			anyFound = GSC_TRUE;
		}
		
		loc = gsc_next_forwards(&it);
	}
	
    if (anyFound) {
        d->n_groups++;
        return newGroup;
    } else {
        return GSC_NO_GROUP;
    }

}

/** Allocates the genotypes with values of a label in a particular range to a new group.
 *
 * Searches through every genotype (or every member of the group, if a specific group is
 * given) to find all genotypes with value in the `whichLabel`th label between `valueLowBound`
 * and `valueHighBound` inclusive, and puts those genotypes into a new group.
 *
 * Returns 0 for invalid parameters, or if no genotypes were found that fit the
 * criteria to be moved to the new group.
 *
 * @shortnamed{split_by_label_range}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group If `group` > 0, then only genotypes with group number `group`
 * AND `whichLabel` value `valueToSplit` will be moved to the new group. Otherwise,
 * all genotypes with `whichLabel` value `valueToSplit` will be moved to the next group.
 * @param whichLabel the label id of the relevant label.
 * @param valueLowBound the minimum value of the label for the genotypes that will be
 * moved to the new group.
 * @param valueHighBound the maximum value of the label for the genotypes that will be
 * moved to the new group.
 * @returns the group number of the new group to which the genotypes with that value
 * for that label were allocated, or 0 if no genotypes that fit the criteria were found.
 */
gsc_GroupNum gsc_split_by_label_range( gsc_SimData* d, const gsc_GroupNum group, const gsc_LabelID whichLabel, const int valueLowBound, const int valueHighBound) {
    int labelIndex;
    if (whichLabel.id == GSC_NO_LABEL.id || (labelIndex = gsc_get_index_of_label(d, whichLabel)) < 0) {
        warning( "Nonexistent label %d\n", whichLabel.id);
        return GSC_NO_GROUP;
    }
    if (valueLowBound > valueHighBound) {
        warning( "Empty range %d to %d: no group created\n", valueLowBound, valueHighBound);
        return GSC_NO_GROUP;
    }

    gsc_GroupNum newGroup = gsc_get_new_group_num(d);
    int anyFound = GSC_FALSE; 
	
    gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d, group);
	gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
		if (gsc_get_label_value(loc,labelIndex) >= valueLowBound &&
				gsc_get_label_value(loc,labelIndex) <= valueHighBound) {
			gsc_set_group(loc,newGroup);
			anyFound = GSC_TRUE;
		}
		
		loc = gsc_next_forwards(&it);
	}

    if (anyFound) {
        d->n_groups++;
        return newGroup;
    } else {
        return GSC_NO_GROUP; // no values with that label
    }
}


/** Split by some quality (generic function)
 *
 *  Somequality: you don't know how many variants of that quality there are,
 *  so you don't initially know how many groups you will have.
 * @see gsc_scaffold_split_by_someallocation
 *
 * Takes a parameter @a somequality_tester that uses a gsc_GenoLocation, somequality_data,
 * the maximum possible number of groups you're splitting into, the current number of new groups
 * that have allocations, and the list of potential group numbers for new groups. It should return
 * a group number (taken from some position in that final parameter) if the genotype
 * at that gsc_GenoLocation is to be added to one of the groups that already has allocations,
 * or return GSC_NO_GROUP if it is to be allocated to one of the groups beyond the [value of the
 * second-last parameter]th (This is so that this caller function can keep allocate new group
 * numbers and keep track of how many subgroups have been found).
 *
 * If the number of variants ends up being larger than @a maxentries_results, further variants
 * will still be allocated to new groups, but will not saved to @a results. The calling function
 * can know this is the case if the returned value is larger than the parameter @a maxentries_results.
 *
 * @return number of groups created. May be less or more than @a maxentries_results.
*/
unsigned int gsc_scaffold_split_by_somequality( gsc_SimData* d, const gsc_GroupNum group_id,  
		void* somequality_data,
        gsc_GroupNum (*somequality_tester)(gsc_GenoLocation, void*, unsigned int, unsigned int, gsc_GroupNum*),
        unsigned int maxentries_results, gsc_GroupNum* results) {
	
	// Access existing groups (to be used to find unused group numbers,
	// and to find maximum number of groups we'd be able to create)
    GSC_CREATE_BUFFER(currentgroups,gsc_GroupNum,d->n_groups);
    GSC_CREATE_BUFFER(currentsizes,unsigned int,d->n_groups);
    int n_groups = gsc_get_existing_group_counts(d, currentgroups, currentsizes);
    int bookmark = 0;
    gsc_GroupNum nextgroup = GSC_NO_GROUP;

    unsigned int splitgroupsize = 0;
    for (int i = 0; i < n_groups; ++i) {
        if (currentgroups[i].num == group_id.num) {
            splitgroupsize = currentsizes[i];
            //GSC_FREE(currentsizes);
            break;
        }
    }
    if (splitgroupsize < 1) {
        return 0;
    }
	
    GSC_DELETE_BUFFER(currentsizes);
	unsigned int subgroupsfound = 0;
    GSC_CREATE_BUFFER(outgroups,gsc_GroupNum,splitgroupsize);
	
	gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,group_id);
    gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
        // Return group number if it should be assigned to an already-extant group. Otherwise return GSC_NO_GROUP and this generic caller function will allocated it one.
		gsc_GroupNum assignedgroup = somequality_tester(loc, somequality_data,splitgroupsize, subgroupsfound, outgroups);
		
        if (assignedgroup.num == GSC_NO_GROUP.num) {
			nextgroup = gsc_get_next_free_group_num(n_groups,currentgroups,&bookmark,nextgroup);
			assignedgroup = nextgroup;
			outgroups[subgroupsfound] = nextgroup;
			subgroupsfound++;	
		}
		
		gsc_set_group(loc,assignedgroup);
		
        loc = gsc_next_forwards(&it);
    }
	
    GSC_DELETE_BUFFER(currentgroups);
    d->n_groups += subgroupsfound - 1;
	
	if (maxentries_results < subgroupsfound) { 
		memcpy(results,outgroups,sizeof(gsc_GroupNum)*maxentries_results);
		fprintf(stderr,"Output vector size is not large enough to hold all created groups: "
                " output list of gsc_GroupNums has been truncated.\n");
	} else {
		memcpy(results,outgroups,sizeof(gsc_GroupNum)*subgroupsfound);
	}
    GSC_DELETE_BUFFER(outgroups);
	return subgroupsfound;
}

static gsc_GroupNum gsc_helper_split_by_quality_halfsibtemplate(gsc_GenoLocation loc, void* datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results, gsc_PedigreeID (*getparent)(gsc_GenoLocation)) {
	gsc_PedigreeID* familyidentities = (gsc_PedigreeID*) datastore;
	
	for (int j = 0; j < groupsfound; ++j) {
		if (getparent(loc).id == familyidentities[j].id) {
            return results[j];
		}
	}
	
	familyidentities[groupsfound] = getparent(loc);
    return GSC_NO_GROUP;
}

static gsc_GroupNum gsc_helper_split_by_quality_halfsib1(gsc_GenoLocation loc, void* datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results) {
	return gsc_helper_split_by_quality_halfsibtemplate(loc,datastore,maxgroups,groupsfound,results,gsc_get_first_parent);
}

static gsc_GroupNum gsc_helper_split_by_quality_halfsib2(gsc_GenoLocation loc, void* datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results) {
	return gsc_helper_split_by_quality_halfsibtemplate(loc,datastore,maxgroups,groupsfound,results,gsc_get_second_parent);
}

/** Split a group into families of half-siblings by shared first or second parent.
 * 
 * Split a group into a set of smaller groups, each containing the
 * genotypes from the original group that share one parent.
 * The shared parent can be either the first or second parent,
 * based on the value of the parameter parent. That is, if parent is 1,
 * within the halfsib families produced, all genotypes will share the
 * same first parent, but may have different second parents.
 * The number of new groups produced depends on the number of unique
 * first/second parents in the set of genotypes in the provided group.
 *
 * Individuals with unknown parent will be grouped together.
 *
 * Stops executing if group is empty or has only one member.
 *
 * If more than maxentries_results groups are created by this function, only 
 * that many results will be saved into the results vector, though all the 
 * one-member groups will be created.
 *
 * @shortnamed{split_into_halfsib_families}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param parent 1 to group together genotypes that share the same first parent,
 * 2 group those with the same second parent. Raises an error
 * if this parameter is not either of those values.
 * @param maxentries_results maximum number of group numbers that can be saved 
 * into the results vector.
 * @param results Pointer to a vector into which to save the identifiers of the newly
 * created family groups. Should have at least enough space for [maxentries_results]
 * identifiers.
 * @return the number of halfsib families identified/number of groups created.
 */
unsigned int gsc_split_into_halfsib_families( gsc_SimData* d, const gsc_GroupNum group_id, const int parent, unsigned int maxentries_results, gsc_GroupNum* results) {
	if (!(parent == 1 || parent == 2)) {
		warning( "Value error: `parent` must be 1 or 2.");
        results = NULL;
		return 0;
	}
	
	//gsc_PedigreeID* familyidentities = gsc_malloc_wrap(sizeof(gsc_PedigreeID)*maxgroups);
    int maxgroups = gsc_get_group_size(d, group_id); // sadinefficient we have to do this
    GSC_CREATE_BUFFER(familyidentities,gsc_PedigreeID,maxgroups);
	
	unsigned int gcount;
	if (parent == 1) {
        gcount = gsc_scaffold_split_by_somequality(d, group_id, (void*)familyidentities,gsc_helper_split_by_quality_halfsib1, maxentries_results, results);
	} else {
        gcount = gsc_scaffold_split_by_somequality(d, group_id, (void*)familyidentities, gsc_helper_split_by_quality_halfsib2, maxentries_results, results);
	}	

    GSC_DELETE_BUFFER(familyidentities);
	return gcount;
}

static gsc_GroupNum gsc_helper_split_by_quality_family(gsc_GenoLocation loc, void* datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results) {
    gsc_PedigreeID** familyidentities = (gsc_PedigreeID**) datastore;
	
	for (int j = 0; j < groupsfound; ++j) {
		if (gsc_get_first_parent(loc).id == familyidentities[0][j].id &&
              gsc_get_second_parent(loc).id == familyidentities[1][j].id) {
            return results[j];
		}
	}
	
	familyidentities[0][groupsfound] = gsc_get_first_parent(loc);
	familyidentities[1][groupsfound] = gsc_get_second_parent(loc);
    return GSC_NO_GROUP;
}

/** Split a group into families by their pedigrees.
 *
 * Split a group into a set of smaller groups, each of which contains the
 * genotypes in the original group that share a particular pair of parents.
 * The number of new groups produced depends on the number of parent-combinations
 * in the set of genotypes in the provided group.
 *
 * Individuals with both parents unknown will be grouped together.
 *
 * If more than maxentries_results groups are created by this function, only 
 * that many results will be saved into the results vector, though all the 
 * one-member groups will be created.
 *
 * Stops executing if group is empty or has only one member.
 *
 * @shortnamed{split_into_families}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param maxentries_results maximum number of group numbers that can be saved 
 * into the results vector.
 * @param results Pointer to a vector into which to save the identifiers of the newly
 * created family groups. Should have at least enough space for [maxentries_results]
 * identifiers.
 * @return the number of families identified/number of family groups created.
 * Also serves as the length of the array in @a results
 */
unsigned int gsc_split_into_families(gsc_SimData* d, const gsc_GroupNum group_id, unsigned int maxentries_results, gsc_GroupNum* results) {
    gsc_PedigreeID* familyidentities[2];
    unsigned int maxgroups = gsc_get_group_size(d, group_id); // sadinefficient we have to do this
    if (maxgroups < 2) {
        return 0;
    }

    GSC_CREATE_BUFFER(p1identity,gsc_PedigreeID,maxgroups);
    GSC_CREATE_BUFFER(p2identity,gsc_PedigreeID,maxgroups);
	familyidentities[0] = p1identity;
	familyidentities[1] = p2identity;
	
    unsigned int out = gsc_scaffold_split_by_somequality(d, group_id, (void*)familyidentities, gsc_helper_split_by_quality_family, maxentries_results, results);

    GSC_DELETE_BUFFER(p1identity);
    GSC_DELETE_BUFFER(p2identity);

    return out;
}

static gsc_GroupNum gsc_helper_split_by_quality_individuate(gsc_GenoLocation loc, void* datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results) {
    return GSC_NO_GROUP;
}

/** Split a group into n one-member groups.
 * 
 * Give every individual in the group a new group number that does not
 * belong to any other existing group (thereby allocating each genotype
 * in the group to a new group of 1).
 *
 * Stops executing if group is empty or has only one member.
 *
 * If more than maxentries_results groups are created by this function, only 
 * that many results will be saved into the results vector, though all the 
 * one-member groups will be created.
 *
 * @shortnamed{split_into_individuals}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param maxentries_results maximum number of group numbers that can be saved 
 * into the results vector.
 * @param results Pointer to a vector into which to save the identifiers of the newly
 * created family groups. Should have at least enough space for [maxentries_results]
 * identifiers.
 * @return the number of groups created (in this case, the same as the size of
 * the original group). Also serves as the length of the array in @a results
 */
unsigned int gsc_split_into_individuals( gsc_SimData* d, const gsc_GroupNum group_id, unsigned int maxentries_results, gsc_GroupNum* results) {
	// **individuate** (verb): to make individuals of.
	// yeah sorry.
    return gsc_scaffold_split_by_somequality(d, group_id, NULL, gsc_helper_split_by_quality_individuate, maxentries_results, results);
}


/** Split a group into two groups of equal size (or size differing only
 * by one, if the original group had an odd number of members) using a
 * random permutation of the group members to determine which goes where.
 *
 * Of the two groups produced, one has the same group number as the original
 * group (parameter group_id) and the other has the return value as its
 * group number.
 *
 * If the original group size was odd, the new group/the return value
 * will have the slightly smaller size.
 *
 * A more general approach to this task: @ref gsc_split_evenly_into_n()
 *
 * An alternate approach to splitting a group in two: @ref gsc_split_randomly_into_two()
 *
 * @shortnamed{split_evenly_into_two}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @returns the group number of the new group to which half the members
 * of the old group have been allocated.
 */
gsc_GroupNum gsc_split_evenly_into_two(gsc_SimData* d, const gsc_GroupNum group_id) {
	// get the shuffle to be our even allocations
    unsigned int size = gsc_get_group_size(d, group_id);
    if (size < 2) {
        if (size < 1) {
            warning("Group %d does not exist\n", group_id.num);
        } else {
            warning("Group %d has only one member so can't be split\n", group_id.num);
        }
        return GSC_NO_GROUP;
    }

    unsigned int even_half = size / 2;
    GSC_CREATE_BUFFER(allocations,unsigned int,size);
    for (unsigned int i = 0; i < size; ++i) {
		allocations[i] = i;
	}
    shuffle_up_to( allocations, size, even_half);

    gsc_GroupNum new_group = gsc_get_new_group_num(d);
    //int anyFound = GSC_FALSE;
	
    gsc_RandomAccessIterator it = gsc_create_randomaccess_iter(d,group_id);
    for (unsigned int i = 0; i < even_half; ++i) {
        gsc_GenoLocation loc = gsc_next_get_nth(&it,allocations[i]);
        if (GSC_IS_VALID_LOCATION(loc)) {
            gsc_set_group(loc,new_group);
            //anyFound = GSC_TRUE;
        }
    }

    GSC_DELETE_BUFFER(allocations);
    gsc_delete_randomaccess_iter(&it);

    //if (anyFound) {
        d->n_groups++;
        return new_group;
    //} else {
    //    return GSC_NO_GROUP;
    //}
}


/** Split by some allocator (generic function)
 *
 *  Allocator: you know how many groups you're splitting into from the beginning.
 *  @see gsc_scaffold_split_by_somequality
 *
 * Takes a parameter @a someallocator that uses a gsc_GenoLocation, gsc_SimData, someallocator_data,
 * the total number of groups you're splitting into, the current number of new groups
 * that have allocations (as a pointer so someallocator can modify it to hold the number
 * of subgroups that have been found, which will become the return value of this function),
 * and the list of group numbers of new groups. It should return
 * a group number (taken from some position in that final parameter) if the genotype
 * at that gsc_GenoLocation is to be added to one of the groups that already has allocations,
 * or return GSC_NO_GROUP if allocation fails. If allocation fails the genotype will remain in
 * the original group.
 *
 * @return number of groups created. May be @a n_outgroups or less.
*/
unsigned int gsc_scaffold_split_by_someallocation( gsc_SimData* d, const gsc_GroupNum group_id, void* someallocator_data,
        gsc_GroupNum (*someallocator)(gsc_GenoLocation, gsc_SimData*, void*, unsigned int, unsigned int*, gsc_GroupNum*),
        unsigned int n_outgroups, gsc_GroupNum* outgroups) {

    // get the n group numbers
    gsc_get_n_new_group_nums(d, n_outgroups, outgroups);

    unsigned int subgroupsfound = 0;
    unsigned int allocationfailures = 0;

    gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d,group_id);
    gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
        gsc_GroupNum assignedgroup = someallocator(loc, d, someallocator_data, n_outgroups, &subgroupsfound, outgroups);
        if (assignedgroup.num != GSC_NO_GROUP.num) {
            gsc_set_group(loc,assignedgroup);
        } else {
            allocationfailures++;
        }

        loc = gsc_next_forwards(&it);
    }

    if (subgroupsfound > 1) {
        d->n_groups += subgroupsfound - 1;
    }
    if (allocationfailures > 0) {
        fprintf(stderr,"While splitting group %i, %i allocations to new groups failed so they remain in the original group.\n",
                group_id.num, allocationfailures);
    }
    return subgroupsfound;

}


static gsc_GroupNum gsc_helper_split_by_allocator_knowncounts(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups,
                                           unsigned int* subgroupsfound, gsc_GroupNum* outgroups) {
	unsigned int* cumulative_counts = (unsigned int*) datastore;
    *subgroupsfound = n_outgroups;
	unsigned int randpos = round(unif_rand() * (cumulative_counts[n_outgroups-1] - 1));
    gsc_GroupNum chosengroup = GSC_NO_GROUP;
    unsigned int j = 0;
    for (; j < n_outgroups; ++j) {
		if (randpos < cumulative_counts[j]) {
            chosengroup = outgroups[j];
			break;
		}
	}
    for (; j < n_outgroups; ++j) {
		cumulative_counts[j]--;
	}
	return chosengroup;
}

/** Split a group into n groups of equal size (or size differing only
 * by one, if n does not perfectly divide the group size.), using a
 * random permutation of the group members to determine which goes where.
 *
 * Of the split groups produced, the first has the same group number as the original
 * group (parameter group_id).
 *
 * A more general approach to this task: @ref gsc_split_into_buckets()
 *
 * @shortnamed{split_evenly_into_n}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 * @return the number of groups created by this process. 
 */
unsigned int gsc_split_evenly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, const int n, gsc_GroupNum* results) {
    if (n <= 1) {
        warning( "Invalid n value: number of fractions into which to split group must be positive.\n");
		return 0;
	}

    int size = gsc_get_group_size(d, group_id); // sadinefficient we have to do this.

    // get the shuffle to be our even allocations
	int each_size = size / n;
	int extra = size % n;
    GSC_CREATE_BUFFER(boxes,int,n);
	for (int i = 0; i < n; ++i) {
		boxes[i] = each_size;
		if (i < extra) {
			boxes[i]++;
		}
		if (i > 0) {
			boxes[i] += boxes[i-1];
		}
	}
	
    unsigned int out;
	if (results == NULL) {
        GSC_CREATE_BUFFER(outgroups,gsc_GroupNum, n);
        out = gsc_scaffold_split_by_someallocation(d, group_id, (void*) boxes, gsc_helper_split_by_allocator_knowncounts, n, outgroups);
        GSC_DELETE_BUFFER(outgroups);
	} else {
        out = gsc_scaffold_split_by_someallocation(d, group_id, (void*) boxes, gsc_helper_split_by_allocator_knowncounts, n, results);
	}
    GSC_DELETE_BUFFER(boxes);
    return out;
}

/** Split a group into n groups of equal size (or size differing only
 * by one, if n does not perfectly divide the group size), using a
 * random permutation of the group members to determine which goes where.
 *
 * Of the split groups produced, the first has the same group number as the original
 * group (parameter group_id).
 *
 * The number of members staying
 * in the old group (group_id) is counts[0]. The number going to the
 * first new group is counts[1], etc.. The number going to the
 * nth group is group_id's group size - sum(counts).
 *
 * The function calculates a random permutation of the group members,
 * then uses cumulative sums to determine to which group the group member
 * is allocated. If the sum of the desired group sizes adds up to more than
 * the number of group members, a warning is raised, and the group numbers for
 * which the cumulative sum of counts is greater than the group size will not
 * be allocated members. That is, the group capacities are filled from first to
 * last, leaving later groups unfilled if there are not enough group members to
 * occupy all capacities.
 *
 * @shortnamed{split_into_buckets}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param counts pointer to an array of length at least n-1 containing the number of
 * members to allocate to each group. The number of members in the last
 * group is group_id's group size - sum(counts).
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 * @return the number of groups created by this process. 
 */
unsigned int gsc_split_into_buckets(gsc_SimData* d, const gsc_GroupNum group_id, const int n, const int* counts, gsc_GroupNum* results) {
    if (n <= 1) {
        warning( "Invalid n value: number of fractions into which to split group must be positive.\n");
        return 0;
    }

    GSC_CREATE_BUFFER(cumulative_counts,int,n);
    cumulative_counts[n-1] = gsc_get_group_size(d, group_id);
    int sum = 0;
    for (int j = 0; j < n - 1; ++j) {
        sum += counts[j];
        cumulative_counts[j] = sum;
    }
    if (cumulative_counts[n-2] > cumulative_counts[n-1]) {
        warning( "Provided capacities are larger than actual group: some buckets will not be filled\n");
    }

    unsigned int gcount;
	if (results == NULL) {
        GSC_CREATE_BUFFER(outgroups,gsc_GroupNum,n);
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, (void*) cumulative_counts, gsc_helper_split_by_allocator_knowncounts, n, outgroups);
        GSC_DELETE_BUFFER(outgroups);
	} else {
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, (void*) cumulative_counts, gsc_helper_split_by_allocator_knowncounts, n, results);
	}

    GSC_DELETE_BUFFER(cumulative_counts);
    return gcount;
}

/** Flip a coin for each member of the group to decide if it should
 * be moved to the new group.
 *
 * There is no guarantee that there will be any genotypes in the new group
 * (if all coin flips were 0) or any genotypes in the old group (if all
 * coin flips were 1). There is no guarantee the two groups will be
 * near the same size.
 *
 * This could be useful for allocating a sex to genotypes.
 *
 * A more general approach to this task: @ref gsc_split_randomly_into_n()
 *
 * An alternate approach to splitting a group in two: @ref gsc_split_evenly_into_two()
 *
 * @shortnamed{split_randomly_into_two}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @returns the group number of the new group to which some members
 * of the old group may have been randomly allocated.
 */
gsc_GroupNum gsc_split_randomly_into_two(gsc_SimData* d, const gsc_GroupNum group_id) {
    gsc_GroupNum outGroup = gsc_get_new_group_num(d);
	int anyFound = GSC_FALSE;
	
    gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d, group_id);
	gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
		anyFound = GSC_TRUE;
		if ((unif_rand() > 0.5)) {
			gsc_set_group(loc,outGroup);
		}
		loc = gsc_next_forwards(&it);
	}

	if (anyFound) {
        d->n_groups++;
		return outGroup;
	} else {
        return GSC_NO_GROUP;
	}
}


static gsc_GroupNum gsc_helper_split_by_allocator_equalprob(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups,
                                             unsigned int* subgroupsfound, gsc_GroupNum* outgroups) {
	int randgroup = round(unif_rand() * (n_outgroups-1));
    if (randgroup < *subgroupsfound) {
        return outgroups[randgroup];
	} else {
        (*subgroupsfound)++;
        return outgroups[*subgroupsfound-1];
	}
}

/** Allocate each member of the group to
 * one of n groups with equal probability.
 *
 * There is no guarantee that all groups will have members.
 * There is no guarantee the groups will be near the same size.
 *
 * Each genotype has equal probability of being allocated to
 * each of n groups. The old group number (group_id) is included
 * as one of these n possible groups.
 *
 * To split by uneven probabilities instead: @ref gsc_split_by_probabilities()
 *
 * @shortnamed{split_randomly_into_n}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 * @return the number of groups created by this process. 
 */
unsigned int gsc_split_randomly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, const int n, gsc_GroupNum* results) {
    if (n <= 1) {
        warning( "Invalid n value: number of fractions in which to split group must be positive.\n");
        return 0;
    }

    unsigned int gcount;
	if (results == NULL) {
        GSC_CREATE_BUFFER(outgroups,gsc_GroupNum,n);
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, NULL, gsc_helper_split_by_allocator_equalprob, n, outgroups);
        GSC_DELETE_BUFFER(outgroups);
	} else {
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, NULL, gsc_helper_split_by_allocator_equalprob, n, results);
	}
    return gcount;
}


static gsc_GroupNum gsc_helper_split_by_allocator_unequalprob(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups,
                                            unsigned int* subgroupsfound, gsc_GroupNum* outgroups) {
	double* cumulative_probs = (double*) datastore;
    *subgroupsfound = n_outgroups;
    double randdraw = unif_rand();
    for (int j = 0; j < n_outgroups; ++j) {
		if (randdraw < cumulative_probs[j]) {
			return outgroups[j];
		}
	}
    // This should not happen if cumulative probs are valid
    return GSC_NO_GROUP;
}

/** Allocate each member of the group to
 * one of n groups with custom probabilities for each group.
 *
 * There is no guarantee that all groups will have members.
 * There is no guarantee the groups will be near the same size.
 *
 * The probability of staying
 * in the old group (group_id) is probs[0]. The probability of going to the
 * first new group is probs[1], etc.. The probability of going to the
 * nth group is 1 minus the sum of all probabilities in probs.
 *
 * The function draws from a uniform distribution, then uses cumulative
 * sums to determine to which group the group member is allocated.
 * If the sum of the probabilities in probs adds up to more than 1,
 * a warning is raised, and the group numbers for which the cumulative
 * sum of probs is greater than 1 have no chance of being allocated members.
 *
 * @shortnamed{split_by_probabilities}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param probs pointer to an array of length n-1 containing the probability
 * of being allocated to each group. The probability of going to the last
 * group is 1 - sum(probs).
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 * @return the number of groups created by this process. 
 */
unsigned int gsc_split_by_probabilities(gsc_SimData* d, const gsc_GroupNum group_id, const int n, const double* probs, gsc_GroupNum* results) {
    if (n <= 1) {
        warning( "Invalid n value: number of fractions in which to split group must be positive.\n");
        return 0;
    }

	// Check the probabilities
    GSC_CREATE_BUFFER(cumulative_probs,double,n);
    cumulative_probs[n-1] = 1.0;
	double sum = 0;
    for (int j = 0; j < n-1; ++j) {
		sum += probs[j];
		cumulative_probs[j] = sum;
		if (cumulative_probs[j] >= 1) {
            warning( "Provided probabilities add up to 1 or more: some buckets will not be filled\n");
            for (; j < n-1; ++j) {
                cumulative_probs[j] = 1;
            }
          //don't bother to calculate more
          break;
		}
	}
	
    unsigned int gcount;
	if (results == NULL) {
        GSC_CREATE_BUFFER(outgroups,gsc_GroupNum,n);
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, (void*) cumulative_probs, gsc_helper_split_by_allocator_unequalprob, n, outgroups);
        GSC_DELETE_BUFFER(outgroups);
	} else {
        gcount = gsc_scaffold_split_by_someallocation(d, group_id, (void*) cumulative_probs, gsc_helper_split_by_allocator_unequalprob, n, results);
	}

    GSC_DELETE_BUFFER(cumulative_probs);
    return gcount;
}


/** Identify group numbers that currently have members.
 *
 *  It saves the extant groups to the @a output vector, and also
 *  realigns @a d->n_groups with the true number of groups.
 *
 *  Frankly this should be deprecated at some point because you can now
 *  call @ref gsc_get_existing_group_counts with the output vector for the counts
 *  set to NULL and get exactly the same result. (in fact that's the current
 *  implementation of this function).
 *
 * @shortnamed{get_existing_groups}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param output Pointer to a location to save the vector of extant groups.
 * It must have at least space for d->n_groups group numbers. Else, if it is
 * NULL, then the group numbers of existing groups are not saved and the only
 * purpose of the function is to return the number of existing groups and update
 * d->n_groups to this number.
 * @returns The number of entries that have been placed in `output`.
 * That is, the number of existing groups.
 */
int gsc_get_existing_groups( gsc_SimData* d, gsc_GroupNum* output) {
    return gsc_get_existing_group_counts(d, output, NULL);
}

/** Identify group numbers that currently have members, and how many members they have
 *
 *  It saves the extant groups to the output vectors, and also
 *  realigns @a d->n_groups with the true number of groups.
 *
 * @shortnamed{get_existing_group_counts}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param out_groups Pointer to a location to save the vector of extant groups.
 * It must have at least space for d->n_groups group numbers. Else, if it is
 * NULL, then the group numbers of existing groups are not saved.
 * @param out_sizes Pointer to a location to save the vector of those groups' sizes.
 * It must have at least space for d->n_groups group sizes. Else, if it is NULL, then
 * the counts of existing groups are not saved.
 * @returns The number of existing groups.
 */
int gsc_get_existing_group_counts( gsc_SimData* d, gsc_GroupNum* out_groups, unsigned int* out_sizes) {
    GSC_CREATE_BUFFER(buckets,unsigned int,d->n_groups+1); // this also creates bucketscap, initalised to d->n_groups+1.
    memset(buckets,0,sizeof(unsigned int)*bucketscap);
    unsigned int filledbuckets = 0;

    gsc_BidirectionalIterator it = gsc_create_bidirectional_iter(d, GSC_NO_GROUP);
    gsc_GenoLocation loc = gsc_next_forwards(&it);
    while (GSC_IS_VALID_LOCATION(loc)) {
        gsc_GroupNum g = gsc_get_group(loc);
        // Unless all group numbers are consecutive starting at 1, the buckets array will need to be resized at some point.
        if (g.num >= bucketscap) {
            unsigned int oldcap = bucketscap;
            unsigned int newbucketcapacity = bucketscap;
            while (g.num >= newbucketcapacity) {
                newbucketcapacity *= 2;
            }
            GSC_STRETCH_BUFFER(buckets,newbucketcapacity);
            if (g.num >= bucketscap) {
                warning("Memory allocation failed. Not all groups found.\n");
                break;
            }
            memset(buckets+oldcap,0,sizeof(unsigned int)*(bucketscap-oldcap));

        }

        buckets[g.num] += 1;
        if (buckets[g.num] == 1) {
            ++filledbuckets;
        }

        loc = gsc_next_forwards(&it);
    }

    // Now save to output and sort.
    unsigned int capacity = filledbuckets;
    if (capacity > d->n_groups) {
        fprintf(stderr,"Found more groups than expected - gsc_SimData.n_groups is outdated somewhere."
                " Trimming output of get_existing_group_ to avoid a crash: not all groups may be shown.\n");
        capacity = d->n_groups;
    }
    unsigned int g_index = 0;
    for (int i = 0; i < bucketscap; ++i) {
        if (buckets[i]) {
            /*if (g_index >= capacity) {
                warning("Found more groups than just a moment ago.");
                --g_index;
                break;
            }*/

            if (out_groups != NULL) {
                out_groups[g_index] = (gsc_GroupNum){.num=i};
            }
            if (out_sizes != NULL) {
                out_sizes[g_index] = buckets[i];
            }
            ++g_index;
        }
    }

    //qsort(*out_groups, g_index, sizeof(gsc_GroupNum), gsc_helper_ascending_gsc_GroupNum_comparer);
    /*for (int i = 0; i < g_index; ++i) {
        (*out_sizes)[i] = buckets[(*out_groups)[i].num];
    }*/
    GSC_DELETE_BUFFER(buckets);
    d->n_groups = g_index;

    return g_index;
}


/** Iterator to get the next currently-free group number.
 *
 * In the vein of @a gsc_get_new_group_num or @a gsc_get_n_new_group_nums, a function
 * that's effectively an iterator for unused group numbers. It's faster than
 * calling either of those multiple times as it can re-use its 'cursor' position
 * to not have to scan the array from the start every time, and it also gets
 * to reuse a set of existing groups collected only once.
 *
 * You will probably want to call @a gsc_get_existing_groups() before this
 * one to get values for the parameters @a n_existing_groups and @a existing_groups.
 * This is a function for internal use, so there won't be much/any bounds or
 * error checking
 *
 * @param n_existing_groups Length of the existing_groups array (eg, return value of
 * @a gsc_get_existing_groups
 * @param existing_groups Pointer to an array of gsc_GroupNums active in simulation (eg, output
 * value of @a gsc_get_existing_groups).
 * @param cursor Index in @a existing_groups that this function has currently checked up to. This
 * value will be updated in the calling function.
 * @param previous Last group number returned by this function, or GSC_NO_GROUP on first call.
 * @return the next sequential currently-unused (according to the memberships in @a existing_groups)
 * group number.
 */
gsc_GroupNum gsc_get_next_free_group_num( const int n_existing_groups, const gsc_GroupNum* existing_groups, int* cursor,  gsc_GroupNum previous) {
    if (existing_groups == NULL) return GSC_NO_GROUP;

    if (*cursor == GSC_UNINIT) {
        *cursor = 0;
    }

    gsc_GroupNum nextgroup = (gsc_GroupNum){.num=previous.num+1};
    // a check here in case previous seems invalid. We need previous so we don't get stuck in a loop
    // of giving the same next 'free' number, but we know what a lower bound on its number should be
    // based on where the cursor is.
    if (*cursor > 0 && nextgroup.num <= existing_groups[(*cursor) - 1].num) {
        nextgroup.num = existing_groups[(*cursor) - 1].num + 1;
    }

    while (*cursor < n_existing_groups) {
        if (nextgroup.num < existing_groups[*cursor].num) {
            break;
        }

        ++(*cursor);
        ++nextgroup.num;
    }
    return nextgroup;
}

/** Function to identify the next sequential integer that does not
 * identify a group that currently has member(s).
 *
 * This calls gsc_get_existing_groups() every time, so for better
 * speed, functions that do repeated group creation, like
 * gsc_split_into_individuals(), are
 * recommended to use gsc_get_n_new_group_nums() (if they know the
 * number of groups they need) or their own implementation, rather
 * than calling this function repeatedly.
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @return the next sequential currently-unused group number.
 */
gsc_GroupNum gsc_get_new_group_num( gsc_SimData* d) {
    // Make sure we get all existing groups
    if (d->m == NULL || (d->m->n_genotypes == 0 && d->m->next == NULL)) {
        return (gsc_GroupNum){.num=1};
    }

    int n_groups = (d->n_groups > 0) ? d->n_groups : 5;
    GSC_CREATE_BUFFER(existing_groups,gsc_GroupNum,n_groups);
    n_groups = gsc_get_existing_groups(d, existing_groups);

	int i = 0;
	int gn = 1;

	while (i < n_groups) {
        if (gn < existing_groups[i].num) {
			break;
		}

		++i;
		++gn;
	}
    GSC_DELETE_BUFFER(existing_groups);
    return (gsc_GroupNum){.num=gn};
}


/** Function to identify the next n sequential integers that do not
 * identify a group that currently has member(s).
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param n the number of group numbers to generate
 * @param result pointer to an array of length at least n where
 * the new group numbers generated can be saved.
 */
void gsc_get_n_new_group_nums( gsc_SimData* d, const int n, gsc_GroupNum* result) {
    // Make sure we get all existing groups
    int n_groups;
    GSC_CREATE_BUFFER(existing_groups,gsc_GroupNum,d->n_groups);
    n_groups = gsc_get_existing_groups(d, existing_groups);

    int existingi = 0;
    int gn = 0;

    // i: current index of `results` (the array of currently empty group numbers)
    // gn: group number being checked against existing_groups. if not in there is added to
    //     the list of results
    // existingi: current index of existing_groups
    for (int i = 0; i < n; ++i) {
        ++gn;
        while (existingi < n_groups) {
            if (gn < existing_groups[existingi].num) {
                break;
            }

            ++existingi;
            ++gn;
        }
        result[i] = (gsc_GroupNum){.num=gn};
    }
    GSC_DELETE_BUFFER(existing_groups);
}

/** Function to identify the next sequential integer that is not
 *  already allocated to a label in the simulation.
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @return the next sequential currently-unused label id, an integer
 * greater than 0.
 */
gsc_LabelID gsc_get_new_label_id( const gsc_SimData* d ) {
    // label_ids must be in sequential order
    gsc_LabelID new = {.id=1};
    int i = 0;

    while (i < d->n_labels) {
        if (new.id < d->label_ids[i].id) {
            break;
        }

        ++i;
        ++(new.id);
    }

    return new;
}

/** Function to identify the next sequential integer that is not
 *  already allocated to a marker effect set ID in the simulation.
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @return the next sequential currently-unused marker effect set id, an integer
 * greater than 0.
 */
gsc_EffectID gsc_get_new_eff_set_id( const gsc_SimData* d ) {
    // label_ids must be in sequential order
    gsc_EffectID new = { .id=1 };
    int i = 0;

    while (i < d->n_eff_sets) {
        if (new.id < d->eff_set_ids[i].id) {
            break;
        }

        ++i;
        ++(new.id);
    }

    return new;
}

/** Function to identify the next sequential integer that is not
 *  already allocated to a map ID in the simulation
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @return the next sequential currently-unused recombination map id, an integer
 * greater than 0.
 */
gsc_MapID gsc_get_new_map_id( const gsc_SimData* d) {
    // map IDs must be in sequential order
    gsc_MapID new = { .id=1 };
    int i = 0;

    while (i < d->genome.n_maps) {
        if (new.id < d->genome.map_ids[i].id) {
            break;
        }

        ++i;
        ++(new.id);
    }

    return new;
}


/** Function to identify the label lookup index of a label identifier.
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param label a label id
 * @return the index in d->label_ids, d->label_defaults, and the
 * ->labels table in gsc_AlleleMatrix where the data for this label
 * is stored, or -1 (GSC_UNINIT)
 * if the label with that id could not be found.
 */
int gsc_get_index_of_label( const gsc_SimData* d, const gsc_LabelID label ) {
    if (d->n_labels <= 0) { return GSC_UNINIT; } // immediate fail
    if (d->n_labels == 1) { return (d->label_ids[0].id == label.id) ? 0 : GSC_UNINIT ; }

    // If there's at least two labels then we binary search.
    int first = 0;
    int last = d->n_labels;
    int mid;

    while (first <= last) {
        mid = (first + last) / 2;

        if (d->label_ids[mid].id == label.id) {
            return mid;
        } else if (d->label_ids[mid].id < label.id) {
            first = mid + 1;
        } else {
            last = mid - 1;
        }

    }

    return GSC_UNINIT;
}

/** Function to identify the lookup index of a marker effect set identifier.
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param eff_set_id a marker effect set id
 * @return the index in d->e where the data for this effect set
 * is stored, or -1 (GSC_UNINIT)
 * if the effect set with that id could not be found.
 */
int gsc_get_index_of_eff_set( const gsc_SimData* d, const gsc_EffectID eff_set_id ) {
    if (d->n_eff_sets <= 0) { return GSC_UNINIT; } // immediate fail
    if (d->n_eff_sets == 1) { return (d->eff_set_ids[0].id == eff_set_id.id) ? 0 : GSC_UNINIT ; }

    // If there's at least two labels then we binary search.
    int first = 0;
    int last = d->n_eff_sets;
    int mid;

    while (first <= last) {
        mid = (first + last) / 2;

        if (d->eff_set_ids[mid].id == eff_set_id.id) {
            return mid;
        } else if (d->eff_set_ids[mid].id < eff_set_id.id) {
            first = mid + 1;
        } else {
            last = mid - 1;
        }

    }

    return GSC_UNINIT;
}

/** Function to identify the lookup index of a recombination map identifier.
 *
 * @param d the simulation containing the map
 * @param map a map id
 * @return the index in g->maps where the information for this map
 * is stored, or -1 (GSC_UNINIT)
 * if the map with that id could not be found.
 */
int gsc_get_index_of_map( const SimData* d, const gsc_MapID map ) {
    if (d->genome.n_maps <= 0) { return GSC_UNINIT; } // immediate fail
    if (d->genome.n_maps == 1) { return (d->genome.map_ids[0].id == map.id) ? 0 : GSC_UNINIT ; }

    // If there's at least two labels then we binary search.
    int first = 0;
    int last = d->genome.n_maps;
    int mid;

    while (first <= last) {
        mid = (first + last) / 2;

        if (d->genome.map_ids[mid].id == map.id) {
            return mid;
        } else if (d->genome.map_ids[mid].id < map.id) {
            first = mid + 1;
        } else {
            last = mid - 1;
        }

    }

    return GSC_UNINIT;
}

//-----------------------------------Data Access-----------------------------------

/** Function to count the number of genotypes that currently belong to
 * the specified group.
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * This goes through and checks every genotype in the gsc_SimData, because
 * there is currently no centralised tracking of groups.
 *
 * @shortnamed{get_group_size}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group whose members need to
 * be counted.
 * @returns the number of genotypes currently belonging to this group.
 */
int gsc_get_group_size( const gsc_SimData* d, const gsc_GroupNum group_id) {
    if (group_id.num == GSC_NO_GROUP.num) {
        return 0; // it is not a group so it does not have a size
    }
    const gsc_AlleleMatrix* m = d->m;
	int size = 0;
	int i;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
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
 * @see gsc_get_group_size()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_genes}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group we want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a char** vector with length at least large enough to fit the whole group.
 * The function will fill this vector with pointers to the allele strings of each member
 * of the group. The vector's contents are only
 * shallow copies that should not be freed.
 * @returns The number of entries of `output` that have been filled. Equal to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_genes( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output) {
    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
	}
	int i, genes_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
                output[genes_i] = m->alleles[i];
				++genes_i;
			}
		}

		if (m->next == NULL) {
            return group_size;
		} else {
			m = m->next;
		}
	}
}

/** Gets a shallow copy of the names of each member of the group.
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_names}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a char** vector with length at least large enough to fit the whole group.
 * The function will fill this vector with pointers to the names of each member
 * of the group. The vector's contents are only
 * shallow copies that should not be freed.
 * @returns The number of entries of `output` that have been filled. Equal to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_names( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output) {
    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }
	int i, names_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
                output[names_i] = m->names[i];
				++names_i;
			}
		}

		if (m->next == NULL) {
            return group_size;
		} else {
			m = m->next;
		}
	}
}

/** Gets the ids of each member of the group.
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_ids}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a vector with length at least large enough to fit the whole group.
 * The function will fill this vector with the ids of each member of the group.
 * @returns The number of entries of `output` that have been filled. Equal to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_ids( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, gsc_PedigreeID *output) {
    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
                output[ids_i] = m->ids[i];
				++ids_i;
			}
		}

		if (m->next == NULL) {
            return group_size;
		} else {
			m = m->next;
		}
	}
}

/** Gets the indexes (0-based, from the start of the linked list in the gsc_SimData)
 * of each member of the group.
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_indexes}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a vector with length at least large enough to fit the whole group.
 * The function will fill this vector with the indexes of each member of the group.
 * @returns The number of entries of `output` that have been filled. Equal to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_indexes(const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, unsigned int* output) {
    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }
	int i, total_i = 0, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i, ++total_i) {
            if (m->groups[i].num == group_id.num) {
                output[ids_i] = total_i;
				++ids_i;
			}
		}

		if (m->next == NULL) {
            return group_size;
		} else {
			m = m->next;
		}
	}
}

/** Gets the breeding values/breeding values/fitnesses of each member of the group
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_bvs}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a vector with length at least large enough to fit the whole group.
 * The function will fill this vector with the breeding values of each member of the group.
 * @returns The number of entries of `output` that have been filled. Equal to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_bvs( const gsc_SimData* d, const gsc_GroupNum group_id, const gsc_EffectID effID, int group_size, double* output) {
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }

    gsc_DecimalMatrix dm_bvs = gsc_calculate_group_bvs(d, group_id, effID );

	for (int i = 0; i < dm_bvs.cols; ++i) {
        output[i] = dm_bvs.matrix[0][i];
	}

	gsc_delete_dmatrix(&dm_bvs);

    return group_size;
}

/** Gets the ids of either the first or second parent of each member
 * of the group.
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_names()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_parent_ids}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param whichParent 1 to get the first parent of each group member, 2 to get the second.
 * Raises an error and returns NULL if this parameter is not either of those values.
 * @param output a vector with length at least large enough to fit the whole group.
 * The function will fill this vector with the ids the chosen parent of each member of the group.
 * @returns GSC_UNINIT if `parent`'s value is incorrect; otherwise the number of entries of `output`
 * that have been filled. This is to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_parent_ids( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, const int whichParent, gsc_PedigreeID* output) {
    if (!(whichParent == 1 || whichParent == 2)) {
		warning( "Value error: `parent` must be 1 or 2.");
        return GSC_UNINIT;
	}
    int parent = whichParent - 1;

    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
                output[ids_i] = m->pedigrees[parent][i];
				++ids_i;
			}
		}

		if (m->next == NULL) {
            return group_size;
		} else {
			m = m->next;
		}
	}
}

/** Gets the names of either the first or second parent of each member
 * of the group. Names may be NULL.
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_pedigrees()
 *
 * @shortnamed{get_group_parent_names}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param whichParent 1 to get the first parent of each group member, 2 to get the second.
 * Raises an error and returns NULL if this parameter is not either of those values.
 * @param output a char** vector with length at least large enough to fit the whole group.
 * The function will fill this vector with pointers to the names of the chosen parent
 * of each member of the group. The vector's contents are only
 * shallow copies that should not be freed.
 * @returns GSC_UNINIT if `parent`'s value is incorrect; otherwise the number of entries of `output`
 * that have been filled. This is to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_parent_names( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, const int whichParent, char** output) {
    if (!(whichParent == 1 || whichParent == 2)) {
        warning( "Value error: `parent` must be 1 or 2.");
        return GSC_UNINIT;
    }
    int parent = whichParent - 1;

    const gsc_AlleleMatrix* m = d->m;
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }

    int i, ids_i = 0;
    while (1) {
        for (i = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
                if (m->pedigrees[parent][i].id != GSC_NO_PEDIGREE.id) {
                    output[ids_i] = gsc_get_name_of_id(d->m, m->pedigrees[parent][i]);
                } else {
                    output[ids_i] = NULL;
                }
                ++ids_i;
            }
        }

        if (m->next == NULL) {
            return group_size;
        } else {
            m = m->next;
        }
    }
}

/** Gets the full pedigree string (as per gsc_save_group_full_pedigree() )
 * of each member of the group.
 *
 * This function is not fast, or particularly memory-efficient. To avoid recursive
 * string-building/concatenation because that's pretty tricky in C, this just
 * calls gsc_save_group_full_pedigree() to a temporary file, then reads that file
 * back in to make a list of strings. Contact the package maintainers if you'd
 * benefit from a faster version of this function, otherwise I probably won't bother
 * improving this.
 *
 * If you just want to have the full pedigrees for further analysis, consider
 * @ref gsc_save_group_full_pedigree(). With the current implementation that one will
 * be faster than this anyway.
 *
 * @shortnamed{get_group_pedigrees}
 *
 * @see gsc_get_group_size()
 * @see gsc_get_group_genes()
 * @see gsc_get_group_names()
 * @see gsc_get_group_ids()
 * @see gsc_get_group_indexes()
 * @see gsc_get_group_bvs()
 * @see gsc_get_group_parent_ids()
 * @see gsc_get_group_parent_names()
 * @see gsc_save_group_full_pedigree()
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1 (GSC_UNINIT). This enables fewer calls to gsc_get_group_size when using multiple
 * group-data-getter functions.
 * @param output a char** vector with length at least large enough to fit the whole group.
 * The function will fill this vector with pointers to strings containing the full pedigree
 * of each member of the group. The pedigree strings are dynamically
 * allocated on the heap, and so should be freed once finished using them. This is because,
 * unlike other data-getter functions, full pedigree data is not stored in the simulation
 * but must be generated specifically to answer this function call.
 * @returns GSC_UNINIT if there are no genotypes in the group or if reading/writing
 * to the temporary file failed; 0 otherwise.
 * @returns GSC_UNINIT if reading/writing to the temporary file failed;
 * otherwise the number of entries of `output` that have been filled. This is to `group_size`
 * if group size was set; equal to the actual size of the group if `group_size` was -1.
 */
int gsc_get_group_pedigrees( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output) {
	char* fname = "gS_gpptmp";
	FILE* fp = fopen(fname, "w");
	gsc_save_group_full_pedigree(fp, d, group_id);
	fclose(fp);

	FILE* fp2;
	if ((fp2 = fopen(fname, "r")) == NULL) {
		warning( "Failed to use temporary file.\n");
        return GSC_UNINIT;
	}

	// Create the list that we will return
    if (group_size <= 0) { // group_size == GSC_UNINIT
        group_size = gsc_get_group_size( d, group_id );
        if (group_size == 0) { return 0; }
    }

	// read one line at a time
	//unsigned int n;
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
        output[i] = gsc_malloc_wrap(sizeof(char) * size,GSC_TRUE);
        while ((nextc = fgetc(fp2)) != '\n' && nextc != EOF) {
            output[i][index] = nextc;
			++index;

			if (index >= size) {
				size *= 2;
                char* temp = realloc(output[i], sizeof(char) * size);
				if (temp == NULL) {
                    GSC_FREE(output[i]);
					warning( "Memory allocation of size %u failed.\n", size);
                    output[i] = NULL;
				} else {
                    output[i] = temp;
				}
			}
		}
        output[i][index] = '\0';
	}

	fclose(fp2);
	remove(fname);

    return group_size;
}

/*---------------------- matrix-operations.c dregs -------------------*/

/** Generates a matrix of c columns, r rows with all 0.
 *
 * @param r the number of rows/first index for the new matrix.
 * @param c the number of columns/second index for the new matrix.
 * @returns a gsc_DecimalMatrix with r rows, c cols, and a matrix of
 * the correct size, with all values zeroed.
 */
gsc_DecimalMatrix gsc_generate_zero_dmatrix(const int r, const int c) {
	gsc_DecimalMatrix zeros;
	zeros.rows = r;
	zeros.cols = c;

    zeros.matrix = gsc_malloc_wrap(sizeof(double*) * r,GSC_TRUE);
	for (int i = 0; i < r; ++i) {
        zeros.matrix[i] = gsc_malloc_wrap(sizeof(double) * c,GSC_TRUE);
		for (int j = 0; j < c; ++j) {
			zeros.matrix[i][j] = 0.0;
		}
	}
	return zeros;
}

/** Multiply a gsc_DecimalMatrix to a vector, and add that product to the first column
 * of a provided gsc_DecimalMatrix.
 *
 * Performs the fused multiply-add operation: `result + a * b` and saves it to `result`.
 *
 * Assumes that the vector `b` has the same number of entries as `a` has columns. That
 * is, assumes the dimensions are valid for multiplication.
 *
 * @param result pointer to the gsc_DecimalMatrix to whose first row the product a*b
 * will be added.
 * @param a pointer to the gsc_DecimalMatrix. Multiply this to b.
 * @param b a vector of doubles, assumed to be the same length as `a` has columns.
 * Multiplied to a.
 * @returns 0 on success, nonzero on failure.
 */
int gsc_add_matrixvector_product_to_dmatrix(gsc_DecimalMatrix* result, const gsc_DecimalMatrix* a, const double* b) {
    /*if (a->cols != b->rows) {
		error( "Dimensions invalid for matrix multiplication.");
	}*/

	//if (result->rows != a->rows || result->cols != b->cols) {
	if (result->cols != a->rows) { //these dimensions make it so result is one heap array, using only first row.
        warning( "Dimensions invalid for adding to result: %d does not fit in %d\n", a->rows, result->cols);
		return 1;
	}

	double cell;

	for (int i = 0; i < result->cols; ++i) {
        cell = 0;
        for (int j = 0; j < a->cols; ++j) {
            // for each cell, we loop through each of the pairs adjacent to it.
            cell += (a->matrix[i][j]) * b[j];
        }

        result->matrix[0][i] += cell;
	}

	return 0;

}

/** Multiply two sets of a gsc_DecimalMatrix and vector, and add both products to
 * the first column of a provided gsc_DecimalMatrix.
 *
 * Performs the double fused multiply-add operation:
 * `result + amat * avec + bmat * bvec`
 * and saves it to `result`.
 *
 * The matrices `amat` and `bmat` must have the same dimensions, and must have
 * the same number of rows as `result` has columns.
 *
 * Assumes that the vectors have the same number of entries as their matrices have
 * columns. That is, assumes the dimensions are valid for multiplication.
 *
 * @param result pointer to the gsc_DecimalMatrix to whose first row both products
 * will be added.
 * @param amat pointer to the first gsc_DecimalMatrix. Multiply this to avec.
 * @param avec a vector of doubles, assumed to be the same length as `amat` has columns.
 * Multiplied to `amat`.
 * @param bmat pointer to the second gsc_DecimalMatrix. Multiply this to bvec.
 * @param bvec a second vector of doubles, assumed to be the same length as
 * `bmat` has columns. Multiplied to `bmat`.
 * @returns 0 on success, nonzero on failure.
 */
int gsc_add_doublematrixvector_product_to_dmatrix(gsc_DecimalMatrix* result, const gsc_DecimalMatrix* amat, const double* avec,
                                              const gsc_DecimalMatrix* bmat, const double* bvec) {

	if (result->cols != amat->rows) { //these dimensions make it so result is one heap array, using only first row.
        warning( "Dimensions invalid for adding to result: %d does not fit in %d\n", amat->rows, result->cols);
		return 1;
	}
	if (result->cols != bmat->rows) {
        warning( "Dimensions invalid for adding to result: %d does not fit in %d\n", bmat->rows, result->cols);
		return 1;
	}
	if (amat->cols != bmat->cols) {
        warning( "Dimensions of the two products are uneven: length %d does not match length %d\n", amat->cols, bmat->cols);
		return 1;
	}

	double cell;

	for (int i = 0; i < result->cols; ++i) {
        cell = 0;
        for (int j = 0; j < amat->cols; ++j) {
            // for each cell, we loop through each of the pairs adjacent to it.
            cell += (amat->matrix[i][j]) * avec[j];
            cell += (bmat->matrix[i][j]) * bvec[j];
        }

        result->matrix[0][i] += cell;
	}

	return 0;
}


/** Deletes a gsc_DecimalMatrix and frees its memory. m will now refer
 * to an empty matrix, with every pointer set to null and dimensions set to 0.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void gsc_delete_dmatrix(gsc_DecimalMatrix* m) {
	if (m->matrix != NULL) {
		for (int i = 0; i < m->rows; i++) {
			if (m->matrix[i] != NULL) {
				GSC_FREE(m->matrix[i]);
			}
		}
		GSC_FREE(m->matrix);
		m->matrix = NULL;
	}
	m->cols = 0;
	m->rows = 0;
}


/*--------------------------------Deleting-----------------------------------*/

/** Deletes all genotypes belonging to a particular group.
 *
 *  This includes all of their details. Persistent ids (used to track
 *  pedigree) will not be re-used.
 *
 * Uses a call to gsc_condense_allele_matrix() to ensure that the gsc_SimData
 * remains valid after deletion.
 *
 * @shortnamed{delete_group}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param group_id the group number of the subset of data to be cleared
 */
void gsc_delete_group(gsc_SimData* d, const gsc_GroupNum group_id) {
	gsc_AlleleMatrix* m = d->m;
	int i, deleted, total_deleted = 0;
	while (1) {

		for (i = 0, deleted = 0; i < m->n_genotypes; ++i) {
            if (m->groups[i].num == group_id.num) {
				// delete data
				if (m->names[i] != NULL) {
					GSC_FREE(m->names[i]);
					m->names[i] = NULL;
				}
				if (m->alleles[i] != NULL) {
					GSC_FREE(m->alleles[i]);
					m->alleles[i] = NULL;
				}
                m->ids[i] = GSC_NO_PEDIGREE;
                m->pedigrees[0][i] = GSC_NO_PEDIGREE;
                m->pedigrees[1][i] = GSC_NO_PEDIGREE;
                m->groups[i] = GSC_NO_GROUP;
				++deleted;
			}
		}
		m->n_genotypes -= deleted;
		total_deleted += deleted;

		if (m->next == NULL) {
			gsc_condense_allele_matrix( d );
			Rprintf("%d genotypes were deleted\n", total_deleted);
            d->n_groups--;
			return;
		} else {
			m = m->next;
		}
	}
}

/** Deletes a particular set of marker effects from memory
*
 * @shortnamed{delete_eff_set}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param effID the identifier of the effect set to be deleted
 */
void gsc_delete_eff_set(gsc_SimData* d, gsc_EffectID effID) {
    int which_ix = gsc_get_index_of_eff_set(d, effID);
    if (which_ix == GSC_UNINIT) {
        warning( "Nonexistent effect set %d\n", effID.id);
        return;
    }

    if (d->n_eff_sets == 1) {
        gsc_delete_effect_matrix(d->e);
        d->n_eff_sets = 0;
        GSC_FREE(d->e);
        GSC_FREE(d->eff_set_ids);
        d->e = NULL;
        d->eff_set_ids = NULL;
    } else {
        d->n_eff_sets--;

        gsc_delete_effect_matrix(d->e + which_ix);
        gsc_EffectMatrix* newE = gsc_malloc_wrap(sizeof(*d->e)*d->n_eff_sets,GSC_FALSE);
        if (newE == NULL) {
            gsc_EffectMatrix cleared = d->e[which_ix];
            for (int i = which_ix; i < d->n_eff_sets-1; ++i) {
                d->e[i] = d->e[i+1];
            }
            d->e[d->n_eff_sets] = cleared;
        } else {
            memcpy(newE, d->e, sizeof(*d->e)*which_ix);
            memcpy(newE + which_ix, d->e + which_ix + 1, sizeof(*d->e)*(d->n_eff_sets - which_ix));
            GSC_FREE(d->e);
            d->e = newE;
        }

        gsc_EffectID* newIDs = gsc_malloc_wrap(sizeof(*d->eff_set_ids)*d->n_eff_sets,GSC_FALSE);
        if (newIDs == NULL) {
            for (int i = which_ix; i < d->n_eff_sets-1; ++i) {
                d->eff_set_ids[i] = d->eff_set_ids[i+1];
            }
            d->eff_set_ids[d->n_eff_sets] = NO_EFFECTSET;
        } else {
            memcpy(newIDs, d->eff_set_ids, sizeof(*d->eff_set_ids)*which_ix);
            memcpy(newIDs + which_ix, d->eff_set_ids + which_ix + 1, sizeof(*d->eff_set_ids)*(d->n_eff_sets - which_ix));
            GSC_FREE(d->eff_set_ids);
            d->eff_set_ids = newIDs;
        }
    }
}

/** Clears memory of this label from the simulation and all its genotypes.
 *
 * @shortnamed{delete_label}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param which_label the label id of the label to be destroyed
 */
void gsc_delete_label(gsc_SimData* d, const gsc_LabelID which_label) {
    int label_ix;
    if (which_label.id == GSC_NO_LABEL.id || (label_ix = gsc_get_index_of_label(d, which_label)) < 0) {
        warning( "Nonexistent label %d\n", which_label.id);
        return;
    }

    if (d->n_labels == 1) {
        // Delete 'em all
        d->n_labels = 0;
        GSC_FREE(d->label_ids);
        d->label_ids = NULL;
        GSC_FREE(d->label_defaults);
        d->label_defaults = NULL;

        gsc_AlleleMatrix* m = d->m;
        do {

            GSC_FREE(m->labels[0]);
            GSC_FREE(m->labels);
            m->labels = NULL;

        } while ((m = m->next) != NULL);

    } else {
        // Reduce the list of labels in the gsc_SimData
        gsc_LabelID* new_label_ids = gsc_malloc_wrap(sizeof(gsc_LabelID) * (d->n_labels - 1),GSC_FALSE);
        if (new_label_ids == NULL) {
            for (int i = label_ix; i < d->n_labels - 1; ++i) {
                d->label_ids[i] = d->label_ids[i+1];
            }
            d->label_ids[d->n_labels] = GSC_NO_LABEL;
        } else {
            memcpy(new_label_ids,d->label_ids,sizeof(*d->label_ids)*label_ix);
            memcpy(new_label_ids + label_ix,d->label_ids + label_ix + 1, sizeof(*d->label_ids)*(d->n_labels - 1 - label_ix));
            GSC_FREE(d->label_ids);
            d->label_ids = new_label_ids;
        }

        int* new_label_defaults = gsc_malloc_wrap(sizeof(int) * (d->n_labels - 1),GSC_FALSE);
        if (new_label_defaults == NULL) {
            for (int i = label_ix; i < d->n_labels - 1; ++i) {
                d->label_defaults[i] = d->label_defaults[i+1];
            }
            // no need to overwrite default
        } else {
            memcpy(new_label_defaults,d->label_defaults,sizeof(*d->label_defaults)*label_ix);
            memcpy(new_label_defaults + label_ix,d->label_defaults + label_ix + 1, sizeof(*d->label_defaults)*(d->n_labels - 1 - label_ix));
            GSC_FREE(d->label_defaults);
            d->label_defaults = new_label_defaults;
        }
        d->n_labels --;

        // Remove the label from the gsc_AlleleMatrix linked list
        gsc_AlleleMatrix* m = d->m;
        do {
            GSC_FREE(m->labels[label_ix]);

            m->n_labels = d->n_labels;
            int** new_label_lookups = gsc_malloc_wrap(sizeof(int*) * (m->n_labels),GSC_FALSE);
            if (new_label_lookups == NULL) {
                for (int i = label_ix; i < m->n_labels; ++i) {
                    m->labels[i] = m->labels[i+1];
                }
                m->labels[m->n_labels + 1] = NULL;
            } else {
                memcpy(new_label_lookups, m->labels, sizeof(*m->labels)*label_ix);
                memcpy(new_label_lookups + label_ix, m->labels + label_ix + 1, sizeof(*m->labels)*(m->n_labels - label_ix));
                GSC_FREE(m->labels);
                m->labels = new_label_lookups;
            }
        } while ((m = m->next) != NULL);
    }
}

/** Deletes and clears the memory of a gsc_KnownGenome object and its children.
 *
 * Children include all recombination maps associated with this KnownGenome object
 *
 * @param g pointer to the structure whose data is to be cleared and memory freed.
 */
void gsc_delete_genome(gsc_KnownGenome* g) {
    if (g->marker_names != NULL) {
        for (int i = 0; i < g->n_markers; i++) {
            if (g->marker_names[i] != NULL) {
                GSC_FREE(g->marker_names[i]);
            }
        }
        GSC_FREE(g->marker_names);
        g->marker_names = NULL;
    }
    if (g->names_alphabetical != NULL) {
        GSC_FREE(g->names_alphabetical);
        g->names_alphabetical = NULL;
    }
    g->n_markers = 0;

    if (g->map_ids != NULL) {
        GSC_FREE(g->map_ids);
        g->map_ids = NULL;
    }

    if (g->maps != NULL) {
        for (int i = 0; i < g->n_maps; ++i) {
            gsc_delete_recombination_map_nointegrity(g->maps + i);
        }
        GSC_FREE(g->maps);
        g->maps = NULL;
    }
    g->n_maps = 0;
}

/** Deletes a particular recombination map from memory
*
 * @shortnamed{delete_recombination_map}
 *
 * @param d the gsc_SimData struct on which to perform the operation
 * @param which_map the ID of the recombination map to be deleted.
 */
void gsc_delete_recombination_map(gsc_SimData* d, const gsc_MapID which_map) {
    int map_ix;
    if (which_map.id == GSC_NO_LABEL.id || (map_ix = gsc_get_index_of_map(d, which_map)) < 0) {
        warning( "Nonexistent recombination map %d\n", which_map.id);
        return;
    }

    if (d->genome.n_maps == 1) {
        GSC_FREE(d->genome.map_ids);
        d->genome.map_ids = NULL;
        gsc_delete_recombination_map_nointegrity(&d->genome.maps[0]);
        GSC_FREE(d->genome.maps);
        d->genome.maps = NULL;
        d->genome.n_maps = 0;
    } else {
        d->genome.n_maps--;
        gsc_delete_recombination_map_nointegrity(&d->genome.maps[map_ix]);
        gsc_RecombinationMap* tmplist = gsc_malloc_wrap(sizeof(*d->genome.maps)*d->genome.n_maps, GSC_FALSE);
        if (tmplist == NULL) {
            gsc_RecombinationMap clearedmap = d->genome.maps[map_ix];
            for (int i = map_ix; i < d->genome.n_maps - 1; ++i) {
                d->genome.maps[i] = d->genome.maps[i+1];
            }
            d->genome.maps[d->genome.n_maps] = clearedmap;
        } else {
            memcpy(tmplist, d->genome.maps, sizeof(*d->genome.maps)*map_ix);
            memcpy(tmplist + map_ix, d->genome.maps + map_ix + 1, sizeof(*d->genome.maps)*(d->genome.n_maps - map_ix));
            GSC_FREE(d->genome.maps);
            d->genome.maps = tmplist;
        }

        gsc_MapID* tmpids = gsc_malloc_wrap(sizeof(*d->genome.map_ids)*d->genome.n_maps, GSC_FALSE);
        if (tmpids == NULL) {
            for (int i = map_ix; i < d->genome.n_maps - 1; ++i) {
                d->genome.map_ids[i] = d->genome.map_ids[i+1];
            }
            d->genome.map_ids[d->genome.n_maps] = NO_MAP;
        } else {
            memcpy(tmpids, d->genome.map_ids, sizeof(*d->genome.map_ids)*map_ix);
            memcpy(tmpids + map_ix, d->genome.map_ids + map_ix + 1, sizeof(*d->genome.map_ids)*(d->genome.n_maps - map_ix));
            GSC_FREE(d->genome.map_ids);
            d->genome.map_ids = tmpids;
        }
    }
}

/** Deletes and clears the memory of a gsc_RecombinationMap struct.
 *
 * "No integrity" because it does not check that the corresponding gsc_MapID is
 * correspondingly removed from simulation memory.
 *
 * For a safer/user-friendlier version of this function, @see gsc_delete_recombination_map
 *
 * @param m pointer to the structure whose data is to be cleared and memory freed.
 */
void gsc_delete_recombination_map_nointegrity(gsc_RecombinationMap* m) {
    if (m->chrs != NULL) {
        for (int i = 0; i < m->n_chr; ++i) {
            switch (m->chrs[i].type) {
            case GSC_LINKAGEGROUP_SIMPLE:
                m->chrs[i].map.simple.expected_n_crossovers = 0;
                m->chrs[i].map.simple.first_marker_index = 0;
                m->chrs[i].map.simple.n_markers = 0;
                if (m->chrs[i].map.simple.dists != NULL) {
                    GSC_FREE(m->chrs[i].map.simple.dists);
                    m->chrs[i].map.simple.dists = NULL;
                }
                break;
            case GSC_LINKAGEGROUP_REORDER:
                m->chrs[i].map.reorder.expected_n_crossovers = 0;
                m->chrs[i].map.reorder.n_markers = 0;
                if (m->chrs[i].map.reorder.dists != NULL) {
                    GSC_FREE(m->chrs[i].map.reorder.dists);
                    m->chrs[i].map.reorder.dists = NULL;
                }
                if (m->chrs[i].map.reorder.marker_indexes != NULL) {
                    GSC_FREE(m->chrs[i].map.reorder.marker_indexes);
                    m->chrs[i].map.reorder.marker_indexes = NULL;
                }
                break;
            }
        }
        GSC_FREE(m->chrs);
        m->chrs = NULL;
    }
}

/** Delete the gsc_AlleleMatrix linked list from m onwards and frees its memory.
 *
 * Freeing its memory includes freeing the gsc_AlleleMatrix, which was allocated on the heap
 * by @ref gsc_create_empty_allelematrix(). All matrices further along in the linked list chain
 * (eg pointed to by m->next or a chain of ->next pointers) will be similarly deleted.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void gsc_delete_allele_matrix(gsc_AlleleMatrix* m) {
	if (m == NULL) {
		return;
	}
	gsc_AlleleMatrix* next;
	do {
		/* free the big data matrix */
        for (int i = 0; i < CONTIG_WIDTH; i++) {
            if (m->alleles[i] != NULL) {
                GSC_FREE(m->alleles[i]);
            }

        }

		// free names
        for (int i = 0; i < CONTIG_WIDTH; i++) {
            if (m->names[i] != NULL) {
                GSC_FREE(m->names[i]);
            }
        }

        // free labels
        for (int i = 0; i < m->n_labels; ++i) {
            if (m->labels[i] != NULL) {
                GSC_FREE(m->labels[i]);
            }
        }
        GSC_FREE(m->labels);

		next = m->next;
		GSC_FREE(m);
	} while ((m = next) != NULL);
}

/** Deletes an gsc_EffectMatrix object and frees its memory.
 *
 * m will now refer
 * to an empty matrix, with every pointer set to null and dimensions set to 0.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void gsc_delete_effect_matrix(gsc_EffectMatrix* m) {
	gsc_delete_dmatrix(&(m->effects));
	if (m->effect_names != NULL) {
		GSC_FREE(m->effect_names);
	}
	m->effect_names = NULL;
}

/** Deletes a gsc_SimData object and frees its memory.
 *
 * @shortnamed{delete_simdata}
 *
 * @param m pointer to the struct whose data is to be cleared and memory freed.
 */
void gsc_delete_simdata(gsc_SimData* m) {
	if (m == NULL) {
		return;
	}

    gsc_delete_genome(&(m->genome));

    if (m->n_eff_sets > 0) {
        GSC_FREE(m->eff_set_ids);
        for (int i = 0; i < m->n_eff_sets; ++i) {
            gsc_delete_effect_matrix(&(m->e[i]));
        }
        GSC_FREE(m->e);
    }

	gsc_delete_allele_matrix(m->m);

    if (m->n_labels > 0) {
        if (m->label_ids != NULL) {
            GSC_FREE(m->label_ids);
        }
        if (m->label_defaults != NULL) {
            GSC_FREE(m->label_defaults);
        }
    }

	if (m != NULL) {
		GSC_FREE(m);
	}
}

/** Delete a gsc_MarkerBlocks struct.
 *
 * Deletes a gsc_MarkerBlocks object and frees its associated memory. b will now refer
 * to an empty struct, with every pointer set to null and number of markers set to 0.
 *
 * @shortnamed{delete_markerblocks}
 *
 * @param b pointer to the struct whose data is to be cleared and memory freed.
 */
void gsc_delete_markerblocks(gsc_MarkerBlocks* b) {
	for (int i = 0; i < b->num_blocks; ++i) {
		GSC_FREE(b->markers_in_block[i]);
	}
	GSC_FREE(b->markers_in_block);
	b->markers_in_block = NULL;
	GSC_FREE(b->num_markers_in_block);
	b->num_markers_in_block = NULL;
	b->num_blocks = 0;

	return;
}


/** Deletes a gsc_BidirectionalIterator object.
 *
 *  A gsc_BidirectionalIterator has no heap memory to free, so calling
 *  this function is mostly unnecessary. The function will set all
 *  values in the struct to uninitialised/null values, and will set
 *  the iterator to think it is both at the start and end of its
 *  sequence, so that any next_* functions will not attempt to
 *  search for group member locations.
 *
 * @shortnamed{delete_bidirectional_iter}
 *
 * @param it pointer to the struct whose data is to be cleared.
 */
void gsc_delete_bidirectional_iter(gsc_BidirectionalIterator* it) {
    it->d = NULL;
    //it->group = GSC_NO_GROUP;
    it->localPos = GSC_UNINIT;
    it->cachedAM = NULL;
    it->cachedAMIndex = GSC_UNINIT;
    it->atEnd = GSC_TRUE;
    it->atStart = GSC_TRUE;
}

/** Deletes a gsc_RandomAccessIterator object and frees its memory.
 *
 *  All values in the struct will be set to uninitialised/null values,
 *  except for groupSize, which will be set to 0 to so that any
 *  next_* functions called on the iterator will not attempt to
 *  search for group member locations.
 *
 * @shortnamed{delete_randomaccess_iter}
 *
 * @param it pointer to the struct whose data is to be cleared and memory freed.
 */
void gsc_delete_randomaccess_iter(gsc_RandomAccessIterator* it) {
    it->d = NULL;
    //it->group = GSC_NO_GROUP;
    if (it->cacheSize > 0) {
        GSC_FREE(it->cache);
    }
    it->cache = NULL;
    it->cacheSize = 0;
    it->largestCached = GSC_UNINIT;
    it->groupSize = 0;
}

/*-------------------------------gsc_SimData loaders-----------------------------*/

/** Open a file for reading with gsc_TableFileReader
 *
 * On successfully opening file, fills the TableFileReader buffer for the first time.
 */
gsc_TableFileReader gsc_tablefilereader_create(const char* filename) {
    FILE* fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        warning( "Failed to open file %s.\n", filename);
    }

    gsc_TableFileReader tfr = { .fp = fp,
                                .buf = { 0 },
                                .buf_fill = 0,
                                .cursor = 0,
                              };

    if (fp != NULL) {
        tfr.buf_fill = fread(tfr.buf,1,sizeof(tfr.buf),fp);
    }
    return tfr;
}

/** Close a gsc_TableFileReader's file pointer.
 */
void gsc_tablefilereader_close(gsc_TableFileReader* tbl) {
    if (tbl->fp != NULL) { fclose(tbl->fp); }
    tbl->fp = NULL;
}

/** Read another buffer's worth of characters from a gsc_TableFileReader's file.
 *
 * @warning This overwrites any characters previously saved in the TableFileReader
 * buffer. The pointers of any @a gsc_TableFileCell read from this table may become invalid
 * if they are shallow copies. If you need to retain any cell values, consider using
 * @a gsc_tablefilecell_deep_copy() before calling this function.
 */
void gsc_helper_tablefilereader_refill_buffer(gsc_TableFileReader* tbl) {
    tbl->cursor = 0;
    if (tbl->fp != NULL) {
        tbl->buf_fill = fread(tbl->buf,1,sizeof(tbl->buf),tbl->fp);
    } else {
        tbl->buf_fill = 0;
    }
}

/** Classify the character under the cursor of a TableFileReader as cell contents or otherwise
 *
 * Does not update @a tbl->cursor, so repeated calls of this same function without updating the
 * cursor in between will return the same result.
 */
enum gsc_TableFileCurrentStatus gsc_helper_tablefilereader_classify_char(gsc_TableFileReader* tbl) {
    if (tbl->buf_fill <= tbl->cursor) {
        if (tbl->buf_fill < sizeof(tbl->buf)) { // last read did not fill the entire buffer
            return GSC_TABLEFILE_ERROR_EOF;
        }
        return GSC_TABLEFILE_ERROR_EOBUF;
    }

    switch (tbl->buf[tbl->cursor]) {
    case '\r': // allow '\r' or '\r\n' as end of lines. also allow '\n' as end of line (see following case)
    case '\n':
        return GSC_TABLEFILE_NEWLINE;
    case '\t':
    case ' ':
    case ',':
        return GSC_TABLEFILE_COLUMNGAP;
    default:
        return GSC_TABLEFILE_CONTENTS;
    }
}

/** Allocate memory to store a deep copy of a gsc_TableFileCell, if previously only a shallow copy.
 *
 *  The deep copy will be a null-terminated string even if the shallow copy was not null-terminated.
 *
 * After this call, the cell is stored in heap memory and will need to be freed once the
 * cell is no longer needed. Schema for doing this:
 * if (!mycell.isCellShallow) { GSC_FREE(mycell.cell); }
 */
void gsc_tablefilecell_deep_copy(gsc_TableFileCell* c) {
    if (c->cell_len > 0 && c->isCellShallow) {
        char* deepcell = gsc_malloc_wrap(sizeof(char)*(c->cell_len+1), GSC_TRUE);
        memcpy(deepcell,c->cell,sizeof(char)*c->cell_len);
        deepcell[c->cell_len] = '\0';
        c->cell = deepcell;
        c->isCellShallow = GSC_FALSE;
    }
}

/** Read forwards in TableFileReader and return the next cell's contents, as
 *  well as how many column gaps and newlines preceeded it.
 *
 * Cells can be of unlimited length, as long as they fit in memory.
 */
gsc_TableFileCell gsc_tablefilereader_get_next_cell(gsc_TableFileReader* tbl) {
    gsc_TableFileCell cur = { .isCellShallow = GSC_TRUE, .cell = NULL, .cell_len = 0,
                              .predCol = 0, .predNewline = 0, .eof = GSC_FALSE };

    GSC_CREATE_BUFFER(tmpcell,char,1);
    unsigned int tmpix = 0;
    unsigned int tblbuf_offset = 0;
    unsigned int tblbuf_len = 0;
    int predCarriageReturn = 0; // for detecting /r/n as a single "newline"
    int warned = GSC_FALSE;

    while (1) {
        enum gsc_TableFileCurrentStatus type = gsc_helper_tablefilereader_classify_char(tbl);
        if (0 < predCarriageReturn) { --predCarriageReturn; } // decremented each time step

        if (0 == cur.cell_len) {
            switch (type) {
            case GSC_TABLEFILE_NEWLINE:
                if (tbl->buf[tbl->cursor] == '\r') {
                    predCarriageReturn = 2; // will have value 1 at next loop iteration, then will fall back to 0
                }
                if (!(predCarriageReturn && tbl->buf[tbl->cursor] == '\n')) {
                    ++cur.predNewline;
                }
                cur.predCol = 0;
                ++tbl->cursor;
                break;

            case GSC_TABLEFILE_COLUMNGAP:
                ++tbl->cursor;
                ++cur.predCol;
                break;

            case GSC_TABLEFILE_ERROR_EOBUF:
                // just refill as we have no contents we need to save yet
                gsc_helper_tablefilereader_refill_buffer(tbl);
                if (0 < predCarriageReturn) { ++predCarriageReturn; } // should not tick down the counter this loop iteration
                break;

            case GSC_TABLEFILE_CONTENTS:
                tblbuf_offset = tbl->cursor; tblbuf_len = 1; // in case we need to make a deep copy later.
                cur.cell = tbl->buf + tbl->cursor;
                ++cur.cell_len;
                ++tbl->cursor;
                break;

            default:
                ++tbl->cursor;
                cur.eof = GSC_TRUE;
                return cur;
            }

        } else { // have found the cell, just need to read the rest of it
            switch (type) {
            case GSC_TABLEFILE_CONTENTS:
                ++tbl->cursor;
                ++tblbuf_len;
                ++cur.cell_len;
                break;

            case GSC_TABLEFILE_ERROR_EOBUF:
                cur.isCellShallow = GSC_FALSE;

                if (!warned && tblbuf_len > 8192) {
                    warned = GSC_TRUE;
                    fprintf(stderr,"Warning: very long cell identified beginning %c%c%c%c%c%c. Column separators may have failed to be recognised\n",
                            tmpcell[0],tmpcell[1],tmpcell[2],tmpcell[3],tmpcell[4],tmpcell[5]);
                }

                GSC_STRETCH_BUFFER(tmpcell,tmpix + tblbuf_len + 1);
                memcpy(tmpcell+tmpix,tbl->buf+tblbuf_offset,sizeof(char)*tblbuf_len);
                tmpix += tblbuf_len;
                tmpcell[tmpix] = '\0';

                tblbuf_offset = 0; tblbuf_len = 0;
                gsc_helper_tablefilereader_refill_buffer(tbl);
                break;

            case GSC_TABLEFILE_ERROR_EOF:
                ++tbl->cursor;
                cur.eof = GSC_TRUE;
                // fall through
            default: // newline or column gap or end of file discovered: save and return.
                if (!cur.isCellShallow) {
                    cur.cell = gsc_malloc_wrap(sizeof(char)*(cur.cell_len + 1),GSC_TRUE);
                    memcpy(cur.cell,tmpcell,sizeof(char)*tmpix);
                    if (0 < tblbuf_len) {
                        memcpy(cur.cell+tmpix,tbl->buf+tblbuf_offset,sizeof(char)*tblbuf_len);
                    }
                    cur.cell[cur.cell_len] = '\0';
                    GSC_DELETE_BUFFER(tmpcell);
                }
                return cur;
            }
        }
    }
}

/** Return whether or not a marker name is present in the tracked markers, and at what index
 *
 * @param target name of the marker that is to be located
 * @param g genome containing list of tracked markers to search within
 * @param out NULL if the output index is of no interest, or a pointer to a place to save a
 * the index of the located marker in the KnownGenome object on success otherwise.
 * @return GSC_TRUE if the marker was located and its index saved to @a outindex, GSC_FALSE if
 * the marker could not be located.
 */
int gsc_get_index_of_genetic_marker(const char* target, gsc_KnownGenome g, unsigned int* out) {
    unsigned int first = 0, last = g.n_markers - 1;
    int index = (first + last) / 2;
    int comparison = strcmp(target,*(g.names_alphabetical[index]));
    while (comparison != 0 && first <= last) {
        if (comparison == 0) {
            if (out != NULL) *out = g.names_alphabetical[index] - g.marker_names;
            return GSC_TRUE;
        } else if (comparison > 0) {
            first = index + 1;
            if (first >= g.n_markers) { return GSC_FALSE; }
        } else {
            if (index == 0) { return GSC_FALSE; }
            last = index - 1;
        }

        // index has been updated, no matter the branch.
        index = (first + last) / 2;
        comparison = strcmp(target, *(g.names_alphabetical[index]));
    }

    if (first > last) {
        return GSC_FALSE;
    }
    if (out != NULL) *out = g.names_alphabetical[index] - g.marker_names;
    return GSC_TRUE;
}

/** Return the next cell from a queue of cells until the queue is exhausted, and thereafter read new
 *  cells from a TableFileReader
 *
 * @param queue pointer to queue containing @a queuesize TableFileCells. The pointer will be updated to point
 * one cell forwards (so that it still points to the 'start' of the queue) if the item returned comes from
 * the queue.
 * @param queuesize number of entries remaining in the queue. Will be updated to be one fewer if the item returned
 * comes from the queue
 */
static gsc_TableFileCell gsc_helper_tablefilereader_get_next_cell_wqueue(gsc_TableFileReader* tf, gsc_TableFileCell** queue, unsigned int* queuesize) {
    gsc_TableFileCell ncell;
    if (*queuesize > 0) {
        ncell = *queue[0];
        /*for (int i = 1; i < *queuesize; ++i) { *queue[i-1] = *queue[i]; }*/
        ++*queue;
        --*queuesize;
    } else {
        ncell = gsc_tablefilereader_get_next_cell(tf);
    }
    return ncell;
}

/** Header row reading and processing for map and effect set files.
 *
 * Task 1: identifies if there is a header (in which case the returned queue
 * will not contain the header cells) or there is no header (in which case
 * the returned queue will contain the first row of cells).
 *
 * The first three cells are identified as a header row if there is a newline
 * after the third cell (and none between the first three cells), and the
 * entries in those three header cells correspond to the strings in the
 * @a canonical_titles parameter, in any order. If the first row does match the
 * requirements of a header row, the ordering of those @a canonical_titles
 * in the header row will be saved to @a col_order.
 *
 * The function returns the number of cells read but not processed. These cells read
 * but not processed (i.e. not any identified header cells) will be saved to
 * @a unprocessedqueue. At maximum 4 cells will be read, so @a unprocessedqueue
 * must be a pointer with a location for space for at least 4 TableFileCells.
 *
 * @param canonical_titles vector of three column titles, which the three potential
 * header cells will be matched against
 * @param col_order if the first row is determined to be a header row, entry i in this
 * vector will be filled with the column number (1, 2, or 3) that corresponds to the
 * i-th title in @a canonical_titles
 * @param unprocessedqueue vector of length at least 4 which will be overwritten with the
 * cells that were read but determined to not be part of the header.
 * @param queuesize space to save number of entries in the queue @a unprocessedqueue after
 * execution of this function
 * @returns GSC_TRUE (1) if the cell has a valid header, GSC_FALSE (0) if the cell has no header
 * and -1 if the cell does not have 3 columns or there is some other formatting error.
 */
static int gsc_helper_parse_3cell_header(gsc_TableFileReader* tf, const char** canonical_titles, int* col_order,
                                          gsc_TableFileCell* unprocessedqueue, unsigned int* queuesize) {
    const int ncells = 3;
    const int INVALID_COLS = -1;

    // assume unprocessedqueue has at least 4 spaces, titles has 3 entries and so does col_order
    unsigned int newest = 0;
    unsigned int onelinefile = GSC_FALSE;
    for (; newest < ncells; ++newest) {
        unprocessedqueue[newest] = gsc_tablefilereader_get_next_cell(tf);
        if ((unprocessedqueue[newest]).eof) {
            if (newest + 1 < ncells) {
                if (!((unprocessedqueue[newest]).isCellShallow)) { GSC_FREE((unprocessedqueue[newest]).cell); }
                *queuesize = newest; // newest does not exist in return group
                return INVALID_COLS;
            } else { // column 3 has eof at the end of it. That's okay.
                onelinefile = GSC_TRUE;
                *queuesize = newest + 1;
            }
        } else if ((unprocessedqueue[newest]).predNewline) {
            *queuesize = newest+1; // newest index included in return group
            return INVALID_COLS;
        }
        gsc_tablefilecell_deep_copy(unprocessedqueue + newest);
    }
     if (!onelinefile) {
        unprocessedqueue[newest] = gsc_tablefilereader_get_next_cell(tf); // four cell
        *queuesize = newest + 1;
     }

    // Check which ordering of the three titles it is.
    for (int i1 = 0; i1 < ncells; ++i1) {
        if (strcmp((unprocessedqueue[i1]).cell,canonical_titles[0]) == 0) {
            for (int inc = 1; inc < ncells; ++inc) {
                int i2 = (i1 + inc) % ncells;
                int i3 = (i1 + (ncells - inc)) % ncells;
                if (strcmp((unprocessedqueue[i2]).cell,canonical_titles[1]) == 0 &&
                        strcmp((unprocessedqueue[i3]).cell,canonical_titles[2]) == 0) {
                    col_order[0] = i1 + 1;
                    col_order[1] = i2 + 1;
                    col_order[2] = i3 + 1;
                    GSC_FREE((unprocessedqueue[0]).cell);
                    GSC_FREE((unprocessedqueue[1]).cell);
                    GSC_FREE((unprocessedqueue[2]).cell);
                    unprocessedqueue[0] = unprocessedqueue[3];
                    *queuesize = 1;
                    return GSC_TRUE;
                }
            }
        }
    }

    return GSC_FALSE;
}

/** Extract the contents of a genetic map file.
 *
 * By default, the file's columns are assumed to be in the order:
 * 1. marker name
 * 2. chromosome/linkage group (name or number)
 * 3. position in cM along that group
 * Extra columns are not permitted. If a header row is included, its fields
 * should be "marker", "chr" and "pos". These names can be reordered if the
 * file presents the columns in a different order to the default order
 * shown above.
 *
 * Marker names may contain any characters except tabs, spaces, commas, and newlines. Chromosome
 * names must be integers (eg 1) or alphanumeric strings (eg 1A). Positions should
 * be able to be parsed as a floating point number by the compiling computer's
 * C standard library (usually allows integers, decimals, and scientific notation,
 * but not fractions).
 *
 * Columns may be space-separated or tab-separated. Column separators need not
 * be consistent throughout the file, and may be repeated. That is, any consecutive
 * sequence of spaces and tabs is interpreted as one column separator. Any single
 * or consecutive pair of '\n' and '\r' characters are
 * inferred as representing a newline. Like with column separators, there is no requirement
 * for consistency across the file in which character or pair of characters represents
 * the end of a row.
 *
 * In the current implementation, it attempts to load the entire file into memory at a
 * time, so the system must have enough memory to hold the file if it is to run this function.
 *
 * @param filename name of the file to parse
 * @param out location to save a heap-allocated vector containing the contents of each line
 * of the file.
 * @return number of lines successfully read from the file
 */
static unsigned int gsc_helper_parse_mapfile(const char* filename, struct gsc_MapfileUnit** out) {
    if (filename == NULL) return 0;

    gsc_TableFileReader tf = gsc_tablefilereader_create(filename);

    unsigned int row = 1;
    unsigned int col = 1;

    gsc_TableFileCell cellsread[4] = { 0 };
    gsc_TableFileCell* cellqueue = cellsread;
    unsigned int queue_size;
    const char* titles[] = { "marker", "chr", "pos"};
    int colnums[] = { 1, 2, 3 };
    int header = gsc_helper_parse_3cell_header(&tf, titles, colnums, cellqueue, &queue_size);
    if (header == GSC_TRUE) {
        Rprintf("(Loading %s) Format: map file with header\n", filename);
    } else if (header == GSC_FALSE) {
        Rprintf("(Loading %s) Format: map file without header\n", filename);
    } else {
        Rprintf("(Loading %s) Failure: Cannot identify the expected 3 columns of the map file\n", filename);
        gsc_tablefilereader_close(&tf);
        return 0;
    }
    int marker_colnum = colnums[0], chr_colnum = colnums[1], pos_colnum = colnums[2];

    int goodrow = (header) ? GSC_FALSE : GSC_TRUE; // discard first row if it's a header, keep if it's not.
    unsigned int goodrow_counter = 0;

    char* marker = NULL;
    unsigned long chr = 0;
    double pos = 0;
    char* conversionflag;

    GSC_CREATE_BUFFER(buffer,struct gsc_MapfileUnit,CONTIG_WIDTH);

    gsc_TableFileCell ncell;
    do {
        ncell = gsc_helper_tablefilereader_get_next_cell_wqueue(&tf, &cellqueue, &queue_size);

        // Update row/col position and save predecessor row
        if (ncell.cell != NULL) {
            if (ncell.predNewline) {
                if (goodrow) { // save predecessor row
                    buffer[goodrow_counter].name = marker;
                    buffer[goodrow_counter].chr  = chr;
                    buffer[goodrow_counter].pos  = pos;

                    ++goodrow_counter;
                    if (goodrow_counter >= buffercap) {
                        GSC_STRETCH_BUFFER(buffer,2*row);
                    }
                    marker = NULL;
                } else if (marker != NULL) {
                    GSC_FREE(marker);
                }
                row += ncell.predNewline;
                goodrow = GSC_TRUE;
                col = 1;
            }
            col += (ncell.predCol > 0) ? GSC_TRUE : GSC_FALSE;

            // Parse this cell
            if (ncell.cell_len == 0) {
                goodrow = GSC_FALSE;
            } if (col == marker_colnum) {
                gsc_tablefilecell_deep_copy(&ncell);
                marker = ncell.cell;
                ncell.isCellShallow = GSC_TRUE; // so it isn't freed.

            } else if (col == chr_colnum) {
                char tmp = ncell.cell[ncell.cell_len]; ncell.cell[ncell.cell_len] = '\0';
                chr = strtoul(ncell.cell,&conversionflag,36);
                ncell.cell[ncell.cell_len] = tmp;
                if (conversionflag != ncell.cell + ncell.cell_len) { // unsuccessful read
                    //warning("Entry at row %i column %i of file %s could not be parsed as an integer or alphanumeric string\n", row, chr_colnum, filename);
                    goodrow = GSC_FALSE;
                }

            } else if (col == pos_colnum) {
                char tmp = ncell.cell[ncell.cell_len]; ncell.cell[ncell.cell_len] = '\0';
                pos = strtod(ncell.cell,&conversionflag);
                ncell.cell[ncell.cell_len] = tmp;
                if (conversionflag != ncell.cell + ncell.cell_len) { // unsuccessful read
                    goodrow = GSC_FALSE;
                    //warning("Entry at row %i column %i of file %s could not be parsed as a numeric value\n", row, pos_colnum, filename);
                }

            } else {
                goodrow = GSC_FALSE;
            }

            // Reset to get next cell.
            if (!ncell.isCellShallow) { GSC_FREE(ncell.cell); }
        }
    } while (!ncell.eof);

    if (col == 3) {
        if (goodrow) { // save predecessor row
            buffer[goodrow_counter].name = marker;
            buffer[goodrow_counter].chr  = chr;
            buffer[goodrow_counter].pos  = pos;

            ++goodrow_counter;
            marker = NULL;
        } else if (marker != NULL) {
            GSC_FREE(marker);
        }
    }
    //row -= ncell.predNewline; // don't count trailing newlines in stats.

    Rprintf("(Loading %s) %u marker(s) with map positions were loaded. Failed to parse %u line(s).\n", filename, (unsigned int) goodrow_counter, (unsigned int) (row - header - goodrow_counter));
    gsc_tablefilereader_close(&tf);

    // Check outputs. We don't delete the buffers because we want to leave them alive with our callers holding the handles.
    GSC_FINALISE_BUFFER(buffer,*out,goodrow_counter);
    return goodrow_counter;
}


/** Discard markers whose names are not present in a gsc_KnownGenome
 *
 * All markers in the genome @a g are retained in the list. All markers not
 * in that genome are freed and deleted from the @a markerlist. The new length of
 * the marker list after deletions is returned.
 *
 * Markers without names to join on are discarded.
 *
 * @param g genome object containing the list of markers to compare markers in @a markerlist against
 * @param n_markers_in_list length of @a markerlist
 * @param markerlist Pointer to a list of name/chromosome/position entries for various markers.
 * Only the names will be read, but for each entry, if its name does not exist in the list of markers
 * to keep, the entry will be deleted and the remainder of the marker list shuffled backwards. Markers
 * without names are automatically discarded.
 * @return new number of markers in markerlist, either @a n_markers_in_list or smaller.
 */
static unsigned int gsc_helper_str_markerlist_leftjoin( gsc_KnownGenome g, unsigned int n_markers_in_list, struct gsc_MapfileUnit** markerlist ) {
    struct gsc_MapfileUnit* rlist = *markerlist;
    unsigned int n_joined = 0;
    /*unsigned int consecutivity_bias; // we cache the index of the last name we found and pre-check whether the next marker
    // in the list is the next marker in the genome. For the case where people organise their genotype file and genetic
    // map file in the same order. Edit: decided this is not likely enough a situation to build this in.*/

    for (unsigned int i = 0; i < n_markers_in_list; ++i) {
        if (rlist[i].name != NULL) {
            unsigned int nameix;
            if (gsc_get_index_of_genetic_marker(rlist[i].name, g, &nameix)) {
                if (n_joined != i) {
                    rlist[n_joined] = rlist[i];
                }
                n_joined++;

            } else { // discard this marker. n_joined lags behind i by one more step.
                GSC_FREE(rlist[i].name);

            }
        }
    }

    return n_joined;
}


/** Save a RecombinationMap to the SimData and allocate it a mapID.
 */
static gsc_MapID gsc_helper_insert_recombmap_into_simdata(gsc_SimData* d, gsc_RecombinationMap map) {
    int newmapindex = 0;
    if (d->genome.n_maps > 0) {
        newmapindex = d->genome.n_maps;

        gsc_MapID* tmpMapIDs = gsc_malloc_wrap(sizeof(gsc_MapID)*(newmapindex+1),GSC_TRUE);
        memcpy(tmpMapIDs,d->genome.map_ids,sizeof(gsc_MapID)*newmapindex);
        GSC_FREE(d->genome.map_ids);
        d->genome.map_ids = tmpMapIDs;

        gsc_RecombinationMap* tmpMaps = gsc_malloc_wrap(sizeof(gsc_RecombinationMap)*(newmapindex+1),GSC_TRUE);
        memcpy(tmpMaps,d->genome.maps,sizeof(gsc_RecombinationMap)*newmapindex);
        GSC_FREE(d->genome.maps);
        d->genome.maps = tmpMaps;

    } else {
        d->genome.map_ids = gsc_malloc_wrap(sizeof(gsc_MapID)*1,GSC_TRUE);
        d->genome.maps = gsc_malloc_wrap(sizeof(gsc_RecombinationMap)*1,GSC_TRUE);
    }
    d->genome.map_ids[newmapindex] = gsc_get_new_map_id(d);
    d->genome.n_maps++;
    d->genome.maps[newmapindex] = map;

    return d->genome.map_ids[newmapindex];
}

/** Save an EffectMatrix to the SimData and allocate it an EffectID.
 */
static gsc_EffectID gsc_helper_insert_eff_set_into_simdata(gsc_SimData* d, gsc_EffectMatrix effset) {
    int neweffsetindex = 0;
    if (d->n_eff_sets > 0) {
        neweffsetindex = d->n_eff_sets;

        gsc_EffectID* tmpIDs = gsc_malloc_wrap(sizeof(gsc_EffectID)*(neweffsetindex+1),GSC_TRUE);
        memcpy(tmpIDs,d->eff_set_ids,sizeof(gsc_EffectID)*neweffsetindex);
        GSC_FREE(d->eff_set_ids);
        d->eff_set_ids = tmpIDs;

        gsc_EffectMatrix* tmpMats = gsc_malloc_wrap(sizeof(gsc_EffectMatrix)*(neweffsetindex+1),GSC_TRUE);
        memcpy(tmpMats,d->e,sizeof(gsc_EffectMatrix)*neweffsetindex);
        GSC_FREE(d->e);
        d->e = tmpMats;

    } else {
        d->eff_set_ids = gsc_malloc_wrap(sizeof(gsc_EffectID)*1,GSC_TRUE);
        d->e = gsc_malloc_wrap(sizeof(gsc_EffectMatrix)*1,GSC_TRUE);
    }
    d->eff_set_ids[neweffsetindex] = gsc_get_new_eff_set_id(d);
    d->n_eff_sets++;
    d->e[neweffsetindex] = effset;

    return d->eff_set_ids[neweffsetindex];
}

/** Sort markerlist by chromosome name, and by position within each chromosome.
 *
 * @param n_markers length of @a markerlist
 * @param markerlist list of MapFileUnits to sort
 * Could return number of unique chromosomes found in the markerlist, but we aren't worried about that level of performance right now.
 */
static void gsc_helper_sort_markerlist(unsigned int n_markers, struct gsc_MapfileUnit* markerlist) {
    if (n_markers < 2) { return; }

    // sort by linkage group
    qsort(markerlist,n_markers,sizeof(*markerlist),gsc_helper_mapfileunit_ascending_chr_comparer);

    // sort each linkage group by pos
    //int n_chr = 1;
    unsigned int chr_start = 0;
    unsigned long current_chr = markerlist[0].chr;

    for (unsigned int i = 1; i < n_markers; ++i) {
        if (markerlist[i].chr != current_chr) { // found end of current chr
            //n_chr++;
            qsort(markerlist + chr_start, i - chr_start,
                    sizeof(*markerlist), gsc_helper_mapfileunit_ascending_d_comparer);

            chr_start = i;
            current_chr = markerlist[i].chr;
        }
    }

    qsort(markerlist + chr_start, n_markers - chr_start,
            sizeof(*markerlist), gsc_helper_mapfileunit_ascending_d_comparer);
    //return n_chr;
}

/** Parse a list of markers/chrs/positions into a gsc_RecombinationMap and save to SimData
 *
 * It partitions the provided list by chr, then sorts by position. Positions are interpreted
 * as positions in cM for the purpose of calculating recombination probabilities between
 * adjacent markers in the sorted list. Markers not present in the SimData's list of known
 * markers are discarded. Markers with no names are discarded.
 *
 * @param d SimData into which to load this RecombinationMap
 * @param n_markers length of @a markerlist
 * @param markerlist @see gsc_helper_parse_mapfile
 * @param n_chr If known, number of unique
 * @return gsc_MapID of the gsc_RecombinationMap that was just loaded into the simulation.
 */
gsc_MapID gsc_create_recombmap_from_markerlist(gsc_SimData* d, unsigned int n_markers, struct gsc_MapfileUnit* markerlist) {
    if (n_markers == 0) return NO_MAP;

    GSC_CREATE_BUFFER(chr_nmembers,unsigned int,40);
    memset(chr_nmembers,0,sizeof(*chr_nmembers)*40);
    chr_nmembers[0] = 1;
    unsigned int n_chr = 1;
    unsigned long current_chr = markerlist[0].chr;
    for (unsigned int i = 1; i < n_markers; ++i) {
        while (i < n_markers && markerlist[i].name == NULL) {
            ++i;
        }
        if (current_chr != markerlist[i].chr) {
            // First of next
            if (n_chr >= chr_nmemberscap) {
                GSC_STRETCH_BUFFER(chr_nmembers,2*n_chr);
                memset(chr_nmembers+n_chr,0,sizeof(*chr_nmembers)*n_chr);
            }
            ++n_chr;
            current_chr = markerlist[i].chr;
            chr_nmembers[n_chr-1] = 1;
        } else {
            ++(chr_nmembers[n_chr-1]);
        }
    }

    gsc_RecombinationMap map = {.n_chr=n_chr, .chrs=gsc_malloc_wrap(sizeof(gsc_LinkageGroup) * n_chr, GSC_TRUE) };

    // Populate the map. Each chr/linkage group may be "Simple" or "Reordered"
    unsigned int could_not_match = 0;
    unsigned int current_marker = 0;
    unsigned int first_marker;
	int n_bad_chr = 0;
    for (int chr_ix = 0; chr_ix < map.n_chr; ++chr_ix) {		
        first_marker = current_marker;
        double chrdist = markerlist[first_marker + chr_nmembers[chr_ix] - 1].pos - markerlist[first_marker].pos;
        double* lgdists = gsc_malloc_wrap(sizeof(double)*(chr_nmembers[chr_ix]),GSC_TRUE);

        char found_first = GSC_FALSE;
		// n_goodmembers == 0 is a guard on firsts_coord_in_genome, but we 
		// still initialise it here (to a value too high to be reasonable)
		// because the compiler can't tell that.
        unsigned int firsts_coord_in_genome = d->genome.n_markers; 
        unsigned int n_goodmembers = 0;
        unsigned int* marker_coords = NULL;

        unsigned int endpt = first_marker + chr_nmembers[chr_ix];
        for (; current_marker < endpt; ++current_marker) { // simple recombination map, if possible
            if (markerlist[current_marker].name == NULL) {
                continue;
            }

            if (!found_first) {
                unsigned int coord;
                if (!gsc_get_index_of_genetic_marker(markerlist[current_marker].name, d->genome, &coord)) {
                    could_not_match++;
                } else {
                    found_first = GSC_TRUE;
                    first_marker = current_marker;
                    firsts_coord_in_genome = coord;
                    lgdists[n_goodmembers] = (markerlist[current_marker].pos - markerlist[first_marker].pos) / chrdist;
                    n_goodmembers++;
                }
            } else if (firsts_coord_in_genome + n_goodmembers < d->genome.n_markers &&
                       strcmp(markerlist[current_marker].name, d->genome.marker_names[firsts_coord_in_genome + n_goodmembers]) == 0) {
                // we are a simple linkage group still so far.
                lgdists[n_goodmembers] = (markerlist[current_marker].pos - markerlist[first_marker].pos) / chrdist;
                n_goodmembers++;
            } else {
                // Just discovered we are a reordered linkage group. Copy over the marker indexes that were as expected.
                marker_coords = gsc_malloc_wrap(sizeof(*marker_coords)*(chr_nmembers[chr_ix]),GSC_TRUE);
                for (unsigned int backfill = 0; backfill < n_goodmembers; ++backfill) {
                    marker_coords[backfill] = firsts_coord_in_genome + backfill;
                }
                break;
            }

        }
        for (; current_marker < endpt; ++current_marker) { // reordered recombination map, if previous failed.
            if (markerlist[current_marker].name == NULL) {
                continue;
            }

            unsigned int coord;
            if (!gsc_get_index_of_genetic_marker(markerlist[current_marker].name, d->genome, &coord)) {
                ++could_not_match;
            } else {
                marker_coords[n_goodmembers] = coord;
                lgdists[n_goodmembers] = (markerlist[current_marker].pos - markerlist[first_marker].pos) / chrdist;
                ++n_goodmembers;
            }
        }

		if (n_goodmembers == 0) { // || firsts_coord_in_genome >= d->genome.n_markers) {
			n_bad_chr++;
        } else if (marker_coords == NULL) {
			int chr_ix_actual = chr_ix-n_bad_chr;
            map.chrs[chr_ix_actual].type = GSC_LINKAGEGROUP_SIMPLE;
            map.chrs[chr_ix_actual].map.simple.expected_n_crossovers = chrdist / 100;
            map.chrs[chr_ix_actual].map.simple.n_markers = n_goodmembers;
            map.chrs[chr_ix_actual].map.simple.first_marker_index = firsts_coord_in_genome;
            map.chrs[chr_ix_actual].map.simple.dists = lgdists;
        } else {
			int chr_ix_actual = chr_ix-n_bad_chr;
            map.chrs[chr_ix_actual].type = GSC_LINKAGEGROUP_REORDER;
            map.chrs[chr_ix_actual].map.reorder.expected_n_crossovers = chrdist / 100;
            map.chrs[chr_ix_actual].map.reorder.n_markers = n_goodmembers;
            map.chrs[chr_ix_actual].map.reorder.marker_indexes = marker_coords;
            map.chrs[chr_ix_actual].map.reorder.dists = lgdists;
        }
    }
	GSC_DELETE_BUFFER(chr_nmembers);
	map.n_chr = map.n_chr-n_bad_chr;
	if (map.n_chr == 0) {
		GSC_FREE(map.chrs);
		return NO_MAP;
	}
    return gsc_helper_insert_recombmap_into_simdata(d, map);
}


/** Create a uniformly-spaced gsc_RecombinationMap from a list of marker names and save to SimData
 *
 * The recombination map produced has one chromosome/linkage group. The markers in this
 * chromosome are equally spaced and their
 * ordering matches the ordering of @a markernames.
 *
 * If @a markernames is NULL, @a n_markers is ignored and the uniformly-spaced
 * recombination map is created from all markers listed in @a d->genome.
 *
 * @param d SimData into which to load this RecombinationMap
 * @param n_markers length of @a markernames
 * @param markernames names of the markers to create this simple recombination map from.
 * @param expected_n_recombinations
 * @return gsc_MapID of the gsc_RecombinationMap that was just loaded into the simulation.
 */
gsc_MapID gsc_create_uniformspaced_recombmap(gsc_SimData* d, unsigned int n_markers, char** markernames, double expected_n_recombinations) {
    gsc_RecombinationMap map = {.n_chr=1, .chrs=gsc_malloc_wrap(sizeof(gsc_LinkageGroup)*1, GSC_TRUE) };

    if (markernames == NULL) {
        if (d->genome.n_markers == 0) return NO_MAP;

        double* lgdists = gsc_malloc_wrap(sizeof(double)*d->genome.n_markers,GSC_TRUE);
        double lgdist = 1./d->genome.n_markers;
        for (unsigned int i = 0; i < d->genome.n_markers; ++i) { lgdists[i] = lgdist; }

        map.chrs[0].type = GSC_LINKAGEGROUP_SIMPLE;
        map.chrs[0].map.simple.expected_n_crossovers = expected_n_recombinations;
        map.chrs[0].map.simple.n_markers = d->genome.n_markers;
        map.chrs[0].map.simple.first_marker_index = 0;
        map.chrs[0].map.simple.dists = lgdists;
    } else {
        if (n_markers == 0) return NO_MAP;

        // markernames could still be simple or reordered compared to the d->genome, so need to check that first.
        char found_first = GSC_FALSE;
        unsigned int could_not_match;
        unsigned int firsts_coord_in_genome = d->genome.n_markers;
        unsigned int chrmarker_ix = 0;

        unsigned int* marker_coords = NULL;
        for (unsigned int i = 0; i < n_markers; ++i) {
            if (!found_first || marker_coords != NULL) {
                // We are first or we are a reordered linkage group. Find what index in the genome the next marker is stored at.
                unsigned int coord;

                if (markernames[i] == NULL) {
                    could_not_match++;
                } else if (!gsc_get_index_of_genetic_marker(markernames[i], d->genome, &coord )) {
                    could_not_match++;
                } else if (!found_first) {
                    found_first = GSC_TRUE;
                    firsts_coord_in_genome = coord;
                    chrmarker_ix++;
                } else { // must be the case that we have marker_coords != NULL and are a reordered linkage group
                    marker_coords[chrmarker_ix] = coord;
                    chrmarker_ix++;
                }

            } else if (firsts_coord_in_genome < d->genome.n_markers && 
					strcmp(markernames[i], d->genome.marker_names[firsts_coord_in_genome + i]) == 0) {
                // are a simple linkage group still so far.
                chrmarker_ix++;

            } else {
                // Just discovered we are a reordered linkage group. Copy over the marker indexes that were as expected.
                marker_coords = gsc_malloc_wrap(sizeof(*marker_coords)*n_markers,GSC_TRUE);
                for (unsigned int backfill = 0; backfill < chrmarker_ix; ++backfill) {
                    marker_coords[backfill] = firsts_coord_in_genome + backfill;
                }

                if (markernames[i] == NULL) {
                    could_not_match++;
                } else if  (!gsc_get_index_of_genetic_marker(markernames[i], d->genome, &(marker_coords[chrmarker_ix]) )) {
                    could_not_match++;
                } else {
                    chrmarker_ix++;
                }
            }
        }

        double* lgdists = gsc_malloc_wrap(sizeof(double)*chrmarker_ix,GSC_TRUE);
        double lgdist = 1./n_markers;
        for (unsigned int i = 0; i < chrmarker_ix; ++i) { lgdists[i] = lgdist; }

        if (marker_coords == NULL) {
            map.chrs[0].type = GSC_LINKAGEGROUP_SIMPLE;
            map.chrs[0].map.simple.expected_n_crossovers = expected_n_recombinations;
            map.chrs[0].map.simple.n_markers = chrmarker_ix;
            map.chrs[0].map.simple.first_marker_index = firsts_coord_in_genome;
            map.chrs[0].map.simple.dists = lgdists;
        } else {
            map.chrs[0].type = GSC_LINKAGEGROUP_REORDER;
            map.chrs[0].map.reorder.expected_n_crossovers = expected_n_recombinations;
            map.chrs[0].map.reorder.n_markers = chrmarker_ix;
            map.chrs[0].map.reorder.marker_indexes = marker_coords;
            map.chrs[0].map.reorder.dists = lgdists;
        }
    }

    return gsc_helper_insert_recombmap_into_simdata(d,map);
}

/** Load a genetic map to a gsc_SimData object.
 *
 * Can be called on an empty SimData, in which case it loads the first genetic map,
 * sets up the list of known markers in the genome, and warns that no genotypes are loaded
 * yet so many simulation functions will not yet be able to run. Can also be called on a
 * SimData with existing map and genotypes, in which case the map file's contents are loaded
 * as another recombination map.
 *
 * The file's format should be:
 *
 * marker	chr	pos
 *
 * [marker name] 	[chr]	[pos]
 *
 * [marker name] 	[chr]	[pos]
 *
 * ...
 *
 * The header line is optional. If no header line is provided, it is assumed that the columns are in
 * the above order (marker, then chromosome, then position). If a header is provided, the columns
 * "marker", "chr", and "pos" can be in any order.
 *
 * Map positions are assumped to be in centimorgans.
 *
 * Extra columns are not permitted. Columns may be space-separated or tab-separated.
 * Column separators need not be consistent throughout the file, and may be made of
 * multiple characters. That is, any consecutive sequence of spaces, tabs, and commas is
 * interpreted as one column separator.
 *
 * Marker names can be of any length and any characters excluding column separator (tab, space, or comma)
 * and newline characters.
 *
 * @shortnamed{load_mapfile}
 *
 * @param d pointer to gsc_SimData to be populated
 * @param filename string name/path of file containing genetic map data.
*/
gsc_MapID gsc_load_mapfile(SimData* d, const char* filename) {
    if (filename == NULL) return NO_MAP;

    struct gsc_MapfileUnit* mapcontents = NULL;
    unsigned int nrows = gsc_helper_parse_mapfile(filename,&mapcontents);
    if (nrows == 0 || mapcontents == NULL) {
        if (mapcontents != NULL) {
            GSC_FREE(mapcontents);
        }
        return NO_MAP;
    }

    char freeMapNames = GSC_TRUE;
    if (d->genome.n_markers > 0) {
        // if genome is already set, leftjoin on those markers.
        unsigned int new_nrows = gsc_helper_str_markerlist_leftjoin(d->genome, nrows, &mapcontents);
        if (new_nrows < nrows) {
            Rprintf("Discarded %u markers when loading map %s because they do not appear in the primary map.\n", (unsigned int) (nrows - new_nrows), filename);
        }
        nrows = new_nrows;
        gsc_helper_sort_markerlist(nrows,mapcontents);
    } else {
        // else set up the list of markers tracked by the simulation
        gsc_helper_sort_markerlist(nrows,mapcontents);
        d->genome = (gsc_KnownGenome){
            .n_markers = nrows,
            .marker_names = gsc_malloc_wrap(sizeof(char**)*nrows,GSC_TRUE),
            .names_alphabetical = gsc_malloc_wrap(sizeof(char**)*nrows,GSC_TRUE),
            .n_maps = 0,
            .map_ids = NULL,
            .maps = NULL
        };
        for (int i = 0; i < d->genome.n_markers; ++i) {
            d->genome.marker_names[i] = mapcontents[i].name;
            d->genome.names_alphabetical[i] = &(d->genome.marker_names[i]);
        }
        qsort(d->genome.names_alphabetical,d->genome.n_markers,sizeof(*d->genome.names_alphabetical),gsc_helper_indirect_alphabetical_str_comparer);

        freeMapNames = GSC_FALSE;
        //printf( "Warning: loading genetic map before loading any founder genotypes. Many simulation operations will not yet run.\n");
    }

    gsc_MapID map = gsc_create_recombmap_from_markerlist(d, nrows, mapcontents);
    if (freeMapNames) {
        for (int i = 0; i < nrows; ++i) {
            GSC_FREE(mapcontents[i].name);
        }
    }
    GSC_FREE(mapcontents);

    return map;
}


/** Populates a gsc_SimData combination with effect values.
 *
 * If the gsc_SimData does not have a list of tracked markers already (i.e. no map
 * data loaded), this function will be unable to match any marker names so
 * will fail to load any marker effects.
 *
 * The file's format should be:
 *
 * marker allele eff
 *
 * [marker] [allele] [effect]
 *
 * [marker] [allele] [effect]
 *
 * ...
 *
 * The header line is optional. If no header line is provided, it is assumed that the columns are in
 * the above order (marker, then allele, then marker effect). If a header is provided, the columns
 * "marker", "allele", and "eff" can be in any order.
 *
 * Extra columns are not permitted. Columns may be space-separated or tab-separated.
 * Column separators need not be consistent throughout the file, and may be made of
 * multiple characters. That is, any consecutive sequence of spaces, tabs, and commas is
 * interpreted as one column separator.
 *
 * Marker names can be of any length and any characters excluding column separator (tab, space, or comma)
 * and newline characters. Rows where the marker name does not match any marker in the
 * SimData's tracked markers will be ignored.
 *
 * @shortnamed{load_effectfile}
 *
 * @param d pointer to gsc_SimData to be populated.
 * @param filename string name/path of file containing effect values.
 * @returns the ID of the set of marker effects just loaded
*/
gsc_EffectID gsc_load_effectfile(gsc_SimData* d, const char* filename) {
    if (filename == NULL) return GSC_NO_EFFECTSET;
    if (d->genome.n_markers == 0) return GSC_NO_EFFECTSET;

    gsc_TableFileReader tf = gsc_tablefilereader_create(filename);

    unsigned int row = 1;
    unsigned int col = 1;

    gsc_TableFileCell cellsread[4] = { 0 };
    gsc_TableFileCell* cellqueue = cellsread;
    const char* titles[] = { "marker", "allele", "eff"};
    int colnums[] = { 1, 2, 3 };
    unsigned int queuesize;
    int header = gsc_helper_parse_3cell_header(&tf, titles, colnums, cellqueue, &queuesize);
    if (header == GSC_TRUE) {
        Rprintf("(Loading %s) Format: effect file with header\n", filename);
    } else if (header == GSC_FALSE) {
        Rprintf("(Loading %s) Format: effect file without header\n", filename);
    } else {
        Rprintf("(Loading %s) Failure: Cannot identify the expected 3 columns of the effect file\n", filename);
        gsc_tablefilereader_close(&tf);
        return NO_EFFECTSET;
    }
    int marker_colnum = colnums[0], allele_colnum = colnums[1], eff_colnum = colnums[2];
	

    int goodrow = (header) ? GSC_FALSE : GSC_TRUE; // discard first row if it's a header, keep if it's not.
    unsigned int goodrow_counter = 0;
    unsigned int allele_counter = 0;

    GSC_CREATE_BUFFER(effset_alleles,char,2);
    GSC_CREATE_BUFFER(effset_rows,double*,2);
    unsigned int n_effset_rows = 0;
	
    unsigned int markerix;
    char allele = '\0';
    unsigned int alleleix = n_effset_rows + 1; // invalid value. 
    double effect = 0;
    char* conversionflag;

    gsc_TableFileCell ncell;

    do {
        ncell = gsc_helper_tablefilereader_get_next_cell_wqueue(&tf, &cellqueue, &queuesize);

        if (ncell.cell != NULL) { // so that we can cope with missing final newline
            // Update row/col position and save predecessor row
            if (ncell.predNewline) {
                if (goodrow && col >= 3) { // save predecessor row.
                    if (alleleix == n_effset_rows) {
                        ++allele_counter;
                        ++n_effset_rows;
                        if (effset_allelescap < n_effset_rows) {
                            GSC_STRETCH_BUFFER(effset_alleles,2*n_effset_rows);
                            GSC_STRETCH_BUFFER(effset_rows,2*n_effset_rows);
                        }
                        effset_alleles[alleleix] = allele;
                        effset_rows[alleleix] = gsc_malloc_wrap(sizeof(double)*d->genome.n_markers,GSC_TRUE);
                        memset(effset_rows[alleleix],0,sizeof(double)*d->genome.n_markers);
                    } 
					if (alleleix < n_effset_rows) { // this is specifically not an "else if" 
						effset_rows[alleleix][markerix] = effect;
						++goodrow_counter;
					}
                }
                row += ncell.predNewline;
                goodrow = GSC_TRUE;
                col = 1;
            }
            col += (ncell.predCol > 0) ? 1 : 0; // multiple column spacers treated as one

            // Parse this cell
            if (ncell.cell_len == 0) {
                goodrow = GSC_FALSE;
            } else if (col == marker_colnum) {
                char tmp = ncell.cell[ncell.cell_len]; ncell.cell[ncell.cell_len] = '\0';
                int validmarker = gsc_get_index_of_genetic_marker(ncell.cell,d->genome,&markerix);
                ncell.cell[ncell.cell_len] = tmp;
                if (!validmarker) {
                    goodrow = GSC_FALSE;
                    //warning("Entry at row %i column %i of file %s does not match the name of a tracked marker\n", row, marker_colnum, filename);
                }

            } else if (col == allele_colnum) {
                if (ncell.cell_len > 1) {
                    goodrow = GSC_FALSE;
                    //warning("Entry at row %i column %i of file %s was too long to represent a single allele\n", row, allele_colnum, filename);
                }
                allele = ncell.cell[0];
                for (alleleix = 0; alleleix < n_effset_rows; ++alleleix) {
                    if (effset_alleles[alleleix] == allele) {
                        break;
                    }
                } // leave this loop with alleleix set to match the allele, or equal to n_effset_rows if it's a new allele.

            } else if (col == eff_colnum) {
                char tmp = ncell.cell[ncell.cell_len]; ncell.cell[ncell.cell_len] = '\0';
                effect = strtod(ncell.cell,&conversionflag);
                ncell.cell[ncell.cell_len] = tmp;
                if (conversionflag != ncell.cell + ncell.cell_len) { // unsuccessful read
                    goodrow = GSC_FALSE;
                    //warning("Entry at row %i column %i of file %s could not be parsed as a numeric value\n", row, eff_colnum, filename);
                }

            } else {
                goodrow = GSC_FALSE;
            }

            // Reset
            if (!ncell.isCellShallow) { GSC_FREE(ncell.cell); }
        }
    } while (!ncell.eof);

    if (col == 3 && goodrow) { // save predecessor row
        ++goodrow_counter;
        if (alleleix == n_effset_rows) {
            ++allele_counter;
            ++n_effset_rows;
            if (effset_allelescap < n_effset_rows) {
                GSC_STRETCH_BUFFER(effset_alleles,2*n_effset_rows);
                GSC_STRETCH_BUFFER(effset_rows,2*n_effset_rows);
            }
            effset_alleles[alleleix] = allele;
            effset_rows[alleleix] = gsc_malloc_wrap(sizeof(double)*d->genome.n_markers,GSC_TRUE);
            memset(effset_rows[alleleix],0,sizeof(double)*d->genome.n_markers);
        }
        effset_rows[alleleix][markerix] = effect;
    }

    Rprintf("(Loading %s) %u effect value(s) spanning %u allele(s) were loaded. Failed to parse %u line(s).\n", filename, (unsigned int) goodrow_counter, (unsigned int) allele_counter, (unsigned int) (row - header - goodrow_counter));
    gsc_tablefilereader_close(&tf);

    if (n_effset_rows > 0) {
        gsc_EffectMatrix effset = { 0 };
        GSC_FINALISE_BUFFER(effset_alleles, effset.effect_names, n_effset_rows);
        GSC_FINALISE_BUFFER(effset_rows, effset.effects.matrix, n_effset_rows);
        effset.effects.rows = n_effset_rows;
        effset.effects.cols = d->genome.n_markers;
        return gsc_helper_insert_eff_set_into_simdata(d, effset);
    } else {
        GSC_DELETE_BUFFER(effset_alleles);
        GSC_DELETE_BUFFER(effset_rows);
        return GSC_NO_EFFECTSET;
    }
}

/** Identify what formatting a genotype matrix is representing alleles as
 *
 * @see gsc_GenotypeFileCellStyle
 *
 * For IUPAC encoding: https://genome.ucsc.edu/goldenPath/help/iupac.html
 */
static enum gsc_GenotypeFileCellStyle gsc_helper_genotype_matrix_identify_cell_style(gsc_TableFileCell c) {
    switch (c.cell_len) {
    case 1:
        switch (c.cell[0]) {
        case '0':
        case '1':
        case '2':
            return GSC_GENOTYPECELLSTYLE_COUNT;
        case 'G': // G
        case 'A': // A
        case 'T': // T
        case 'C': // C
        case 'R': // G/A
        case 'Y': // T/C
        case 'M': // A/C
        case 'K': // G/T
        case 'S': // G/C
        case 'W': // A/T
        case 'N': // any
            return GSC_GENOTYPECELLSTYLE_ENCODED;
        default:
            break;
        }
        break;
    case 2:
        if (c.cell[0] == 'm') { // m[numeric] case, which is probably a marker not an allele pair
            switch (c.cell[1]) {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                return GSC_GENOTYPECELLSTYLE_UNKNOWN;
            default:
                break;
            }
        }
        return GSC_GENOTYPECELLSTYLE_PAIR;
    case 3:
        if (c.cell[1] == '/') {
            return GSC_GENOTYPECELLSTYLE_SLASHPAIR;
        }
        break;
    default:
        break;
    }
    return GSC_GENOTYPECELLSTYLE_UNKNOWN;
}

/** Parse a string and save it as the alleles of a genotype at a particular location and genetic marker.
 *
 * If the style in which the allele string is encoded does not denote phase,
 * linkage phase is randomly chosen.
 */
static void gsc_helper_genotypecell_to_allelematrix(GenoLocation loc, unsigned int markerix,
                                                    enum gsc_GenotypeFileCellStyle style, char* cell, gsc_SimData* forrng) {
    char* pos = loc.localAM->alleles[loc.localPos] + 2*markerix;
    int phase = 0;
    switch (style) {
    case GSC_GENOTYPECELLSTYLE_PAIR:
        pos[0] = cell[0];
        pos[1] = cell[1];
        break;
    case GSC_GENOTYPECELLSTYLE_SLASHPAIR:
        pos[0] = cell[0];
        pos[1] = cell[2];
        break;
    case GSC_GENOTYPECELLSTYLE_COUNT:
        switch (cell[0]) {
        case '0':
            pos[0] = 'T';
            pos[1] = 'T';
            break;
        case '1':
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'A';
            pos[1-phase] = 'T';
            break;
        case '2':
            pos[0] = 'A';
            pos[1] = 'A';
            break;
        }
        break;
    case GSC_GENOTYPECELLSTYLE_ENCODED:
        switch (cell[0]) {
        case 'G': // G
            pos[0] = 'G';
            pos[1] = 'G';
            break;
        case 'A': // A
            pos[0] = 'A';
            pos[1] = 'A';
            break;
        case 'T': // T
            pos[0] = 'T';
            pos[1] = 'T';
            break;
        case 'C': // C
            pos[0] = 'C';
            pos[1] = 'C';
            break;
        case 'R': // G/A
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'G';
            pos[1-phase] = 'A';
            break;
        case 'Y': // T/C
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'T';
            pos[1-phase] = 'C';
            break;
        case 'M': // A/C
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'A';
            pos[1-phase] = 'C';
            break;
        case 'K': // G/T
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'G';
            pos[1-phase] = 'T';
            break;
        case 'S': // G/C
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'G';
            pos[1-phase] = 'C';
            break;
        case 'W': // A/T
            phase = (unif_rand() > 0.5);
            pos[phase]   = 'A';
            pos[1-phase] = 'T';
            break;
        default:
            break;
        }
        break;
    default: break;
    }
}

/** Create a new gsc_EmptyListNavigator, including an empty AlleleMatrix suitable for inserting into
 *  the provided SimData once the emptylist is finalised.
 */
static struct gsc_EmptyListNavigator gsc_create_emptylistnavigator(gsc_SimData* d, gsc_GroupNum allocation_group) {
    struct gsc_EmptyListNavigator me = { .d=d, .localPos = 0, .alloctogroup = allocation_group, .currentid = d->current_id };
    me.firstAM = gsc_create_empty_allelematrix(me.d->genome.n_markers, me.d->n_labels, me.d->label_defaults, 0);
    me.localAM = me.firstAM;
    return me;
}

/** Reset the cursor of a gsc_EmptyListNavigator to the first genotype
 *
 * @return GenoLocation for the first sequential genotype in the emptylist
 */
static gsc_GenoLocation gsc_emptylistnavigator_get_first(struct gsc_EmptyListNavigator* it) {
    it->localAM = it->firstAM;
    it->localPos = 0;
    if (1 > it->localAM->n_genotypes) {
        it->localAM->n_genotypes = 1;
        it->localAM->alleles[0] = gsc_malloc_wrap(sizeof(char) * (it->localAM->n_markers<<1),GSC_TRUE);
        memset(it->localAM->alleles[0], 0, sizeof(char) * (it->localAM->n_markers<<1));
        it->localAM->names[0] = NULL;
        it->localAM->groups[0] = it->alloctogroup;
        ++(it->currentid.id);
        it->localAM->ids[0] = it->currentid;
    }
    return (gsc_GenoLocation){.localAM=it->localAM, .localPos =it->localPos};
}

/** Get the next sequential genotype in an gsc_EmptyListNavigator
 *
 * @return GenoLocation for the next sequential genotype in the emptylist
 */
static gsc_GenoLocation gsc_emptylistnavigator_get_next(struct gsc_EmptyListNavigator* it) {
    if (CONTIG_WIDTH - 1 == it->localPos) {
        if (NULL == it->localAM->next) {
            gsc_AlleleMatrix* next = gsc_create_empty_allelematrix(it->d->genome.n_markers, it->d->n_labels, it->d->label_defaults, 0);
            it->localAM->next = next;
            it->localAM = next;
            it->localPos = 0;
        } else {
            it->localAM = it->localAM->next;
            it->localPos = 0;
        }
    } else {
        ++(it->localPos);
    }

    if (it->localAM->n_genotypes <= it->localPos) {
        if (1 < it->localPos - it->localAM->n_genotypes) {
            warning("EmptyListNavigator invalid\n");
            return INVALID_GENO_LOCATION;
        }
        ++(it->localAM->n_genotypes);

        it->localAM->alleles[it->localPos] = gsc_malloc_wrap(sizeof(char) * (it->localAM->n_markers<<1),GSC_TRUE);
        memset(it->localAM->alleles[it->localPos], 0, sizeof(char) * (it->localAM->n_markers<<1));
        it->localAM->names[it->localPos] = NULL;
        it->localAM->groups[it->localPos] = it->alloctogroup;
        ++(it->currentid.id);
        it->localAM->ids[it->localPos] = it->currentid;
    }

    return (gsc_GenoLocation){.localAM=it->localAM, .localPos =it->localPos};
}

/** Push emptylist edited genotypes into the SimData
 *
 *  This is the expected final step in use of gsc_EmptyListNavigator
 */
static void gsc_emptylistnavigator_finaliselist(struct gsc_EmptyListNavigator* it) {
    if (NULL == it->d->m) {
        it->d->m = it->firstAM;
    } else {
        gsc_AlleleMatrix* listend = it->d->m;
        while (NULL != listend->next) {
            listend = listend->next;
        }
        listend->next = it->firstAM;
        gsc_condense_allele_matrix(it->d);
    }
    it->d->current_id = it->currentid;
}

/** Determine whether a genotype matrix is row- or column-oriented
 *
 * The decision of whether a genotype matrix is storing markers as columns or rows is based on the following factors:
 * - If no names from the column headers of the file cannot be matched to marker names from the SimData 
 * (either because no markers are tracked by the SimData, because there is no column header in this matrix,
 * or simply because none of the names match names tracked by the SimData), it is assumed that rows in the 
 * genotype matrix represent markers.
 * - If any name from the column header of the genotype file can be matched to the name of a genetic marker tracked 
 * by the SimData, then it is assumed that columns in the genotype matrix represent markers. 
 * 
 * genomicSimulation requires marker names for file loading, so if markers are represented by columns in the matrix,
 * the format specifier returned will also specify that the matrix does have a header row.
 *
 * @param d SimData object the genotype file is assumed to be loaded to. Required because determining 
 * matrix orientation requires checking marker names against the markers tracked in the SimData.
 * @param cellqueue array of cells read from the matrix. Requires at least the first row of the matrix, preferably
 * the first row of the matrix plus the first two cells (header and first body cell) of the second row of the matrix.
 * @param firstrowlen number of cells in the first row of the matrix. These cells can be accessed at indices 0-{firstrowlen-1}
 * in @a cellqueue.
 * @param queuelen number of cells saved sequentially in @a cellqueue. It is assumed that queuelen >= firstrowlen.
 * @param format established details about the formatting of the genotype matrix
 * @param filenameforlog The name of the matrix file or other identifier for this matrix that should be used 
 * when printing out logging information.
 * @return a format specification for the genotype matrix, with details about the orientation of the matrix added if
 * those were not already supplied.
 */
static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_orientation(const SimData* d, 
        const gsc_TableFileCell* cellqueue, const unsigned int firstrowlen, const unsigned int queuelen,
        struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog) {
            
    if (format.markers_as_rows == GSC_TRUE || format.markers_as_rows == GSC_FALSE) {
        // pass
    } else if (d->genome.n_maps == 0) {
        // If there is no genetic map, we cannot check the row/column headers to see if any of them match the marker names..
        // Default to markers being rows

        Rprintf("(Loading %s) Format axis: genetic markers are -rows-, founder lines are |columns| (by assumption when no genetic map is loaded)\n", filenameforlog);
        Rprintf("(Loading %s) No genetic map is loaded, will invent a map with equal spacing of these genetic markers (1cM apart)\n", filenameforlog);
        format.markers_as_rows = GSC_TRUE;

    } else if (format.has_header == GSC_FALSE) {
        Rprintf("(Loading %s) Format axis: genetic markers are -rows-, founder lines are |columns| "
                "(by assumption when matrix has no header row)\n", filenameforlog);
        format.markers_as_rows = GSC_TRUE;

    } else {
        // Note: by here, either the user has told us there is a header row, or we get to detect whether there is one. So will investigate it by comparing names to what's in our map
        // taken from older function gsc_helper_genotypefile_matrix_check_markers_are_rows
        int firstsafeheaderindex = -1;
        if (firstrowlen > 1) {
            firstsafeheaderindex = 1;
        } else if (firstrowlen == 1 && queuelen > firstrowlen + 1) { // second row has more than one cell read.
            firstsafeheaderindex = 0; // assume there's no corner cell
            format.has_header = GSC_TRUE;
        }

        if (firstsafeheaderindex >= 0) {
            // Don't check the "first" cell. It might be a corner cell between the two headers, whose value should be ignored
            // Check the next cell in the first row.
            if (gsc_get_index_of_genetic_marker(cellqueue[firstsafeheaderindex].cell, d->genome, NULL)) {
                Rprintf("(Loading %s) Format axis: genetic markers are |columns|, founder lines are -rows-\n", filenameforlog);
                format.markers_as_rows = GSC_FALSE;
                format.has_header = GSC_TRUE;
                return format;
            }

            // If that wasn't a match, check the first row header, if it exists:
            if (queuelen > firstrowlen && !cellqueue[firstrowlen].eof &&
                    gsc_get_index_of_genetic_marker(cellqueue[firstrowlen].cell, d->genome, NULL)) {
                Rprintf("(Loading %s) Format axis: genetic markers are -rows-, founder lines are |columns|\n", filenameforlog);
                format.markers_as_rows = GSC_TRUE;
                return format;
            }

            // Check remaining column headers
            for (unsigned int i = firstsafeheaderindex + 1; i < firstrowlen; ++i) {
                if (gsc_get_index_of_genetic_marker(cellqueue[i].cell, d->genome, NULL)) {
                    Rprintf("(Loading %s) Format axis: genetic markers are |columns|, founder lines are -rows-\n", filenameforlog);
                    format.markers_as_rows = GSC_FALSE;
                    format.has_header = GSC_TRUE;
                    return format;
                }
            }

        }
        Rprintf("(Loading %s) Format axis: genetic markers are -rows-, founder lines are |columns| (by default file format)\n", filenameforlog);
        format.markers_as_rows = GSC_TRUE;
    }
    return format;

}

/** Determine the style in which alleles are stored in a genotype matrix 
 *
 * @param cellqueue array of cells read from the matrix. It should contain at least one body cell, if the matrix has any boy cells: that is, at least two cells if @a format.has_header is false, or an entire first row plus two cells if @a format.has_header is unknown or true.
 * @param firstrowlen number of cells in the first row of the matrix.
 * @param queuelen number of cells saved sequentially in @a cellqueue. It is assumed that queuelen >= firstrowlen.
 * @param format established details about the formatting of the genotype matrix
 * @param filenameforlog The name of the matrix file or other identifier for this matrix that should be used 
 * when printing out logging information. 
 * @return a format specification for the genotype matrix, with details about the style of
 * allele encoding used in the body of the matrix added if those were not already supplied.
 */
static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_cellstyle(const gsc_TableFileCell* cellqueue, 
        const unsigned int firstrowlen, const unsigned int queuelen, struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog) {

    int style_detected = GSC_FALSE;
    int single_col_file = GSC_FALSE;
    // 1. Detect format if not yet provided.
    if (format.cell_style == GSC_GENOTYPECELLSTYLE_UNKNOWN) {
        style_detected = GSC_TRUE;

        if (firstrowlen == queuelen || cellqueue[firstrowlen].eof) { // There is only one row. Short-circuiting necessary
            // if there is also only one column, we have no body cells to detect the style of
            if (firstrowlen > 1) {
                // Detection path for a single-line file. If it has a header, then this value might end up ignored
                format.cell_style = gsc_helper_genotype_matrix_identify_cell_style(cellqueue[1]);
            } else {
                single_col_file = GSC_TRUE; // one-cell file. needs the warning.
            }
        } else { // there is more than one row
            // If there is only one column, there are no body cells with style to detect
            if (firstrowlen + 1 < queuelen && cellqueue[firstrowlen+1].predNewline < 1) {
                // Detection path. There exists a second cell on the second line that we can read
                format.cell_style = gsc_helper_genotype_matrix_identify_cell_style(cellqueue[firstrowlen+1]);
            } else {
                single_col_file = GSC_TRUE;
            }
        }
    }
    
    // 2. Print cell style detection logs 
    if (style_detected) {
        switch(format.cell_style) {
            case GSC_GENOTYPECELLSTYLE_PAIR:      Rprintf("(Loading %s) Allele format: phased allele pairs\n", filenameforlog); break;
            case GSC_GENOTYPECELLSTYLE_SLASHPAIR: Rprintf("(Loading %s) Allele format: phased allele pairs (slash-separated)\n", filenameforlog); break;
            case GSC_GENOTYPECELLSTYLE_COUNT:     Rprintf("(Loading %s) Allele format: reference allele counts (phase will be randomised)\n", filenameforlog); break;
            case GSC_GENOTYPECELLSTYLE_ENCODED:   Rprintf("(Loading %s) Allele format: IUPAC encoded pair (phase will be randomised)\n", filenameforlog); break;
            case GSC_GENOTYPECELLSTYLE_UNKNOWN:
                if (single_col_file || firstrowlen == queuelen ||
                        (firstrowlen + 1 == queuelen && cellqueue[firstrowlen].eof && cellqueue[firstrowlen].cell_len == 0)) {
                    Rprintf("(Loading %s) Warning: empty genotype matrix. No genotypes will be loaded.\n", filenameforlog);
                } else {
                    fprintf(stderr,"(Loading %s) Failure: Unable to determine the formatting of pairs of alleles."
                    " Check genomicSimulation manual for accepted allele pair encodings\n", filenameforlog);
                }
        }
    }
    
    return format;
}

/** Determine whether a genotype matrix has a header row or not
 *
 * It is assumed that, during the detection process, @see gsc_helper_genotypefile_matrix_detect_orientation 
 * and @see gsc_helper_genotypefile_matrix_detect_cellstyle
 * will be called before this function. Doing so in the opposite order risks that the detection
 * functions assume a file has an invalid configuration even though there is another 
 * potential configuration that would be valid to load.
 *
 * This function does not guarantee that the returned format has a known (true or false)
 * has_header field. has_header may still be unknown.
 *
 * The decision as to whether the genotype matrix has a header row is based on:
 * - If the file is a single column, then it has no header row. 
 * - If any cell in the first row does not match @a format.cell_style, then the first
 * row is assumed to be a header row. 
 * - If all cells in the first row match @a format.cell_style, then there is assumed to
 * be no header row. 
 *
 * @param cellqueue array of cells read from the matrix. It should contain at least two cells, and ideally the whole of the first row of the matrix.
 * @param firstrowlen number of cells in the first row of the matrix.
 * @param queuelen number of cells saved sequentially in @a cellqueue. It is assumed that queuelen >= firstrowlen.
 * @param format established details about the formatting of the genotype matrix
 * @param filenameforlog The name of the matrix file or other identifier for this matrix that should be used 
 * when printing out logging information. 
 * @return a format specification for the genotype matrix, with the field defining whether
 * there is a header row in the matrix added, if that were not already supplied.
 */
static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_header(const gsc_TableFileCell* cellqueue, 
        const unsigned int firstrowlen, const unsigned int queuelen, struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog) {
    // Validity check: if genetic markers are columns, header row is mandatory
    if (format.has_header == GSC_FALSE && format.markers_as_rows == GSC_FALSE) {
        Rprintf("(Loading %s) Failure: genetic markers cannot be represented by columns when matrix has no header row\n", filenameforlog);
        format.has_header = GSC_UNINIT;
        return format;
    }

    // Detect header if we need to detect it.
    if (format.has_header != GSC_FALSE && format.has_header != GSC_TRUE) {
        if (firstrowlen == 1) { 
            // we could have a single-column file (no header assumed), or 
            // we could be a two-column file with no corner cell (must have a header)
            if (queuelen > 2) {
                if (cellqueue[2].eof || cellqueue[2].predNewline) {
                    format.has_header = GSC_FALSE; // single column file
                } else {
                    format.has_header = GSC_TRUE; 
                }
            } // else can't draw any conclusions.
            
        } else if (format.cell_style != GSC_GENOTYPECELLSTYLE_UNKNOWN) {
            // Idea: if we find a cell in the first row that doesn't match the expected cell style, then that first row is probably a header
            format.has_header = GSC_FALSE;
            for (unsigned int i = 1; i < firstrowlen; ++i) { // ignore first cell in row, it could be a corner cell or row header
                if (gsc_helper_genotype_matrix_identify_cell_style(cellqueue[i]) != format.cell_style) {
                    format.has_header = GSC_TRUE;
                    break;
                }
            }
        } // else don't know how to detect.

        switch (format.has_header) {
            case GSC_FALSE: Rprintf("(Loading %s) Format: genotype matrix without header row\n", filenameforlog); break;
            case GSC_TRUE: Rprintf("(Loading %s) Format: genotype matrix with header row\n", filenameforlog); break;
            default: warning("(Loading %s) Failure: Unable to determine whether file has header row\n", filenameforlog); break;
        }
    }

    return format;
}

/** Determine whether a genotype matrix has a corner cell or not 
 *
 * A corner cell is a non-blank value in the first-row/first-column of the matrix,
 * if the matrix has both a header row and header column.
 *
 * @param ncellsfirstrow number of cells in the first row of the matrix.
 * @param ncellssecondrow number of cells in the second row of the matrix.
 * @param secondrowheaderisempty truthy if the value in the first column of the second 
 * row of the matrix is blank, falsy if there is a non-blank value in the first cell 
 * of that second row of the matrix. (note: genomicSimulation's TableFileReader
 * takes any blank space (including spaces at the beginning of a row) as cell dividers.
 * That interface is the assumed input for calculating these parameters).
 * @returns GSC_TRUE if there is a corner cell, GSC_FALSE if there is not a corner cell,
 * and GSC_UNINIT if this was unable to be determined.
 */
static int gsc_helper_genotypefile_matrix_detect_cornercell_presence(const unsigned int ncellsfirstrow, const unsigned int ncellssecondrow, const int secondrowheaderisempty) {
    if (ncellssecondrow == ncellsfirstrow + 1) {
        return GSC_FALSE;
    } else if (ncellssecondrow == ncellsfirstrow) {
        if (secondrowheaderisempty) {
            return GSC_FALSE; //genotype name is simply empty, making the second row look one column shorter than reality
        } else {
            return GSC_TRUE;
        }
    } else if (ncellssecondrow == ncellsfirstrow - 1 && secondrowheaderisempty) {
        return GSC_TRUE; // genotype name on row 2 is empty but corner cell is not
    } else {
        return GSC_UNINIT;
    }
}

/** Give genomicSimulation hints on the format of a genotype matrix file to be loaded.
 *
 * Sometimes genomicSimulation's automatic file formatting detection may misinterpret
 * the formatting of a genotype matrix (eg assuming markers are columns, when they are
 * actually rows of the matrix; assuming there is no header row even though there is one;
 * being unable to determine that the body of the matrix are alternate allele counts because 
 * there are some confusingly-placed "NA"s). For particularly large files, the file formatting
 * detection process might slow down file imports or require more memory. 
 *
 * To bypass part or all of the formatting detection steps when importing a genotype matrix
 * file, this function can be used to provide the final parameter for
 * @see gsc_load_genotypefile or @see gsc_load_data_files.
 * 
 * @param has_header GSC_TRUE if the genotype matrix to be imported definitely has a header
 * row, GSC_FALSE if the genotype matrix has no header row, or some other value (eg GSC_UNINIT)
 * to not bypass the header detection steps of the import process. 
 * @param markers_as_rows GSC_TRUE if each row in the genotype matrix represents a genetic
 * marker, GSC_FALSE if each column of the genotype matrix represents a genetic marker, or 
 * some other value (eg GSC_UNINIT) to not bypass the orientation detection steps of the import
 * process.
 * @param cell_style The style in which the alleles of a candidate at a marker are encoded 
 * in the body cells of the genotype matrix. Use GSC_GENOTYPECELLSTYLE_UNKNOWN to not 
 * bypass the cell style detection step of the import process.
 * @returns a structure to pass as the final parameter of @see gsc_load_genotypefile or @see gsc_load_data_files
 */
gsc_FileFormatSpec gsc_define_matrix_format_details(const int has_header, const int markers_as_rows, const enum gsc_GenotypeFileCellStyle cell_style) {
    return (FileFormatSpec){.filetype=GSC_GENOTYPEFILE_MATRIX,
                            .spec={(struct gsc_GenotypeFile_MatrixFormat){.cell_style=cell_style,
                                                                        .has_header=has_header,
                                                                        .markers_as_rows=markers_as_rows}}};
}

/** Loads a genotype file, with or without existing genome model in the SimData
 *
 * First, it reads the first row to determine if it is a header or not. This,
 * and early information from the second row, can be used to determine:
 * - the number of columns in the file
 * - whether markers are stored as rows and parent lines as columns, or vice-versa
 * - how the pair of alleles for each SNP is encoded in the body of the matrix @see gsc_GenotypeFileCellStyle
 *
 * Rows with more entries than the first will have their excess entries ignored.
 * Rows with fewer entries than the first will also provide a warning that there is an
 * incorrect number of cells in that row, and will simply have default null '\0' alleles
 * in that marker/parent line.
 *
 * @return group number for the set of genotypes loaded from the file
 */
static gsc_GroupNum gsc_load_genotypefile_matrix(gsc_SimData* d, const char* filename, const gsc_FileFormatSpec format) {
    if (filename == NULL) return NO_GROUP;
    if (format.filetype != GSC_GENOTYPEFILE_MATRIX && format.filetype != GSC_GENOTYPEFILE_UNKNOWN) {
        warning("Non-genotype-matrix format specification provided to genotype matrix file loader function\n");
        return NO_GROUP;
    }

    // Part 1: Detect file formatting details
    struct gsc_GenotypeFile_MatrixFormat format_detected = 
        { .has_header = GSC_UNINIT, .markers_as_rows = GSC_UNINIT, .cell_style = GSC_GENOTYPECELLSTYLE_UNKNOWN };
    if (format.filetype == GSC_GENOTYPEFILE_MATRIX) {
        format_detected = format.spec.matrix;
    }
    unsigned int queuesize = 0;

    gsc_TableFileReader tbl = gsc_tablefilereader_create(filename);
    // Read one row + 2 cells (if possible)
    GSC_CREATE_BUFFER(cellsread,gsc_TableFileCell,100);
    unsigned int ncellsread = 0;
    do {
        cellsread[ncellsread] = gsc_tablefilereader_get_next_cell(&tbl);
        gsc_tablefilecell_deep_copy(&cellsread[ncellsread]);
        ++ncellsread;
        if (ncellsread >= cellsreadcap) {
            GSC_STRETCH_BUFFER(cellsread,2*ncellsread);
        }
    } while (!cellsread[ncellsread-1].eof && (ncellsread <= 1 || !cellsread[ncellsread-1].predNewline));
    unsigned int ncellsfirstrow = (cellsread[ncellsread-1].eof && cellsread[ncellsread-1].cell_len > 0) ? ncellsread : ncellsread - 1;
    if (!cellsread[ncellsread-1].eof) { // read one more cell if possible
        cellsread[ncellsread] = gsc_tablefilereader_get_next_cell(&tbl);
        gsc_tablefilecell_deep_copy(&cellsread[ncellsread]);
        ++ncellsread;
        if (ncellsread >= cellsreadcap) {
            GSC_STRETCH_BUFFER(cellsread,2*ncellsread);
        }
    }
    queuesize = ncellsread; // so that we know how many to free if we failure_exit
    if (ncellsread <= 1) { // file is an EOF only
        goto failure_exit;
    }
    //int is_onecol_file = cellsread[ncellsfirstrow + 1].predNewline > 0 || ncellsread == 2; // ncellsread == 2 means we read one cell, then an EOF
    int is_onerow_file = ncellsread == ncellsfirstrow || cellsread[ncellsfirstrow].eof; // short-circuiting essential!

    format_detected = gsc_helper_genotypefile_matrix_detect_orientation(d, cellsread, ncellsfirstrow, ncellsread, format_detected, filename);
    format_detected = gsc_helper_genotypefile_matrix_detect_cellstyle(cellsread, ncellsfirstrow, ncellsread, format_detected, filename);
    format_detected = gsc_helper_genotypefile_matrix_detect_header(cellsread, ncellsfirstrow, ncellsread, format_detected, filename);
    if ((format_detected.has_header != GSC_FALSE && format_detected.has_header != GSC_TRUE) || 
        (format_detected.markers_as_rows != GSC_FALSE && format_detected.markers_as_rows != GSC_TRUE) ||
        format_detected.cell_style == GSC_GENOTYPECELLSTYLE_UNKNOWN) {
            goto failure_exit;
    }

    int format_has_corner_cell = GSC_UNINIT;
    // If markers as columns, we do need to know how many cells are in the second row in order to detect a corner cell
    if (!format_detected.markers_as_rows && !is_onerow_file) {
        // Read rest of second row
        while (!cellsread[ncellsread-1].eof && !cellsread[ncellsread-1].predNewline) {
            cellsread[ncellsread] = gsc_tablefilereader_get_next_cell(&tbl);
            gsc_tablefilecell_deep_copy(&cellsread[ncellsread]);
            ++ncellsread;
            if (ncellsread >= cellsreadcap) {
                GSC_STRETCH_BUFFER(cellsread,2*ncellsread);
            }
        }
        // Detect corner cell
        queuesize = ncellsread; // so that we know how many to free if we failure_exit
        unsigned int ncellssecondrow = ncellsread - ncellsfirstrow - 1;
        format_has_corner_cell = gsc_helper_genotypefile_matrix_detect_cornercell_presence(ncellsfirstrow, ncellssecondrow, cellsread[ncellsfirstrow].predCol > 0);
        if (format_has_corner_cell == GSC_UNINIT) {
            warning( "(Loading %s) Failure: Header row length and second row length do not align\n", filename);
            goto failure_exit;
        }
    }

    // Create the queue of cells to parse (exclude header from this queue, because it needs to be dealt with differently)
    gsc_TableFileCell* cellqueue = cellsread;
    //queuesize = ncellsread; (already done above)
    if (format_detected.has_header) {
        cellqueue = cellsread + ncellsfirstrow;
        queuesize = ncellsread - ncellsfirstrow;
    }

    // PART 2: Create uniform-spaced map, if we have no map currently
    int build_map_from_rows = GSC_FALSE;
    if (d->genome.n_markers == 0) {
        if (format_detected.markers_as_rows) {
            build_map_from_rows = GSC_TRUE;
            // We're going to have to do an independent read of the file to extract these. Will be a bit slower.
            gsc_TableFileReader tbl2 = gsc_tablefilereader_create(filename);
            gsc_TableFileCell cell = gsc_tablefilereader_get_next_cell(&tbl2);
            unsigned int nmarkersread = format_detected.has_header ? 0 : 1;
            do {
             cell = gsc_tablefilereader_get_next_cell(&tbl2);
                if (cell.predNewline) { ++nmarkersread; }
            } while (!cell.eof);
            gsc_tablefilereader_close(&tbl2);
            if (cell.predNewline) { // there's a newline before eof, so no real actual last row
                --nmarkersread;
            }

            d->genome.n_markers = nmarkersread;
            d->genome.marker_names = gsc_malloc_wrap(sizeof(*d->genome.marker_names)*d->genome.n_markers, GSC_TRUE);
            d->genome.names_alphabetical = gsc_malloc_wrap(sizeof(*d->genome.names_alphabetical)*d->genome.n_markers, GSC_TRUE);
            gsc_create_uniformspaced_recombmap(d,0,NULL,d->genome.n_markers);
            
        } else { // markers as columns
            if (!format_detected.has_header) { // you should not be able to get here. // assert(format_detected.has_header == GSC_TRUE);
                warning( "(Loading %s) Failure: Genotype matrix with markers as columns but no header row is an unsupported file type (there is no way to tell which column is which marker)\n", filename);
                goto failure_exit;
            }
            
            unsigned int i = format_has_corner_cell ? 1 : 0; // starting index for iterating through names
            d->genome.n_markers = ncellsfirstrow - i;
            d->genome.marker_names = gsc_malloc_wrap(sizeof(*d->genome.marker_names)*d->genome.n_markers, GSC_TRUE);
            d->genome.names_alphabetical = gsc_malloc_wrap(sizeof(*d->genome.names_alphabetical)*d->genome.n_markers, GSC_TRUE);
            for (unsigned int j = 0; j < d->genome.n_markers; ++i, ++j) {
                //gsc_tablefilecell_deep_copy(&cellqueue[i]); // already deep copied
                d->genome.marker_names[j] = cellsread[i].cell;
                cellsread[i].isCellShallow = GSC_TRUE; // prevent deletion
                d->genome.names_alphabetical[j] = &d->genome.marker_names[j];
            }
            qsort(d->genome.names_alphabetical,d->genome.n_markers,sizeof(*d->genome.names_alphabetical),gsc_helper_indirect_alphabetical_str_comparer);
            gsc_create_uniformspaced_recombmap(d,0,NULL,d->genome.n_markers); // create based on the markers we've saved in 'genome'
        }
    }

    // PART 3: Parse file into an AlleleMatrix
    
    gsc_GroupNum group = gsc_get_new_group_num(d);
    struct gsc_EmptyListNavigator it = gsc_create_emptylistnavigator(d, group);
    unsigned int nvalidmarker = 0;
    unsigned int n_cols = 0;
    if (format_detected.markers_as_rows) {

        gsc_GenoLocation loc;
        gsc_TableFileCell ncell;
        n_cols = (format_detected.has_header) ? ncellsfirstrow + 1 : ncellsfirstrow; // assume first row has no corner cell for now
        int first = GSC_TRUE; 
        int have_valid_marker = GSC_FALSE; unsigned int markerix;
        unsigned int column = 0;
        unsigned int row = 0;
        do {
            ncell = gsc_helper_tablefilereader_get_next_cell_wqueue(&tbl,&cellqueue,&queuesize);

            if (ncell.cell != NULL) {
                if (ncell.predNewline || first) {
                    char tmp = ncell.cell[ncell.cell_len]; ncell.cell[ncell.cell_len] = '\0';
                    
                    if (build_map_from_rows) {
                        have_valid_marker = GSC_TRUE;
                        if (first) {
                            markerix = 0;
                        } else {
                            markerix++;
                        }
                        gsc_tablefilecell_deep_copy(&ncell);
                        d->genome.marker_names[markerix] = ncell.cell;
                        ncell.isCellShallow = GSC_TRUE; // prevent deletion
                    } else {
                        have_valid_marker = gsc_get_index_of_genetic_marker(ncell.cell, d->genome, &markerix);
                    }
                    
                    nvalidmarker += have_valid_marker;
                    ncell.cell[ncell.cell_len] = tmp;
                    // Then, after reading first row, detect what our expected row length is, if defaults don't suit.
                    if (row == 1 && format_detected.has_header) { 
                        if (column + 1 != ncellsfirstrow && column + 1 != ncellsfirstrow + 1) {
                            warning( "(Loading %s) Failure: Header row length and second row length do not align\n", filename);
                            goto failure_exit;
                        } else {
                            n_cols = column + 1;
                        }
                    }
                    first = GSC_FALSE;
                    column = 0;
                    ++row;

                } else if (ncell.predCol) { // any number of column spacers treated as one column gap when reading a genotype matrix
                    ++column;
                    if (have_valid_marker && column < n_cols) {
                        loc = (1 == column) ? gsc_emptylistnavigator_get_first(&it) : gsc_emptylistnavigator_get_next(&it);
                        gsc_helper_genotypecell_to_allelematrix(loc,markerix,format_detected.cell_style,ncell.cell,d);
                    } // Note we ignore all extra cells in all rows
                }
            }

            if (!ncell.isCellShallow) { GSC_FREE(ncell.cell); }
        } while (!ncell.eof);
        if (row == 1 && format_detected.has_header) {
            if (column + 1 != ncellsfirstrow && column + 1 != ncellsfirstrow + 1) {
                warning( "(Loading %s) Failure: Header row length and second row length do not align\n", filename);
                goto failure_exit;
            } else {
                n_cols = column + 1;
            }
        }

        // Then save the genotype names
        if (format_detected.has_header) {
            format_has_corner_cell = gsc_helper_genotypefile_matrix_detect_cornercell_presence(ncellsfirstrow, n_cols, cellsread[ncellsfirstrow].predCol > 0);
            unsigned int i = format_has_corner_cell ? 1 : 0;
            gsc_GenoLocation loc;
            for (unsigned int j = 0; i < ncellsfirstrow; ++i, ++j) {
                loc = (j == 0) ? gsc_emptylistnavigator_get_first(&it) : gsc_emptylistnavigator_get_next(&it);
                // assert(!cellsread[i].isShallowCopy);
                gsc_set_name(loc,cellsread[i].cell); // using names here so no need to free them. Since they're in cellsread 
                cellsread[i].isCellShallow = GSC_TRUE; // prevent deletion
            }
        }

        // Then finalise the map, if we're creating one:
        if (build_map_from_rows) {
            for (unsigned int i = 0; i < d->genome.n_markers; ++i) {
                d->genome.names_alphabetical[i] = &d->genome.marker_names[i];
            }
            qsort(d->genome.names_alphabetical,d->genome.n_markers,sizeof(*d->genome.names_alphabetical),gsc_helper_indirect_alphabetical_str_comparer);
        }
        
    } else { // markers as columns
        if (!format_detected.has_header) { // you should not be able to get here.
            warning( "(Loading %s) Failure: Genotype matrix with markers as columns but no header row is an unsupported file type (there is no way to tell which column is which marker)\n", filename);
            goto failure_exit;
        }

        // Identify the marker corresponding to each column
        unsigned int i = format_has_corner_cell ? 1 : 0;
        unsigned int n_col = ncellsfirstrow + (1-i);
        unsigned int* markerixs = gsc_malloc_wrap(sizeof(*markerixs)*ncellsfirstrow,GSC_TRUE);
        for (unsigned int j = 0; i < ncellsfirstrow; ++i, ++j) {
            markerixs[j] = d->genome.n_markers;
            nvalidmarker += gsc_get_index_of_genetic_marker(cellsread[i].cell, d->genome, &markerixs[j]);
        }

        // Read the table
        int first = GSC_TRUE; 
        unsigned int row = 0;
        unsigned int column = 0; // we count column numbers from 1 for the first body cell. sorry for the inconsistency with the branch of the if statement above.
        gsc_GenoLocation loc = gsc_emptylistnavigator_get_first(&it);
        gsc_TableFileCell ncell;
        do {
            ncell = gsc_helper_tablefilereader_get_next_cell_wqueue(&tbl,&cellqueue,&queuesize);

            if (ncell.cell != NULL) {
                if (ncell.predNewline) {
                    loc = (first) ? gsc_emptylistnavigator_get_first(&it) : gsc_emptylistnavigator_get_next(&it);
                    first = GSC_FALSE;
                    
                    ++row;
                    column = 0;
                    if (ncell.predCol) { // missing name.
                        gsc_set_name(loc,NULL);
                    } else {
                        gsc_tablefilecell_deep_copy(&ncell);
                        gsc_set_name(loc,ncell.cell);
                        ncell.isCellShallow = GSC_TRUE; // so it does not get deleted
                    }
                }

                if (ncell.predCol) {
                    ++column;
                    if (column < n_col && markerixs[column-1] < d->genome.n_markers) {
                        gsc_helper_genotypecell_to_allelematrix(loc,markerixs[column-1],format_detected.cell_style,ncell.cell,d);
                    }
                }
            }

            if (!ncell.isCellShallow) { GSC_FREE(ncell.cell); } 
        } while (!ncell.eof);

        GSC_FREE(markerixs);
        
    }

    // PART 4: Tidy and clean and exit
    unsigned int ngenos = 0;
    AlleleMatrix* tmpam = it.firstAM;
    do {
        ngenos += tmpam->n_genotypes;
    } while ((tmpam = tmpam->next) != NULL);
    Rprintf("(Loading %s) %u genotype(s) of %u marker(s) were loaded.\n", filename, (unsigned int) ngenos, (unsigned int) nvalidmarker);
	if (ngenos == 0) {
		gsc_delete_allele_matrix(it.firstAM);
		goto failure_exit;
	}
	gsc_emptylistnavigator_finaliselist(&it);
    ++d->n_groups;

    // ... cleaning up the header row
    if (format_detected.has_header) {
        for (int j = 0; j < ncellsfirstrow; ++j) { 
            if (!cellsread[j].isCellShallow) { GSC_FREE(cellsread[j].cell); }
        }
    }
    GSC_DELETE_BUFFER(cellsread);
    gsc_tablefilereader_close(&tbl);
    return group;

    failure_exit:
        // Clean up structures and return, having loaded no genotypes
        // ... cleaning up unprocessed cells in the queue
        for (int i = 1; i <= queuesize; ++i) { 
            if (!cellsread[ncellsread-i].isCellShallow) {
                GSC_FREE(cellsread[ncellsread-i].cell); 
                cellsread[ncellsread-i].isCellShallow = GSC_TRUE;
            }
        }
        // ... cleaning up the header row
        if (format_detected.has_header) {
            for (int j = 0; j < ncellsfirstrow; ++j) { 
                if (!cellsread[j].isCellShallow) { GSC_FREE(cellsread[j].cell); }
            }
        }
        GSC_DELETE_BUFFER(cellsread);
        gsc_tablefilereader_close(&tbl);
        return NO_GROUP;
}

/** Load a set of genotypes to a gsc_SimData object.
 *
 * Can be called on an empty SimData, in which case it invents a genetic map
 * comprising of the markers present in the genotype file, all on a single
 * chromosome and spaced 1cM apart, or on a SimData with existing maps or genotypes,
 * in which case markers whose names do not appear in the SimData's list of tracked
 * markers are discarded.
 *
 * @see enum gsc_GenotypeFileType
 * for the types of files that can be loaded
 *
 * @shortnamed{load_genotypefile}
 *
 * @param d pointer to gsc_SimData to be populated
 * @param filename string name/path of file containing genotype data.
 * @param format details on the format of the file (option, set to GSC_DETECT_FILE_FORMAT to auto-detect format).
*/
gsc_GroupNum gsc_load_genotypefile(SimData* d, const char* filename, const gsc_FileFormatSpec format) {
    return gsc_load_data_files(d,filename,NULL,NULL,format).group;
}

/** Populates a gsc_SimData object with marker allele data, a genetic map, and
 * (optionally) marker effect values.
 *
 * This is the suggested first function to call to set up a genomicSimulation simulation.
 *
 * It will attempt to use file extensions to determine how to parse the files
 *
 * @shortnamed{load_data_files}
 *
 * @param d pointer to gsc_SimData to be populated
 * @param data_file string containing name/path of file containing SNP marker
 * allele data.
 * @param map_file string name/path of file containing genetic map data (optional, set to NULL if not needed).
 * @param effect_file string name/path of file containing effect values (optional, set to NULL if not wanted).
 * @param format details on the format of the file (option, set to GSC_DETECT_FILE_FORMAT to auto-detect format).
 * @returns a @a gsc_MultiIDSet entry, containing the group number of the founding group, the map id of the recombination
 * map loaded, and the effect set id of the loaded effect file if applicable.
*/
struct gsc_MultiIDSet gsc_load_data_files(gsc_SimData* d, const char* genotype_file,
                                             const char* map_file, const char* effect_file, const gsc_FileFormatSpec format) {
    // Parse file suffix for file type, if it was not already provided
    enum gsc_GenotypeFileType type = format.filetype;

    if (type == GSC_GENOTYPEFILE_UNKNOWN) {
        type = GSC_GENOTYPEFILE_MATRIX;
        char* suffix = strrchr(genotype_file,'.');
        if (suffix != NULL) {
            if (strcmp(suffix,".bed") == 0) {
                type = GSC_GENOTYPEFILE_BED;
            } else if (strcmp(suffix,".ped") == 0) {
                type = GSC_GENOTYPEFILE_PED;
            } else if (strcmp(suffix,".vcf") == 0) {
                type = GSC_GENOTYPEFILE_VCF;
            }
        }
    }

    struct gsc_MultiIDSet out = { .group=NO_GROUP, .map=NO_MAP, .effSet=NO_EFFECTSET };

    switch (type) {
    case GSC_GENOTYPEFILE_BED:
        //if (detectedtype) { Rprintf("Will attempt to parse %s as a plink .bed file\n", filename); }
        warning("plink .bed file parsing not yet implemented\n");
        break;
    case GSC_GENOTYPEFILE_PED:
        warning("plink .ped file parsing not yet implemented\n");
        break;
    case GSC_GENOTYPEFILE_VCF:
        warning("vcf file parsing not yet implemented\n");
        break;
    default:
        //Rprintf("(Loading files) Will treat %s as a genotype matrix (see genomicSimulation's default input file types)\n", genotype_file);
        out.map = gsc_load_mapfile(d, map_file);
        out.group = gsc_load_genotypefile_matrix(d, genotype_file, format);
        out.effSet = gsc_load_effectfile(d, effect_file);
    }

    return out;
}

/*--------------------------Recombination counts-----------------------------*/

/** Identify markers in the genotype of `offspring` where recombination from its parents
 * occured. This function is a little lower-level (see the kinds of parameters required) and
 * so a wrapper like gsc_calculate_recombinations_from_file is suggested for end users.
 * @see gsc_calculate_recombinations_from_file()
 *
 * The function reads start to end along each chromosome. At each marker, it checks if
 * the alleles the offspring has could only have come from one parent/there is known parentage
 * of that allele. If that is the case, it saves the provided id number of the source parent
 * to the matching position in the result vector. If it is not the case, its behaviour depends
 * on the `certain` parameter.
 *
 * Parents do not have to be directly identified as parents by the pedigree functionality of
 * this library. A sample usage is performing a cross then multiple generations of selfing,
 * then comparing the final inbred to the original two lines of the cross.
 *
 * @param d pointer to the gsc_SimData struct whose genetic map matches the provided genotypes.
 * @param mapid ID of the map from which to calculate potential historical crossovers, or NO_MAP
 * to use the first-loaded/primary map by default
 * @param parent1 a character vector containing one parent's alleles at each marker in the
 * gsc_SimData.
 * @param p1num an integer that will be used to identify areas of the genome that come
 * from the first parent in the returned vector.
 * @param parent2 a character vector containing the other parent's alleles at each marker in the
 * gsc_SimData.
 * @param p2num an integer that will be used to identify areas of the genome that come
 * from the second parent in the returned vector.
 * @param offspring a character vector containing the alleles at each marker in the
 * gsc_SimData of the genotype whose likely recombinations we want to identify.
 * @param certain a boolean. If TRUE, markers where the parent of origin cannot be identified
 * will be set to 0, if FALSE, the value will be set to the id of the parent that provided
 * the most recently identified allele in that chromosome.
 * @returns a heap vector of length `d->n_markers` containing the id of the parent of origin
 * at each marker in the `offspring` genotype.
*/
int* gsc_calculate_min_recombinations_fw1(gsc_SimData* d, gsc_MapID mapid, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int certain) {
    if (d->genome.n_maps < 1) {
        warning("Need at least one recombination map loaded to estimate recombinations\n");
        return NULL;
    }
    int mapix = 0;
    if (mapid.id != NO_MAP.id) { mapix = gsc_get_index_of_map(d, mapid); }
    if (mapix >= d->genome.n_maps) {
        warning("We don't have that recombination maps loaded\n");
        return NULL;
    }
    gsc_RecombinationMap map = d->genome.maps[mapix];

    int* origins = gsc_malloc_wrap(sizeof(int) * d->genome.n_markers,GSC_TRUE);
    memset(origins,0,sizeof(*origins)*d->genome.n_markers);
    int p1match, p2match;
    int previous = 0;


    for (int chr = 0; chr < map.n_chr; ++chr) {
        R_CheckUserInterrupt();

        switch (map.chrs[chr].type) {
        case GSC_LINKAGEGROUP_SIMPLE:
            for (int i = 0; i < map.chrs[chr].map.simple.n_markers; ++i) {
                p1match = gsc_has_same_alleles(parent1, offspring, i);
                p2match = gsc_has_same_alleles(parent2, offspring, i);
                if (p1match && !p2match) {
                    origins[map.chrs[chr].map.simple.first_marker_index + i] = p1num;
                    previous = p1num;
                } else if (p2match && !p1match) {
                    origins[map.chrs[chr].map.simple.first_marker_index + i] = p2num;
                    previous = p2num;
                } else {
                    if (certain) {
                        origins[map.chrs[chr].map.simple.first_marker_index + i] = 0;
                    } else {
                        origins[map.chrs[chr].map.simple.first_marker_index + i] = previous;
                    }
                }
            }
            break;

        case GSC_LINKAGEGROUP_REORDER:
            for (int i = 0; i < map.chrs[chr].map.reorder.n_markers; ++i) {
                p1match = gsc_has_same_alleles(parent1, offspring, i);
                p2match = gsc_has_same_alleles(parent2, offspring, i);
                if (p1match && !p2match) {
                    origins[map.chrs[chr].map.reorder.marker_indexes[i]] = p1num;
                    previous = p1num;
                } else if (p2match && !p1match) {
                    origins[map.chrs[chr].map.reorder.marker_indexes[i]] = p2num;
                    previous = p2num;
                } else {
                    if (certain) {
                        origins[map.chrs[chr].map.reorder.marker_indexes[i]] = 0;
                    } else {
                        origins[map.chrs[chr].map.reorder.marker_indexes[i]] = previous;
                    }
                }
            }
            break;
        }

    }
	return origins;
}

/** Identify markers in the genotype of `offspring` where recombination from its parents
 * occured, as judged by the marker itself and a short window around it.
 * This function is a little lower-level (see the kinds of parameters required) and
 * so a wrapper like gsc_calculate_recombinations_from_file is suggested for end users.
 * @see gsc_calculate_recombinations_from_file()
 *
 * The function reads start to end along each chromosome. At each marker, it checks if
 * the alleles the offspring has in the window centered at that marker could have come
 * from one parent but could not have come from the other/there is known parentage
 * of that allele. If that is the case, it saves the provided id number of the source parent
 * to the matching position in the result vector. If it is not the case, its behaviour depends
 * on the `certain` parameter.
 *
 * Parents do not have to be directly identified as parents by the pedigree functionality of
 * this library. A sample usage is performing a cross then multiple generations of selfing,
 * then comparing the final inbred to the original two lines of the cross.
 *
 * Behaviour when the window size is not an odd integer has not been tested.
 *
 * @param d pointer to the gsc_SimData struct whose genetic map matches the provided genotypes.
 * @param mapid ID of the map from which to calculate potential historical crossovers, or NO_MAP
 * to use the first-loaded/primary map by default
 * @param parent1 a character vector containing one parent's alleles at each marker in the
 * gsc_SimData.
 * @param p1num an integer that will be used to identify areas of the genome that come
 * from the first parent in the returned vector.
 * @param parent2 a character vector containing the other parent's alleles at each marker in the
 * gsc_SimData.
 * @param p2num an integer that will be used to identify areas of the genome that come
 * from the second parent in the returned vector.
 * @param offspring a character vector containing the alleles at each marker in the
 * gsc_SimData of the genotype whose likely recombinations we want to identify.
 * @param window_size an odd integer representing the number of markers to check for known parentage
 * around each marker
 * @param certain a boolean. If TRUE, markers where the parent of origin cannot be identified
 * will be set to 0, if FALSE, the value will be set to the id of the parent that provided
 * the most recently identified allele in that chromosome.
 * @returns a heap vector of length `d->n_markers` containing the id of the parent of origin
 * at each marker in the `offspring` genotype.
*/
int* gsc_calculate_min_recombinations_fwn(gsc_SimData* d, gsc_MapID mapid, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int window_size, int certain) {
    if (d->genome.n_maps < 1) {
        warning("Need at least one recombination map loaded to estimate recombinations\n");
        return NULL;
    }
    int mapix = 0;
    if (mapid.id != NO_MAP.id) { mapix = gsc_get_index_of_map(d, mapid); }
    if (mapix >= d->genome.n_maps) {
        warning("We don't have that recombination maps loaded\n");
        return NULL;
    }
    gsc_RecombinationMap map = d->genome.maps[mapix];


    int* origins = gsc_malloc_wrap(sizeof(int) * d->genome.n_markers,GSC_TRUE);
    memset(origins,0,sizeof(*origins)*d->genome.n_markers);
    int p1match, p2match;
    int previous = 0, window_range = (window_size - 1)/2, i;

    for (int chr = 0; chr < map.n_chr; ++chr) {
        R_CheckUserInterrupt();

        switch (map.chrs[chr].type) {
        case GSC_LINKAGEGROUP_SIMPLE:
            for (i = 0; i < window_range; ++i) {
                origins[map.chrs[chr].map.simple.first_marker_index + i] = 0;
            }
            for (; i < map.chrs[chr].map.simple.n_markers - window_range; ++i) {
                p1match = gsc_has_same_alleles_window(parent1, offspring, i, window_size);
                p2match = gsc_has_same_alleles_window(parent2, offspring, i, window_size);
                if (p1match && !p2match) {
                    origins[map.chrs[chr].map.simple.first_marker_index + i] = p1num;
                    previous = p1num;
                } else if (p2match && !p1match) {
                    origins[map.chrs[chr].map.simple.first_marker_index + i] = p2num;
                    previous = p2num;
                } else {
                    if (certain) {
                        origins[map.chrs[chr].map.simple.first_marker_index + i] = 0;
                    } else {
                        origins[map.chrs[chr].map.simple.first_marker_index + i] = previous;
                    }
                }
            }
            for (; i < map.chrs[chr].map.simple.n_markers; ++i) {
                origins[map.chrs[chr].map.simple.first_marker_index + i] = 0;
            }
            break;

        case GSC_LINKAGEGROUP_REORDER:
            for (i = 0; i < window_range; ++i) {
                origins[map.chrs[chr].map.reorder.marker_indexes[i]] = 0;
            }
            for (; i < map.chrs[chr].map.reorder.n_markers - window_range; ++i) {
                p1match = gsc_has_same_alleles_window(parent1, offspring, i, window_size);
                p2match = gsc_has_same_alleles_window(parent2, offspring, i, window_size);
                if (p1match && !p2match) {
                    origins[map.chrs[chr].map.reorder.marker_indexes[i]] = p1num;
                    previous = p1num;
                } else if (p2match && !p1match) {
                    origins[map.chrs[chr].map.reorder.marker_indexes[i]] = p2num;
                    previous = p2num;
                } else {
                    if (certain) {
                        origins[map.chrs[chr].map.reorder.marker_indexes[i]] = 0;
                    } else {
                        origins[map.chrs[chr].map.reorder.marker_indexes[i]] = previous;
                    }
                }
            }
            for (; i < map.chrs[chr].map.reorder.n_markers; ++i) {
                origins[map.chrs[chr].map.reorder.marker_indexes[i]] = 0;
            }
            break;
        }

    }
	return origins;
}

/** Provides guesses as to the location of recombination events that led to the
 * creation of certain genotypes from certain other genotypes.
 *
 * The input file (which pairs up which targets and their parents the calculation
 * should be carried out on) should have format:
 *
 * [target name]	[parent1name]	[parent2name]
 *
 * [target name]	[parent1name]	[parent2name]
 *
 * ...
 *
 * The tab-separated output file produced by this function will have format:
 *
 * 	[marker 1 name]	[marker 2 name]...
 *
 * [target name]	[tab-separated recombination vector, containing the index at
 * each marker of the parent the function guesses the target's alleles came from, or
 * 0 if this is unknow]
 *
 * ...
 *
 * Parents do not have to be directly identified as parents by the pedigree functionality of
 * this library. A sample usage is performing a cross then multiple generations of selfing,
 * then comparing the final inbred to the original two lines of the cross.
 *
 * @param d pointer to the gsc_SimData struct containing the genotypes and map under consideration.
 * @param input_file string containing the name of the file with the pairs of parents
 * and offsprings of which to calculate recombinations
 * @param output_file string containing the filename to which to save the results.
 * @param window_len an odd integer representing the number of markers to check for known parentage
 * around each marker
 * @param certain TRUE to fill locations where parentage is unknown with 0, FALSE
 * to fill locations where parentage is unknown with the most recent known parent
 * @returns 0 on success.
 */
int gsc_calculate_recombinations_from_file(gsc_SimData* d, const char* input_file, const char* output_file,
		int window_len, int certain) {
	struct gsc_TableSize t = gsc_get_file_dimensions(input_file, '\t');
	//open file
	FILE* fp;
	if ((fp = fopen(input_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", input_file);
	}
	FILE* fpo;
	if ((fpo = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}

	// print header.
    for (int j = 0; j < d->genome.n_markers; ++j) {
        fprintf(fpo, "\t%s", d->genome.marker_names[j]);
	}

	int combin_i[3];
	char* combin_genes[3];
	char buffer[3][50];
	int* r;
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s %s \n", buffer[0], buffer[1], buffer[2]);
		combin_i[0] = gsc_get_index_of_name(d->m, buffer[0]);
		combin_i[1] = gsc_get_index_of_name(d->m, buffer[1]);
		combin_i[2] = gsc_get_index_of_name(d->m, buffer[2]);
		combin_genes[0] = gsc_get_genes_of_index(d->m, combin_i[0]);
		combin_genes[1] = gsc_get_genes_of_index(d->m, combin_i[1]);
		combin_genes[2] = gsc_get_genes_of_index(d->m, combin_i[2]);

		if (window_len == 1) {
            r = gsc_calculate_min_recombinations_fw1(d, NO_MAP, combin_genes[1],
                    gsc_get_id_of_index(d->m, combin_i[1]).id, combin_genes[2],
                    gsc_get_id_of_index(d->m, combin_i[2]).id, combin_genes[0], certain);
		} else {
            r = gsc_calculate_min_recombinations_fwn(d, NO_MAP, combin_genes[1],
                    gsc_get_id_of_index(d->m, combin_i[1]).id, combin_genes[2],
                    gsc_get_id_of_index(d->m, combin_i[2]).id, combin_genes[0], window_len, certain);
		}
		fprintf(fpo, "\n%s", buffer[0]);
        for (int j = 0; j < d->genome.n_markers; ++j) {
			fprintf(fpo, "\t%d", r[j]);
		}
		GSC_FREE(r);
	}

	fclose(fp);
	fwrite("\n", sizeof(char), 1, fpo);
	fflush(fpo);
	fclose(fpo);
	return 0;
}


/*--------------------------------Crossing-----------------------------------*/

/** Fills a char* with the simulated result of meiosis (reduction and
 * recombination) from the marker alleles of a given parent.
 *
 * It generates the number of crossover events in each chromosome by drawing
 * from a Poisson distribution with a certain expected number of crossovers,
 * where the expected number of crossovers is saved in the chosen recombination
 * map.
 *
 * It generates the positions of those crossover events from a uniform distribution.
 *
 * It picks a random one of the gametes generated by picking randomly which column
 * of the parent's alleles to start with. When crossover events occur it starts reading
 * the other column.
 *
 * @shortnamed{generate_gamete}
 *
 * @param d pointer to the gsc_SimData object containing map positions for the markers
 * that make up the rows of `parent_table`. sort_markers() and locate_chromosomes()
 * should have been called previously.
 * @param parent_genome the char* containing the parent's genome as a character string
 * made up of sequential pairs of alleles for each marker in d->markers.
 * @param output the char* to which to save the gamete. It saves the alleles every second
 * character, starting at 0, so that calling gsc_generate_gamete(..., offspring_genome) &
 * gsc_generate_gamete(..., offspring_genome + 1) can be used to generate both halves of its genome.
 * @param map_index index of the recombination map in `d` that should
 * be used to choose positions and probabilities of recombination. @see gsc_get_index_of_map
*/
void gsc_generate_gamete(gsc_SimData* d, const char* parent_genome, char* output, const unsigned int map_index) {
	// assumes rand is already seeded
	if (parent_genome == NULL) {
        warning( "Could not generate this gamete: no parent provided\n");
		return;
	}
    if (map_index >= d->genome.n_maps) {
        warning( "Could not generate this gamete: invalid map provided\n");
        return;
    }
    gsc_RecombinationMap map = d->genome.maps[map_index];

	// treat each chromosome separately.
    GSC_CREATE_BUFFER(crossover_where, double, 100);
    for (int chr = 0; chr < d->genome.maps[map_index].n_chr; ++chr) {

        // Task 1: How many crossovers
        int num_crossovers;
        switch (map.chrs[chr].type) {
            case GSC_LINKAGEGROUP_SIMPLE:
                num_crossovers = Rf_rpois( map.chrs[chr].map.simple.expected_n_crossovers);
                break;
            case GSC_LINKAGEGROUP_REORDER:
                num_crossovers = Rf_rpois( map.chrs[chr].map.reorder.expected_n_crossovers);
                break;
            default:
                warning( "Linkage group type of group with index %i of map with index %ui is corrupted.\n", chr, map_index);
				num_crossovers = 0;
        }

        // Task 2: Find positions of all crossovers
        if (num_crossovers > crossover_wherecap) {
            GSC_STRETCH_BUFFER(crossover_where,num_crossovers);
        }
        for (int i = 0; i < num_crossovers; ++i) {
            crossover_where[i] = unif_rand();
        }
        if (num_crossovers > 1) {
            qsort(crossover_where, num_crossovers, sizeof(double), gsc_helper_ascending_double_comparer);
        }

        // Task 3: Read off the gamete that those crossovers produce.
        int which = (unif_rand() > 0.5); // if this is 0, we start with the left haplotype
        int up_to_crossover = 0;
        switch (map.chrs[chr].type) {
            case GSC_LINKAGEGROUP_SIMPLE:
                for (int i = 0; i < map.chrs[chr].map.simple.n_markers; ++i) {
                    // is it time to invert which parent haplotype we're reading?
                    while (up_to_crossover < num_crossovers &&
                           map.chrs[chr].map.simple.dists[i] > crossover_where[up_to_crossover]) {
                        which = 1 - which;
                        up_to_crossover++;
                    }
                    output[2*(i + map.chrs[chr].map.simple.first_marker_index)] =
                            parent_genome[2*(i + map.chrs[chr].map.simple.first_marker_index) + which];
                }
                break;
            case GSC_LINKAGEGROUP_REORDER:
                for (int i = 0; i < map.chrs[chr].map.reorder.n_markers; ++i) {
                    // is it time to invert which parent haplotype we're reading?
                    while (up_to_crossover < num_crossovers &&
                           map.chrs[chr].map.reorder.dists[i] > crossover_where[up_to_crossover]) {
                        which = 1 - which;
                        up_to_crossover++;
                    }
                    output[2*map.chrs[chr].map.reorder.marker_indexes[i]] =
                            parent_genome[2*map.chrs[chr].map.reorder.marker_indexes[i] + which];
                }
                break;
            default:
                break;
        }
    }
    GSC_DELETE_BUFFER(crossover_where);
}

/** Get the alleles of the outcome of producing a doubled haploid from
 * a gamete from a given parent.
 *
 * One gamete is generated, then doubled. The output will be perfectly
 * homozygous.
 *
 * @see gsc_generate_gamete(), on which this is based.
 *
 * @shortnamed{generate_doubled_haploid}
 *
 * @param d pointer to the gsc_SimData object that includes genetic map data
 * needed to simulate meiosis and the value of n_markers
 * @param parent_genome a 2x(n_markers) array of characters containing the
 * alleles of the first parent
 * @param output a 2x(n_marker) array of chars which will be overwritten
 * with the offspring genome.
 * @param map_index index of the recombination map in `d` that should
 * be used to choose positions and probabilities of recombination. @see gsc_get_index_of_map
*/
void gsc_generate_doubled_haploid(gsc_SimData* d, const char* parent_genome, char* output, const unsigned int map_index) {
    /* For cache reasons it'll be better to copy-paste gsc_generate_gamete with
     * one extra line added to the inner loop, than to generate a single gamete
     * and then scan over `output` again to copy it. */

	// assumes rand is already seeded
	if (parent_genome == NULL) {
		warning( "Could not make this doubled haploid\n");
		return;
	}
    if (map_index >= d->genome.n_maps) {
        warning( "Could not generate this gamete: invalid map provided\n");
        return;
    }
    gsc_RecombinationMap map = d->genome.maps[map_index];

    // treat each chromosome separately.
    GSC_CREATE_BUFFER(crossover_where, double, 100);
    for (int chr = 0; chr < d->genome.maps[map_index].n_chr; ++chr) {

        // Task 1: How many crossovers
        int num_crossovers;
        switch (map.chrs[chr].type) {
            case GSC_LINKAGEGROUP_SIMPLE:
                num_crossovers = Rf_rpois( map.chrs[chr].map.simple.expected_n_crossovers);
                break;
            case GSC_LINKAGEGROUP_REORDER:
                num_crossovers = Rf_rpois( map.chrs[chr].map.reorder.expected_n_crossovers);
                break;
            default:
                warning( "Linkage group type of group with index %i of map with index %ui is corrupted.\n", chr, map_index);
				num_crossovers = 0;
        }

        // Task 2: Find positions of all crossovers
        if (num_crossovers > crossover_wherecap) {
            GSC_STRETCH_BUFFER(crossover_where,num_crossovers);
        }
        for (int i = 0; i < num_crossovers; ++i) {
            crossover_where[i] = unif_rand();
        }
        if (num_crossovers > 1) {
            qsort(crossover_where, num_crossovers, sizeof(double), gsc_helper_ascending_double_comparer);
        }

        // Task 3: Read off the gamete that those crossovers produce.
        int which = (unif_rand() > 0.5); // if this is 0, we start with the left haplotype
        int up_to_crossover = 0;
        switch (map.chrs[chr].type) {
            case GSC_LINKAGEGROUP_SIMPLE:
                for (int i = 0; i < map.chrs[chr].map.simple.n_markers; ++i) {
                    // is it time to invert which parent haplotype we're reading?
                    while (up_to_crossover < num_crossovers &&
                           map.chrs[chr].map.simple.dists[i] > crossover_where[up_to_crossover]) {
                        which = 1 - which;
                        up_to_crossover++;
                    }
                    int pos = i + map.chrs[chr].map.simple.first_marker_index;
                    output[2*pos] = parent_genome[2*pos + which];
                    output[2*pos + 1] = output[2*pos];  // haploid doubling happens here
                }
                break;
            case GSC_LINKAGEGROUP_REORDER:
                for (int i = 0; i < map.chrs[chr].map.reorder.n_markers; ++i) {
                    // is it time to invert which parent haplotype we're reading?
                    while (up_to_crossover < num_crossovers &&
                           map.chrs[chr].map.reorder.dists[i] > crossover_where[up_to_crossover]) {
                        which = 1 - which;
                        up_to_crossover++;
                    }
                    int pos = map.chrs[chr].map.reorder.marker_indexes[i];
                    output[2*pos] = parent_genome[2*pos + which];
                    output[2*pos + 1] = output[2*pos]; // haploid doubling happens here
                }
                break;
            default:
                break;
        }
    }
    GSC_DELETE_BUFFER(crossover_where);
}


/** Get an identical copy of a given genotype.
*
 * @shortnamed{generate_clone}
 *
 * @param d pointer to the gsc_SimData object that includes genetic map data
 * needed to simulate meiosis and the value of n_markers
 * @param parent_genome a 2x(n_markers) array of characters containing the
 * alleles of the first parent
 * @param output a 2x(n_marker) array of chars which will be overwritten
 * with the offspring genome.
*/
void gsc_generate_clone(gsc_SimData* d, const char* parent_genome, char* output) {
    for (int j = 0; j < d->genome.n_markers; ++j) {
        output[2*j] = parent_genome[2*j];
        output[2*j + 1] = parent_genome[2*j + 1];
    }
    return;
}


/** Opens file for writing save-as-you-go pedigrees in accordance with gsc_GenOptions */
static FILE* gsc_helper_genoptions_save_pedigrees_setup(const gsc_GenOptions g) {
	FILE* fp = NULL;
	if (g.will_save_pedigree_to_file) { 					
        char tmpname_p[NAME_LENGTH];
        if (g.filename_prefix != NULL) {
            strncpy(tmpname_p, g.filename_prefix,
					sizeof(char)*(NAME_LENGTH-13));
        } else {
            strcpy(tmpname_p, "out");
        }
		strcat(tmpname_p, "-pedigree.txt");
		fp = fopen(tmpname_p, "w");	
	}
	return fp;
}
/** Opens file for writing save-as-you-go breeding values in accordance with gsc_GenOptions.
 *
 * @param effIndexp a location to save the index of the chosen effect set. 
 * It will be set to GSC_UNINIT if the effect set ID in @a g was invalid.
 */
static FILE* gsc_helper_genoptions_save_bvs_setup(const gsc_SimData* d, const gsc_GenOptions g, int* effIndexp) {
	FILE* fe = NULL;
    if (g.will_save_bvs_to_file.id != GSC_NO_EFFECTSET.id) {
        *effIndexp = gsc_get_index_of_eff_set(d,g.will_save_bvs_to_file);
        if (*effIndexp != GSC_UNINIT) {
            char tmpname_b[NAME_LENGTH];
            if (g.filename_prefix != NULL) {
                strncpy(tmpname_b, g.filename_prefix, 
						sizeof(char)*(NAME_LENGTH-7));
            } else {
                strcpy(tmpname_b, "out");
            }
            strcat(tmpname_b, "-bv.txt");
            fe = fopen(tmpname_b, "w");
        }
    }
	return fe;
}
/** Opens file for writing save-as-you-go genotypes in accordance with gsc_GenOptions */
static FILE* gsc_helper_genoptions_save_genotypes_setup(const gsc_SimData* d, const gsc_GenOptions g) {
	FILE* fg = NULL;
	if (g.will_save_alleles_to_file) {
        char tmpname_g[NAME_LENGTH];
        if (g.filename_prefix != NULL) {
            strncpy(tmpname_g, g.filename_prefix,
					sizeof(char)*(NAME_LENGTH-13));
        } else {
            strcpy(tmpname_g, "out");
        }
		strcat(tmpname_g, "-genotype.txt");
		fg = fopen(tmpname_g, "w");
        gsc_save_names_header(fg,d->genome.n_markers,(const char**)d->genome.marker_names);
	}
	return fg;
}
/** save-as-you-go (pedigrees)
 *
 * Saves pedigrees in gsc_AlleleMatrix @a tosave to a file. Performs 
 * null-checks on file pointer @a fp, and fails silently if checks do. 
 */
static void gsc_helper_genoptions_save_pedigrees(FILE* fp, gsc_SimData* d, gsc_AlleleMatrix* tosave) {
	if (fp) {
        gsc_save_allelematrix_full_pedigree( fp, tosave, d);
	}
}
/** save-as-you-go (breeding values)
 * 
 * Calculates and saves breeding values of genotypes in an gsc_AlleleMatrix 
 * to a file. Performs null-checks on file pointer and @a effIndex, 
 * and fails silently if checks do. 
 */
static void gsc_helper_genoptions_save_bvs(FILE* fe, gsc_EffectMatrix* effMatrices, int effIndex, gsc_AlleleMatrix* tosave) {
	if (fe && effIndex != GSC_UNINIT) {
        gsc_DecimalMatrix eff = gsc_calculate_bvs( tosave, effMatrices + effIndex );
        gsc_save_manual_bvs( fe, &eff, tosave->ids, (const char**) tosave->names);
		gsc_delete_dmatrix( &eff);
	}
}
/** save-as-you-go (genotypes/alleles)
 *
 * Saves genotypes in an gsc_AlleleMatrix to a file. Performs null-checks 
 * on file pointer and fails silently if checks do. 
 */
static void gsc_helper_genoptions_save_genotypes(FILE* fg, gsc_AlleleMatrix* tosave) {
	if (fg) {
		gsc_save_allele_matrix( fg, tosave );
	}
}
/** Apply gsc_GenOptions naming scheme and gsc_PedigreeID allocation to a single 
 * gsc_AlleleMatrix.
 */
static void gsc_helper_genoptions_give_names_and_ids(gsc_AlleleMatrix* am, gsc_SimData* d, const gsc_GenOptions g) {
	if (g.will_name_offspring) {
		gsc_set_names(am, g.offspring_name_prefix, d->current_id.id, 0);
	}
	if (g.will_allocate_ids) {
        for (int j = 0; j < am->n_genotypes; ++j) {
			++(d->current_id.id);
			am->ids[j] = d->current_id;
		}
	}
}
/** Make new genotypes (generic function)
 *
 * Applies all settings from gsc_GenOptions @a g. 
 *
 * Takes two parameter functions: @a parentChooser and @a offspringGenerator, described below.
 * 
 * @a parentChooser should function as a 'generator' of parents 
 * for the new genotypes. Every time the function is called, it 
 * should return a truthy value and save the gsc_GenoLocations of its 
 * chosen parents in its final parameter, or return a falsy value 
 * if there are no more offspring to be created. It has access to 
 * three informational parameters to help it choose parents: 
 * in order: @a parentIterator, passed by the calling function, 
 * @a datastore, passed by the calling function, and a pointer to 
 * a counter representing the number of times @a parentChooser 
 * has been called by this function prior to the current call. This 
 * generic function trusts that @a parentChooser will not return a 
 * truthy value if setting the parent choice gsc_GenoLocations to invalid
 * values (of a type that cannot be handled by the corresponding 
 * @a offspringGenerator).
 *
 * @a offspringGenerator should make the call that generates alleles
 * for a new genotype and saves them to a given gsc_GenoLocation. (It may 
 * also make any other necessary specific modifications to other data
 * about the new genotype). It will be called @a g.family_size times
 * for each truthy return value of @a parentChooser. As parameters,
 * it receives a pointer to the gsc_SimData, and @a datastore as passed 
 * by the calling function and potentially modified by @a parentChooser,
 * the gsc_GenoLocations of up to two parents, as set by @a parentChooser, 
 * and the location to which it is to save the new genotype.
 *
 * @return group number of new group created, or GSC_NO_GROUP if no group
 * was created or the new group was empty.
 */
gsc_GroupNum gsc_scaffold_make_new_genotypes(gsc_SimData* d, const gsc_GenOptions g,
		void* parentIterator, void* datastore, 
        int (*parentChooser)(void*, void*, unsigned int*, gsc_ParentChoice[static 2]),
        void (*offspringGenerator)(gsc_SimData*, void*, gsc_ParentChoice[static 2], gsc_GenoLocation) ) {
    if (g.family_size < 1 || d == NULL ||
		parentChooser == NULL || offspringGenerator == NULL) {
        return GSC_NO_GROUP;
	}
	
	// create the buffer we'll use to save the output crosses before they're printed.
    gsc_AlleleMatrix* offspring = gsc_create_empty_allelematrix(d->genome.n_markers, d->n_labels, d->label_defaults, CONTIG_WIDTH);
	unsigned int fullness = 0;
	unsigned int counter = 0;
    gsc_ParentChoice parents[2] = { { GSC_INVALID_GENO_LOCATION, GSC_UNINIT }, { GSC_INVALID_GENO_LOCATION, GSC_UNINIT } };

	gsc_AlleleMatrix* last = NULL;
    gsc_GroupNum output_group = GSC_NO_GROUP;
	if (g.will_save_to_simdata) {
		last = d->m; // for saving to simdata
		while (last->next != NULL) {
			last = last->next;
		}
		output_group = gsc_get_new_group_num( d );
	}

	// open the output files, if applicable
	FILE* fp = gsc_helper_genoptions_save_pedigrees_setup(g);
	int effIndex = GSC_UNINIT;
	FILE* fe = gsc_helper_genoptions_save_bvs_setup(d,g,&effIndex);
	FILE* fg = gsc_helper_genoptions_save_genotypes_setup(d,g);
	
	GetRNGstate();
	// loop through each combination
    while (parentChooser(parentIterator, datastore, &counter, parents)) {
		++counter;
		for (int f = 0; f < g.family_size; ++f, ++fullness) {
			R_CheckUserInterrupt();

			// when offspring buffer is full, save these outcomes to the file.
			if (fullness >= CONTIG_WIDTH) {
                offspring->n_genotypes = CONTIG_WIDTH;
				gsc_helper_genoptions_give_names_and_ids(offspring, d, g);
				gsc_helper_genoptions_save_pedigrees(fp, d, offspring);
                gsc_helper_genoptions_save_bvs(fe, d->e, effIndex, offspring);
				gsc_helper_genoptions_save_genotypes(fg, offspring);
				
				if (g.will_save_to_simdata) {
                    last->next = offspring;
					last = last->next;
                    offspring = gsc_create_empty_allelematrix(d->genome.n_markers, d->n_labels, d->label_defaults, CONTIG_WIDTH);
				}
				fullness = 0; //reset the count and start refilling the matrix
			}
			
			// do the cross.
            gsc_GenoLocation offspringPos = { .localAM=offspring, .localPos=fullness };
            offspringGenerator(d, datastore, parents, offspringPos);
			offspring->groups[fullness] = output_group;
			if (g.will_track_pedigree) {
                offspring->pedigrees[0][fullness] = gsc_get_id(parents[0].loc);
                offspring->pedigrees[1][fullness] = gsc_get_id(parents[1].loc);
			}
		}
	}
	PutRNGstate();

    offspring->n_genotypes = fullness;
	gsc_helper_genoptions_give_names_and_ids(offspring, d, g);
	gsc_helper_genoptions_save_pedigrees(fp, d, offspring);
    gsc_helper_genoptions_save_bvs(fe, d->e, effIndex, offspring);
	gsc_helper_genoptions_save_genotypes(fg, offspring);
	
	if (fp) fclose(fp);
	if (fe) fclose(fe);
	if (fg) fclose(fg);
	
	if (counter > 0 && g.will_save_to_simdata) {
        last->next = offspring;
        d->n_groups++;
		gsc_condense_allele_matrix( d );
		return output_group;
	} else {
        gsc_delete_allele_matrix( offspring );
        return GSC_NO_GROUP;
	}
}

/** parentChooser function parameter for gsc_make_random_crosses.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_random_crosses
 *
 * Chooses parents randomly without overruning parent usage caps 
 * and without permitting selfing. Guarantees @a parents contains 
 * two valid gsc_GenoLocations at the time it returns a truthy value.
 */
static int gsc_helper_parentchooser_cross_randomly(void* parentIterator, void* datastore, unsigned int* counter,
                                                   gsc_ParentChoice parents[static 2]) {
	gsc_RandomAccessIterator* it = (gsc_RandomAccessIterator*) parentIterator;
	unsigned int* datastoreint = (unsigned int*) datastore;
	unsigned int n_crosses = datastoreint[0];
	unsigned int cap = datastoreint[1];
	unsigned int nparents = datastoreint[2];
    //unsigned int mapindex = datastoreint[3];
    unsigned int* parent_uses = datastoreint + 4;
	unsigned int parentixs[2] = { 0 };

    if (*counter < n_crosses && (cap == 0 || (*counter) < cap * nparents)) {
		// get parents, randomly. Must not be identical or already been used too many times.
		parentixs[0] = gsc_randomdraw_replacementrules(it[0].d, nparents, cap, parent_uses, GSC_UNINIT);
		parentixs[1] = gsc_randomdraw_replacementrules(it[0].d, nparents, cap, parent_uses, parentixs[0]);
		
		if (cap > 0) { 
			parent_uses[parentixs[0]] += 1; 
			parent_uses[parentixs[1]] += 1;
		}
		
		// Neither of these should fail, if nparents is good. 
        parents[0].loc = gsc_next_get_nth(parentIterator, parentixs[0]);
        parents[1].loc = gsc_next_get_nth(parentIterator, parentixs[1]);
        // Reiterate map. Might save us a read to not bother checking their values first.
        parents[0].mapindex = datastoreint[3];
        parents[1].mapindex = datastoreint[3];
		// This will cut short gsc_scaffold_make_new_genotypes execution if either parent is invalid.
        return GSC_IS_VALID_LOCATION(parents[0].loc) && GSC_IS_VALID_LOCATION(parents[1].loc);
	} else {
		return GSC_FALSE;
	}
}

/** offspringGenerator function parameter for gsc_make_random_crosses.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_random_crosses
 *
 * Generates new alleles from two separate parents. Does not check 
 * they are valid parents.
 */
static void gsc_helper_make_offspring_cross(gsc_SimData* d, void* datastore, gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere) {
	// (silly name)
    gsc_generate_gamete(d, gsc_get_alleles(parents[0].loc), (gsc_get_alleles(putHere)  ), parents[0].mapindex);
    gsc_generate_gamete(d, gsc_get_alleles(parents[1].loc), (gsc_get_alleles(putHere)+1), parents[1].mapindex);
}

/** Check input parameters of random crossing functions. (helper function)
 *
 * @return 0 if the parameters fail these checks, number of members of 
 * @a from_group otherwise.
 */
static int gsc_helper_random_cross_checks(gsc_SimData* d, const gsc_GroupNum from_group, const int n_crosses, const int cap) {
	int g_size = gsc_get_group_size( d, from_group); // might be a better way to do this using the iterator.
    if (g_size < 1) {
        warning("Group %d does not exist.\n", from_group.num);
        return 0;
	}
	
	if (n_crosses < 1) {
        warning("Invalid n_crosses value provided: n_crosses must be greater than 0.\n");
        return 0;
    }

    if (cap < 0) {
        warning("Invalid cap value provided: cap can't be negative.\n");
        return 0;
    }
    if (cap > 0 && cap*g_size < n_crosses) {
        warning("Invalid cap value provided: cap of %d uses on %d parents too small to make %d crosses.\n", cap, g_size, n_crosses);
        return 0;
    }

	return g_size;
}	

/** Performs random crosses among members of a group. 
 * 
 * The group must have at least two members. Selfing is not considered a form of 
 * random crossing. The resulting genotypes are allocated to a new group.
 *
 * Preferences in gsc_GenOptions are applied to this cross. The family_size parameter
 * in gsc_GenOptions allows you to repeat each particular randomly-chosen cross a
 * certain number of times.
 *
 * Parents are drawn uniformly from the group when picking which crosses to make.
 *
 * @shortnamed{make_random_crosses}
 *
 * @param d pointer to the gsc_SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param from_group group number from which to draw the parents.
 * @param cap If set, the maximum number of times each member of from_group can be
 * used as the parent of a cross. Set to 0 for no restriction on the number of offspring
 * produced by a given member of from_group
 * @param n_crosses number of random pairs of parents to cross.
 * @param which_map recombination map to use to generate gametes from members of @a from_group. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
gsc_GroupNum gsc_make_random_crosses(gsc_SimData* d, const gsc_GroupNum from_group, const int n_crosses, const int cap,
                                     const gsc_MapID which_map, const gsc_GenOptions g) {
    int g_size = gsc_helper_random_cross_checks(d, from_group, n_crosses*2, cap);
	if (g_size == 0) {
        return GSC_NO_GROUP;
    } else if (g_size == 1) {
        warning("Group %d must contain multiple individuals to be able to perform random crossing\n", from_group.num);
        return GSC_NO_GROUP;
    }
    if (d->genome.n_maps < 1) {
        warning("Crossing requires at least one recombination map loaded\n");
        return GSC_NO_GROUP;
    }
    unsigned int map_index = 0;
    if (which_map.id != NO_MAP.id) { map_index = gsc_get_index_of_map(d, which_map); }
    if (map_index < 0) {
        warning("Could not find recombination map with identifier %i.\n", which_map.id);
        return GSC_NO_GROUP;
    }
	
    unsigned int* datastore = NULL; // cap = 0 means unlimited uses. Otherwise we need to track number of times each is used.
    if (cap > 0) {
        datastore = gsc_malloc_wrap(sizeof(unsigned int)*(g_size+4),GSC_TRUE);
        memset(datastore + 4,0,sizeof(unsigned int)*g_size);
    } else {
        datastore = gsc_malloc_wrap(sizeof(unsigned int)*4,GSC_TRUE);
	}
	datastore[0] = n_crosses;
	datastore[1] = cap;
	datastore[2] = g_size;
    datastore[3] = map_index;

	gsc_RandomAccessIterator parentit = gsc_create_randomaccess_iter( d, from_group);
	
    gsc_GroupNum offspring = gsc_scaffold_make_new_genotypes(d, g, (void*) &parentit, (void*) datastore,
                                                             gsc_helper_parentchooser_cross_randomly,
                                                             gsc_helper_make_offspring_cross );

	gsc_delete_randomaccess_iter(&parentit);
	GSC_FREE(datastore);
	return offspring;
}

/** Randomly pick a number in a range, optionally with a cap on how many times a 
 * number can be picked, and optionally required to be different to the last pick.
 *
 * Used in random crossing functions. 
 *
 * Draws a random integer from the range [0, max). All numbers have equal probability
 * of being drawn. The number will not be the  same as @a noCollision 
 *(set @a noCollision to max or greater to make all numbers in the range possible results),
 * and the number will fulfil @a member_uses[{number}] < @a cap, which can be used 
 * for selection without replacement or with only a certain number of replacements.
 *
 * @param d gsc_SimData, only used as the source of the random number generator (in genomicSimulationC version). 
 * @param max upper bound (non-inclusive) of the range to draw from.
 * @param cap maximum number of uses (as tracked by @a member_uses) of each number 
 * in the range. If, for a given number "num" in the range, @a member_uses["num"] is 
 * greater than or equal to @a cap, then the draw "num" will be discarded and the 
 * @param member_uses array of length @a max. See @a cap.
 * @param noCollision this integer cannot be the return value.
 * @return Random integer from the range 0 (inclusive) to @a max (exclusive) that fulfils 
 * the @a cap and @a noCollision conditions, or GSC_UNINIT if input parameters made it 
 * impossible to draw any number.
 */
unsigned int gsc_randomdraw_replacementrules(gsc_SimData* d, unsigned int max, unsigned int cap, unsigned int* member_uses, unsigned int noCollision) {
    if (max < 1 || (max == 1 && noCollision == 0)) {
		return GSC_UNINIT;
	}
	unsigned int parentix = 0;
	if (cap > 0) {  // n uses of each parent is capped at a number cap.
		do {
            parentix = round(unif_rand() * (max - 1));
        } while (parentix == noCollision || member_uses[parentix] >= cap);
	} else { // no cap on usage of each parent.
		do {
            parentix = round(unif_rand() * (max - 1));
        } while (parentix == noCollision);
	}
	return parentix;
}

/** parentChooser function parameter for gsc_make_random_crosses_between.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_random_crosses_between
 *
 * Chooses parents randomly without overruning parent usage caps 
 * and without permitting selfing. First and second parent come 
 * from different groups and may have different usage caps.
 * Guarantees @a parents contains 
 * two valid gsc_GenoLocations at the time it returns a truthy value.
 */
static int gsc_helper_parentchooser_cross_randomly_between(void* parentIterator, void* datastore, unsigned int* counter,
                                                           gsc_ParentChoice parents[static 2]) {
	// caller function should guarantee that nparents is not 1. How would you make a nonselfed cross then?
	gsc_RandomAccessIterator* it = (gsc_RandomAccessIterator*) parentIterator;
	unsigned int* datastoreint = (unsigned int*) datastore;
    unsigned int n_crosses = datastoreint[0];
    unsigned int* caps = datastoreint + 1; // caps[0], caps[1]
    unsigned int* nparents = datastoreint + 3; //nparents[0], nparents[1]
    unsigned int* map_index = datastoreint + 5; // map_index[0], map_index[1]
    unsigned int* parent_uses = datastoreint + 7; // parent_uses (nparents[0] + nparents[1] long)
	unsigned int parentixs[2] = { 0 };
	
    if (*counter < n_crosses && (caps[0] == 0 || (*counter) < caps[0] * nparents[0])
                             && (caps[1] == 0 || (*counter) < caps[1] * nparents[1])) {
		// get parents, randomly. Must not be identical or already been used too many times.
		parentixs[0] = gsc_randomdraw_replacementrules(it[0].d, nparents[0], caps[0], parent_uses, GSC_UNINIT);
		parentixs[1] = gsc_randomdraw_replacementrules(it[1].d, nparents[1], caps[1],  parent_uses + caps[0]*nparents[0], GSC_UNINIT);
			
		if (caps[0] > 0) { 
			parent_uses[parentixs[0]] += 1; 
		} 
		if (caps[1] > 0) {
			parent_uses[caps[0]*nparents[0] + parentixs[1]] += 1;
		}

        parents[0].loc = gsc_next_get_nth(it+0, parentixs[0]);
        parents[1].loc = gsc_next_get_nth(it+1, parentixs[1]);
        parents[0].mapindex = map_index[0];
        parents[1].mapindex = map_index[1];
        return GSC_IS_VALID_LOCATION(parents[0].loc) && GSC_IS_VALID_LOCATION(parents[1].loc);
	}
	return GSC_FALSE;
}

/** Performs random crosses where the first parent comes from one group and the second from
 *  another. 
 * 
 * The group must have at least two members.
 * The resulting genotypes are allocated to a new group. 
 *
 * Preferences in gsc_GenOptions are applied to this cross. The family_size parameter
 * in gsc_GenOptions allows you to repeat each particular randomly-chosen cross a
 * certain number of times.
 *
 * Parents are drawn uniformly from the group when picking which crosses to make.
 *
 * Parameters set_parent_gp1 and set_parent_gp2 are deprecated and removed!
 * Use gsc_make_group_from and gsc_combine_groups to temporarily move an individual to their own
 * group if you wish to cross randomly from a group to an individual. 
 *
 * @shortnamed{make_random_crosses_between}
 *
 * @param d pointer to the gsc_SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param group1 group number from which to draw the first parent.
 * @param group2 group number from which to draw the second parent.
 * @param n_crosses number of random pairs of parents to cross.
 * @param cap1 If set, the maximum number of times each member of group1 can be
 * used as the parent of a cross. Set to 0 for no restriction on the number of offspring
 * produced by a given member of group1
 * @param cap2 If set, the maximum number of times each member of group2 can be
 * used as the parent of a cross. Set to 0 for no restriction on the number of offspring
 * produced by a given member of group2
 * @param map1 recombination map to use to generate gametes from members of @a group1. If
 * NO_MAP, uses first/primary map by default
 * @param map2 recombination map to use to generate gametes from members of @a group2. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
gsc_GroupNum gsc_make_random_crosses_between(gsc_SimData*d, const gsc_GroupNum group1, const gsc_GroupNum group2, const int n_crosses,
                                             const int cap1, const int cap2, const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g) {
    int group1_size = gsc_helper_random_cross_checks(d, group1, n_crosses, cap1);
    int group2_size = gsc_helper_random_cross_checks(d, group2, n_crosses, cap2);
	if (group1_size == 0 || group2_size == 0) {
        return GSC_NO_GROUP;
	}
    if (d->genome.n_maps < 1) {
        warning("Crossing requires at least one recombination map loaded\n");
        return GSC_NO_GROUP;
    }
    unsigned int map1_index = 0;
    if (map1.id != NO_MAP.id) { map1_index = gsc_get_index_of_map(d, map1); }
    if (map1_index == GSC_UNINIT) {  //|| map_index > d->genome.n_maps) {
        warning("Could not find recombination map with identifier %i.\n", map1.id);
        return GSC_NO_GROUP;
    }
    unsigned int map2_index = 0;
    if (map2.id != NO_MAP.id) { map2_index = gsc_get_index_of_map(d, map2); }
    if (map2_index == GSC_UNINIT) {  //|| map_index > d->genome.n_maps) {
        warning("Could not find recombination map with identifier %i.\n", map2.id);
        return GSC_NO_GROUP;
    }

	unsigned int* datastore = NULL;
	int hascap1 = cap1 > 0 ? 1 : 0;
	int hascap2 = cap2 > 0 ? 1 : 0;
	if (hascap1 || hascap2) {
        datastore = gsc_malloc_wrap(sizeof(unsigned int)*(hascap1*group1_size + hascap2*group2_size + 7),GSC_TRUE);
        memset(datastore + 7,0,sizeof(unsigned int)*(hascap1*group1_size + hascap2*group2_size));
    } else {
        datastore = gsc_malloc_wrap(sizeof(unsigned int)*7,GSC_TRUE);
	}
    datastore[0] = n_crosses;
    datastore[1] = cap1;
    datastore[2] = cap2;
    datastore[3] = group1_size;
    datastore[4] = group2_size;
    datastore[5] = map1_index;
    datastore[6] = map2_index;
	
	gsc_RandomAccessIterator parentit[2] = { gsc_create_randomaccess_iter( d, group1 ),
										 gsc_create_randomaccess_iter( d, group2 ),
										};
	
    gsc_GroupNum offspring = gsc_scaffold_make_new_genotypes(d, g, (void*) parentit, (void*) datastore,
                                                             gsc_helper_parentchooser_cross_randomly_between,
                                                             gsc_helper_make_offspring_cross );

	gsc_delete_randomaccess_iter(&parentit[0]);
	gsc_delete_randomaccess_iter(&parentit[1]);
	GSC_FREE(datastore);
	return offspring;
	
}

/** parentChooser function parameter for gsc_make_targeted_crosses.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_targeted_crosses
 *
 * Locates the next pair of parents in the provided sequence of combinations. 
 * Guarantees @a parents contains 
 * two valid gsc_GenoLocations at the time it returns a truthy value.
 */
static int gsc_helper_parentchooser_cross_targeted(void* parentIterator, void* datastore, unsigned int* counter,
                                                   gsc_ParentChoice parents[static 2]) {
	gsc_RandomAccessIterator* it = (gsc_RandomAccessIterator*) parentIterator;
    int** datastoreintp = (int**) datastore;
    int ncrosses = **datastoreintp;
    if (*counter >= ncrosses) {
        return GSC_FALSE;
    }

    int** crosses = datastoreintp+1;
    int* maps = datastoreintp[3];
	
    parents[0].loc = gsc_next_get_nth(it, crosses[0][*counter]);
    parents[1].loc = gsc_next_get_nth(it, crosses[1][*counter]);
    parents[0].mapindex = maps[0];
    parents[1].mapindex = maps[1];

	// This will cut short gsc_scaffold_make_new_genotypes execution if either parent is invalid.
    return GSC_IS_VALID_LOCATION(parents[0].loc) && GSC_IS_VALID_LOCATION(parents[1].loc);
}

/** Performs the crosses of pairs of parents whose indexes are provided in an
 * array. 
 *
 * The resulting genotypes are allocated to a new group.
 *
 * Preferences in gsc_GenOptions are applied to this cross. The family_size parameter
 * in gsc_GenOptions allows you to repeat each particular cross a
 * certain number of times.
 *
 * Previously had a parameter combinations[2][n_combinations] instead of firstParents
 * and secondParents. This was changed to lower the boilerplate needs of calling this
 * function: now there is no need to create a 2-wide int* to hold the two separate parent
 * vectors if they already exist.
 *
 * @shortnamed{make_targeted_crosses}
 *
 * @param d pointer to the gsc_SimData object that includes genetic map data and
 * allele data needed to simulate crossing.
 * @param n_combinations the number of pairs to cross.
 * @param firstParents a vector of indexes of parents to be the first parent of each cross.
 * The vector must have at least [n_combinations] entries.
 * @param secondParents a vector of indexes of parents to be the second parent of each
 * cross. The vector must have at least [n_combinations] entries.
 * firstParents[0] is crossed to secondParents[0], firstParents[1] is crossed to 
 * secondParents[1], and so forth.
 * @param map1 recombination map to use to generate gametes from @a firstParents. If
 * NO_MAP, uses first/primary map by default
 * @param map2 recombination map to use to generate gametes from @a secondParents. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
gsc_GroupNum gsc_make_targeted_crosses(gsc_SimData* d, const int n_combinations, const int* firstParents, const int* secondParents,
                                       const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g) {
	if (n_combinations < 1) {
        warning("Invalid n_combinations value provided: n_combinations must be greater than 0.\n");
        return GSC_NO_GROUP;
	}
    if (d->genome.n_maps < 1) {
        warning("Crossing requires at least one recombination map loaded\n");
        return GSC_NO_GROUP;
    }
    int mapindexes[2] = { 0, 0 };
    if (map1.id != NO_MAP.id) { mapindexes[0] = gsc_get_index_of_map(d, map1); }
    if (mapindexes[0] == GSC_UNINIT) {  //|| map_index > d->genome.n_maps) {
        warning("Could not find recombination map with identifier %i.\n", map1.id);
        return GSC_NO_GROUP;
    }
    if (map2.id != NO_MAP.id) { mapindexes[1] = gsc_get_index_of_map(d, map2); }
    if (mapindexes[1] == GSC_UNINIT) {  //|| map_index > d->genome.n_maps) {
        warning("Could not find recombination map with identifier %i.\n", map2.id);
        return GSC_NO_GROUP;
    }
	
    const int* datastore[4] = { &n_combinations, firstParents, secondParents, mapindexes } ;
    gsc_RandomAccessIterator parentit = gsc_create_randomaccess_iter( d, GSC_NO_GROUP );
	
    gsc_GroupNum offspring = gsc_scaffold_make_new_genotypes(d, g, (void*) &parentit, (void*) datastore,
                                                             gsc_helper_parentchooser_cross_targeted,
                                                             gsc_helper_make_offspring_cross );

	gsc_delete_randomaccess_iter(&parentit);
	return offspring;
}

/** parentChooser function parameter for gsc_self_n_times.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_self_n_times
 * @see gsc_make_doubled_haploids
 *
 * Locates the next group member. Guarantees @a parents contains 
 * two valid gsc_GenoLocations (which are the same) at the time it returns a truthy value.
 */
static int gsc_helper_parentchooser_selfing(void* parentIterator, void* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]) {
	gsc_BidirectionalIterator* it = (gsc_BidirectionalIterator*) parentIterator;
	
    parents[0].loc = gsc_next_forwards(it);
    parents[0].mapindex = *((unsigned int*) datastore);
    parents[1] = parents[0];

    return GSC_IS_VALID_LOCATION(parents[0].loc);
}

/** offspringGenerator function parameter for gsc_self_n_times.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_self_n_times
 *
 * Datastore should be a pointer to an unsigned int containing 
 * the number of generations to self for.
 *
 * Generates new alleles by selfing from one parent. Only 
 * looks at the first parent and does not check it is valid.
 */
static void gsc_helper_make_offspring_self_n_times(gsc_SimData* d, void* datastore, gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere) {
    unsigned int n = ((unsigned int*) datastore)[1];
	
	// error checking parents are the same is not done.
	// error checking n >= 1 is not done.
	
    char* tmpparent = gsc_get_alleles(parents[0].loc);
    unsigned int map = parents[0].mapindex;
    GSC_CREATE_BUFFER(tmpchild,char,d->genome.n_markers<<1);
	char* output = gsc_get_alleles(putHere);
	int n_oddness = n % 2; 
	for (int i = 0; i < n; ++i) {
		if (i % 2 == n_oddness) {
            gsc_generate_gamete(d, tmpparent, tmpchild, map);
            gsc_generate_gamete(d, tmpparent, tmpchild+1, map);
			tmpparent = tmpchild;
		} else {
            gsc_generate_gamete(d, tmpparent, output, map);
            gsc_generate_gamete(d, tmpparent, output+1, map);
			tmpparent = output;
		}
	}
    GSC_DELETE_BUFFER(tmpchild);
}

/** Selfs each member of a group for a certain number of generations.
 *
 * The resulting genotypes are allocated to a new group.
 *
 * Only the genotype after all n generations is saved. Intermediate steps will be lost.
 *
 * Preferences in gsc_GenOptions are applied to this operation. The family_size parameter
 * in gsc_GenOptions allows you to generate multiple selfed offspring from each member
 * of the group. These multiple selfed offspring all originate from independent
 * processes of selfing the parent with itself then selfing that intermediate offspring
 * with itself for the remaining n-1 generations.
 *
 * @shortnamed{self_n_times}
 *
 * @param d pointer to the gsc_SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param n number of generations of selfing simulation to carry out.
 * @param group the genotypes on which to perform these n generations of selfing.
 * @param which_map recombination map to use to generate gametes from members of @a group. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
gsc_GroupNum gsc_self_n_times(gsc_SimData* d, const unsigned int n, const gsc_GroupNum group, const gsc_MapID which_map, const gsc_GenOptions g) {
	/*int group_size = gsc_get_group_size( d, group);
	if (group_size < 1) {
        warning("Group %d does not exist.\n", group.num);
        return GSC_NO_GROUP;
	}*/
	if (n < 1) {
        warning("Invalid n value provided: Number of generations must be greater than 0.\n");
        return GSC_NO_GROUP;
	}
    if (d->genome.n_maps < 1) {
        warning("Selfing requires at least one recombination map loaded\n");
        return GSC_NO_GROUP;
    }
    unsigned int map_index = 0;
    if (which_map.id != NO_MAP.id) { map_index = gsc_get_index_of_map(d, which_map); }
    if (map_index == GSC_UNINIT) {  //|| map_index > d->genome.n_maps) {
        warning("Could not find recombination map with identifier %i.\n", which_map.id);
        return GSC_NO_GROUP;
    }
    unsigned int datastore[2] = { map_index, n };
	
	gsc_BidirectionalIterator parentit = gsc_create_bidirectional_iter( d, group );
	
    return gsc_scaffold_make_new_genotypes(d, g, (void*) &parentit, (void*) datastore,
                                           gsc_helper_parentchooser_selfing,
                                           gsc_helper_make_offspring_self_n_times );
}

/** offspringGenerator function parameter for gsc_make_doubled_haploids.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_doubled_haploids
 *
 * Generates new alleles by making a doubled haploid. Only 
 * looks at the first parent and does not check it is valid.
 */
static void gsc_helper_make_offspring_doubled_haploids(gsc_SimData* d, void* datastore, gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere) {
    gsc_generate_doubled_haploid(d, gsc_get_alleles(parents[0].loc), gsc_get_alleles(putHere), parents[0].mapindex);
}

/** Creates a doubled haploid from each member of a group.
 *
 * The resulting genotypes are allocated to a new group.
 *
 * Preferences in gsc_GenOptions are applied to this operation. The family_size parameter
 * in gsc_GenOptions allows you to generate multiple doubled haploid offspring from each member
 * of the group. These multiple offspring all originate from independent
 * processes of generating a gamete then doubling its alleles.
 *
 * @shortnamed{make_doubled_haploids}
 *
 * @param d pointer to the gsc_SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param group the genotypes on which to perform the operation.
 * @param which_map recombination map to use to generate gametes from members of @a group. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
gsc_GroupNum gsc_make_doubled_haploids(gsc_SimData* d, const gsc_GroupNum group, const gsc_MapID which_map, const gsc_GenOptions g) {
	/*int group_size = gsc_get_group_size( d, group);
	if (group_size < 1) {
        warning("Group %d does not exist.\n", group.num);
        return GSC_NO_GROUP;
	}*/
    if (d->genome.n_maps < 1) {
        warning("Crossing requires at least one recombination map loaded\n");
        return GSC_NO_GROUP;
    }
    unsigned int map_index = 0;
    if (which_map.id != NO_MAP.id) { map_index = gsc_get_index_of_map(d, which_map); }
    if (map_index == GSC_UNINIT) {
        warning("Could not find recombination map with identifier %i.\n", which_map.id);
        return GSC_NO_GROUP;
    }
	
	gsc_BidirectionalIterator parentit = gsc_create_bidirectional_iter( d, group );
	
    return gsc_scaffold_make_new_genotypes(d, g, (void*) &parentit, (void*) &map_index,
                                           gsc_helper_parentchooser_selfing,
                                           gsc_helper_make_offspring_doubled_haploids );
}

/** parentChooser function parameter for gsc_make_clones.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_clones
 *
 * Very similar to @ref gsc_helper_parentchooser_selfing, but also accesses the datastore
 * to see if it needs to inherit the name from its parent.
 *
 * Locates the next group member. Guarantees @a parents contains 
 * two valid gsc_GenoLocations (which are the same) at the time it returns a truthy value.
 */
static int gsc_helper_parentchooser_cloning(void* parentIterator, void* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]) {
	gsc_BidirectionalIterator* it = (gsc_BidirectionalIterator*) parentIterator;
	
    parents[0].loc = gsc_next_forwards(it);
	parents[1] = parents[0];

    if (GSC_IS_VALID_LOCATION(parents[0].loc)) {
		char inherit_names = *( ((char**) datastore)[0] );
		if (inherit_names) {
            ((char**) datastore)[1] = gsc_get_name(parents[0].loc);
		}
		return GSC_TRUE;
	} else {
		return GSC_FALSE;
	}
}

/** offspringGenerator function parameter for gsc_make_clones.
 *
 * @see gsc_scaffold_make_new_genotypes
 * @see gsc_make_clones
 *
 * Generates new alleles by making a clone of its parent. Only 
 * looks at the first parent and does not check it is valid.
 */
static void gsc_helper_make_offspring_clones(gsc_SimData* d, void* datastore, gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere) {
	char inherit_names = *( ((char**) datastore)[0] );
	if (inherit_names) {
        char* tmpname = gsc_malloc_wrap(sizeof(char)*(strlen(((char**) datastore)[1]) + 1),GSC_TRUE);
        strcpy(tmpname, ((char**) datastore)[1]);
        gsc_set_name(putHere,tmpname);
	}
	
    gsc_generate_clone(d, gsc_get_alleles(parents[0].loc), gsc_get_alleles(putHere));
}

/** Creates an identical copy of each member of a group.
 *
 * The resulting genotypes are allocated to a new group.
 *
 * Preferences in gsc_GenOptions are applied to this operation. The family_size parameter
 * in gsc_GenOptions allows you to generate multiple cloned offspring from each member
 * of the group.
 *
 * If pedigree tracking and ID allocation are active in gsc_GenOptions, clones are given
 * individual IDs and are children of their single progenitor parent. If the inherit_names
 * parameter is 1/truthy, it overrides whatever naming settings are present in gsc_GenOptions
 * in favour of giving each clone the exact name of the individual it was cloned from.
 *
 * Clones currently keep the default value for every custom label.
 *
 * @shortnamed{make_clones}
 *
 * @param d pointer to the gsc_SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param group the genotypes on which to perform the operation.
 * @param g options for the genotypes created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
gsc_GroupNum gsc_make_clones(gsc_SimData* d, const gsc_GroupNum group, const int inherit_names, gsc_GenOptions g) {
    /*int group_size = gsc_get_group_size( d, group);
    if (group_size < 1) {
        warning("Group %d does not exist.\n", group.num);
        return GSC_NO_GROUP;
    }*/
	
	char do_inherit_names = inherit_names;
	char* nameinherit[2] = {&do_inherit_names, NULL };

    gsc_BidirectionalIterator parentit = gsc_create_bidirectional_iter( d, group );
	
    return gsc_scaffold_make_new_genotypes(d, g, (void*) &parentit, (void*) nameinherit,
                                           gsc_helper_parentchooser_cloning,
                                           gsc_helper_make_offspring_clones );
}


/** Perform crosses between all pairs of parents
 * in the group `from_group` and allocates the resulting offspring to a new group.
 *
 * If the group has n members, there will be $n * (n-1) / 2$ offspring produced.
 *
 * Preferences in gsc_GenOptions are applied to this cross.
 *
 * @shortnamed{make_all_unidirectional_crosses}
 *
 * @param d pointer to the gsc_SimData object containing or markers and parent alleles
 * @param from_group group number from which to do all these crosses.
 * @param which_map recombination map to use to generate gametes from members of @a from_group. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the gsc_AlleleMatrix created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
gsc_GroupNum gsc_make_all_unidirectional_crosses(gsc_SimData* d, const gsc_GroupNum from_group, const gsc_MapID mapID, const gsc_GenOptions g) {
    int group_size = gsc_get_group_size( d, from_group );
	if (group_size < 2) {
		if (group_size == 1) {
            warning("Group %d does not have enough members to perform crosses\n", from_group.num);
		} else {
            warning("Group %d does not exist.\n", from_group.num);
		}
        return GSC_NO_GROUP;
	}
    GSC_CREATE_BUFFER(group_indexes,unsigned int,group_size);
    gsc_get_group_indexes( d, from_group, group_size, group_indexes );

	// number of crosses = number of entries in upper triangle of matrix
	//    = half of (n entries in matrix - length of diagonal)
	// 	  = half of (lmatrix * lmatrix - lmatrix);
	int n_crosses = group_size * (group_size - 1) / 2; //* g.family_size;

    GSC_CREATE_BUFFER(combos0,int,n_crosses);
    GSC_CREATE_BUFFER(combos1,int,n_crosses);
    int* combinations[2] = {combos0, combos1};
	int cross_index = 0;
	for (int i = 0; i < group_size; ++i) {
		for (int j = i + 1; j < group_size; ++j) {
			combinations[0][cross_index] = group_indexes[i];
			combinations[1][cross_index] = group_indexes[j];

			++cross_index;
		}
	}

    GSC_DELETE_BUFFER(group_indexes);
    gsc_GroupNum out = gsc_make_targeted_crosses(d, n_crosses, combinations[0], combinations[1], mapID, mapID, g);
    GSC_DELETE_BUFFER(combos0);
    GSC_DELETE_BUFFER(combos1);
    return out;
}

/** Find the top m percent of a group and perform random crosses between those
 * top individuals. The resulting offspring are allocated to a new group.
 *
 * Preferences in gsc_GenOptions are applied to this cross.
 *
 * @shortnamed{make_n_crosses_from_top_m_percent}
 *
 * @param d pointer to the gsc_SimData object containing or markers and parent alleles
 * @param n number of random crosses to make.
 * @param m percentage (as a decimal, i.e. 0.2 for 20percent). Take the best [this portion]
 * of the group for performing the random crosses.
 * @param group group number from which to identify top parents and do crosses
 * @param effID Identifier of the set of marker effects to be used to calculate the breeding values
 * with which to rank top parents
 * @param which_map recombination map to use to generate gametes from members of @a group. If
 * NO_MAP, uses first/primary map by default
 * @param g options for the gsc_AlleleMatrix created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
gsc_GroupNum gsc_make_n_crosses_from_top_m_percent(gsc_SimData* d, const int n, const int m, const gsc_GroupNum group,
                                                   const gsc_MapID mapID, const gsc_EffectID effID, const gsc_GenOptions g) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (n < 1) {
        warning("Invalid n value provided: Number of crosses must be greater than 0.\n");
        return GSC_NO_GROUP;
    }
    if (effIndex == GSC_UNINIT) {
        warning("Cannot find marker effect set with id %d\n", effID.id);
        return GSC_NO_GROUP;
    }
    if (m < 1 || m > 100) {
        warning("Invalid m value provided: Percent to select must be between 1 and 100.\n");
        return GSC_NO_GROUP;
    }

	// move the top m% to a new group
	int group_size = gsc_get_group_size(d, group);
    if (group_size < 1) {
        warning("Group %d does not exist.\n", group.num);
        return GSC_NO_GROUP;
    }

    int n_top_group = group_size * m / 100;
	Rprintf("There are %d lines in the top %d%%\n", n_top_group, m);

    gsc_GroupNum topgroup = gsc_split_by_bv(d, group, effID, n_top_group, GSC_FALSE);

	// do the random crosses
    gsc_GroupNum gp = gsc_make_random_crosses(d, topgroup, 0, n, mapID, g);

	// unconvert from a group
    gsc_GroupNum to_combine[] = {group, topgroup};
	gsc_combine_groups(d, 2, to_combine);

	return gp;
}

/** Perform crosses between pairs of parents identified by name in a file
 * and allocate the resulting offspring to a new group.
 *
 * The input file should have format:
 *
 * [parent1 name] [parent2 name]
 *
 * [parent1 name] [parent2 name]
 *
 * ...
 *
 * where each row represents a separate cross to carry out.
 *
 * Preferences in gsc_GenOptions are applied to this cross.
 *
 * @shortnamed{make_crosses_from_file}
 *
 * @param d pointer to the gsc_SimData object containing or markers and parent alleles
 * @param input_file file instructing which crosses to perform
 * @param map1 recombination map to use to generate gametes from the first parent in the
 * crosses in the file. If NO_MAP, uses first/primary map by default
 * @param map2 recombination map to use to generate gametes from the second parent in the
 * crosses in the file. If NO_MAP, uses first/primary map by default
 * @param g options for the gsc_AlleleMatrix created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
gsc_GroupNum gsc_make_crosses_from_file(gsc_SimData* d, const char* input_file, const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g) {
	struct gsc_TableSize t = gsc_get_file_dimensions(input_file, '\t');
	if (t.num_rows < 1) {
		warning( "No crosses exist in that file\n");
        return GSC_NO_GROUP;
	}

	//open file
	FILE* fp;
	if ((fp = fopen(input_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", input_file);
	}

    GSC_CREATE_BUFFER(combos0,int,t.num_rows);
    GSC_CREATE_BUFFER(combos1,int,t.num_rows);
    int* combinations[2] = {combos0,combos1};
    char buffer[2][NAME_LENGTH];
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s \n", buffer[0], buffer[1]);
		combinations[0][i] = gsc_get_index_of_name(d->m, buffer[0]);
		combinations[1][i] = gsc_get_index_of_name(d->m, buffer[1]);
	}

	fclose(fp);
    gsc_GroupNum out = gsc_make_targeted_crosses(d, t.num_rows, combinations[0], combinations[1], map1, map2, g);
    GSC_DELETE_BUFFER(combos0);
    GSC_DELETE_BUFFER(combos1);
    return out;
}

/** Perform crosses between previously-generated offspring of pairs of parents
 * identified by name in a file. The resulting offspring are allocated to a new
 * group.
 *
 * The input file should have format:
 *
 * [gparent1 name] [gparent2 name] [gparent3 name] [gparent4 name]
 *
 * [gparent1 name] [gparent2 name] [gparent3 name] [gparent4 name]
 *
 * ...
 *
 * where each row represents a separate cross to carry out.
 *
 * Results of the cross [gparent1 name] x [gparent2 name] and
 * [gparent3 name] x [gparent4 name] must have already been generated
 * with pedigree tracking options turned on. Messages will be printed
 * to stderr if such offspring cannot be found.
 *
 * Preferences in gsc_GenOptions are applied to this cross.
 *
 * @shortnamed{make_double_crosses_from_file}
 *
 * @param d pointer to the gsc_SimData object containing or markers and parent alleles
 * @param input_file file instructing which crosses to perform
 * @param map1 recombination map to use to generate gametes from the first parent in each of the
 * crosses in the file. If NO_MAP, uses first/primary map by default
 * @param map2 recombination map to use to generate gametes from the second parent in each of the
 * crosses in the file. If NO_MAP, uses first/primary map by default
 * @param g options for the gsc_AlleleMatrix created. @see gsc_GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
gsc_GroupNum gsc_make_double_crosses_from_file(gsc_SimData* d, const char* input_file, const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g) {
	struct gsc_TableSize t = gsc_get_file_dimensions(input_file, '\t');
	if (t.num_rows < 1) {
		warning( "No crosses exist in that file\n");
        return GSC_NO_GROUP;
	}

	//open file
	FILE* fp;
	if ((fp = fopen(input_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", input_file);
	}

    GSC_CREATE_BUFFER(combos0,int,t.num_rows);
    GSC_CREATE_BUFFER(combos1,int,t.num_rows);
    int* combinations[2] = {combos0,combos1};
    char buffer[4][NAME_LENGTH];
    const char* to_buffer[] = {buffer[0], buffer[1], buffer[2], buffer[3]};
    gsc_PedigreeID g0_id[4];
	int f1_i[2];
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s %s %s \n", buffer[0], buffer[1], buffer[2], buffer[3]);
		gsc_get_ids_of_names(d->m, 4, to_buffer, g0_id);
        if (g0_id[0].id == GSC_NO_PEDIGREE.id || g0_id[1].id == GSC_NO_PEDIGREE.id || g0_id[2].id == GSC_NO_PEDIGREE.id || g0_id[3].id == GSC_NO_PEDIGREE.id) {
			warning( "Could not go ahead with the line %d cross - g0 names not in records\n", i);
            combinations[0][i] = GSC_UNINIT;
            combinations[1][i] = GSC_UNINIT;
			continue;
		}

		// identify two parents
		f1_i[0] = gsc_get_index_of_child(d->m, g0_id[0], g0_id[1]);
		f1_i[1] = gsc_get_index_of_child(d->m, g0_id[2], g0_id[3]);
		if (f1_i[0] < 0 || f1_i[1] < 0) {
			// try different permutations of the four grandparents.
            f1_i[0] = gsc_get_index_of_child(d->m, g0_id[0], g0_id[2]);
            f1_i[1] = gsc_get_index_of_child(d->m, g0_id[1], g0_id[3]);
			if (f1_i[0] < 0 || f1_i[1] < 0) {
				f1_i[0] = gsc_get_index_of_child(d->m, g0_id[0], g0_id[3]);
				f1_i[1] = gsc_get_index_of_child(d->m, g0_id[1], g0_id[2]);
				if (f1_i[0] < 0 || f1_i[1] < 0) {
					warning( "Could not go ahead with the line %d cross - f1 children do not exist for this quartet\n", i);
                    combinations[0][i] = GSC_UNINIT;
                    combinations[1][i] = GSC_UNINIT;
					continue;
				}
			}
		}

		//add them to a combinations list
		combinations[0][i] = f1_i[0];
		combinations[1][i] = f1_i[1];

	}

	fclose(fp);
    gsc_GroupNum out = gsc_make_targeted_crosses(d, t.num_rows, combinations[0], combinations[1], map1, map2, g);
    GSC_DELETE_BUFFER(combos0);
    GSC_DELETE_BUFFER(combos1);
    return out;
}


/*--------------------------------Fitness------------------------------------*/

/** Takes the `top_n` individuals in the group with the best breeding values/fitnesses
 * and puts them in a new group. The new group number is returned.
 *
 * @shortnamed{split_by_bv}
 *
 * @param d pointer to the gsc_SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param group group number from which to split the top individuals.
 * @param effID Identifier of the set of marker effects to use to calculate BV.
 * @param top_n The number of individuals to put in the new group.
 * @param lowIsBest boolean, if TRUE the `top_n` with the lowest breeding value
 * will be selected, if false the `top_n` with the highest breeding value are.
 * @returns the group number of the newly-created split-off group
 */
gsc_GroupNum gsc_split_by_bv(gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID, const int top_n, const int lowIsBest) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT || d->e[effIndex].effects.rows < 1 || d->m == NULL) {
        error( "Either effect matrix or allele matrix does not exist\n");
    }

    unsigned int group_size = gsc_get_group_size( d, group );
    GSC_CREATE_BUFFER(group_indexes,unsigned int,group_size);
    gsc_get_group_indexes( d, group, group_size, group_indexes );
	
	if (group_size <= top_n) {
		// well we'll just have to move em all
        gsc_GroupNum migration = gsc_make_group_from(d, group_size, group_indexes);
		return migration;
	}
	
	// This should be ordered the same as the indexes
    gsc_DecimalMatrix fits = gsc_calculate_group_bvs( d, group, effID ); // 1 by group_size matrix
	
	// get an array of pointers to those fitnesses
    GSC_CREATE_BUFFER(p_fits,double*,fits.cols);
	for (int i = 0; i < fits.cols; i++) {
		p_fits[i] = &(fits.matrix[0][i]);
	}

	// sort descending
	if (lowIsBest) {
        qsort(p_fits, fits.cols, sizeof(double*), gsc_helper_ascending_pdouble_comparer);
	} else {
        qsort(p_fits, fits.cols, sizeof(double*), gsc_helper_descending_pdouble_comparer);
	}

	// save the indexes of the best n
    GSC_CREATE_BUFFER(top_individuals,unsigned int,top_n);
	for (int i = 0; i < top_n; i++) {
        top_individuals[i] = group_indexes[p_fits[i] - fits.matrix[0]];
	}
	gsc_delete_dmatrix(&fits);
    GSC_DELETE_BUFFER(p_fits);
    GSC_DELETE_BUFFER(group_indexes);

	// send those n to a new group
    gsc_GroupNum out = gsc_make_group_from(d, top_n, top_individuals);
    GSC_DELETE_BUFFER(top_individuals);
    return out;
}

/** Calculates the fitness metric/breeding value for each genotype in the gsc_AlleleMatrix
* in a certain group, and returns the results in a gsc_DecimalMatrix.
*
* The breeding value is calculated for each genotype by taking the sum of (number of
* copies of this allele at this marker times this allele's effect at this marker)
* for each marker for each different allele. To do
* this, the vector containing each allele's effect values is multiplied by a matrix
* containing the counts of that allele at each marker.
*
* The function exits with error code 1 if no marker effect file is loaded.
*
 * @shortnamed{calculate_group_bvs}
*
* @param d pointer to the gsc_SimData object to which the groups and individuals belong.
* It must have a marker effect file loaded to successfully run this function.
* @param group calculate breeding values for each genotype in the group with this group number.
* @param effID Identifier of the marker effect set to be used to calculate these breeding values
* @returns A gsc_DecimalMatrix containing the score for each individual in the group.
*/
gsc_DecimalMatrix gsc_calculate_group_bvs(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID) {
	// check that both of the items to be multiplied exist.
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT || d->e[effIndex].effects.rows < 1 || d->m == NULL) {
		error( "Either effect matrix or allele matrix does not exist\n");
	}
    gsc_EffectMatrix e = d->e[effIndex];

    int group_size = gsc_get_group_size( d, group );
	gsc_DecimalMatrix sum = gsc_generate_zero_dmatrix(1, group_size);
    if (group_size < 1) {
        warning("Group %d does not exist.\n", group.num);
        return sum;
    }
    gsc_DecimalMatrix counts = gsc_generate_zero_dmatrix(group_size, d->genome.n_markers);
    gsc_DecimalMatrix counts2 = gsc_generate_zero_dmatrix(group_size, d->genome.n_markers);

    int i = 0; // highest allele index

    for (; i < e.effects.rows - 1; i += 2) {
		// get the allele counts in counts
        gsc_calculate_group_count_matrix_pair(d, group, e.effect_names[i], &counts, e.effect_names[i+1], &counts2);

		// multiply counts with effects and add to bv sum
        gsc_add_doublematrixvector_product_to_dmatrix(&sum, &counts, e.effects.matrix[i], &counts2, e.effects.matrix[i+1]);
	}
    if (i < e.effects.rows) { // deal with the last odd-numbered allele
        gsc_calculate_group_count_matrix(d, group, e.effect_names[i], &counts);
        gsc_add_matrixvector_product_to_dmatrix(&sum, &counts, e.effects.matrix[i]);
	}

	gsc_delete_dmatrix(&counts);
	gsc_delete_dmatrix(&counts2);

	return sum;
}

/** Calculates the fitness metric/breeding value for each genotype in the gsc_AlleleMatrix,
 * and returns the results in a gsc_DecimalMatrix struct.
 *
 * The breeding value is calculated for each genotype by taking the sum of (number of
 * copies of this allele at this marker times this allele's effect at this marker)
 * for each marker for each different allele. To do
 * this, the vector containing each allele's effect values is multiplied by a matrix
 * containing the counts of that allele at each marker.
 *
 * The function exits with error code 1 if no marker effect file is loaded.
 *
 * @shortnamed{calculate_bvs}
 *
 * @param m pointer to the gsc_AlleleMatrix object to which the genotypes belong.
 * @param e pointer to the gsc_EffectMatrix that effect values have been loaded into.
 * @returns A gsc_DecimalMatrix containing the score for each individual in the group.
 */
gsc_DecimalMatrix gsc_calculate_bvs( const gsc_AlleleMatrix* m, const gsc_EffectMatrix* e) {
	// check that both of the items to be multiplied exist.
    if (e->effects.rows < 1 || m == NULL) {
		error( "Either effect matrix or allele matrix does not exist\n");
	}

	gsc_DecimalMatrix sum = gsc_generate_zero_dmatrix(1, m->n_genotypes);
	gsc_DecimalMatrix counts = gsc_generate_zero_dmatrix(m->n_genotypes, m->n_markers);
	gsc_DecimalMatrix counts2 = gsc_generate_zero_dmatrix(m->n_genotypes, m->n_markers);

    int i = 0; // highest allele index

	for (; i < e->effects.rows - 1; i += 2) {
		// get the allele counts in counts
		gsc_calculate_count_matrix_pair(m, e->effect_names[i], &counts, e->effect_names[i+1], &counts2);

		// multiply counts with effects and add to bv sum
		gsc_add_doublematrixvector_product_to_dmatrix(&sum, &counts, e->effects.matrix[i], &counts2, e->effects.matrix[i+1]);
	}
	if (i < e->effects.rows) { // deal with the last odd-numbered allele
		gsc_calculate_count_matrix(m, e->effect_names[i], &counts);
		gsc_add_matrixvector_product_to_dmatrix(&sum, &counts, e->effects.matrix[i]);
	}

	gsc_delete_dmatrix(&counts);
	gsc_delete_dmatrix(&counts2);

	return sum;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in a group.
 * Returns the result in a pre-created gsc_DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @shortnamed{calculate_group_count_matrix}
 *
 * @param d pointer to the gsc_SimData.
 * @param group group number for which to calculate counts of the allele, or 0
 * to count alleles for all genotypes in the simulation
 * @param allele the single-character allele to be counting.
 * @param counts pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int gsc_calculate_group_count_matrix( const gsc_SimData* d, const gsc_GroupNum group, const char allele, gsc_DecimalMatrix* counts) {
    if (group.num == GSC_NO_GROUP.num) {
        return GSC_UNINIT; // @@
    } else {
        int groupSize = gsc_get_group_size(d, group);
        if (counts->rows < groupSize || counts->cols < d->genome.n_markers) {
            fprintf(stderr, "`counts` is the wrong size to be filled: needs %u by %d but is %d by %d\n",
                    (unsigned int) d->genome.n_markers, groupSize, counts->rows, counts->cols);
            return 1;
        }

        GSC_CREATE_BUFFER(genes,char*,groupSize);
        gsc_get_group_genes(d, group, groupSize, genes);

        for (int i = 0; i < groupSize; ++i) {
            R_CheckUserInterrupt();
            for (int j = 0; j < d->genome.n_markers; ++j) {
                int cell_sum = 0;
                if (genes[i] != NULL) {
                    if (genes[i][2*j] == allele)     cell_sum += 1;
                    if (genes[i][2*j + 1] == allele) cell_sum += 1;
                }
                counts->matrix[i][j] = cell_sum;
            }
        }
        GSC_DELETE_BUFFER(genes);
        return 0;
    }
}

/** Calculates the number of times at each marker that two particular alleles appear
 * for each genotype in a group.
 * Returns the result in two pre-created gsc_DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * This function is the same as gsc_calculate_group_count_matrix(), but for
 * two alleles at a time, for the purpose of loop unrolling in breeding value
 * calculations.
 *
 * @shortnamed{calculate_group_count_matrix_pair}
 *
 * @param d pointer to the gsc_SimData.
 * @param group group number for which to calculate counts of the allele, or 0 to
 * count alleles for all genotypes in the simulation
 * @param allele the first single-character allele to be counting.
 * @param counts pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @param allele2 the second single-character allele to be counting.
 * @param counts2 pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele2` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int gsc_calculate_group_count_matrix_pair( const gsc_SimData* d, const gsc_GroupNum group, const char allele, gsc_DecimalMatrix* counts, const char allele2, gsc_DecimalMatrix* counts2) {
    if (group.num == GSC_NO_GROUP.num) {
        return GSC_UNINIT; //@@
    } else {
        int groupSize = gsc_get_group_size(d, group);
        if (counts->rows < groupSize || counts->cols < d->genome.n_markers) {
            fprintf(stderr, "`counts` is the wrong size to be filled: needs %u by %d but is %d by %d\n",
                    (unsigned int) d->genome.n_markers, groupSize, counts->rows, counts->cols);
            return 1;
        }
        if (counts2->rows < groupSize || counts2->cols < d->genome.n_markers) {
            fprintf(stderr, "`counts2` is the wrong size to be filled: needs %u by %d but is %d by %d\n",
                    (unsigned int) d->genome.n_markers, groupSize, counts2->rows, counts2->cols);
            return 1;
        }

        GSC_CREATE_BUFFER(genes,char*,groupSize);
        gsc_get_group_genes(d, group, groupSize, genes);

        for (int i = 0; i < groupSize; ++i) {
            R_CheckUserInterrupt();
            if (genes[i] == NULL) {
                continue;
            }

            for (int j = 0; j < d->genome.n_markers; ++j) {
                int cell_sum = 0;
                int cell_sum2 = 0;
                if      (genes[i][2*j] == allele)      { ++cell_sum; }
                else if (genes[i][2*j] == allele2)     { ++cell_sum2;}
                if      (genes[i][2*j + 1] == allele)  { ++cell_sum; }
                else if (genes[i][2*j + 1] == allele2) { ++cell_sum2;}
                counts->matrix[i][j] = cell_sum;
                counts2->matrix[i][j] = cell_sum2;
            }
        }
        GSC_DELETE_BUFFER(genes);
        return 0;
    }
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in a particular gsc_AlleleMatrix (does not go through the linked list).
 * Returns the result in a pre-created gsc_DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @shortnamed{calculate_count_matrix}
 *
 * @param m pointer to the gsc_AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the single-character allele to be counting.
 * @param counts pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the gsc_AlleleMatrix.
 * @returns 0 on success, nonzero on failure.
 * */
int gsc_calculate_count_matrix( const gsc_AlleleMatrix* m, const char allele, gsc_DecimalMatrix* counts) {
	//gsc_DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, m->n_genotypes);
    if (counts->rows < m->n_genotypes || counts->cols < m->n_markers) {
        fprintf(stderr, "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                m->n_markers, m->n_genotypes, counts->rows, counts->cols);
        return 1;
    }

    for (int i = 0; i < m->n_genotypes; ++i) {
        R_CheckUserInterrupt();
        if (m->alleles[i] == NULL) {
            continue;
        }

        for (int j = 0; j < m->n_markers; ++j) {
            int cell_sum = 0;
            if (m->alleles[i][2*j] == allele)     cell_sum += 1;
            if (m->alleles[i][2*j + 1] == allele) cell_sum += 1;
            counts->matrix[i][j] = cell_sum;
        }
    }
    return 0;
}

/** Calculates the number of times at each marker that two particular alleles appear
 * for each genotype in an gsc_AlleleMatrix (does not go through the linked list).
 * Returns the result in two pre-created gsc_DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * This function is the same as gsc_calculate_count_matrix(), but for
 * two alleles at a time, for the purpose of loop unrolling in breeding value
 * calculations.
 *
 * @shortnamed{calculate_count_matrix_pair}
 *
 * @param m pointer to the gsc_AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the first single-character allele to be counting.
 * @param counts pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @param allele2 the second single-character allele to be counting.
 * @param counts2 pointer to the gsc_DecimalMatrix into which to put the number
 * of `allele2` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int gsc_calculate_count_matrix_pair( const gsc_AlleleMatrix* m, const char allele, gsc_DecimalMatrix* counts, const char allele2, gsc_DecimalMatrix* counts2) {
	if (counts->rows < m->n_genotypes || counts->cols < m->n_markers) {
        fprintf(stderr, "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                m->n_markers, m->n_genotypes, counts->rows, counts->cols);
        return 1;
	}
    if (counts2->rows < m->n_genotypes || counts2->cols < m->n_markers) {
        fprintf(stderr, "`counts2` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                m->n_markers, m->n_genotypes, counts2->rows, counts2->cols);
        return 1;
	}

	for (int i = 0; i < m->n_genotypes; ++i) {
		R_CheckUserInterrupt();
		if (m->alleles[i] == NULL) {
			continue;
		}

		for (int j = 0; j < m->n_markers; ++j) {
			int cell_sum = 0;
			int cell_sum2 = 0;
			if      (m->alleles[i][2*j] == allele)      { ++cell_sum; }
			else if (m->alleles[i][2*j] == allele2)     { ++cell_sum2;}
			if      (m->alleles[i][2*j + 1] == allele)  { ++cell_sum; }
			else if (m->alleles[i][2*j + 1] == allele2) { ++cell_sum2;}
			counts->matrix[i][j] = cell_sum;
			counts2->matrix[i][j] = cell_sum2;
		}
	}
	return 0;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in the linked list starting at the given gsc_AlleleMatrix.
 * Returns the result as a gsc_DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @shortnamed{calculate_full_count_matrix}
 *
 * @param m pointer to the gsc_AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the single-character allele to be counting.
 * @returns gsc_DecimalMatrix containing counts for every genotype in the gsc_AlleleMatrix
 * linked list.
 * */
gsc_DecimalMatrix gsc_calculate_full_count_matrix( const gsc_AlleleMatrix* m, const char allele) {
    const gsc_AlleleMatrix* currentm = m;
	int size = 0;
	do {
		size += currentm->n_genotypes;
	} while ((currentm = currentm->next) != NULL);

	gsc_DecimalMatrix counts = gsc_generate_zero_dmatrix(size, m->n_markers);

	currentm = m;
	int currenti = 0;
	for (int i = 0; i < size; ++i, ++currenti) {
		R_CheckUserInterrupt();
		if (currenti >= currentm->n_genotypes) {
			currenti = 0;
			currentm = currentm->next;
		}

        for (int j = 0; j < m->n_markers; ++j) {
			int cell_sum = 0;
			if (currentm->alleles[i][2*j] == allele)     cell_sum += 1;
			if (currentm->alleles[i][2*j + 1] == allele) cell_sum += 1;
			counts.matrix[i][j] = cell_sum;
		}
	}
	return counts;
}

/** Divide the genotype into blocks where each block contains all markers within
 * a 1/n length section of each chromosome in the map, and return the resulting
 * blocks in a struct.
 *
 * Chromosomes where there are no markers tracked are excluded/no blocks are created
 * in those chromosomes.
 *
 * Chromosomes where only one marker is tracked put the marker in the first block, and
 * have all other blocks empty.
 *
 * Empty blocks have blocks.num_markers_in_block[index] equal to 0 and contain a null pointer
 * at blocks.markers_in_block[index].
 *
 * Remember to call the gsc_MarkerBlocks destructor gsc_delete_markerblocks() on the returned
 * struct.
 *
 * @shortnamed{create_evenlength_blocks_each_chr}
 *
 * @param d pointer to the gsc_SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param mapid ID of the recombination map whose markers and chromosomes will be divided
 * into blocks, or NO_MAP to use the first-loaded/primary map by default.
 * @param mapid recombination map to use to determine block lengths and marker allocations to
 * blocks. If this value is NO_MAP, uses first/primary map.
 * @param n number of blocks into which to split each chromosome.
 * @returns a struct containing the markers identified as belonging to each block.
 */
gsc_MarkerBlocks gsc_create_evenlength_blocks_each_chr(const gsc_SimData* d, const gsc_MapID mapid, const int n) {
    gsc_MarkerBlocks blocks;
    blocks.num_blocks = 0;

    if (d->genome.n_maps < 1) {
        warning("Creating blocks by chromosome length requires at least one recombination map loaded\n");
        return blocks;
    }
    int mapix = 0;
    if (mapid.id != NO_MAP.id) { mapix = gsc_get_index_of_map(d, mapid); }
    if (mapix >= d->genome.n_maps) {
        warning("We don't have that recombination maps loaded. Using default map\n");
        mapix = 0;
    }
    gsc_RecombinationMap map = d->genome.maps[mapix];

    if (n < 1) {
        warning("Invalid n value: number of blocks must be positive.\n");
        return blocks;
    }
    if (map.n_chr < 1) {
        warning("Map has no chromosomes, so it cannot be divided into blocks.\n");
    }

    blocks.num_blocks = n * map.n_chr;
    blocks.num_markers_in_block = gsc_malloc_wrap(sizeof(*blocks.num_markers_in_block) * blocks.num_blocks,GSC_TRUE);
    blocks.markers_in_block = gsc_malloc_wrap(sizeof(*blocks.markers_in_block) * blocks.num_blocks,GSC_TRUE);
    for (int i = 0; i < blocks.num_blocks; ++i) {
        blocks.num_markers_in_block[i] = 0;
        blocks.markers_in_block[i] = NULL;
    }

    GSC_CREATE_BUFFER(temp_markers_in_block, int, 128);
    int bi = 0;

    for (int chr = 0; chr < map.n_chr; ++chr) {
        int current_block_filling = 0; //counter of how many blocks we have in this chr so far
        double chrpos = 0;
        bi = 0;

        // loop through each marker in this chromosome
        switch (map.chrs[chr].type) {
        case GSC_LINKAGEGROUP_SIMPLE:
            if (map.chrs[chr].map.simple.n_markers == 1) {
                int b = chr*n + 0;
                blocks.markers_in_block[b] = gsc_malloc_wrap(sizeof(*blocks.markers_in_block[b]), GSC_TRUE);
                blocks.markers_in_block[b][0] = map.chrs[chr].map.simple.first_marker_index;
                ++(blocks.num_markers_in_block[b]);
            } else {
                for (int i = 0; i < map.chrs[chr].map.simple.n_markers; ++i) {
                    R_CheckUserInterrupt();
                    chrpos += map.chrs[chr].map.simple.dists[i];
                    while (current_block_filling < n - 1 && chrpos > current_block_filling / n) {
                        int b = chr*n + current_block_filling;
                        if (blocks.num_markers_in_block[b] > 0) {
                            unsigned int bcapacity = sizeof(*blocks.markers_in_block[b])*blocks.num_markers_in_block[b];
                            blocks.markers_in_block[b] = gsc_malloc_wrap(bcapacity, GSC_TRUE);
                            memcpy(blocks.markers_in_block[b],temp_markers_in_block,bcapacity);
                        }

                        ++current_block_filling;
                        bi = 0;
                    }

                    // save marker
                    if (bi >= temp_markers_in_blockcap) {
                        GSC_STRETCH_BUFFER(temp_markers_in_block,2*bi);
                    }
                    temp_markers_in_block[bi] = map.chrs[chr].map.simple.first_marker_index + i;
                    ++(blocks.num_markers_in_block[chr*n + current_block_filling]);
                    ++bi;

                }

                int b = chr*n + current_block_filling;
                if (blocks.num_markers_in_block[b] > 0) {
                    unsigned int bcapacity = sizeof(*blocks.markers_in_block[b])*blocks.num_markers_in_block[b];
                    blocks.markers_in_block[b] = gsc_malloc_wrap(bcapacity, GSC_TRUE);
                    memcpy(blocks.markers_in_block[b],temp_markers_in_block,bcapacity);
                }
            }
            break;

        case GSC_LINKAGEGROUP_REORDER:
            if (map.chrs[chr].map.simple.n_markers == 1) {
                int b = chr*n + 0;
                blocks.markers_in_block[b] = gsc_malloc_wrap(sizeof(*blocks.markers_in_block[b]), GSC_TRUE);
                blocks.markers_in_block[b][0] = map.chrs[chr].map.reorder.marker_indexes[0];
                ++(blocks.num_markers_in_block[b]);
            } else {
                for (int i = 0; i < map.chrs[chr].map.reorder.n_markers; ++i) {
                    R_CheckUserInterrupt();
                    chrpos += map.chrs[chr].map.reorder.dists[i];
                    while (current_block_filling < n - 1 && chrpos > current_block_filling / n) {
                        int b = chr*n + current_block_filling;
                        if (blocks.num_markers_in_block[b] > 0) {
                            unsigned int bcapacity = sizeof(*blocks.markers_in_block[b])*blocks.num_markers_in_block[b];
                            blocks.markers_in_block[b] = gsc_malloc_wrap(bcapacity, GSC_TRUE);
                            memcpy(blocks.markers_in_block[b],temp_markers_in_block,bcapacity);
                        }

                        ++current_block_filling;
                        bi = 0;
                    }

                    // save marker
                    if (bi >= temp_markers_in_blockcap) {
                        GSC_STRETCH_BUFFER(temp_markers_in_block,2*bi);
                    }
                    temp_markers_in_block[bi] = map.chrs[chr].map.reorder.marker_indexes[i];
                    ++(blocks.num_markers_in_block[chr*n + current_block_filling]);
                    ++bi;

                }

                int b = chr*n + current_block_filling;
                if (blocks.num_markers_in_block[b] > 0) {
                    unsigned int bcapacity = sizeof(*blocks.markers_in_block[b])*blocks.num_markers_in_block[b];
                    blocks.markers_in_block[b] = gsc_malloc_wrap(bcapacity, GSC_TRUE);
                    memcpy(blocks.markers_in_block[b],temp_markers_in_block,bcapacity);
                }
            }
            break;
        }

	}

   GSC_DELETE_BUFFER(temp_markers_in_block);

	return blocks;
}

/** Given a file containing definitions of blocks of markers, process that file
 * and return a struct containing the definitions of those blocks.
 *
 * The block file is designed after the output from a call to the R SelectionTools
 * package's `st.def.hblocks` function. It should have the format (tab-separated):
 *
 * Chrom	Pos	Name	Class	Markers
 *
 * [ignored]	[ignored]	[ignored]	[ignored]	[semicolon];[separated];[list]
 * ;[of];[marker];[names];[belonging];[to];[this];[block];
 *
 * ...
 *
 * @shortnamed{load_blocks}
 *
 * @param d pointer to the gsc_SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param block_file string containing filename of the file with blocks
 * @returns a struct containing the markers identified as belonging to each block
 * according to their definitions in the file.
 */
gsc_MarkerBlocks gsc_load_blocks(const gsc_SimData* d, const char* block_file) {
	struct gsc_TableSize ts = gsc_get_file_dimensions(block_file, '\t');

	gsc_MarkerBlocks blocks;
	blocks.num_blocks = ts.num_rows - 1;
    blocks.num_markers_in_block = gsc_malloc_wrap(sizeof(int) * blocks.num_blocks,GSC_TRUE);
    blocks.markers_in_block = gsc_malloc_wrap(sizeof(int*) * blocks.num_blocks,GSC_TRUE);

	FILE* infile;
	if ((infile = fopen(block_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", block_file);
		//return blocks;
	}

    int bufferlen = d->genome.n_markers;
    GSC_CREATE_BUFFER(markername,char,CONTIG_WIDTH);
    GSC_CREATE_BUFFER(markerbuffer,int,bufferlen);
	int bi = 0; // block number

	// Ignore the first line
	fscanf(infile, "%*[^\n]\n");

	// Loop through rows of the file (each row corresponds to a block)
	while (fscanf(infile, "%*d %*f %*s %*s ") != EOF) {
	//for (int bi = 0; bi < n_blocks; ++bi) {

		// Indexes in play:
		//		bi: index in the blocks struct's arrays of the current block/line in the file
		//		ni: number of characters so far in the name of the next marker being read from the file
		//		mi: number of markers that have so far been read from the file for this block
		blocks.num_markers_in_block[bi] = 0;
		int c, ni = 0, mi = 0;

        memset(markerbuffer, 0, sizeof(int) * bufferlen);
		while ((c = fgetc(infile)) != EOF && c !='\n') {
			if (c == ';') {
				markername[ni] = '\0';

				// identify the index of this marker and save it in the temporary marker buffer `markerbuffer`
                int markerindex = gsc_get_from_unordered_str_list(markername, d->genome.n_markers, (const char**) d->genome.marker_names);
				if (markerindex >= 0) {
					++(blocks.num_markers_in_block[bi]);
					markerbuffer[mi] = markerindex;
					++mi;
				}

				ni = 0;
			} else {
				markername[ni] = c;
				++ni;
			}
		}

		// copy the markers belonging to this block into the struct
        blocks.markers_in_block[bi] = gsc_malloc_wrap(sizeof(int) * mi,GSC_TRUE);
		for (int i = 0; i < mi; ++i) {
			blocks.markers_in_block[bi][i] = markerbuffer[i];
		}

		++bi;
	}

    GSC_DELETE_BUFFER(markerbuffer);
    GSC_DELETE_BUFFER(markername);
	fclose(infile);
	return blocks;
}

/** Given a set of blocks of markers in a file, for each genotype in a group,
 * calculate the local fitness metric/breeding value for the first allele at
 * each marker in the block, and the local fitness metric/breeding value
 * for the second allele at each marker in the block, then save
 * the result to a file. This gives block effects/local BVs for each haplotype of each
 * individual in the group.
 *
 * The output file will have format:
 *
 * [genotype name]_1 [effect for first block, first allele]
 * [effect for second block, second allele] ...
 *
 * [genotype name]_2 [effect for first block, second allele] [effect for second block,
 * second allele] ...
 *
 * [genotype 2 name]_1 ...
 *
 *
 * @shortnamed{calculate_group_local_bvs}
 *
 * @param d pointer to the gsc_SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param b struct containing the blocks to use
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @param output_file string containing the filename of the file to which output
 * block effects/local BVs will be saved.
 * @param group group number from which to split the top individuals.
 */
void gsc_calculate_group_local_bvs(const gsc_SimData* d, const gsc_MarkerBlocks b, const gsc_EffectID effID, const char* output_file, const gsc_GroupNum group) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return;
    }
    gsc_EffectMatrix e = d->e[effIndex];

	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}

    GSC_CREATE_BUFFER(buffer,char,CONTIG_WIDTH);

	int gsize = gsc_get_group_size(d, group);
    GSC_CREATE_BUFFER(ggenos,char*,gsize);
    GSC_CREATE_BUFFER(gnames,char*,gsize);
    gsc_get_group_genes(d, group, gsize, ggenos);
    gsc_get_group_names(d, group, gsize, gnames);
    GSC_CREATE_BUFFER(gids,gsc_PedigreeID,gsize);
    gsc_get_group_ids(d,group,gsize,gids);

	double beffect;

	// for each group member
	for (int i = 0; i < gsize; ++i) {
		// for each block
        if (gnames[i] != NULL) {
            sprintf(buffer, "%s_1", gnames[i]);
        } else {
            sprintf(buffer, "%d_1", gids[i].id);
        }
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);

		// for each block
		for (int j = 0; j < b.num_blocks; ++j) {
			beffect = 0;

			// calculate the local BV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
                for (int q = 0; q < e.effects.rows; ++q) {
                    if (ggenos[i][2 * b.markers_in_block[j][k]] == e.effect_names[q]) {
                        beffect += e.effects.matrix[q][b.markers_in_block[j][k]];
					}
				}
			}

			// print the local BV
			fprintf(outfile, " %lf", beffect);
			fflush(outfile);
		}

        if (gnames[i] != NULL) {
            sprintf(buffer, "\n%s_2", gnames[i]);
        } else {
            sprintf(buffer, "\n%d_2", gids[i].id);
        }
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);

		// for each block for the second haplotype
		for (int j = 0; j < b.num_blocks; ++j) {
			beffect = 0;
			// calculate the local BV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
                for (int q = 0; q < e.effects.rows; ++q) {
                    if (ggenos[i][2 * b.markers_in_block[j][k] + 1] == e.effect_names[q]) {
                        beffect += e.effects.matrix[q][b.markers_in_block[j][k]];
					}
				}
			}

			// print the local BV
			fprintf(outfile, " %lf", beffect);
			fflush(outfile);
		}
		fwrite("\n", sizeof(char), 1, outfile);
	}

    GSC_DELETE_BUFFER(ggenos);
    GSC_DELETE_BUFFER(gnames);
    GSC_DELETE_BUFFER(gids);
    GSC_DELETE_BUFFER(buffer);
	fflush(outfile);
	fclose(outfile);
}

/** Given a set of blocks of markers in a file, for each genotype saved,
 * calculate the local BV for the first allele at each marker in the block, and
 * the local BV for the second allele at each marker in the block, then save
 * the result to a file. This gives block effects for each haplotype of each
 * individual currently saved to the gsc_SimData.
 *
 * The output file will have format:
 *
 * [genotype name]_1 [effect for first block, first allele]
 * [effect for second block, second allele] ...
 *
 * [genotype name]_2 [effect for first block, second allele] [effect for second block,
 * second allele] ...
 *
 * [genotype 2 name]_1 ...
 *
 *
 * @shortnamed{calculate_local_bvs}
 *
 * @param d pointer to the gsc_SimData object to which individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param b struct containing the blocks to use
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @param output_file string containing the filename of the file to which output
 * block effects/local BV will be saved.
 */
void gsc_calculate_local_bvs(const gsc_SimData* d, const gsc_MarkerBlocks b, const gsc_EffectID effID, const char* output_file) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return;
    }
    gsc_EffectMatrix e = d->e[effIndex];

	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}

    GSC_CREATE_BUFFER(buffer,char,CONTIG_WIDTH);

	int gsize = 0;
	gsc_AlleleMatrix* m = d->m;
	do {
		gsize += m->n_genotypes;
	} while ((m = m->next) != NULL);

	double beffect;

	// for each group member
	m = d->m;
	int total_i = 0;
	do {
		for (int i = 0; i < m->n_genotypes; ++i, ++total_i) {
            // for each group member
            if (m->names[i] != NULL) {
                sprintf(buffer, "%s_1",  m->names[i]);
            } else {
                sprintf(buffer, "%d_1", m->ids[i].id);
            }
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);

			// for each block
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;

				// calculate the local BV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
                    for (int q = 0; q < e.effects.rows; ++q) {
                        if (m->alleles[i][2 * b.markers_in_block[j][k]] == e.effect_names[q]) {
                            beffect += e.effects.matrix[q][b.markers_in_block[j][k]];
						}
					}
				}

				// print the local BV
				fprintf(outfile, " %lf", beffect);
				fflush(outfile);
			}

            if (m->names[i] != NULL) {
                sprintf(buffer, "\n%s_2",  m->names[i]);
            } else {
                sprintf(buffer, "\n%d_2", m->ids[i].id);
            }
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);

			// for each block for the second haplotype
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;
				// calculate the local BV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
                    for (int q = 0; q < e.effects.rows; ++q) {
                        if (m->alleles[i][2 * b.markers_in_block[j][k] + 1] == e.effect_names[q]) {
                            beffect += e.effects.matrix[q][b.markers_in_block[j][k]];
						}
					}
				}

				// print the local BV
				fprintf(outfile, " %lf", beffect);
				fflush(outfile);
			}
			fwrite("\n", sizeof(char), 1, outfile);
		}
	} while ((m = m->next) != NULL);

    GSC_DELETE_BUFFER(buffer);
	fflush(outfile);
	fclose(outfile);
}

/** Takes a look at the currently-loaded effect values and creates a string
 * containing the allele with the highest effect value for each marker as ordered
 * in the gsc_SimData. Returns this as a null-terminated heap array of characters of
 * length d->n_markers + one null byte.
 *
 * The return value should be freed when usage is finished!
 *
 * The gsc_SimData must be initialised with marker effects for this function to succeed.
 *
 * @shortnamed{calculate_optimal_haplotype}
 *
 * @param d pointer to the gsc_SimData containing markers and marker effects.
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @returns a heap array filled with a null-terminated string containing the highest
 * scoring allele at each marker.
 */
char* gsc_calculate_optimal_haplotype(const gsc_SimData* d, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
		return NULL;
	}
    gsc_EffectMatrix e = d->e[effIndex];

    char* optimal = gsc_malloc_wrap(sizeof(char)* (d->genome.n_markers + 1),GSC_TRUE);

    for (int i = 0; i < d->genome.n_markers; ++i) {
        char best_allele = e.effect_names[0];
        double best_score = e.effects.matrix[0][i];
        for (int a = 1; a < e.effects.rows; ++a) {
            if (e.effects.matrix[a][i] > best_score) {
                best_score = e.effects.matrix[a][i];
                best_allele = e.effect_names[a];
			}
		}
		optimal[i] = best_allele;
	}
    optimal[d->genome.n_markers] = '\0';
	return optimal;
}


/** Calculates the highest-breeding-value haplotype that can be created from the
 *  alleles present in a given group.
 *
 * The return value should be freed when usage is finished!
 *
 * The gsc_SimData must be initialised with marker effects for this function to succeed.
 *
 * @shortnamed{calculate_optimal_possible_haplotype}
 *
 * @param d pointer to the gsc_SimData containing markers and marker effects.
 * @param group group number to look at the available alleles in
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @returns a heap array filled with a null-terminated string containing the highest
 * scoring allele at each marker.
 */
char* gsc_calculate_optimal_possible_haplotype(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return NULL;
    }
    gsc_EffectMatrix e = d->e[effIndex];
    // assumes no alleles in the matrix are spaces.

    int gsize = gsc_get_group_size(d, group);
    GSC_CREATE_BUFFER(ggenes,char*,gsize);
    gsc_get_group_genes(d, group, gsize, ggenes);

    char* optimal = gsc_malloc_wrap(sizeof(char)* (d->genome.n_markers + 1),GSC_TRUE);

    // for each locus
    for (int j = 0; j < d->genome.n_markers; ++j) {
        char best_allele = '\0';
        double best_score;
        for (int i = 0; i < gsize; ++i) {

            // If the allele is different to the previous best (guaranteed if best_allele is not initialised)
            if (ggenes[i][2*j] != best_allele) {
                // Find it and see if it scores better.
                for (int a = 0; a < e.effects.rows; ++a) {

                    if (e.effect_names[a] == ggenes[i][2*j] &&
                            (best_allele == '\0' || e.effects.matrix[a][j] > best_score)) { // if it scores better than current best

                        best_allele = ggenes[i][2*j];
                        best_score = e.effects.matrix[a][j];

                        break;
                    }

                }
            }

            // Repeat for second allele of the group member at that locus
            if (ggenes[i][2*j + 1] != best_allele) {
                // Find it and see if it scores better.
                for (int a = 0; a < e.effects.rows; ++a) {

                    if (e.effect_names[a] == ggenes[i][2*j + 1] &&
                            (best_allele == '\0' || e.effects.matrix[a][j] > best_score)) { // if it scores better than current best

                        best_allele = ggenes[i][2*j + 1];
                        best_score = e.effects.matrix[a][j];

                        break;
                    }

                }
            }
        }
        optimal[j] = best_allele;
    }

    GSC_DELETE_BUFFER(ggenes);
    optimal[d->genome.n_markers] = '\0';
    return optimal;
}


/** Takes a look at the currently-loaded effect values and returns the highest possible
 * breeding value any (diploid) genotype could have using those effect values.
 *
 * The gsc_SimData must be initialised with marker effects for this function to succeed.
 *
 * @shortnamed{calculate_optimal_bv}
 *
 * @param d pointer to the gsc_SimData containing markers and marker effects.
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @returns the fitness metric/breeding value of the best/ideal genotype.
 */
double gsc_calculate_optimal_bv(const gsc_SimData* d, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return 0;
    }
    gsc_EffectMatrix e = d->e[effIndex];

	double best_gebv = 0;

    for (int i = 0; i < d->genome.n_markers; ++i) {
		// Find the allele with the highest effect
        double best_score = e.effects.matrix[0][i];
        for (int a = 1; a < e.effects.rows; ++a) {
            if (e.effects.matrix[a][i] > best_score) {
                best_score = e.effects.matrix[a][i];
			}
		}

		// add that highest allele to the score twice over
		best_gebv += (2*best_score);
	}

	return best_gebv;
}

/** Calculates the breeding value of the highest breeding-value genotype that can be
 *  created from the alleles present in a given group.
 *
 *  The highest value genotype is completely homozygous with the same alleles as
 *  the haplotype from gsc_calculate_optimal_possible_haplotype(), because of the additive
 *  model of trait effects.
 *
 * The gsc_SimData must be initialised with marker effects for this function to succeed.
 *
 * @shortnamed{calculate_optimal_possible_bv}
 *
 * @param d pointer to the gsc_SimData containing markers and marker effects.
 * @param  group number to look at the available alleles in
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @returns the fitness metric/breeding value of the best genotype
 */
double gsc_calculate_optimal_possible_bv(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return 0;
    }
    gsc_EffectMatrix e = d->e[effIndex];

    // assumes no alleles in the matrix are spaces.

    int gsize = gsc_get_group_size(d, group);
	if (gsize == 0) {
		warning("Nonexistent group %d\n", group.num);
        return 0;
	}
    GSC_CREATE_BUFFER(ggenes,char*,gsize);
    gsc_get_group_genes(d, group, gsize, ggenes);

    double total_score = 0;
    char best_allele;
    double best_score;

    // for each locus
    for (int j = 0; j < d->genome.n_markers; ++j) {
        best_allele = '\0';
		best_score = 0;
        for (int i = 0; i < gsize; ++i) {

            // If the allele is different to the previous best (guaranteed if best_allele is not initialised)
            if (ggenes[i][2*j] != best_allele) {
                // Find it and see if it scores better.
                for (int a = 0; a < e.effects.rows; ++a) {

                    if (e.effect_names[a] == ggenes[i][2*j] &&
                            (best_allele == '\0' || e.effects.matrix[a][j] > best_score)) { // if it scores better than current best

                        best_allele = ggenes[i][2*j];
                        best_score = e.effects.matrix[a][j];

                        break;
                    }

                }
            }

            // Repeat for second allele of the group member at that locus
            if (ggenes[i][2*j + 1] != best_allele) {
                // Find it and see if it scores better.
                for (int a = 0; a < e.effects.rows; ++a) {

                    if (e.effect_names[a] == ggenes[i][2*j + 1] &&
                            (best_allele == '\0' || e.effects.matrix[a][j] > best_score)) { // if it scores better than current best

                        best_allele = ggenes[i][2*j + 1];
                        best_score = e.effects.matrix[a][j];

                        break;
                    }

                }
            }
        }
        total_score += (2*best_score);
    }

    GSC_DELETE_BUFFER(ggenes);
    return total_score;
}

/** Takes a look at the currently-loaded effect values and returns the lowest possible
 * breeding value any (diploid) genotype could score using those effect values.
 *
 * The gsc_SimData must be initialised with marker effects for this function to succeed.
 *
 * @shortnamed{calculate_minimal_bv}
 *
 * @param d pointer to the gsc_SimData containing markers and marker effects.
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 * @returns the fitness metric/breeding value of the worst genotype.
 */
double gsc_calculate_minimal_bv(const gsc_SimData* d, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return 0;
    }
    gsc_EffectMatrix e = d->e[effIndex];

	double worst_gebv = 0;
	double worst_score;

    for (int i = 0; i < d->genome.n_markers; ++i) {
		// Find the allele with the highest effect
        worst_score = e.effects.matrix[0][i];
        for (int a = 1; a < e.effects.rows; ++a) {
            if (e.effects.matrix[a][i] < worst_score) {
                worst_score = e.effects.matrix[a][i];
			}
		}

		// add that highest allele to the score twice over
		worst_gebv += (2*worst_score);
	}

	return worst_gebv;
}

/*--------------------------------Printing-----------------------------------*/

/** Prints the markers contained in a set of blocks to a file. Column separators are tabs.
 *
 * The printing format is:
 *
 * Chrom	Pos	Name	Class	Markers
 *
 * 0	0	b0	b	m1;m2;m3;m4;
 *
 * 0	0	b0	b	m7;m9;
 *
 * ...
 *
 * where m1, m2, m3, m4 are the names of the markers in the first block and
 * m7 and m9 are the names of the markers in the second block. Note that chromosome
 * numbers and positions are not actually specified, just replaced with 0. Only the last
 * column is meaningful.
 *
 * @shortnamed{save_markerblocks}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData whose data we print
 * @param b gsc_MarkerBlocks struct containing the groupings of markers to print.
*/
void gsc_save_markerblocks(FILE* f, const gsc_SimData* d, const gsc_MarkerBlocks b) {
	const char header[] = "Chrom\tPos\tName\tClass\tMarkers\n";
	fwrite(header, sizeof(char)*strlen(header), 1, f);

	// for the moment we do not name or give locations of different blocks
	const char unspecified[] = "0\t0\tb0\tb\t";
	const int unspeci_length = strlen(unspecified);

	for (int i = 0; i < b.num_blocks; ++i) {
		fwrite(unspecified, sizeof(char)*unspeci_length, 1, f);

		for (int j = 0; j < b.num_markers_in_block[i]; ++j) {
			int k = b.markers_in_block[i][j];

            fwrite(d->genome.marker_names[k], sizeof(char)*strlen(d->genome.marker_names[k]), 1, f);
            fputc(';',f);
		}

		fwrite("\n", sizeof(char), 1, f);
	}

	fflush(f);
	return;

}

/** Prints a list of names as a tab-separated row to a file.
 *
 * Used, for example, for printing marker names as a header
 * in a genotype output file.
 *
 * @shortnamed{save_names_header}
 *
 * @param f file pointer opened for writing to put the output
 * @param n number of names to print
 * @param names list of strings/names to print. Must contain at least [n] names.
*/
void gsc_save_names_header(FILE* f, unsigned int n, const char** names) {
    if (names == NULL) return;
    for (int i = 0; i < n; ++i) {
        fwrite("\t", sizeof(char), 1, f);
        if (names[i] != NULL) { // assume all-or-nothing with marker names
            fwrite(names[i], sizeof(char), strlen(names[i]), f);
        }
    }
    fwrite("\n", sizeof(char), 1, f);
}

/** Prints all the genotype data saved in the linked list of AlleleMatrices
 * starting with `m` to a file. Uses the following format:
 *
 * [id]OR[name]	[allele pairs for each marker]
 *
 * [id]OR[name]	[allele pairs for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved
 *
 * @shortnamed{save_allele_matrix}
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the gsc_AlleleMatrix whose data we print
*/
void gsc_save_allele_matrix(FILE* f, const gsc_AlleleMatrix* m) {
	do {
		for (int i = 0; i < m->n_genotypes; ++i) {
			// print the name or ID of the individual.
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			} else {
				//fwrite(group_contents + i, sizeof(int), 1, f);
                fprintf(f, "%d", m->ids[i].id);
			}


			for (int j = 0; j < m->n_markers; ++j) {
				//fprintf(f, "\t%c%c", m->alleles[j][2*i], m->alleles[j][2*i + 1]);
				fwrite("\t", sizeof(char), 1, f);
				fwrite(m->alleles[i] + 2*j, sizeof(char), 1, f);
				fwrite(m->alleles[i] + 2*j + 1, sizeof(char), 1, f);
			}
			///fprintf(f, "\n");
			fwrite("\n", sizeof(char), 1, f);
		}
	} while ((m = m->next) != NULL);

	fflush(f);

}

/** Prints all the gene data saved in the linked list starting with `m` to the
 * file. Uses the following format:
 *
 * 		[id]OR[name]	[id]OR[name] ...
 *
 * [marker name]	[allele pairs for each marker]
 *
 * [marker name]	[allele pairs for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved
 *
 * @shortnamed{save_transposed_allele_matrix}
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the gsc_AlleleMatrix whose data we print
 * @param markers array of strings that correspond to names of the markers.
*/
void gsc_save_transposed_allele_matrix(FILE* f, const gsc_AlleleMatrix* m, const char** markers) {
	// Count number of genotypes in the AM
    const gsc_AlleleMatrix* currentm = m; // current matrix
	int tn_genotypes = 0;
	do {
		tn_genotypes += currentm->n_genotypes;
	} while ((currentm = currentm->next) != NULL);

	currentm = m;

	/* Print header */
	for (int i = 0, currenti = 0; i < tn_genotypes; ++i, ++currenti) {
		if (currenti >= currentm->n_genotypes) {
			currenti = 0;
			currentm = currentm->next;
		}
		if (currentm->names[currenti] != NULL) { // assume all-or-nothing with marker names
			//fprintf(f, "\t%s", markers[i]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(currentm->names[currenti], sizeof(char), strlen(currentm->names[currenti]), f);
		} else {
            fprintf(f, "\t%d", currentm->ids[currenti].id);
		}

	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

	for (int j = 0; j < m->n_markers; ++j) {
		if (markers != NULL && markers[j] != NULL) {
			fwrite(markers[j], sizeof(char), strlen(markers[j]), f);
		}

		currentm = m;

		for (int i = 0, currenti = 0; i < tn_genotypes; ++i, ++currenti) {
			if (currenti >= currentm->n_genotypes) {
				currenti = 0;
				currentm = currentm->next;
			}

			fwrite("\t", sizeof(char), 1, f);
			fwrite(currentm->alleles[currenti] + 2*j, sizeof(char), 1, f);
			fwrite(currentm->alleles[currenti] + 2*j + 1, sizeof(char), 1, f);
		}

		fwrite("\n", sizeof(char), 1, f);
	}

	fflush(f);
}

/** Prints the genotypes of each individual in a given group to a file, with
 * the following format.
 *
 * 		[marker name]	[marker name]
 *
 * [id]OR[name]	[allele pairs for each marker]
 *
 * [id]OR[name]	[allele pairs for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved
 *
 * @shortnamed{save_group_genotype}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the genotypes of the group and
 * the marker names.
 * @param group_id group number of the group of individuals whose genotypes to print.
*/
void gsc_save_group_genotypes(FILE* f, gsc_SimData* d, gsc_GroupNum group_id) {
	/* Get the stuff we'll be printing. */
	int group_size = gsc_get_group_size( d, group_id);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group_id.num);
        return;
    }
    GSC_CREATE_BUFFER(names,char*,group_size);
    GSC_CREATE_BUFFER(ids,gsc_PedigreeID,group_size);
    gsc_get_group_names( d, group_id, group_size, names );
    gsc_get_group_ids( d, group_id, group_size, ids );

	/* Print header */
	//fwrite(&group_id, sizeof(int), 1, f);
    fprintf(f, "%d", group_id.num);
    if (d->genome.marker_names != NULL) {
        for (int i = 0; i < d->genome.n_markers; ++i) {
			// assume all-or-nothing with marker names
			//fprintf(f, "\t%s", markers[i]);
			fwrite("\t", sizeof(char), 1, f);
            fwrite(d->genome.marker_names[i], sizeof(char), strlen(d->genome.marker_names[i]), f);
		}
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

    GSC_CREATE_BUFFER(alleles,char*,group_size);
    gsc_get_group_genes( d, group_id, group_size, alleles );

	/* Print the body */
	for (int i = 0; i < group_size; ++i) {
		// print the name or ID of the individual.
		if (names[i] != NULL) {
			fwrite(names[i], sizeof(char), strlen(names[i]), f);
		} else {
			//fwrite(group_contents + i, sizeof(int), 1, f);
            fprintf(f, "%d", ids[i].id);
		}

        for (int j = 0; j < d->genome.n_markers; ++j) {
			//fprintf(f, "\t%c%c", m->alleles[j][2*i], m->alleles[j][2*i + 1]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(alleles[i] + 2*j, sizeof(char), 1, f);
			fwrite(alleles[i] + 2*j + 1, sizeof(char), 1, f);
		}
		///fprintf(f, "\n");
		fwrite("\n", sizeof(char), 1, f);
	}
	fflush(f);
    GSC_DELETE_BUFFER(alleles);
    GSC_DELETE_BUFFER(names);
    GSC_DELETE_BUFFER(ids);
}

/** Prints the genotypes of each individual in a given group to a file, with
 * the following format.
 *
* 		[id]OR[name]	[id]OR[name] ...
 *
 * [marker name]	[allele pairs for each marker]
 *
 * [marker name]	[allele pairs for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved
 *
 * @shortnamed{save_transposed_group_genotype}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the genotypes of the group and
 * the marker names.
 * @param group_id group number of the group of individuals whose genotypes to print.
*/
void gsc_save_transposed_group_genotypes(FILE* f, const gsc_SimData* d, const gsc_GroupNum group_id) {
	/* Get the stuff we'll be printing. */
	int group_size = gsc_get_group_size( d, group_id);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group_id.num);
        return;
    }

    GSC_CREATE_BUFFER(names,char*,group_size);
    GSC_CREATE_BUFFER(ids,gsc_PedigreeID,group_size);
    gsc_get_group_names( d, group_id, group_size, names );
    gsc_get_group_ids( d, group_id, group_size, ids );

	/* Print header */
    fprintf(f, "%d", group_id.num);
	for (int i = 0; i < group_size; ++i) {
		fwrite("\t", sizeof(char), 1, f);
        if (names[i] != NULL) {
            fwrite(names[i], sizeof(char), strlen(names[i]), f);
        } else {
            //fwrite(group_contents + i, sizeof(int), 1, f);
            fprintf(f, "%d", ids[i].id);
        }
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

    GSC_DELETE_BUFFER(names);
    GSC_DELETE_BUFFER(ids);
    GSC_CREATE_BUFFER(alleles,char*,group_size);
    gsc_get_group_genes( d, group_id, group_size, alleles );

	/* Print the body */
    for (int i = 0; i < d->genome.n_markers; ++i) {
		// print the name or ID of the individual.
        if (d->genome.marker_names != NULL && d->genome.marker_names[i] != NULL) {
            fwrite(d->genome.marker_names[i], sizeof(char), strlen(d->genome.marker_names[i]), f);
		}

		for (int j = 0; j < group_size; ++j) {
			fwrite("\t", sizeof(char), 1, f);
			fwrite(alleles[j] + 2*i, sizeof(char), 1, f);
			fwrite(alleles[j] + 2*i + 1, sizeof(char), 1, f);
		}
		///fprintf(f, "\n");
		fwrite("\n", sizeof(char), 1, f);
	}
	fflush(f);
    GSC_DELETE_BUFFER(alleles);
}


/** Print the number of copies of a particular allele at each marker of each genotype
 * in the gsc_SimData to a file. The following tab-separated format is used:
 *
 * 		[id]OR[name]	[id]OR[name] ...
 *
 * [marker name]	[allele count for each marker]
 *
 * [marker name]	[allele count for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved.
 *
 * @shortnamed{save_count_matrix}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param allele the allele character to count
 */
void gsc_save_count_matrix(FILE* f, const gsc_SimData* d, const char allele) {
	gsc_DecimalMatrix counts = gsc_calculate_full_count_matrix(d->m, allele);

	gsc_AlleleMatrix* currentm = d->m;
	// print the header
	for (int i = 0, currenti = 0; i < counts.rows; ++i, ++currenti) {
		if (currenti >= currentm->n_genotypes) {
			currenti = 0;
			currentm = currentm->next;
		}
		fwrite("\t", sizeof(char), 1, f);
		if (currentm->names[currenti] != NULL) {
			fwrite(currentm->names[currenti], sizeof(char), strlen(currentm->names[currenti]), f);
        } else {
            fprintf(f, "%d", currentm->ids[currenti].id);
        }
	}

	fwrite("\n", sizeof(char), 1, f);

	// Print the body
    for (int i = 0; i < d->genome.n_markers; ++i) { // loop through markers
        if (d->genome.marker_names != NULL && d->genome.marker_names[i] != NULL) {
            fwrite(d->genome.marker_names[i], sizeof(char), strlen(d->genome.marker_names[i]), f);
		}

		for (int j = 0; j < counts.rows; ++j) { // loop through genotypes
			// print the matrix entries
            fwrite("\t", sizeof(char), 1, f);
            fprintf(f, "%d", (int) counts.matrix[j][i]);
		}
		//print the newline
		if (i + 1 < counts.rows) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}

	gsc_delete_dmatrix(&counts);
	fflush(f);
}

/** Print the number of copies of a particular allele at each marker of each genotype
 * in a group to a file. The following tab-separated format is used:
 *
 * 		[id]OR[name]	[id]OR[name] ...
 *
 * [marker name]	[allele count for each marker]
 *
 * [marker name]	[allele count for each marker]
 *
 * ...
 *
 * ID will be printed if the genotype does not have a name saved.
 *
 * @shortnamed{save_group_count_matrix}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * allele count of.
 * @param allele the allele character to count
 */
void gsc_save_group_count_matrix(FILE* f, const gsc_SimData* d, const char allele, const gsc_GroupNum group) {
	unsigned int group_size = gsc_get_group_size( d, group);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group.num);
        return;
    }

    GSC_CREATE_BUFFER(group_names,char*,group_size);
    gsc_get_group_names( d, group, group_size, group_names );
    GSC_CREATE_BUFFER(group_ids,gsc_PedigreeID,group_size);
    gsc_get_group_ids(d,group,group_size, group_ids);

    fprintf(f, "%d", group.num);
	// print the header
	for (int i = 0; i < group_size; ++i) {
		fwrite("\t", sizeof(char), 1, f);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
        } else {
            fprintf(f, "%d", group_ids[i].id);
        }
	}
    GSC_DELETE_BUFFER(group_names);
    GSC_DELETE_BUFFER(group_ids);
    gsc_DecimalMatrix counts = gsc_generate_zero_dmatrix(group_size,d->genome.n_markers);
    gsc_calculate_group_count_matrix(d,group,allele,&counts);

	fwrite("\n", sizeof(char), 1, f);

	// Print the body
    for (int i = 0; i < d->genome.n_markers; ++i) { // loop through markers
        if (d->genome.marker_names != NULL && d->genome.marker_names[i] != NULL) {
            fwrite(d->genome.marker_names[i], sizeof(char), strlen(d->genome.marker_names[i]), f);
		}

		for (int j = 0; j < group_size; ++j) { // loop through genotypes
			// print the matrix entries
            fwrite("\t", sizeof(char), 1, f);
            fprintf(f, "%d", (int) counts.matrix[j][i]);
		}
		//print the newline
        if (i + 1 < counts.rows) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}

	gsc_delete_dmatrix(&counts);
	fflush(f);
}

/** Print the parents of each genotype in the gsc_SimData to a file. The following
 * tab-separated format is used:
 *
 * [genotype name]	[parent 1 name]	[parent 2 name]
 *
 * [genotype name]	[parent 1 name]	[parent 2 name]
 *
 * ...
 *
 * The parents are identified by the two ids saved in the genotype's
 * pedigrees field in the gsc_AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id is 0 (which means the parent is unknown) no
 * name or id is printed for that parent.
 *
 * @shortnamed{save_one_step_pedigree}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the genotypes and their pedigrees
 */
void gsc_save_one_step_pedigree(FILE* f, const gsc_SimData* d) {
	char* name;
	gsc_AlleleMatrix* m = d->m;

	do {
		for (int i = 0; i < m->n_genotypes; ++i) {
			/*Group member name*/
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			} else {
                fprintf(f, "%d", m->ids[i].id);
			}
			fwrite("\t", sizeof(char), 1, f);


            /* Parent 1 */
            if (m->pedigrees[0][i].id != GSC_NO_PEDIGREE.id) {
                name = gsc_get_name_of_id( d->m, m->pedigrees[0][i] );
                if (name != NULL) {
                    fwrite(name, sizeof(char), strlen(name), f);
                } else {
                    fprintf(f, "%d", m->pedigrees[0][i].id);
                }
            }
            fwrite("\t", sizeof(char), 1, f);

            /* Parent 2 */
            if (m->pedigrees[1][i].id != GSC_NO_PEDIGREE.id) {
                name = gsc_get_name_of_id( d->m, m->pedigrees[1][i]);
                if (name != NULL) {
                    fwrite(name, sizeof(char), strlen(name), f);
                } else {
                    fprintf(f, "%d", m->pedigrees[1][i].id);
                }
            }

			fwrite("\n", sizeof(char), 1, f);
		}
	} while ((m = m->next) != NULL);
	fflush(f);
}

/** Print the parents of each genotype in a group to a file. The following
 * tab-separated format is used:
 *
 * [group member name]	[parent 1 name]	[parent 2 name]
 *
 * [group member name]	[parent 1 name]	[parent 2 name]
 *
 * ...
 *
 * The parents are identified by the two ids saved in the genotype's
 * pedigrees field in the gsc_AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0 (which means the parent is unknown) no
 * name or id is printed for that parent.
 *
 * @shortnamed{save_group_one_step_pedigree}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * immediate parents of.
 */
void gsc_save_group_one_step_pedigree(FILE* f, const gsc_SimData* d, const gsc_GroupNum group) {
	int group_size = gsc_get_group_size( d, group);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group.num);
        return;
    }

    GSC_CREATE_BUFFER(group_contents,gsc_PedigreeID,group_size);
    gsc_get_group_ids( d, group, group_size, group_contents );
    GSC_CREATE_BUFFER(group_names,char*,group_size);
    gsc_get_group_names( d, group, group_size, group_names );
    GSC_CREATE_BUFFER(parent1s,gsc_PedigreeID,group_size);
    GSC_CREATE_BUFFER(parent2s,gsc_PedigreeID,group_size);
    gsc_get_group_parent_ids(d, group, group_size, 1, parent1s );
    gsc_get_group_parent_ids(d, group, group_size, 2, parent2s );
	char* name;

	for (int i = 0; i < group_size; i++) {
		/*Group member name*/
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		} else {
			//fwrite(group_contents + i, sizeof(int), 1, f);
            fprintf(f, "%d", group_contents[i].id);
		}
		fwrite("\t", sizeof(char), 1, f);

        // Prints both parents, even if they're the same one.
        /* Parent 1 */
        if (parent1s[i].id != GSC_NO_PEDIGREE.id) {
            name = gsc_get_name_of_id( d->m, parent1s[i]);
            if (name != NULL) {
                fwrite(name, sizeof(char), strlen(name), f);
            } else {
                fprintf(f, "%d", parent1s[i].id);
            }
        }
        fwrite("\t", sizeof(char), 1, f);

        /* Parent 2 */
        if (parent2s[i].id != GSC_NO_PEDIGREE.id) {
            name = gsc_get_name_of_id( d->m, parent2s[i]);
            if (name != NULL) {
                fwrite(name, sizeof(char), strlen(name), f);
            } else {
                fprintf(f, "%d", parent2s[i].id);
            }
        }

		fwrite("\n", sizeof(char), 1, f);
	}
	fflush(f);
    GSC_DELETE_BUFFER(parent1s);
    GSC_DELETE_BUFFER(parent2s);
    GSC_DELETE_BUFFER(group_contents);
    GSC_DELETE_BUFFER(group_names);
}

/** Print the full known pedigree of each genotype in the gsc_SimData
 * to a file. The following
 * tab-separated format is used:
 *
 * [id]	[name]=([parent 1 pedigree],[parent 2 pedigree])
 *
 * [id]	[name]=([parent pedigree])
 *
 * ...
 *
 * Note that this pedigree is recursively constructed, so if a genotype's
 * parents are known, [parent pedigree] is replaced with a string of format
 * name=([parent1 pedigree],[parent2 pedigree]) alike. If the two parents of
 * a genotype are the same individual (i.e. it was produced by selfing, the
 * comma and second parent pedigree are ommitted, as in the second line sample
 * in the format above.
 *
 * The parents of a genotype are identified by the two ids saved in the
 * genotype's pedigrees field in the gsc_AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * @shortnamed{save_full_pedigree}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing all genotypes to print.
 */
void gsc_save_full_pedigree(FILE* f, const gsc_SimData* d) {
	const char newline[] = "\n";

	gsc_AlleleMatrix* m = d->m;

	do {
		for (int i = 0; i < m->n_genotypes; ++i) {
			/*Group member name*/
            fprintf(f, "%d\t", m->ids[i].id);
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			}

            if (m->pedigrees[0][i].id != GSC_NO_PEDIGREE.id || m->pedigrees[1][i].id != GSC_NO_PEDIGREE.id) {
				gsc_save_parents_of(f, d->m, m->pedigrees[0][i], m->pedigrees[1][i]);
			}
			fwrite(newline, sizeof(char), 1, f);
		}
	} while ((m = m->next) != NULL);
	fflush(f);
}

/** Print the full known pedigree of each genotype in a group to a file. The following
 * tab-separated format is used:
 *
 * [id]	[name]=([parent 1 pedigree],[parent 2 pedigree])
 *
 * [id]	[name]=([parent pedigree])
 *
 * ...
 *
 * Note that this pedigree is recursively constructed, so if a genotype's
 * parents are known, [parent pedigree] is replaced with a string of format
 * name=([parent1 pedigree],[parent2 pedigree]) alike. If the two parents of
 * a genotype are the same individual (i.e. it was produced by selfing, the
 * comma and second parent pedigree are ommitted, as in the second line sample
 * in the format above.
 *
 * The parents of a genotype are identified by the two ids saved in the
 * genotype's pedigrees field in the gsc_AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * @shortnamed{save_group_full_pedigree}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * pedigree of.
 */
void gsc_save_group_full_pedigree(FILE* f, const gsc_SimData* d, const gsc_GroupNum group) {
	int group_size = gsc_get_group_size( d, group);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group.num);
        return;
    }

    GSC_CREATE_BUFFER(group_contents,gsc_PedigreeID,group_size);
    gsc_get_group_ids( d, group, group_size, group_contents );
    GSC_CREATE_BUFFER(group_names,char*,group_size);
    gsc_get_group_names( d, group, group_size, group_names );
	const char newline[] = "\n";
    GSC_CREATE_BUFFER(parent1s,gsc_PedigreeID,group_size);
    GSC_CREATE_BUFFER(parent2s,gsc_PedigreeID,group_size);
    gsc_get_group_parent_ids(d, group, group_size, 1, parent1s );
    gsc_get_group_parent_ids(d, group, group_size, 2, parent2s );

	for (int i = 0; i < group_size; i++) {
		/*Group member name*/
        fprintf(f, "%d\t", group_contents[i].id);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}

        if (parent1s[i].id != GSC_NO_PEDIGREE.id || parent2s[i].id != GSC_NO_PEDIGREE.id) {
            gsc_save_parents_of(f, d->m, parent1s[i], parent2s[i]);
		}
		fwrite(newline, sizeof(char), 1, f);
	}
	fflush(f);
    GSC_DELETE_BUFFER(parent1s);
    GSC_DELETE_BUFFER(parent2s);
    GSC_DELETE_BUFFER(group_contents);
    GSC_DELETE_BUFFER(group_names);
}

/** Print the full known pedigree of each genotype in a single gsc_AlleleMatrix
 * to a file. The following
 * tab-separated format is used:
 *
 * [id]	[name]=([parent 1 pedigree],[parent 2 pedigree])
 *
 * [id]	[name]=([parent pedigree])
 *
 * ...
 *
 * Note that this pedigree is recursively constructed, so if a genotype's
 * parents are known, [parent pedigree] is replaced with a string of format
 * name=([parent1 pedigree],[parent2 pedigree]) alike. If the two parents of
 * a genotype are the same individual (i.e. it was produced by selfing, the
 * comma and second parent pedigree are ommitted, as in the second line sample
 * in the format above.
 *
 * The parents of a genotype are identified by the two ids saved in the
 * genotype's pedigrees field in the gsc_AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * Note this does not follow through the linked list of gsc_AlleleMatrix.
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the gsc_AlleleMatrix containing the genotypes to print
 * @param parents pointer to an gsc_AlleleMatrix that heads the linked list
 * containing the parents and other ancestry of the given genotypes.
 */
void gsc_save_allelematrix_full_pedigree(FILE* f, const gsc_AlleleMatrix* m, const gsc_SimData* parents) {
	const char newline[] = "\n";

	for (int i = 0; i < m->n_genotypes; ++i) {
		/*Group member name*/
        fprintf(f, "%d\t", m->ids[i].id);
		if (m->names[i] != NULL) {
			fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
		}

        if (m->pedigrees[0][i].id != GSC_NO_PEDIGREE.id || m->pedigrees[1][i].id != GSC_NO_PEDIGREE.id) {
            gsc_save_parents_of(f, parents->m, m->pedigrees[0][i], m->pedigrees[1][i]);
        }
		fwrite(newline, sizeof(char), 1, f);
	}
	fflush(f);
}

/** Recursively save the parents of a particular id to a file.
 *
 * It saves using the following format:
 *
 * - no characters are saved if the parents of the id are unknown/0
 *
 * - if the id has one parent, repeated twice (it was produced by selfing),
 * print "=(parentname)", where parentname is the name of the parent or its
 * id if it does not have one, followed by whatever is printed by a call
 * to this function on the parent's id.
 *
 * - if the id has two separate parents, print "=(parent1name,parent2name)",
 * where parent1name and parent2name are the name2 of the two parents or their
 * ids if they does not have names, each name immediately followed by whatever
 * is printed by a call to this function on the corresponding parent's id.
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to an gsc_AlleleMatrix that heads the linked list
 * containing the parents and other ancestry of the given id.
 * @param p1 the session-unique id of the first parent to be saved.
 * @param p2 the session-unique id of the second parent to be saved.
 */
void gsc_save_parents_of(FILE* f, const gsc_AlleleMatrix* m, gsc_PedigreeID p1, gsc_PedigreeID p2) {
    gsc_PedigreeID pedigree[2];

	// open brackets
	fwrite("=(", sizeof(char), 2, f);
	char* name;

	// enables us to print only the known parent if one is unknown
    if (p1.id == GSC_NO_PEDIGREE.id || p2.id == GSC_NO_PEDIGREE.id) {
        p1.id = (p1.id >= p2.id) ? p1.id : p2.id; //max of the two
        p2.id = p1.id;
	}

    if (p1.id == p2.id) {
        if (p1.id != GSC_NO_PEDIGREE.id) { //print nothing if both are unknown.
			// Selfed parent
            name = gsc_get_name_of_id( m, p1);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
            } else if (p1.id != GSC_NO_PEDIGREE.id) {
                fprintf(f, "%d", p1.id);
				//fwrite(pedigree, sizeof(int), 1, f);
			}

            if (gsc_get_parents_of_id(m, p1, pedigree) == 0) {
				gsc_save_parents_of(f, m, pedigree[0], pedigree[1]);
			}
		}
	} else {
		// Parent 1
		name = gsc_get_name_of_id( m, p1);
		if (name != NULL) {
			fwrite(name, sizeof(char), strlen(name), f);
        } else if (p1.id != GSC_NO_PEDIGREE.id) {
            fprintf(f, "%d", p1.id);
			//fwrite(pedigree, sizeof(int), 1, f);
		}
		if (gsc_get_parents_of_id(m, p1, pedigree) == 0) {
			gsc_save_parents_of(f, m, pedigree[0], pedigree[1]);
		}

		// separator
		fwrite(",", sizeof(char), 1, f);

		// Parent 2
		name = gsc_get_name_of_id( m, p2);
		if (name != NULL) {
			fwrite(name, sizeof(char), strlen(name), f);
        } else if (p2.id != GSC_NO_PEDIGREE.id) {
            fprintf(f, "%d", p2.id);
			//fwrite(pedigree + 1, sizeof(int), 1, f);
		}

		if (gsc_get_parents_of_id(m, p2, pedigree) == 0) {
			gsc_save_parents_of(f, m, pedigree[0], pedigree[1]);
		}

	}

	// close brackets
	fwrite(")", sizeof(char), 1, f);
}

/** Print the breeding value of each genotype in the gsc_SimData to a file. The following
 * tab-separated format is used:
 *
 * [id]	[name]	[BV]
 *
 * [id]	[name]	[BV]
 *
 * ...
 *
 * The gsc_SimData must have loaded marker effects for this function to succeed.
 *
 * @shortnamed{save_bvs}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 */
void gsc_save_bvs(FILE* f, const gsc_SimData* d, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return;
    }

	gsc_AlleleMatrix* am = d->m;
	const char newline[] = "\n";
	const char tab[] = "\t";
	gsc_DecimalMatrix effects;

	do {
        effects = gsc_calculate_bvs(am, d->e + effIndex);
		for (int i = 0; i < effects.cols; ++i) {
			/*Group member name*/
			//fwrite(group_contents + i, sizeof(int), 1, f);
            fprintf(f, "%d", am->ids[i].id);
			fwrite(tab, sizeof(char), 1, f);
			if (am->names[i] != NULL) {
				fwrite(am->names[i], sizeof(char), strlen(am->names[i]), f);
			}
			fwrite(tab, sizeof(char), 1, f);
			//fwrite(effects.matrix[0], sizeof(float), 1, f);
			fprintf(f, "%f", effects.matrix[0][i]);
			fwrite(newline, sizeof(char), 1, f);
		}
		gsc_delete_dmatrix(&effects);
	} while ((am = am->next) != NULL);
	fflush(f);
}

/** Print the breeding value of each genotype in a group to a file. The following
 * tab-separated format is used:
 *
 * [group member id]	[group member name]	[BV]
 *
 * [group member id]	[group member name]	[BV]
 *
 * ...
 *
 * The gsc_SimData must have loaded marker effects for this function to succeed.
 *
 * @shortnamed{save_group_bvs}
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the gsc_SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * breeding values of.
 * @param effID Identifier of the marker effect set to be used to calculate these breeding values
 */
void gsc_save_group_bvs(FILE* f, const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID) {
    const int effIndex = gsc_get_index_of_eff_set(d, effID);
    if (effIndex == GSC_UNINIT) {
        warning("Nonexistent effect set with id %d\n", effID.id);
        return;
    }

	int group_size = gsc_get_group_size( d, group);
    if (group_size < 1) {
        warning("Group %d does not exist: no data saved.\n", group.num);
        return;
    }
    GSC_CREATE_BUFFER(group_contents,gsc_PedigreeID,group_size);
    gsc_get_group_ids( d, group, group_size, group_contents );
    GSC_CREATE_BUFFER(group_names,char*,group_size);
    gsc_get_group_names( d, group, group_size, group_names );
    gsc_DecimalMatrix effects = gsc_calculate_group_bvs(d, group, effID);
	const char newline[] = "\n";
	const char tab[] = "\t";

	for (int i = 0; i < group_size; ++i) {
		/*Group member name*/
		//fwrite(group_contents + i, sizeof(int), 1, f);
        fprintf(f, "%d", group_contents[i].id);
		fwrite(tab, sizeof(char), 1, f);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}
		fwrite(tab, sizeof(char), 1, f);
		//fwrite(effects.matrix[0], sizeof(float), 1, f);
		fprintf(f, "%f", effects.matrix[0][i]);
		fwrite(newline, sizeof(char), 1, f);
	}

	gsc_delete_dmatrix(&effects);
    GSC_DELETE_BUFFER(group_contents);
    GSC_DELETE_BUFFER(group_names);
	fflush(f);
}

/** Print the provided breeding values of each provided name and id to a file,
 * with the same format as a regular call to `gsc_save_bvs` or `gsc_save_group_bvs`.
 * The following tab-separated format is used:
 *
 * [id]	[name]	[BV]
 *
 * [id]	[name]	[BV]
 *
 * ...
 *
 * The function assumes the array of ids, array of names, and columns of the
 * gsc_DecimalMatrix are ordered the same, so a single index produces the corresponding
 * value from each. It is used by as-you-go savers within crossing functions. For
 * general use consider `gsc_save_bvs` or `gsc_save_group_bvs` instead.
 *
 * @param f file pointer opened for writing to put the output
 * @param e pointer to the gsc_DecimalMatrix containing the BVs in the first row.
 * @param ids array of ids to print alongside the BVs.
 * @param names array of names to print alongside the BVs.
 */
void gsc_save_manual_bvs(FILE* f, const gsc_DecimalMatrix* e, const gsc_PedigreeID* ids, const char** names) {
	char sep[] = "\t";
	char newline[] = "\n";

	for (int i = 0; i < e->cols; ++i) {
		//fwrite(ids + i, sizeof(int), 1, f);
        fprintf(f, "%d", ids[i].id);
		fwrite(sep, sizeof(char), 1, f);
		if (names != NULL && names[i] != NULL) {
			fwrite(names[i], sizeof(char), strlen(names[i]), f);
		}
		fwrite(sep, sizeof(char), 1, f);
		//fwrite(e->matrix[i], sizeof(double), 1, f);
		fprintf(f, "%f", e->matrix[0][i]);

		//print the newline
		fwrite(newline, sizeof(char), 1, f);
	}
	fflush(f);
}

#endif
