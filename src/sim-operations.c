#include "sim-operations.h"

/** Options parameter to run SimData functions in their bare-bones form.*/
const GenOptions BASIC_OPT = {
	.will_name_offspring = FALSE,
	.offspring_name_prefix = NULL,
	.family_size = 1,
	.will_track_pedigree = FALSE,
	.will_allocate_ids = TRUE,
	.filename_prefix = NULL,
	.will_save_pedigree_to_file = FALSE,
	.will_save_bvs_to_file = FALSE,
	.will_save_alleles_to_file = FALSE,
	.will_save_to_simdata = TRUE
};

/* Creator for an empty SimData object */
/*const SimData EMPTY_SIMDATA = {
	.n_markers = 0,
	.markers = NULL,
	.map = {.n_chr = 0, .chr_ends = NULL, .positions = NULL},
	.m = NULL,
	.e = {.effects = {.matrix = NULL, .rows = 0, .cols = 0}, .effect_names = NULL},
	.current_id = 0
};*/

int RAND_HALFPOINT = (RAND_MAX/2);

/** Replace calls to malloc direct with this function, which errors and exits
 * with status 2 if memory allocation fails.
 *
 * @param size size of space to be allocated. Treat like the parameter of a
 * regular malloc call.
 * @returns pointer to the allocated space.
 */
void* get_malloc(size_t size) {
    if (size == 0) {
        warning( "0 memory allocation requested. The maintainer isn't sure why that would happen.\n");
        return NULL;
    }
	void* v = malloc(size);
	if (v == NULL) {
		error( "Memory allocation failed. Exiting.\n");
	}
	return v;
}


/** Creator for an empty AlleleMatrix object of a given size. Includes memory
 * allocation for `n_genotypes` worth of `.alleles`.
 *
 * @param n_markers number of rows/markers to create
 * @param n_genotypes number of individuals to create. This includes filling the first
 * n_genotypes entries of .alleles with heap char* of length n_markers, so that the
 * alleles for these can be added without further memory allocation.
 * @returns pointer to the empty created AlleleMatrix
 */
AlleleMatrix* create_empty_allelematrix(int n_markers, int n_genotypes) {
	AlleleMatrix* m = get_malloc(sizeof(AlleleMatrix));

	m->n_genotypes = n_genotypes;
	m->n_markers = n_markers;
	//m->alleles = get_malloc(sizeof(char*) * CONTIG_WIDTH);
	for (int i = 0; i < n_genotypes; ++i) {
		m->alleles[i] = get_malloc(sizeof(char) * (n_markers<<1));
		memset(m->alleles[i], 0, sizeof(char) * (n_markers<<1));
		//m->ids[i] = 0;
	}
	memset(m->alleles + n_genotypes, 0, sizeof(char*) * (CONTIG_WIDTH - n_genotypes)); // setting the pointers to NULL

	memset(m->ids, 0, sizeof(unsigned int) * CONTIG_WIDTH);
	memset(m->pedigrees[0], 0, sizeof(unsigned int) * CONTIG_WIDTH);
	memset(m->pedigrees[1], 0, sizeof(unsigned int) * CONTIG_WIDTH);
	memset(m->groups, 0, sizeof(unsigned int) * CONTIG_WIDTH);
	memset(m->names, 0, sizeof(char*) * CONTIG_WIDTH); // setting the pointers to NULL

	m->next = NULL;

	return m;
}

/** Creator for an empty SimData object on the heap. This is the main struct
 * that will contain/manage simulation data.
 *
 * @returns pointer to the empty created SimData
 */
SimData* create_empty_simdata() {
	SimData* d = get_malloc(sizeof(SimData));
	d->n_markers = 0;
	d->markers = NULL;
	d->map.n_chr = 0;
	d->map.chr_ends = NULL;
	d->map.chr_lengths = NULL;
	d->map.positions = NULL;
	d->m = NULL;
	d->e.effects.rows = 0;
	d->e.effects.cols = 0;
	d->e.effects.matrix = NULL;
	d->e.effect_names = NULL;
	d->current_id = 0;
	return d;
}



/*------------------------Supporter Functions--------------------------------*/

/** Allocate lifetime-unique ids to each genotype in the range of whole
 * SimData indexes `from_index` to `to_index` inclusive. Not intended to
 * be called by an end user.
 *
 * @param d the SimData struct on which to perform actions
 * @param from_index the starting 0-based index of the range to allocate ids
 * @param to_index the last 0-based index in the range to allocate ids
 */
void set_ids(SimData* d, int from_index, int to_index) {
	if (to_index < from_index) {
		warning( "Bad range for setting ids\n");
		return;
	}

	AlleleMatrix* m = d->m;
	int i, total_i = 0;
	// find the AlleleMatrix from_index belongs to
	while (total_i + m->n_genotypes <= from_index) {
		if (m->next == NULL) {
			warning( "That index does not exist.");
			return;
		} else {
			total_i += m->n_genotypes;
			m = m->next;
		}
	}
	total_i = from_index;

	// check for overflow:
	if (d->current_id > UINT_MAX - to_index + from_index) {
		warning( "IDs will overflow in the process of this calculation. Past this point id-based functions (anything to do with pedigrees) have undefined behaviour.\n");
	}

	// allocate the ids
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i, ++total_i) {
			++ d->current_id;
			m->ids[i] = d->current_id;
			if (total_i >= to_index) {
				return;
			}
		}

		if (m->next == NULL) {
			warning( "Could not set all ids.");
			return;
		} else {
			m = m->next;
		}
	}
}

/** Opens a table file and reads the number of columns and rows
 * (including headers) separated by `sep` into a TableSize struct that is
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
 * @returns TableSize struct with .num_columns and .num_rows filled. These
 * counts include header rows/columns and exclude blank rows.
 */
struct TableSize get_file_dimensions(const char* filename, char sep) {
	struct TableSize details;
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
	int has_length = FALSE;
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
			has_length = FALSE;

		} else if (c == sep) {
			sep_count += 1;
		} else if (has_length == FALSE) {
			has_length = TRUE;
		}
		c = fgetc(fp);
	}
	if (has_length) {
		details.num_rows += 1; // for the last row before EOF
	}

	fclose(fp);
	return details;
}

/** Returns the located index in an array of integers where the integer
 * is `target`. Returns -1 if no match was found.
 * @see get_from_unordered_str_list()
 *
 * The list is assumed to be sorted in ascending order.
 *
 * It uses the bisection method to search the list.
 *
 * @param target the integer to be located
 * @param list the array of integers to search
 * @param list_len length of the array of integers to search
 * @returns Index in `list` where we find the same integer as
 * `target`, or -1 if no match is found.
 */
int get_from_ordered_uint_list(unsigned int target, unsigned int* list, unsigned int list_len) {
	int first = 0, last = list_len - 1;
	int index = (first + last) / 2;
	while (list[index] != target && first <= last) {
		if (list[index] < target) {
			first = index + 1;
		} else {
			last = index - 1;
		}
		index = (first + last) / 2;
	}

	if (first > last) {
		warning( "Search failed: are you sure your list is sorted?\n");
	}
	return index;
}

/** Returns the first located index in an array of strings where the string
 * is the same as the string `target`. Returns -1 if no match was found.
 * @see get_from_ordered_uint_list()
 *
 * The list of strings is not assumed to be sorted.
 *
 * @param target the string to be located
 * @param list the array of strings to search
 * @param list_len length of the array of strings to search
 * @returns Index in `list` where we find the same string as
 * `target`, or -1 if no match is found.
 */
int get_from_unordered_str_list(char* target, char** list, int list_len) {
	for (int i = 0; i < list_len; ++i) {
		if (strcmp(list[i], target) == 0) {
			return i;
		}
	}
	return -1; // did not find a match.
}

/**Produce a random ordering of the first n elements in an array of integers
 * using a (partial) Fisher-Yates shuffle.
 *
 * Modified from https://benpfaff.org/writings/clc/shuffle.html
 *
 * @param sequence the array
 * @param total_n length of the array
 * @param n_to_shuffle After calling this function, the first n_to_shuffle
 * integers in the array will be randomly ordered by a Fischer-Yates shuffle
 */
void shuffle_up_to(int* sequence, size_t total_n, size_t n_to_shuffle) {
	if (n_to_shuffle > 1) {
        size_t maxi = total_n > n_to_shuffle ? n_to_shuffle - 1 : total_n - 1;
		size_t i;
        for (i = 0; i < maxi; ++i) {
			// items before i are already shuffled
			size_t j = i + rand() / (RAND_MAX / (total_n - i) + 1);

			// add the next chosen value to the end of the shuffle
			int t = sequence[j];
			sequence[j] = sequence[i];
			sequence[i] = t;
		}
	}
}

/** Fills the designated section of the `.names` array in an
 * AlleleMatrix with the pattern `prefix`index. This function is not intended
 * to be called by an end user.
 *
 * In future this function could be expanded to allow for different naming formats.
 *
 * The index is padded with zeros depending on the size of
 * `a->n_genotypes`.
 *
 * @param a pointer to the AlleleMatrix whose `.names` to modify
 * @param prefix the prefix to add to the suffix to make the new genotype name
 * @param suffix suffixes start at this value and increment for each additional name
 * @param from_index the new names are added to this index and all those following it in this
 * AlleleMatrix.
*/
void set_names(AlleleMatrix* a, char* prefix, int suffix, int from_index) {
	char sname[NAME_LENGTH];
	char format[NAME_LENGTH];
	if (prefix == NULL) {
		// make it an empty string instead, so it is not displayed as (null)
		prefix = "";
	}
	// use sname to save the number of digits to pad by:
	sprintf(sname, "%%0%dd", get_integer_digits(a->n_genotypes - from_index));  // Creates: %0[n]d
	sprintf(format, "%s%s", prefix, sname);

	++suffix;
	for (int i = from_index; i < a->n_genotypes; ++i) {
		// clear name if it's pre-existing
		if (a->names[i] != NULL) {
			free(a->names[i]);
		}

		// save new name
		sprintf(sname, format, suffix);
		a->names[i] = get_malloc(sizeof(char) * (strlen(sname) + 1));
		strcpy(a->names[i], sname);

		++suffix;
	}
}

/** Count and return the number of digits in `i`.
 *
 * @param i the integer whose digits are to be counted.
 * @returns the number of digits to print `i`
 */
int get_integer_digits(int i) {
	int digits = 0;
	while (i != 0) {
		i = i / 10;
		digits ++;
	}
	return digits;
}


/** Comparator function for qsort. Used to compare an array of MarkerPosition**
 * @see sort_markers()
 *
 * Use this on an array of pointers to values in SimData.map.positions.
 *
 * Sorts lower chromosome numbers before higher chromosome numbers. Within
 * chromosomes, sorts lower positions before higher ones.
 *
 * If chromosome number is 0, the MarkerPosition is considered uninitialised.
 * All unitialised values are moved to the end. The ordering of uninitialised
 * values amongst themselves is undefined.
 *
 * Reference: https://stackoverflow.com/questions/32948281/c-sort-two-arrays-the-same-way
 */
int _simdata_pos_compare(const void *pp0, const void *pp1) {
	MarkerPosition marker0_pos = **(MarkerPosition**)pp0;
	MarkerPosition marker1_pos = **(MarkerPosition**)pp1;
	// first, move any empties towards the end.
	if (marker0_pos.chromosome == 0) {
		if (marker1_pos.chromosome == 0) {
			return 0;
		} else {
			return 1; // move towards end
		}
	} else if (marker1_pos.chromosome == 0) {
		return -1; // move towards end
	} else if (marker0_pos.chromosome == marker1_pos.chromosome) { // on same chromosome
		if (marker0_pos.position < marker1_pos.position) {
			return -1;
		} else if (marker0_pos.position > marker1_pos.position) {
			return 1;
		} else {
			return 0; // are completely equal. Might both be empty
		}
	} else if (marker0_pos.chromosome < marker1_pos.chromosome) {
		return -1;
	} else {//if (marker0_pos.chromosome > marker1_pos.chromosome) {
		return 1;
	}
}

/** Comparator function for qsort. Used to compare an array of doubles* to sort
 * them in descending order of the doubles they point to.
 * @see split_by_bv()
 *
 * Sorts higher numbers before lower numbers. If they are equal, their
 * order after comparison is undefined.
 */
int _descending_double_comparer(const void* pp0, const void* pp1) {
	double d0 = **(double **)pp0;
	double d1 = **(double **)pp1;
	if (d0 > d1) {
		return -1;
	} else {
		return (d0 < d1); // 0 if equal, 1 if d0 is smaller
	}
}

/** Comparator function for qsort. Used to compare an array of doubles* to sort
 * them in ascending order of the doubles they point to.
 * @see split_by_bv()
 *
 * Sorts lower numbers before higher numbers. If they are equal, their
 * order after comparison is undefined.
 */
int _ascending_double_comparer(const void* pp0, const void* pp1) {
	double d0 = **(double **)pp0;
	double d1 = **(double **)pp1;
	if (d0 < d1) {
		return -1;
	} else {
		return (d0 > d1); // 0 if equal, 1 if d0 is smaller
	}
}

/** Comparator function for qsort. Used to compare an array of floats to sort
 * them in ascending order.
 * @see generate_gamete()
 *
 * Sorts lower numbers before higher numbers. If floats are equal, their
 * order after comparison is undefined.
 */
int _ascending_float_comparer(const void* p0, const void* p1) {
	float f0 = *(float *)p0;
	float f1 = *(float *)p1;
	if (f0 < f1) {
		return -1;
	} else {
		return (f0 > f1); // 0 if equal, 1 if f0 is greater
	}
}

/** Comparator function for qsort. Used to compare an array of integers to sort
 * them in ascending order.
 * @see get_existing_groups()
 *
 * Sorts lower numbers before higher numbers. If floats are equal, their
 * order after comparison is undefined.
 */
int _ascending_int_comparer(const void* p0, const void* p1) {
	int f0 = *(int *)p0;
	int f1 = *(int *)p1;
	if (f0 < f1) {
		return -1;
	} else {
		return (f0 > f1); // 0 if equal, 1 if f0 is greater
	}
}

/** Comparator function for qsort. Used to compare an array of integers* to sort
 * the integers pointed to in ascending order.
 * @see get_existing_group_counts()
 *
 * Sorts lower numbers before higher numbers. If floats are equal, their
 * order after comparison is undefined.
 */
int _ascending_int_dcomparer(const void* pp0, const void* pp1) {
	int f0 = **(int **)pp0;
	int f1 = **(int **)pp1;
	if (f0 < f1) {
		return -1;
	} else {
		return (f0 > f1); // 0 if equal, 1 if f0 is greater
	}
}

/** A function to tidy the internal storage of genotypes after addition
 * or deletion of genotypes in the SimData. Not intended to be called by an
 * end user - functions which require it should be calling it already.
 *
 * Ideally, we want all AlleleMatrix structs in the SimData's linked list
 * to have no gaps. That is, if there are more than CONTIG_WIDTH genotypes, the
 * first AM should be full/contain CONTIG_WIDTH genotypes, and so forth, and the
 * AM at the end should have its n genotypes having indexes < n (so all at
 * the start of the AM with no gaps in-between).
 *
 * This function achieves the cleanup by using two pointers: a filler out
 * the front that identifies a genotype that needs to be shifted back/that
 * occurs after a gap, and a filler that identifies each gap and copies
 * the genotype at the filler back into it.
 *
 * @param d The SimData struct on which to operate.
 */
void condense_allele_matrix( SimData* d) {
	int checker, filler;
	AlleleMatrix* checker_m = d->m, *filler_m = d->m;

	// find the first empty space with filler
	while (1) {
		if (filler_m->n_genotypes < CONTIG_WIDTH) {
			for (filler = 0; filler < CONTIG_WIDTH; ++filler) {
				// an individual is considered to not exist if it has no genome.
				if (filler_m->alleles[filler] == NULL) {
					break; // escape for loop
				}
			}
			// assume we've found one, since n_genotypes < CONTIG_WIDTH.
			break; // escape while loop
		}

		if (filler_m->next == NULL) {
			// No empty spaces
			return;
		} else {
			filler_m = filler_m->next;
		}
	}

	// loop through all genotypes with checker, shifting them back when we find them.
	while (1) {
		for (checker = 0; checker < CONTIG_WIDTH; ++checker) {
			if (checker_m->alleles[checker] == NULL) {
				// check our filler has a substitute
				while (filler_m->alleles[filler] == NULL) {
					++filler;
					if (filler >= CONTIG_WIDTH) {
						// move to the next AM
						if (filler_m->next != NULL) {
							filler_m = filler_m->next;
							filler = 0;
						} else {
							return;
						}
					}
				}

				// put in the substitute
				checker_m->alleles[checker] = filler_m->alleles[filler];
				filler_m->alleles[filler] = NULL;
				checker_m->names[checker] = filler_m->names[filler];
				filler_m->names[filler] = NULL;
				checker_m->ids[checker] = filler_m->ids[filler];
				filler_m->ids[filler] = -1;
				checker_m->pedigrees[0][checker] = filler_m->pedigrees[0][filler];
				checker_m->pedigrees[1][checker] = filler_m->pedigrees[1][filler];
				filler_m->pedigrees[0][filler] = 0;
				filler_m->pedigrees[1][filler] = 0;
				checker_m->groups[checker] = filler_m->groups[filler];
				filler_m->groups[filler] = 0;
				if (checker_m != filler_m) {
					++ checker_m->n_genotypes;
					-- filler_m->n_genotypes;

					if (filler_m->n_genotypes == 0) {
						// find the previous matrix in the linked list
						AlleleMatrix* previous = checker_m;
						while (previous->next != filler_m) {
							previous = previous->next;
							if (previous == NULL) {
								error( "Error: filler got somewhere checker can't go.\n");
							}
						}
						// delete this one and cut it out of the list
						previous->next = filler_m->next;
						filler_m->next = NULL; // so we don't delete the rest of the list
						delete_allele_matrix(filler_m);
						// move our filler values on
						if (previous->next == NULL) {
							return;
						} else {
							filler_m = previous->next;
							filler = 0;
						}

					}
				}
			}
		}

		// We're done with the AM, move on to the next one.
		if (checker_m->next == NULL) {
			error( "During allele matrix condensing, checker got ahead of filler somehow\n");
		} else {
			checker_m = checker_m->next;
		}
	}
}



/*----------------------------------Locators---------------------------------*/

/** Returns the name of the genotype with a given id.
 *
 * The function uses a bisection search on the AlleleMatrix where it should be
 * found. Searching by id is therefore fast.
 * @see get_from_ordered_uint_list()
 *
 * This function assumes that ids are never reshuffled in the SimData. This is true
 * as long as condense_allele_matrix's process of moving genotypes while retaining
 * their order is the only genotype-rearranging function in use.
 * This function's algorithm will need to
 * be revisited if different genotype-rearranging processes are implemented.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose name is sought
 * @returns the name of the genotype that has id `id`, as a copy of the pointer
 * to the heap memory where the name is saved (so *don't* free the pointer returned
 * from this function). Returns NULL if the ID does not exist. Note the return value
 * might also be NULL if the genotype of this ID has no name.
 */
char* get_name_of_id( AlleleMatrix* start, unsigned int id) {
	if (id <= 0) {
		warning( "Invalid negative ID %d\n", id);
		return NULL;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;

	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_genotypes != 0 && id >= m->ids[0] && id <= m->ids[m->n_genotypes - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_genotypes);

			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
					return NULL;
				} else {
					m = m->next;
				}
			}

			return m->names[index];

		}

		if (m->next == NULL) {
			warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
			return NULL;
		} else {
			m = m->next;
		}
	}
}

/** Returns the alleles at each marker of the genotype with a given id.
 *
 * The function uses a bisection search on the AlleleMatrix where it should be
 * found. Searching by id is therefore fast.
 * @see get_from_ordered_uint_list()
 *
 * This function assumes that ids are never reshuffled in the SimData. This is true
 * as long as condense_allele_matrix's process of moving genotypes while retaining
 * their order is the only genotype-rearranging function in use.
 * This function's algorithm will need to
 * be revisited if different genotype-rearranging processes are implemented.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose alleles are sought
 * @returns the alleles of the genotype that has id `id`, as a copy of the pointer
 * to the heap memory where the genotype is saved (so *don't* free the pointer returned
 * from this function). It points to a sequence of characters, ordered according to
 * the markers in the SimData to which the AlleleMatrix belongs. Returns NULL if the
 * ID does not exist.
 */
char* get_genes_of_id ( AlleleMatrix* start, unsigned int id) {
	if (id <= 0) {
		warning( "Invalid negative ID %d\n", id);
		return NULL;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;

	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_genotypes != 0 && id >= m->ids[0] && id <= m->ids[m->n_genotypes - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_genotypes);

			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
					return NULL;
				} else {
					m = m->next;
					continue;
				}
			}

			return m->alleles[index];

		}

		if (m->next == NULL) {
			warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
			return NULL;
		} else {
			m = m->next;
		}
	}
}

/** Saves the ids of the parents of a genotype with a particular id to
 * the output array `output`.
 *
 * The function uses a bisection search on the AlleleMatrix where it should be
 * found. Searching by id is therefore fast.
 * @see get_from_ordered_uint_list()
 *
 * This function assumes that ids are never reshuffled in the SimData. This is true
 * as long as condense_allele_matrix's process of moving genotypes while retaining
 * their order is the only genotype-rearranging function in use.
 * This function's algorithm will need to
 * be revisited if different genotype-rearranging processes are implemented.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose parents are sought
 * @param output An array which the calling function can access where this function
 * will put its results.
 * @returns 0 when the id is successfully identified and at least one parent's
 * id is known, 1 if neither parent is known, and 2 if the ID passed in does
 * not exist. The ids of both parents if at least one parent is
 * known/nonzero are saved to the array `output`.
 */
int get_parents_of_id( AlleleMatrix* start, unsigned int id, unsigned int output[2]) {
	if (id <= 0) {
		return 1;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;
	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_genotypes != 0 && id >= m->ids[0] && id <= m->ids[m->n_genotypes - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_genotypes);

			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
					return 2;
				} else {
					m = m->next;
				}
			}

			if (m->pedigrees[0][index] > 0 || m->pedigrees[1][index] > 0) {
				output[0] = m->pedigrees[0][index];
				output[1] = m->pedigrees[1][index];
				return 0;
			}
			return 1; // if neither parent's id is known

		}

		if (m->next == NULL) {
			warning( "Could not find the ID %d: did you prematurely delete this genotype?\n", id);
			return 2;
		} else {
			m = m->next;
		}
	}
}

/** Search for genotypes with certain names in a linked list of AlleleMatrix and
 * save the ids of those names. Exits if any name cannot be found.
 *
 * This function must check every name in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param n_names the length of the array of names which are being sought.
 * @param names an array of names whose ids are being sought.
 * @param output pointer to an array of length at least `n_names` which can
 * be accessed by the calling function. The ids of each name are saved to corresponding
 * indexes in the array this pointer points to.
 */
void get_ids_of_names( AlleleMatrix* start, int n_names, char* names[n_names], unsigned int* output) {
	int found;
	//int ids = malloc(sizeof(int) * n_names);
	AlleleMatrix* m;
	int i, j;

	for (i = 0; i < n_names; ++i) {
		found = FALSE;
		output[i] = -1;
		m = start;
		while (1) {
			// try to identify the name in this AM
			for (j = 0; j < m->n_genotypes; ++j) {
				if (strcmp(m->names[j], names[i]) == 0) {
					found = TRUE;
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
 * list of AlleleMatrix, and return its id. Exits if such a child cannot be found.
 *
 * This function must check every genotype in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the child is assumed to be found.
 * @param parent1id one of the parents of the genotype must have this id.
 * @param parent2id the other parent of the genotype must have this id.
 * @returns the id of the first sequentially located genotype whose parents match the
 * two parent ids provided
 */
unsigned int get_id_of_child( AlleleMatrix* start, unsigned int parent1id, unsigned int parent2id) {
	AlleleMatrix* m = start;
	int j;

	while (1) {
		// try to identify the child in this AM
		for (j = 0; j < m->n_genotypes; ++j) {
			if ((parent1id == m->pedigrees[0][j] && parent2id == m->pedigrees[1][j]) ||
					(parent1id == m->pedigrees[1][j] && parent2id == m->pedigrees[0][j])) {
				return m->ids[j];
			}
		}

		if ((m = m->next) == NULL) {
			warning( "Didn't find the child of %d & %d\n", parent1id, parent2id);
			return 0;
		}
	}
}

/** Search for a genotype with parentage matching two given parent ids in a linked
 * list of AlleleMatrix, and return its index. Exits if such a child cannot be found.
 *
 * This function must check every genotype in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the child is assumed to be found.
 * @param parent1id one of the parents of the genotype must have this id.
 * @param parent2id the other parent of the genotype must have this id.
 * @returns the index (0-based, starting at the start of `start`) of the first sequentially
 * located genotype whose parents match the two parent ids provided
 */
int get_index_of_child( AlleleMatrix* start, unsigned int parent1id, unsigned int parent2id) {
	AlleleMatrix* m = start;
	int j, total_j = 0;

	while (1) {
		// try to identify the child in this AM
		for (j = 0; j < m->n_genotypes; ++j, ++total_j) {
			if ((parent1id == m->pedigrees[0][j] && parent2id == m->pedigrees[1][j]) ||
					(parent1id == m->pedigrees[1][j] && parent2id == m->pedigrees[0][j])) {
				return total_j;
			}
		}

		if ((m = m->next) == NULL) {
			error( "Didn't find the child of %d & %d\n", parent1id, parent2id);
		}
	}
}

/** Search for a genotype with a particular name in a linked
 * list of AlleleMatrix, and return its index. Exits if such a child cannot be found.
 *
 * This function must check every genotype in the linked list for matches, so will be
 * relatively slow.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype is assumed to be found.
 * @param name a string to match to the name of the target
 * @returns the index (0-based, starting at the start of `start`) of the first sequentially
 * located genotype whose name is the same as the provided name.
 */
int get_index_of_name( AlleleMatrix* start, char* name) {
	AlleleMatrix* m = start;
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
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * to locate the id at a particular index.
 * @param index the index of the target. It is assumed to start at 0 at the start
 * of the matrix in `start` and be incremented up until the last entry in the matrix
 * of the last entry in the linked list.
 * @returns the lifetime-unique id of the genotype found at that index.
 */
unsigned int get_id_of_index( AlleleMatrix* start, int index) {
	AlleleMatrix* m = start;
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
			return 0;
		}
	}
}

/** Get the alleles of a genotype by its index. The index is assumed to be 0-based,
 * starting at the first entry of `start` and continuing through the linked list
 * to which `start` belongs. Exits if the linked list does not contain that index.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * to locate the id at a particular index.
 * @param index the index of the target. It is assumed to start at 0 at the start
 * of the matrix in `start` and be incremented up until the last entry in the matrix
 * of the last entry in the linked list.
 * @returns the alleles of the genotype that has idnex `index`, as a copy of the pointer
 * to the heap memory where the genotype is saved (so *don't* free the pointer returned
 * from this function). It points to a sequence of characters, ordered according to
 * the markers in the SimData to which the AlleleMatrix belongs.
 */
char* get_genes_of_index( AlleleMatrix* start, int index) {
	if (index < 0) {
		warning( "Invalid negative index %d\n", index);
		return NULL;
	}
	if (start == NULL) {
		error( "Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;
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
			// for each genotype, check all group numbers.
			for (i = 0; i < m->n_genotypes; ++i) {
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
			// for each genotype, check all group numbers.
			for (i = 0; i < m->n_genotypes; ++i) {

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
		while (indexes_to_split[i] >= total_i + m->n_genotypes) {
			if (m->next == NULL) {
				warning( "Only found %d out of %d indexes\n", i, n);
				return new_group;
			}
			total_i += m->n_genotypes;
			m = m->next;
		}

		m->groups[indexes_to_split[i] - total_i] = new_group;
	}
	return new_group;
}

/** Give every individual in the group a new group number that does not
 * belong to any other existing group (thereby allocating each genotype
 * in the group to a new group of 1).
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store the
 * identifiers of all groups created. The array's length should be the number
 * of members of group_id.
 */
void split_into_individuals( SimData* d, int group_id, int* results) {
	// get pre-existing numbers
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	// have another variable for the next id we can't allocate so we can still free the original.
	int level = 0;

	// looping through all entries
	AlleleMatrix* m = d->m;
	int n_found = 0; //number of group members reallocated.
	int next_id = 0;
	while (1) {
		// check all lines to see if this one belongs to the group.
		for (int i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// change it to a new unique group
				// first, find the next unused group;
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
				if (results != NULL) {
					results[n_found] = next_id;
					++n_found;
				}
				//++next_id;
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
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param results NULL if the caller does not care to know the identifiers of the
 * family groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store the
 * identifiers of all groups created. For safety, unless you know how many
 * groups will be created, the array's length should be the number
 * of members of group_id.
 */
void split_into_families(SimData* d, int group_id, int* results) {
	// get pre-existing numbers
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	// have another variable for the next id we can't allocate so we can still free the original.
	int level = 0;

	// maximum number of families is the number of members in the group
	// but this could be greater than our maximum allocation, so we
	// will have to go through in batches
	int families_found = 0;
	unsigned int family_groups[CONTIG_WIDTH];
	unsigned int family_identities[2][CONTIG_WIDTH];


	// looping through all entries
	AlleleMatrix* m = d->m;
	AlleleMatrix* bookmarkm = NULL;
	int n_found = 0; // number of groups created, for saving results
	int next_id = 0;
	while (1) {
		// check all genotypes to see if this one belongs to the group.
		for (int i = 0; i < m->n_genotypes; ++i) {
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
					if (families_found < CONTIG_WIDTH) {
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

						if (results != NULL) {
							results[n_found] = next_id;
							++n_found;
						}

					} else if (bookmarkm == NULL) {
						// if we are here, then there are more families than we
						// can track at a time. Save the place and come back for the
						// extra families later.
						bookmarkm = m;
					} // else we have a bookmark already and so will be coming back here
				}
			}
		}

		if (m->next == NULL) {
			free(existing_groups);
			if (bookmarkm == NULL) {
				// we have processed all members of the original group
				return;
			} else {
				// we need to go back and process remaining members of the group
				m = bookmarkm;
				families_found = 0;
				existing_groups = get_existing_groups( d, &n_groups);
				level = 0;
				next_id = 0;
			}
		} else {
			m = m->next;
		}
	}
}

/** Split a group into a set of smaller groups, each containing the
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
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param parent 1 to group together genotypes that share the same first parent,
 * 2 group those with the same second parent. Raises an error
 * if this parameter is not either of those values.
 * @param results NULL if the caller does not care to know the identifiers of the
 * family groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store the
 * identifiers of all groups created. For safety, unless you know how many
 * groups will be created, the array's length should be the number
 * of members of group_id.
 */
void split_into_halfsib_families( SimData* d, int group_id, int parent, int* results) {
	if (!(parent == 1 || parent == 2)) {
		warning( "Value error: `parent` must be 1 or 2.");
		return;
	}
	--parent;

	// get pre-existing numbers
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	// have another variable for the next id we can't allocate so we can still free the original.
	int level = 0;

	int families_found = 0;
	unsigned int family_groups[CONTIG_WIDTH];
	unsigned int family_identities[CONTIG_WIDTH];

	// looping through all entries
	AlleleMatrix* m = d->m;
	AlleleMatrix* bookmarkm = NULL;
	int n_found = 0; // number of families, for saving to results
	int next_id = 0;
	while (1) {
		// check all genotypes to see if this one belongs to the group.
		for (int i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// First, see if it is a member of a family we already know.
				for (int j = 0; j < families_found; ++j) {
					if (m->pedigrees[parent][i] == family_identities[j]) {
						m->groups[i] = family_groups[j];
						break;
					}
				}

				// if the group number has not been updated in the above loop
				// (so we don't know this family yet)
				if (m->groups[i] == group_id) {
					if (families_found < CONTIG_WIDTH) {
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
						family_identities[families_found] = m->pedigrees[parent][i];
						++families_found;

						if (results != NULL) {
							results[n_found] = next_id;
							++n_found;
						}

					} else if (bookmarkm == NULL) {
						// if we are here, then there are more families than we
						// can track at a time. Save the place and come back for the
						// extra families later.
						bookmarkm = m;
					} // else we have a bookmark already and so will be coming back here
				}
			}
		}

		if (m->next == NULL) {
			free(existing_groups);
			if (bookmarkm == NULL) {
				// we have processed all members of the original group
				return;
			} else {
				// we need to go back and process remaining members of the group
				m = bookmarkm;
				families_found = 0;
				existing_groups = get_existing_groups( d, &n_groups);
				level = 0;
				next_id = 0;
			}
		} else {
			m = m->next;
		}
	}
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
 * A more general approach to this task: split_evenly_into_n()
 *
 * An alternate approach to splitting a group in two: split_randomly_into_two()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @returns the group number of the new group to which half the members
 * of the old group have been allocated.
 */
int split_evenly_into_two(SimData* d, int group_id) {
	// get the shuffle to be our even allocations
	int size = get_group_size(d, group_id);
	int even_half = size / 2;
	int allocations[size];
	for (int i = 0; i < size; ++i) {
		allocations[i] = i;
	}
	shuffle_up_to(allocations, size, even_half);

	int new_group = get_new_group_num(d);
	AlleleMatrix* m = d->m;
	int groupi = 0;
	int i;

	// loop through group members
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// See if it should be allocated to the new group
				for (int j = 0; j < even_half; ++j) {
					if (allocations[j] == groupi) {
						m->groups[i] = new_group;
						break;
					}
				}
				++groupi;
			}
		}

		if (m->next == NULL) {
			return new_group;
		} else {
			m = m->next;
		}
	}
}

/** Split a group into n groups of equal size (or size differing only
 * by one, if n does not perfectly divide the group size.), using a
 * random permutation of the group members to determine which goes where.
 *
 * Of the split groups produced, the first has the same group number as the original
 * group (parameter group_id).
 *
 * A more general approach to this task: split_by_specific_counts_into_n()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 */
void split_evenly_into_n(SimData* d, int group_id, int n, int* results) {
    if (n <= 1) {
		warning( "Cannot distribute between %d groups\n", n);
		return;
	}

	// get the shuffle to be our even allocations
	int size = get_group_size(d, group_id);

	int each_size = size / n;
	int extra = size % n;
	int boxes[n];
	for (int i = 0; i < n; ++i) {
		boxes[i] = each_size;
		if (i < extra) {
			boxes[i] ++;
		}
	}

	// no code simplification just from having equal sizes is possible so:
	split_by_specific_counts_into_n(d, group_id, n, boxes, results);

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
 * @param d the SimData struct on which to perform the operation
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
 */
void split_by_specific_counts_into_n(SimData* d, int group_id, int n, int* counts, int* results) {
    if (n <= 1) {
		warning( "Cannot distribute between %d groups\n", n);
		return;
	}

	int size = get_group_size(d, group_id);

	int cumulative_counts[n-1];
	int sum = 0;
	for (int j = 0; j < n - 1; ++j) {
		sum += counts[j];
		cumulative_counts[j] = sum;
		if (cumulative_counts[j] >= size) {
            warning( "Provided capacities are larger than actual group: some buckets will not be filled\n");
			//don't bother to calculate more
			break;
		}
	}

	int allocations[size];
	for (int i = 0; i < size; ++i) {
		allocations[i] = i;
	}
    shuffle_up_to(allocations, size, cumulative_counts[n-2]);

	int new_group[n-1];
	get_n_new_group_nums(d,n-1, new_group);

	if (results != NULL) {
		results[0] = group_id;
		memcpy(results + 1, new_group, sizeof(int) * (n-1));
	}

	AlleleMatrix* m = d->m;
	int groupi = 0;
	int i;

	// loop through group members
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// See where it was shuffled
				for (int j = 0; j < size; ++j) {
					if (allocations[j] == groupi) {
						// allocate it to its correct group
						for (int k = n-2; k >=0; --k) {
                            if (j >= cumulative_counts[k]) {
                                m->groups[i] = new_group[k];
                                break;
                            }
						}
                        // else is already in the first group, where it belongs
						break;
					}
				}
				++groupi;
			}
		}

		if (m->next == NULL) {
			return;
		} else {
			m = m->next;
		}
	}

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
 * A more general approach to this task: split_randomly_into_n()
 *
 * An alternate approach to splitting a group in two: split_evenly_into_two()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @returns the group number of the new group to which some members
 * of the old group may have been randomly allocated.
 */
int split_randomly_into_two(SimData* d, int group_id) {
	int new_group = get_new_group_num(d);

	AlleleMatrix* m = d->m;
	int i;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id && (unif_rand() > 0.5)) {
				m->groups[i] = new_group;
			}
		}

		if (m->next == NULL) {
			return new_group;
		} else {
			m = m->next;
		}
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
 * To split by uneven probabilities instead: split_by_probabilities_into_n()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group to be split
 * @param n the number of groups among which to randomly distribute group
 * members.
 * @param results NULL if the caller does not care to know the identifiers of the
 * groups created, or a pointer to an array to which these identifiers
 * should be saved. It is assumed that the array is long enough to store
 * n identifiers.
 */
void split_randomly_into_n(SimData* d, int group_id, int n, int* results) {
    if (n <= 1) {
		warning( "Cannot distribute between %d groups\n", n);
		return;
	}

	// get the n group numbers
	int new_groups[n-1];
	get_n_new_group_nums(d, n-1, new_groups);

	// save the results
	if (results != NULL) {
		results[0] = group_id;
		memcpy(results + 1, new_groups, sizeof(int) * (n-1));
	}

	AlleleMatrix* m = d->m;
	int i, randgroup;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			randgroup = round(unif_rand() * (n-1));
			if (m->groups[i] == group_id && randgroup) {
				m->groups[i] = new_groups[randgroup - 1];
			}
		}

		if (m->next == NULL) {
			return;
		} else {
			m = m->next;
		}
	}
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
 * @param d the SimData struct on which to perform the operation
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
 */
void split_by_probabilities_into_n(SimData* d, int group_id, int n, double* probs, int* results) {
    if (n <= 1) {
		warning( "Cannot distribute between %d groups\n", n);
		return;
	}

	// Check the probabilities
	double cumulative_probs[n-1];
	double sum = 0;
	for (int j = 0; j < n - 1; ++j) {
		sum += probs[j];
		cumulative_probs[j] = sum;
		if (cumulative_probs[j] >= 1) {
            warning( "Provided probabilities add up to 1 or more: some buckets will not be filled\n");
			//don't bother to calculate more
			break;
		}
	}

	// get the n group numbers
	int new_groups[n-1];
	get_n_new_group_nums(d, n-1, new_groups);

	// save the results
	if (results != NULL) {
		results[0] = group_id;
		memcpy(results + 1, new_groups, sizeof(int) * (n-1));
	}

	AlleleMatrix* m = d->m;
	int i;
	double randdraw;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// Allocate to a random group based on probabilities
				randdraw = unif_rand();
				for (int j = n-2; j >= 0; --j) {
                    if (randdraw > cumulative_probs[j]) {
                        m->groups[i] = new_groups[j];
                        break;
                    }
				}
			}
		}

		if (m->next == NULL) {
			return;
		} else {
			m = m->next;
		}
	}
}

/** Identify every group number that currently has members.
 *
 * @param d the SimData struct on which to perform the operation
 * @param n_groups a pointer to a location that can be accessed by the
 * calling function. The number of groups/length of the returned array
 * will be saved here.
 * @returns a heap array of length [value saved at n_groups after this
 * function has run] containing the group numbers that appear in the
 * group memberships in `d`, sorted in ascending order.
 */
int* get_existing_groups( SimData* d, int* n_groups) {
	int eg_size = 50;
	int* existing_groups = malloc(sizeof(int)* eg_size);
	*n_groups = 0;

	AlleleMatrix* m = d->m;
	int i, j;
	int new; // is the group number at this index a new one or not
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			new = 1;
			for (j = 0; j < *n_groups; ++j) {
				if (m->groups[i] == existing_groups[j]) {
					new = 0;
					break;
				}
			}

			if (new) {
				++ (*n_groups);
				existing_groups[*n_groups - 1] = m->groups[i];
			}
		}

		if (m->next == NULL) {
			qsort(existing_groups, *n_groups, sizeof(int), _ascending_int_comparer);
			return existing_groups;
		} else {
			m = m->next;
		}

		if (*n_groups >= eg_size - 1) {
			eg_size *= 1.5;
			int* temp = realloc(existing_groups, eg_size * sizeof(int));
			if (temp == NULL) {
				free(existing_groups);
				error( "Can't get all those groups.\n");
			} else {
				existing_groups = temp;
			}
		}
	}
}

/** Identify every group number that currently has members and
 * the number of genotypes currently allocated to that group.
 *
 * @param d the SimData struct on which to perform the operation
 * @param n_groups a pointer to a location that can be accessed by the
 * calling function. The number of groups/length of the returned array
 * will be saved here.
 * @returns a heap array of length [2][value saved at n_groups after this
 * function has run]. The array at the first index contains the group numbers,
 * and the array at the second index contains the number of members of that group,
 * ordered correspondingly. The arrays are sorted in ascending order by group number.
 * All three components of the returned value (2-long first-index array containing
 * the pointers to the toehr 2, and the two data-containing arrays) are allocated
 * from the heap so should be freed by the calling function when it finishes with
 * them.
 */
int** get_existing_group_counts( SimData* d, int* n_groups) {
	int eg_size = 50;
	int** returnval = get_malloc(sizeof(int*) * 2);
	int* existing_groups = get_malloc(sizeof(int)* eg_size); // group ids
	int* existing_group_counts = get_malloc(sizeof(int) * eg_size); // number in that group
	*n_groups = 0;

	AlleleMatrix* m = d->m;
	int i, j;
	int new; // is the group number at this index a new one or not
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			new = 1;
			for (j = 0; j < *n_groups; ++j) {
				if (m->groups[i] == existing_groups[j]) {
					existing_group_counts[j] += 1;
					new = 0;
					break;
				}
			}

			if (new) {
				++ (*n_groups);
				existing_groups[*n_groups - 1] = m->groups[i];
				existing_group_counts[*n_groups - 1] = 1;

			}
		}

		if (m->next == NULL) {
			if (*n_groups > 0) {
				int* sorting[*n_groups];
				for (int i = 0; i < *n_groups; i++) {
					sorting[i] = &(existing_groups[i]);
				}

				qsort(sorting, *n_groups, sizeof(int*), _ascending_int_dcomparer);
				int location_in_old;
				int* location_origin = existing_groups;
				int* sorted_groups = get_malloc(sizeof(int) * (*n_groups));
				int* sorted_counts = get_malloc(sizeof(int) * (*n_groups));
				for (int i = 0; i < *n_groups; ++i) {
					location_in_old = sorting[i] - location_origin;
					sorted_groups[i] = existing_groups[location_in_old];
					sorted_counts[i] = existing_group_counts[location_in_old];
				}

				free(existing_groups);
				free(existing_group_counts);
				returnval[0] = sorted_groups;
				returnval[1] = sorted_counts;
				return returnval;
			} else {
				return NULL;
			}
		} else {
			m = m->next;
		}

		if (*n_groups >= eg_size - 1) {
			eg_size *= 1.5;
			int* temp = realloc(existing_groups, eg_size * sizeof(int));
			if (temp == NULL) {
				free(existing_groups);
				error( "Can't get all those groups.\n");
			} else {
				existing_groups = temp;
			}
			temp = realloc(existing_group_counts, eg_size * sizeof(int));
			if (temp == NULL) {
				free(existing_group_counts);
				error( "Can't get all those groups.\n");
			} else {
				existing_group_counts = temp;
			}
		}
	}
}

/** Function to identify the next sequential integer that does not
 * identify a group that currently has member(s).
 *
 * This calls get_existing_groups() every time, so for better
 * speed, functions that do repeated group creation, like
 * split_into_individuals(), are
 * recommended to use get_n_new_group_nums() (if they know the
 * number of groups they need) or their own implementation, rather
 * than calling this function repeatedly.
 *
 * @param d the SimData struct on which to perform the operation
 * @return the next sequential currently-unused group number.
 */
int get_new_group_num( SimData* d) {
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	int i = 0;
	int gn = 1;

	while (i < n_groups) {
		if (gn < existing_groups[i]) {
			break;
		}

		++i;
		++gn;
	}
	free(existing_groups);
	return gn;
}

/** Function to identify the next n sequential integers that do not
 * identify a group that currently has member(s).
 *
 * @param d the SimData struct on which to perform the operation
 * @param n the number of group numbers to generate
 * @param result pointer to an array of length at least n where
 * the new group numbers generated can be saved.
 */
void get_n_new_group_nums( SimData* d, int n, int* result) {
	int n_groups = 0;
	int* existing_groups = get_existing_groups( d, &n_groups);
	int existingi = 0;
	int gn = 0;

	for (int i = 0; i < n; ++i) {
		++gn;
		while (existingi < n_groups) {
			if (gn < existing_groups[i]) {
				break;
			}

			++existingi;
			++gn;
		}
		result[i] = gn;
	}
	free(existing_groups);
}


//-----------------------------------Data Access-----------------------------------

/** Function to count the number of genotypes that currently belong to
 * the specified group.
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_bvs()
 * @see get_group_parent_ids()
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
		for (i = 0; i < m->n_genotypes; ++i) {
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
 * @see get_group_parent_ids()
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
 * shallow copies that should not be freed. Returns a null pointer if the group is empty.
 */
char** get_group_genes( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
	}
    char** genes = get_malloc(sizeof(char*) * group_size);
	int i, genes_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
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
 * @see get_group_parent_ids()
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
 * shallow copies that should not be freed. Returns a null pointer if the group is empty.
 */
char** get_group_names( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    char** names = get_malloc(sizeof(char*) * group_size);
	int i, names_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				names[names_i] = m->names[i];
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
 * @see get_group_parent_ids()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the ids of each member of the group.
 * The vector itself is on the heap and should be freed. Returns a null pointer if the group is empty.
 */
unsigned int* get_group_ids( SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    unsigned int* gids = get_malloc(sizeof(unsigned int) * group_size);
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
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
 * @see get_group_parent_ids()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the indexes of each member of the group.
 * The vector itself is on the heap and should be freed. Returns a null pointer if the group is empty.
 */
int* get_group_indexes(SimData* d, int group_id, int group_size) {
	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    int* gis = get_malloc(sizeof(int) * group_size);
	int i, total_i = 0, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i, ++total_i) {
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

/** Gets the breeding values/breeding values/fitnesses of each member of the group
 * @see get_group_size()
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 * @see get_group_parent_ids()
 * @see get_group_parent_names()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the breeding values of each member of the group.
 * The vector itself is on the heap and should be freed. Returns a null pointer if the group is empty.
 */
double* get_group_bvs( SimData* d, int group_id, int group_size) {
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    double* bvs = get_malloc(sizeof(double) * group_size);

	DecimalMatrix dm_bvs = calculate_group_bvs(d, (unsigned int) group_id);

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
 * Raises an error and returns NULL if this parameter is not either of those values.
 * @returns a vector containing the id of either the first or second parent of
 * each member of the group.
 * The vector itself is on the heap and should be freed. Returns a null pointer if the group is empty.
 */
unsigned int* get_group_parent_ids( SimData* d, int group_id, int group_size, int parent) {
	if (!(parent == 1 || parent == 2)) {
		warning( "Value error: `parent` must be 1 or 2.");
		return NULL;
	}
	--parent;

	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    unsigned int* pids = get_malloc(sizeof(unsigned int) * group_size);
	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
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
 * @see get_group_parent_ids()
 * @see get_group_pedigrees()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @param parent 1 to get the first parent of each group member, 2 to get the second.
 * Raises an error and returns NULL if this parameter is not either of those values.
 * @returns  a vector containing pointers to the names either the first or second
 * parent of each member of the group. The vector itself is on the heap and should be
 * freed, but its contents are only shallow copies that should not be freed. Returns a null pointer if the group is empty.
 */
char** get_group_parent_names( SimData* d, int group_id, int group_size, int parent) {
	if (!(parent == 1 || parent == 2)) {
		warning( "Value error: `parent` must be 1 or 2.");
		return NULL;
	}
	--parent;

	AlleleMatrix* m = d->m;
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
    }
    char** pnames = get_malloc(sizeof(char*) * group_size);

	int i, ids_i = 0;
	while (1) {
		for (i = 0; i < m->n_genotypes; ++i) {
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
 * @see get_group_parent_ids()
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
 * specifically to answer this function call. If it failed to write to and
 * read from the temporary file, then it returns NULL; Also returns a null pointer if the group is empty.
 */
char** get_group_pedigrees( SimData* d, int group_id, int group_size) {
	char* fname = "gS_gpptmp";
	FILE* fp = fopen(fname, "w");
	save_group_full_pedigree(fp, d, group_id);
	fclose(fp);

	FILE* fp2;
	if ((fp2 = fopen(fname, "r")) == NULL) {
		warning( "Failed to use temporary file.\n");
		return NULL;
	}

	// Create the list that we will return
    if (group_size <= 0) {
        group_size = get_group_size( d, group_id );
        if (group_size == 0) { warning("Group does not exist\n"); return 0; }
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
					warning( "Memory allocation of size %u failed.\n", size);
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


/*------------------------------------------------------------------------*/
/*----------------------Nicked from matrix-operations.c-------------------*/

/** Generates a matrix of c columns, r rows with all 0.
 *
 * @param r the number of rows/first index for the new matrix.
 * @param c the number of columns/second index for the new matrix.
 * @returns a DecimalMatrix with r rows, c cols, and a matrix of
 * the correct size, with all values zeroed.
 */
DecimalMatrix generate_zero_dmatrix(int r, int c) {
	DecimalMatrix zeros;
	zeros.rows = r;
	zeros.cols = c;

	zeros.matrix = malloc(sizeof(double*) * r);
	for (int i = 0; i < r; ++i) {
		zeros.matrix[i] = malloc(sizeof(double) * c);
		for (int j = 0; j < c; ++j) {
			zeros.matrix[i][j] = 0.0;
		}
	}
	return zeros;
}

/** Multiply a DecimalMatrix to a vector, and add that product to the first column
 * of a provided DecimalMatrix.
 *
 * Performs the fused multiply-add operation: `result + a * b` and saves it to `result`.
 *
 * Assumes that the vector `b` has the same number of entries as `a` has columns. That
 * is, assumes the dimensions are valid for multiplication.
 *
 * @param result pointer to the DecimalMatrix to whose first row the product a*b
 * will be added.
 * @param a pointer to the DecimalMatrix. Multiply this to b.
 * @param b a vector of doubles, assumed to be the same length as `a` has columns.
 * Multiplied to a.
 * @returns 0 on success, nonzero on failure.
 */
int add_matrixvector_product_to_dmatrix(DecimalMatrix* result, DecimalMatrix* a, double* b) {
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

/** Multiply two sets of a DecimalMatrix and vector, and add both products to
 * the first column of a provided DecimalMatrix.
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
 * @param result pointer to the DecimalMatrix to whose first row both products
 * will be added.
 * @param amat pointer to the first DecimalMatrix. Multiply this to avec.
 * @param avec a vector of doubles, assumed to be the same length as `amat` has columns.
 * Multiplied to `amat`.
 * @param bmat pointer to the second DecimalMatrix. Multiply this to bvec.
 * @param bvec a second vector of doubles, assumed to be the same length as
 * `bmat` has columns. Multiplied to `bmat`.
 * @returns 0 on success, nonzero on failure.
 */
int add_doublematrixvector_product_to_dmatrix(DecimalMatrix* result, DecimalMatrix* amat, double* avec, DecimalMatrix* bmat, double* bvec) {

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


/** Deletes a DecimalMatrix and frees its memory. m will now refer
 * to an empty matrix, with every pointer set to null and dimensions set to 0.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void delete_dmatrix(DecimalMatrix* m) {
	if (m->matrix != NULL) {
		for (int i = 0; i < m->rows; i++) {
			if (m->matrix[i] != NULL) {
				free(m->matrix[i]);
			}
		}
		free(m->matrix);
		m->matrix = NULL;
	}
	m->cols = 0;
	m->rows = 0;
}


/*--------------------------------Deleting-----------------------------------*/

/** Clears memory of all genotypes and associated information belonging
 * to a particular group.
 *
 * Uses a call to condense_allele_matrix() to ensure that the SimData
 * remains valid after deletion.
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the subset of data to be cleared
 */
void delete_group(SimData* d, int group_id) {
	AlleleMatrix* m = d->m;
	int i, deleted, total_deleted = 0;
	while (1) {

		for (i = 0, deleted = 0; i < m->n_genotypes; ++i) {
			if (m->groups[i] == group_id) {
				// delete data
				if (m->names[i] != NULL) {
					free(m->names[i]);
					m->names[i] = NULL;
				}
				if (m->alleles[i] != NULL) {
					free(m->alleles[i]);
					m->alleles[i] = NULL;
				}
				m->ids[i] = 0;
				m->pedigrees[0][i] = 0;
				m->pedigrees[1][i] = 0;
				m->groups[i] = 0;
				++deleted;
			}
		}
		m->n_genotypes -= deleted;
		total_deleted += deleted;

		if (m->next == NULL) {
			condense_allele_matrix( d );
			Rprintf("%d genotypes were deleted\n", total_deleted);
			return;
		} else {
			m = m->next;
		}
	}
}

/** Deletes a GeneticMap object and frees its memory. m will now refer
 * to an empty matrix, with every pointer set to null and dimensions set to 0.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void delete_genmap(GeneticMap* m) {
	m->n_chr = 0;
	if (m->chr_ends != NULL) {
		free(m->chr_ends);
	}
	m->chr_ends = NULL;
	if (m->chr_lengths != NULL) {
		free(m->chr_lengths);
	}
	m->chr_lengths = NULL;
	if (m->positions != NULL) {
		free(m->positions);
	}
	m->positions = NULL;
}

/** Deletes the full AlleleMatrix object and frees its memory. m will now refer
 * to an empty matrix, with pointers set to null and dimensions set to 0.
 *
 * freeing its memory includes freeing the AlleleMatrix, which was allocated on the heap
 * by @see create_empty_allelematrix()
 *
 * all subsequent matrixes in the linked list are also deleted.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void delete_allele_matrix(AlleleMatrix* m) {
	if (m == NULL) {
		return;
	}
	AlleleMatrix* next;
	do {
		/* free the big data matrix */
        for (int i = 0; i < m->n_genotypes; i++) {
            if (m->alleles[i] != NULL) {
                free(m->alleles[i]);
            }

        }

		// free names
        for (int i = 0; i < m->n_genotypes; i++) {
            if (m->names[i] != NULL) {
                free(m->names[i]);
            }
        }

		next = m->next;
		free(m);
	} while ((m = next) != NULL);
}

/** Deletes an EffectMatrix object and frees its memory. m will now refer
 * to an empty matrix, with every pointer set to null and dimensions set to 0.
 *
 * @param m pointer to the matrix whose data is to be cleared and memory freed.
 */
void delete_effect_matrix(EffectMatrix* m) {
	delete_dmatrix(&(m->effects));
	if (m->effect_names != NULL) {
		free(m->effect_names);
	}
	m->effect_names = NULL;
}

/** Deletes a SimData object and frees its memory.
 *
 * @param m pointer to the struct whose data is to be cleared and memory freed.
 */
void delete_simdata(SimData* m) {
	if (m == NULL) {
		return;
	}
	// free markers
	if (m->markers != NULL) {
		for (int i = 0; i < m->n_markers; i++) {
			if (m->markers[i] != NULL) {
				free(m->markers[i]);
			}
		}
		free(m->markers);
		m->markers = NULL;
	}
	//m->n_markers = 0;

	// free genetic map and effects
	delete_genmap(&(m->map));
	delete_effect_matrix(&(m->e));

	// free tables of alleles across generations
	delete_allele_matrix(m->m);

	//m->current_id = 0;
	free(m);
}

/** Deletes a MarkerBlocks object and frees its associated memory. b will now refer
 * to an empty struct, with every pointer set to null and number of markers set to 0.
 *
 * @param b pointer to the struct whose data is to be cleared and memory freed.
 */
void delete_markerblocks(MarkerBlocks* b) {
	for (int i = 0; i < b->num_blocks; ++i) {
		free(b->markers_in_block[i]);
	}
	free(b->markers_in_block);
	b->markers_in_block = NULL;
	free(b->num_markers_in_block);
	b->num_markers_in_block = NULL;
	b->num_blocks = 0;

	return;
}

/*-------------------------------SimData loaders-----------------------------*/

/** Populates a SimData combination with marker allele data.
 * @see load_transposed_encoded_genes_to_simdata()
 * Assumes it is starting from a clean/empty SimData.
 *
 * Given a file with the following format:
 *
 * name [line] [line] [line] ... [line]
 *
 * [marker] [SNP pair] [SNP pair] [SNP pair] ... [SNP pair]
 *
 * [marker] [SNP pair] [SNP pair] [SNP pair] ... [SNP pair]
 *
 * ...
 *
 * Where [line] is a code for a line, [marker] is a code for a marker, and
 * [SNP pair] is eg TT, TA.
 *
 * Note: this function should be called first when populating a SimData object -
 * it clears everything in the SimData. This is because all the data in SimData
 * is based on what markers exist in the loaded marker allele file.
 *
 * An output message stating the number of genotypes and number of markers loaded
 * is printed to stdout.
 *
 * @param d pointer to SimData to be populated
 * @param filename string containing name/path of file containing SNP marker
 * allele data.
 * @returns the group number of the loaded genotypes. All genotypes are loaded into
 * the same group.
*/
int load_transposed_genes_to_simdata(SimData* d, const char* filename) {
	struct TableSize t = get_file_dimensions(filename, '\t');
    char cell[4] = "\t%s";
    if (t.num_columns == 1) {
        t = get_file_dimensions(filename, ' ');
        cell[0] = ' ';
    }
    if (t.num_columns == 1) {
        warning( "Only found one column in file %s. File may be using an unsupported separator.\n", filename);
    }

	FILE* fp;
	const int gp = 1;
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.

	// discard the column title that is not a line
	char word[NAME_LENGTH];

	fscanf(fp, "%s", word);

	// now we want to read the header columns.
	// There are num_columns-1 of these because of the 'name' entry
	// this will also create our unique ids
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < CONTIG_WIDTH) {
		current_am = create_empty_allelematrix(t.num_rows - 1, n_to_go);
		d->m = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(t.num_rows - 1, CONTIG_WIDTH);
		d->m = current_am;
		n_to_go -= CONTIG_WIDTH;
		while (n_to_go) {
			if (n_to_go < CONTIG_WIDTH) {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, CONTIG_WIDTH);
				current_am = current_am->next;
				n_to_go -= CONTIG_WIDTH;
			}
		}
		current_am = d->m;
	}

	// load in the genotype names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, cell, word);

		if (i_am >= current_am->n_genotypes) {
			i_am = 0;
			current_am = current_am->next;
		}

		if (current_am == NULL) {
			warning( "Something went wrong during setup\n");
			// will occur if there's some bug in create_empty_allelematrix again
		}

		current_am->names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->names[i_am], word);

	}

	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");

	// set the ids for the genotypes we loaded
	set_ids(d, 0, t.num_columns - 2);

	// get space to put marker names and data we gathered
	d->n_markers = t.num_rows - 1;
	d->markers = get_malloc(sizeof(char*) * (t.num_rows-1));
	//memset(d->markers, '\0', sizeof(char*) * (t.num_rows-1));

	// now read the rest of the table.
	char word2[NAME_LENGTH];
	int badRows = 0;
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.

		// get the row name, store in markers
		fscanf(fp, "%s", word);
		d->markers[j] = get_malloc(sizeof(char) * strlen(word) + 1);
		strncpy(d->markers[j], word, strlen(word) + 1);

		current_am = d->m;
		//d->m->alleles[j] = get_malloc(sizeof(char) * d->m[0].n_genotypes * 2);
		for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
			// looping through the remaining columns in this row.
			fscanf(fp, cell, word2);

			// save the two alleles.
			if (strlen(word2) != 2) {
				++badRows;
				//warning("This file is invalid, but nothing will be done about it.\n");
			}

			if (i_am >= CONTIG_WIDTH) {//current_am->n_genotypes) {
				i_am = 0;
				current_am = current_am->next;
			}

			current_am->alleles[i_am][2*j] = word2[0];
			current_am->alleles[i_am][2*j + 1] = word2[1];
			current_am->groups[i_am] = gp;
		}
	}

	fclose(fp);
	Rprintf("%d genotypes of %d markers were loaded. %d pairs of alleles could not be loaded\n", (t.num_columns - 1), (t.num_rows - 1), badRows);
	return gp;
}

/** Populates a SimData combination with marker allele data.
 * @see load_transposed_genes_to_simdata()
 * Assumes it is starting from a clean/empty SimData.
 *
 * Given a file with the following format:
 *
 * name [line] [line] [line] ... [line]
 *
 * [marker] [encoded] [encoded] [encoded] ... [encoded]
 *
 * [marker] [encoded] [encoded] [encoded] ... [encoded]
 *
 * ...
 *
 * Where [line] is a code for a line, [marker] is a code for a marker, and
 * [encoded] is the standard IUPAC encoding for a particular pair. Because this simulation
 * tracks phase, and this encoding does not, the phase at heterozygous markers is
 * chosen randomly.
 *
 * Code => Alleles key:
 * A => AA    ; C => CC    ; G => GG    ; T => TT   ;
 * R => AG    ; Y => CT    ; S => CG    ; W => AT   ; K => GT   ; M => AC
 *
 * Note: this function should be called first when populating a SimData object -
 * it clears everything in the SimData. This is because all the data in SimData
 * is based on what markers exist in the loaded marker allele file.
 *
 * An output message stating the number of genotypes and number of markers loaded
 * is printed to stdout.
 *
 * @param d pointer to SimData to be populated
 * @param filename string containing name/path of file containing SNP marker
 * allele data.
 * @returns the group number of the loaded genotypes. All genotypes are loaded into
 * the same group.
*/
int load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename) {
	struct TableSize t = get_file_dimensions(filename, '\t');
    char cell[4] = "\t%s";
    if (t.num_columns == 1) {
        t = get_file_dimensions(filename, ' ');
        cell[0] = ' ';
    }
    if (t.num_columns == 1) {
        warning( "Only found one column in file %s. File may be using an unsupported separator.\n", filename);
    }

	FILE* fp;
	const int gp = 1;
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.

	// discard the column title that is not a line
	char word[NAME_LENGTH];

	fscanf(fp, "%s", word);

	// now we want to read the header columns.
	// There are num_columns-1 of these because of the 'name' entry
	// this will also create our unique ids
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < CONTIG_WIDTH) {
		current_am = create_empty_allelematrix(t.num_rows - 1, n_to_go);
		d->m = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(t.num_rows - 1, CONTIG_WIDTH);
		d->m = current_am;
		n_to_go -= CONTIG_WIDTH;
		while (n_to_go) {
			if (n_to_go < CONTIG_WIDTH) {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, CONTIG_WIDTH);
				current_am = current_am->next;
				n_to_go -= CONTIG_WIDTH;
			}
		}
		current_am = d->m;
	}


	// load in the genotypes' names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, cell, word);

		if (i_am >= current_am->n_genotypes) {
			i_am = 0;
			current_am = current_am->next;
		}

		current_am->names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->names[i_am], word);
	}

	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");

	// set the ids for the genotypes we loaded
	set_ids(d, 0, t.num_columns - 2);

	// get space to put marker names and data we gathered
	d->n_markers = t.num_rows - 1;
	d->markers = get_malloc(sizeof(char*) * (t.num_rows-1));
	//memset(d->markers, '\0', sizeof(char*) * (t.num_rows-1));

	// now read the rest of the table.
	GetRNGstate();
	char c, decoded[2];
	int r;
	cell[2] = 'c';
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.

		// get the row name, store in markers
		fscanf(fp, "%s", word);
		d->markers[j] = get_malloc(sizeof(char) * strlen(word) + 1);
		strncpy(d->markers[j], word, strlen(word) + 1);

		current_am = d->m;
		//d->m->alleles[j] = get_malloc(sizeof(char) * d->m[0].n_genotypes * 2);
		for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
			// looping through the remaining columns in this row.
			fscanf(fp, cell, &c);

			if (i_am >= current_am->n_genotypes) {
				i_am = 0;
				current_am = current_am->next;
			}

			// if it's a homozygous code, just copy directly over.
			if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
				current_am->alleles[i_am][(j<<1)] = c;
				current_am->alleles[i_am][(j<<1) + 1] = c;
			} else {
				// choose a random order for the two alleles.
				r = (unif_rand() > 0.5);
				// identify the two alleles
				switch (c) {
					case 'R':
						decoded[0] = 'A'; decoded[1] = 'G'; break;
					case 'Y':
						decoded[0] = 'C'; decoded[1] = 'T'; break;
					case 'S':
						decoded[0] = 'C'; decoded[1] = 'G'; break;
					case 'W':
						decoded[0] = 'A'; decoded[1] = 'T'; break;
					case 'K':
						decoded[0] = 'A'; decoded[1] = 'T'; break;
					case 'M':
						decoded[0] = 'A'; decoded[1] = 'T'; break;
					default:
						current_am->alleles[i_am][(j<<1)] = 0; current_am->alleles[i_am][(j<<1) + 1] = 0;

						continue;
				}

				current_am->alleles[i_am][(j<<1)] = decoded[r];
				current_am->alleles[i_am][(j<<1) + 1] = decoded[1-r];
				current_am->groups[i_am] = gp;
			}
		}
	}
	PutRNGstate();
	fclose(fp);
	return gp;
}

/** Appends genotype data from a file to an existing SimData
 * @see load_transposed_genes_to_simdata()
 *
 * Given a file with the following format:
 *
 * name [line] [line] [line] ... [line]
 *
 * [marker] [SNP pair] [SNP pair] [SNP pair] ... [SNP pair]
 *
 * [marker] [SNP pair] [SNP pair] [SNP pair] ... [SNP pair]
 *
 * ...
 *
 * Where [line] is a code for a line, [marker] is a code for a marker, and
 * [SNP pair] is eg TT, TA.
 *
 * If a given marker does not exist in the SimData's set of markers, it is ignored.
 * for the purposes of loading. No markers can be added to a SimData after the creation
 * step.
 *
 * An output message stating the number of genotypes and number of markers loaded
 * is printed to stdout.
 *
 * @param d pointer to SimData to be populated
 * @param filename string containing name/path of file containing SNP marker
 * allele data.
 * @returns the group number of the loaded genotypes. All genotypes are loaded into
 * the same group.
*/
int load_more_transposed_genes_to_simdata(SimData* d, const char* filename) {
	struct TableSize t = get_file_dimensions(filename, '\t');
    char cell[4] = "\t%s";
    if (t.num_columns == 1) {
        t = get_file_dimensions(filename, ' ');
        cell[0] = ' ';
    }
    if (t.num_columns == 1) {
        warning( "Only found one column in file %s. File may be using an unsupported separator.\n", filename);
    }

	FILE* fp;
	int gp = get_new_group_num(d);
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.

	// discard the column title that is not a line
	char word[NAME_LENGTH];

	fscanf(fp, "%s", word);

	// now we want to read the header columns.
	// There are num_columns-1 of these because of the 'name' entry
	// this will also create our unique ids

	// find the end of the AM chain so far
	AlleleMatrix* last_am = d->m;
	int last_n_genotypes = last_am->n_genotypes;
	while (last_am->next != NULL) {
		last_am = last_am->next;
		last_n_genotypes += last_am->n_genotypes;
	}

	// Create new AMs that will be populated from the file.
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < CONTIG_WIDTH) {
		current_am = create_empty_allelematrix(d->n_markers, n_to_go);
		last_am->next = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
		last_am->next = current_am;
		n_to_go -= CONTIG_WIDTH;
		while (n_to_go) {
			if (n_to_go <= CONTIG_WIDTH) {
				current_am->next = create_empty_allelematrix(d->n_markers, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
				current_am = current_am->next;
				n_to_go -= CONTIG_WIDTH;
			}
		}
		current_am = last_am->next;
	}

	// set the ids for the genotypes we loaded
	set_ids(d, last_n_genotypes, last_n_genotypes + t.num_columns - 2);

	// load in the genotypes' names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, cell, word);

		if (i_am >= current_am->n_genotypes) {
			i_am = 0;
			current_am = current_am->next;
		}

		current_am->names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->names[i_am], word);
	}

	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");

	// get space to put marker names and data we gathered
	// now read the rest of the table.
	char word2[NAME_LENGTH];
	int markeri;
	current_am = last_am->next;
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.

		// get the row name, store in markers
		fscanf(fp, "%s", word);
		markeri = get_from_unordered_str_list(word, d->markers, d->n_markers);

		current_am = last_am->next;

		if (markeri >= 0) {
			for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
				// looping through the remaining columns in this row.
				fscanf(fp, cell, word2);

				// save the two alleles.
				if (strlen(word2) != 2) {
					warning( "This file is invalid, but nothing will be done about it.\n");
				}

				if (i_am >= current_am->n_genotypes) {
					i_am = 0;
					current_am = current_am->next;
				}

				//strncpy(d->m->alleles[i] + (2*j), word2, 2);
				current_am->alleles[i_am][2*markeri] = word2[0];
				current_am->alleles[i_am][2*markeri + 1] = word2[1];
				current_am->groups[i_am] = gp;
			}
		} else {
			warning( "Could not find the marker %s\n", word);
		}
	}
	condense_allele_matrix(d);

	Rprintf("%d genotypes were loaded.\n", t.num_columns - 1);
	fclose(fp);
	return gp;
}

/** Populates a SimData combination with data from a genetic map. Map positions must be in cM.
 *
 * Note: this function should be called second when populating a SimData object,
 * after populating it with marker allele data. This is because this function
 * loads the genmap data corresponding to markers used in the allele file, then
 * deletes markers in SimData that do not have positions for ease of simulation.
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
 * The function assumes the maximum line length is 99 characters.
 * It also assumes that there is only one mapping per marker in the file.
 *
 * @param d pointer to SimData to be populated
 * @param filename string name/path of file containing genetic map data.
*/
void load_genmap_to_simdata(SimData* d, const char* filename) {
	// open our file.
	FILE* fp;
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	// ignore the first line of the file
	fscanf(fp, "%*[^\n]\n");

	int bufferlen = 100; // assume no line is over 100 characters long
	char buffer[bufferlen];
	char marker_name[bufferlen]; // for scanning name from line
	int chr; // for scanning chromosome from line
	float pos; // for scanning postion value from line
	int location; // used for location of marker in m->marker_names
	int positions_loaded = 0;

	if (d->map.positions != NULL) {
		delete_genmap(&(d->map));
	}

	d->map.positions = calloc(sizeof(MarkerPosition) * d->n_markers, sizeof(MarkerPosition));

	// loop through rows of the file (until we've got all our positions)
	while (fgets(buffer, bufferlen, fp) != NULL && (positions_loaded < d->n_markers)) {
		R_CheckUserInterrupt();
		sscanf(buffer, "%s %d %f\n", marker_name, &chr, &pos);

		if ((location = get_from_unordered_str_list( marker_name, d->markers, d->n_markers)) >= 0) {
			// the marker is in our list, so save its position
			d->map.positions[location].chromosome = chr;
			d->map.positions[location].position = pos;
		}
	}

	// count number of markers that don't have positions loaded.
	int n_nopos = 0;
	for (int i = 0; i < d->n_markers; i++) {
		if (d->map.positions[i].chromosome == 0) {
			n_nopos += 1;
		}
	}

	Rprintf("%d markers with map positions. %d markers remain unmapped.\n",
	d->n_markers - n_nopos, n_nopos);

	fclose(fp);

	//Order the markers and positions, eliminating markers with no positions
	if (n_nopos > 0) {
		get_sorted_markers(d, d->n_markers - n_nopos);
		get_chromosome_locations(d);
	}
}

/** Takes a SimData object, and sorts its markers, the rows of its parent gen
 * AlleleMatrix (because they are ordered by the markers), and its genetic map
 * so that the markers are ordered by chromosome number then position.
 *
 * Markers that do not have a position in d->map.positions are deleted from
 * all those three lists. This is done by malloc-ing new memory, copying
 * data over in the new sorted order, and freeing the old memory.
 *
 * Note: only the starting generation in the AlleleMatrix list is reordered, and
 * all additional AlleleMatrix objects are deleted.
 *
 * @param d pointer to the SimData to have its markers, genetic map, and allele
 * matrix sorted. The SimData pointed to by d will be modified by this function.
 * @param actual_n_markers If previously calculated, include the number of
 * markers in `d->markers` that have a position loaded into `d->map.positions`.
 * If this has not been calculated yet, make the value -1 and this function will
 * calculate it. The value is calculated as the number of positions in
 * `d->map.positions` that have a chromosome number of 0.
*/
void get_sorted_markers(SimData* d, int actual_n_markers) {
	MarkerPosition* sortable[d->n_markers];
	for (int i = 0; i < d->n_markers; i++) {
		sortable[i] = &(d->map.positions[i]);
	}

	// if this was not pre-calculated do it now.
	if (actual_n_markers < 0) {
		actual_n_markers = d->n_markers;
		for (int i = 0; i < d->n_markers; i++) {
			if (d->map.positions[i].chromosome == 0) {
				actual_n_markers -= 1;
			}
		}
	}

	/* Sort the pointers */
	qsort(sortable, d->n_markers, sizeof(sortable[0]), _simdata_pos_compare);
	int location_in_old;

	R_CheckUserInterrupt();

	if (d->markers != NULL) {
		char** new_markers = get_malloc(sizeof(char*) * actual_n_markers);
		for (int i = 0; i < actual_n_markers; ++i) {
			location_in_old = sortable[i] - d->map.positions;
			new_markers[i] = d->markers[location_in_old]; // shallow copy
		}

		free(d->markers);
		d->markers = new_markers;
	}

	char* temp;
    if (d->m != NULL) {
		//temp = get_malloc(sizeof(char) * ((actual_n_markers * 2)));
		AlleleMatrix* am = d->m;

		do {
			for (int i = 0; i < am->n_genotypes; ++i) {
				R_CheckUserInterrupt();

				//strncpy(temp, am->alleles[i], sizeof(char) * ((am->n_markers * 2)));
				temp = get_malloc(sizeof(char) * ((actual_n_markers * 2)));

				for (int j = 0; j < actual_n_markers; ++j) {
					location_in_old = sortable[j] - d->map.positions;
					temp[2*j] = am->alleles[i][2*location_in_old];
					temp[2*j + 1] = am->alleles[i][2*location_in_old + 1];

				}
				free(am->alleles[i]);
				am->alleles[i] = temp;
			}
			am->n_markers = actual_n_markers;
		} while ((am = am->next) != NULL);
		//free(temp);
	}

	if (d->e.effects.matrix != NULL) {
		// Don't need to update row names, just matrix.
		DecimalMatrix new_eff = generate_zero_dmatrix(d->e.effects.rows, actual_n_markers);
		for (int i = 0; i < actual_n_markers; ++i) {
			R_CheckUserInterrupt();
			location_in_old = sortable[i] - d->map.positions;
			for (int j = 0; j < d->e.effects.rows; ++j) {
				new_eff.matrix[j][i] = d->e.effects.matrix[j][location_in_old];
			}
		}

		delete_dmatrix(&(d->e.effects));
		d->e.effects = new_eff;
	}

	if (d->map.positions != NULL) {
		MarkerPosition* new_map = get_malloc(sizeof(MarkerPosition) * actual_n_markers);
		for (int i = 0; i < actual_n_markers; ++i) {
			R_CheckUserInterrupt();
			location_in_old = sortable[i] - d->map.positions;
			new_map[i].chromosome = d->map.positions[location_in_old].chromosome;
			new_map[i].position = d->map.positions[location_in_old].position;
		}

		delete_genmap(&(d->map));
		d->map.positions = new_map;
	}
	d->n_markers = actual_n_markers;

}

/** Updates the chr_ends, n_chr and chr_lengths fields in SimData.map.
 *
 * This should only be run on a SimData that has already been ordered.
 * @see sort_markers()
 *
 * The function loops over all MarkerPosition in SimData.map.positions
 * twice, for flexible and minimal memory usage rather than maximum speed.
 *
 * Chromosome lengths are intialised to 0 for chromosomes that contain
 * no markers, and to a flat (meaningless) value of 1 for chromosomes that
 * contain exactly one marker.
 *
 * @param d pointer to the SimData object for which the fields under `map`
 * need to be initialised or updated.
 */
void get_chromosome_locations(SimData *d) {
	// count the chromosomes
	int highest_chr_found = 0;

	for (int i = 0; i < d->n_markers; i++) {
		if (d->map.positions[i].chromosome > highest_chr_found) {
			highest_chr_found = d->map.positions[i].chromosome;
		}
	}
	d->map.n_chr = highest_chr_found;

	if (d->map.chr_ends != NULL) {
		free(d->map.chr_ends);
	}
	if (d->map.chr_lengths != NULL) {
		free(d->map.chr_lengths);
	}

	// identify the start/end points of all chromosomes
	d->map.chr_ends = get_malloc(sizeof(int) * (highest_chr_found + 1));
	d->map.chr_lengths = get_malloc(sizeof(float) * (highest_chr_found));

	highest_chr_found = 0;
	for (int i = 0; i < d->n_markers; i++) {
		R_CheckUserInterrupt();
		if (d->map.positions[i].chromosome == highest_chr_found + 1) {
			highest_chr_found = d->map.positions[i].chromosome;
			d->map.chr_ends[highest_chr_found - 1] = i;
		} else if (d->map.positions[i].chromosome > highest_chr_found) {
			// deal with chromosomes that have no markers
			for (int j = highest_chr_found; j < d->map.positions[i].chromosome; j++) {
				d->map.chr_ends[j] = i;
			}

			highest_chr_found = d->map.positions[i].chromosome;
		}
	}
	// and add on the end index
	d->map.chr_ends[d->map.n_chr] = d->n_markers;

	// calculate lengths
	for (int i = 0; i < d->map.n_chr; i++) {
		if (d->map.chr_ends[i+1] - 1 > d->map.chr_ends[i]) { //more than 1 marker in chr
			d->map.chr_lengths[i] = d->map.positions[d->map.chr_ends[i+1] - 1].position
					- d->map.positions[d->map.chr_ends[i]].position;

		} else if (d->map.chr_ends[i+1] - 1 == d->map.chr_ends[i]) { // exactly 1 marker in chr
			d->map.chr_lengths[i] = 1; //pretty much arbitrary, won't affect crossing anyway
		} else { // no markers tracked at this chr
			d->map.chr_lengths[i] = 0;
		}
	}

	//srand(time(NULL)); //seed the random generator so we're ready to start crossing.
}

/** Populates a SimData combination with effect values. The SimData must already
 * have its allele data and map data loaded (so that it has an ordered `markers`
 * list and no markers that will not be used for simulation.
 *
 * It loads in the file as rows of effects for each allele that appears in the
 * allele data out of 'A', 'C', 'G', 'T'.
 *
 * The file should have format:
 *
 * [marker] [allele] [effect]
 *
 * [marker] [allele] [effect]
 *
 * ...
 *
 * The function assumes the maximum line length is 99 characters.
 * It also assumes that the array ref_alleles is the same
 * length as m's marker_names vector.
 *
 * @param d pointer to SimData to be populated.
 * @param filename string name/path of file containing effect values.
*/
void load_effects_to_simdata(SimData* d, const char* filename) {
	// open our file.
	FILE* fp;
	if ((fp = fopen(filename, "r")) == NULL) {
		warning( "Failed to open file %s.\n", filename);
		return;
	}

	int bufferlen = 100; // assume no line is over 100 characters long
	char buffer[bufferlen];
	char marker_name[bufferlen]; // for scanning name from line
	char allele; // for scanning allele from line
	double effect; // for scanning effect value from line
	int location; // used for location of marker in m->marker_names

	int n_loaded = 0;
	int n_allele = 0; // count the different alleles we're tracking
	const int MAX_SYMBOLS = 25;
	char alleles_loaded[MAX_SYMBOLS + 1];
	memset(alleles_loaded, '\0', MAX_SYMBOLS + 1);
	double* effects_loaded[MAX_SYMBOLS];
	memset(effects_loaded, 0, MAX_SYMBOLS);

	if (d->e.effects.matrix != NULL) {
		delete_dmatrix(&(d->e.effects));
	}
	if (d->e.effect_names != NULL) {
		free(d->e.effect_names);
	}

	// loop through rows of the file
	//for (int i = 0; i < (t.num_rows - 1); i++) {
	while (fgets(buffer, bufferlen, fp) != NULL) {
		R_CheckUserInterrupt();
		//fgets(buffer, bufferlen, fp);
		sscanf(buffer, "%s %c %lf\n", marker_name, &allele, &effect);

		if ((location = get_from_unordered_str_list(  marker_name, d->markers, d->n_markers)) >= 0) {
			// if the marker exists (at index `location`) find the allele index
			int symbol_index;
			char* symbol_location = strchr(alleles_loaded, allele);
			if (symbol_location == NULL) {
				symbol_index = n_allele;
				++n_allele;
				alleles_loaded[symbol_index] = allele;
			} else {
				symbol_index = symbol_location - alleles_loaded; // difference between the pointers
			}

			// now the marker is in our list and the allele value is valid
			// so save the effect value in effects_loaded.
			if (effects_loaded[symbol_index] == NULL) {
				effects_loaded[symbol_index] = calloc(d->n_markers, sizeof(double) * d->n_markers);
			}
			effects_loaded[symbol_index][location] = effect;
			n_loaded += 1;
		}
	}

	d->e.effects.matrix = get_malloc(sizeof(double*) * n_allele);
	d->e.effects.rows = n_allele;
	d->e.effects.cols = d->n_markers;
	d->e.effect_names = get_malloc(sizeof(char) * (n_allele + 1));

	// loop again to save values now we have enough memory.
	for (int i = 0; i < n_allele; i++) {
		d->e.effect_names[i] = alleles_loaded[i];
		d->e.effects.matrix[i] = effects_loaded[i];
	}
	d->e.effect_names[n_allele] = '\0'; // string terminator

	// integer division is intended here.
	Rprintf("%d effect values spanning %d alleles loaded.\n", n_loaded, n_allele);

	fclose(fp);
	return;
}

/** Populates a SimData combination from scratch with marker allele data, a genetic map, and
 * effect values.
 *
 * Note: this function shows the order that files need to be loaded into SimData,
 * i.e. Allele data first (so we know what markers we care about), then genetic
 * map data (so we know what markers we can use in simulation, and rearrange them
 * to be ordered), then effects data (to be saved in the correct order according
 * to newly ordered markers).
 *
 * @param d pointer to SimData to be populated
 * @param data_file string containing name/path of file containing SNP marker
 * allele data.
 * @param map_file string name/path of file containing genetic map data.
 * @param effect_file string name/path of file containing effect values.
*/
int load_all_simdata(SimData* d, const char* data_file, const char* map_file, const char* effect_file) {
	delete_simdata(d); // make this empty.
	d = create_empty_simdata();
	int gp = load_transposed_genes_to_simdata(d, data_file);

	load_genmap_to_simdata(d, map_file);

	if (effect_file != NULL) {
		load_effects_to_simdata(d, effect_file);
	}

	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	return gp;
}


/*--------------------------Recombination counts-----------------------------*/

/** Identify markers in the genotype of `offspring` where recombination from its parents
 * occured. This function is a little lower-level (see the kinds of parameters required) and
 * so a wrapper like calculate_recombinations_from_file is suggested for end users.
 * @see calculate_recombinations_from_file()
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
 * @param d pointer to the SimData struct whose genetic map matches the provided genotypes.
 * @param parent1 a character vector containing one parent's alleles at each marker in the
 * SimData.
 * @param p1num an integer that will be used to identify areas of the genome that come
 * from the first parent in the returned vector.
 * @param parent2 a character vector containing the other parent's alleles at each marker in the
 * SimData.
 * @param p2num an integer that will be used to identify areas of the genome that come
 * from the second parent in the returned vector.
 * @param offspring a character vector containing the alleles at each marker in the
 * SimData of the genotype whose likely recombinations we want to identify.
 * @param certain a boolean. If TRUE, markers where the parent of origin cannot be identified
 * will be set to 0, if FALSE, the value will be set to the id of the parent that provided
 * the most recently identified allele in that chromosome.
 * @returns a heap vector of length `d->n_markers` containing the id of the parent of origin
 * at each marker in the `offspring` genotype.
*/
int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int certain) {
	int* origins = malloc(sizeof(int) * d->n_markers);
	int p1match, p2match;
	int previous = 0;

	// treat each chromosome separately.
	for (int chr = 1; chr <= d->map.n_chr; ++chr) {
		R_CheckUserInterrupt();

		previous = 0;
		for (int i = d->map.chr_ends[chr - 1]; i < d->map.chr_ends[chr]; ++i) {
			p1match = has_same_alleles(parent1, offspring, i);
			p2match = has_same_alleles(parent2, offspring, i);
			if (p1match && !p2match) {
				origins[i] = p1num;
				previous = p1num;
			} else if (p2match && !p1match) {
				origins[i] = p2num;
				previous = p2num;
			} else {
				if (certain) {
					origins[i] = 0;
				} else {
					origins[i] = previous;
				}
			}
		}
	}
	return origins;
}

/** Identify markers in the genotype of `offspring` where recombination from its parents
 * occured, as judged by the marker itself and a short window around it.
 * This function is a little lower-level (see the kinds of parameters required) and
 * so a wrapper like calculate_recombinations_from_file is suggested for end users.
 * @see calculate_recombinations_from_file()
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
 * @param d pointer to the SimData struct whose genetic map matches the provided genotypes.
 * @param parent1 a character vector containing one parent's alleles at each marker in the
 * SimData.
 * @param p1num an integer that will be used to identify areas of the genome that come
 * from the first parent in the returned vector.
 * @param parent2 a character vector containing the other parent's alleles at each marker in the
 * SimData.
 * @param p2num an integer that will be used to identify areas of the genome that come
 * from the second parent in the returned vector.
 * @param offspring a character vector containing the alleles at each marker in the
 * SimData of the genotype whose likely recombinations we want to identify.
 * @param window_size an odd integer representing the number of markers to check for known parentage
 * around each marker
 * @param certain a boolean. If TRUE, markers where the parent of origin cannot be identified
 * will be set to 0, if FALSE, the value will be set to the id of the parent that provided
 * the most recently identified allele in that chromosome.
 * @returns a heap vector of length `d->n_markers` containing the id of the parent of origin
 * at each marker in the `offspring` genotype.
*/
int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int window_size, int certain) {
	int* origins = malloc(sizeof(int) * d->n_markers);
	int p1match, p2match;
	int previous = 0, window_range = (window_size - 1)/2, i;
	int lookable_bounds[2];

	// treat each chromosome separately.
	for (int chr = 1; chr <= d->map.n_chr; ++chr) {
		R_CheckUserInterrupt();

		previous = 0;
		lookable_bounds[0] = d->map.chr_ends[chr - 1] + window_range;
		lookable_bounds[1] = d->map.chr_ends[chr] - window_range;

		for (i = d->map.chr_ends[chr - 1]; i < lookable_bounds[0]; ++i) {
			origins[i] = 0;
		}
		for (; i < lookable_bounds[1]; ++i) {

			p1match = has_same_alleles_window(parent1, offspring, i, window_size);
			p2match = has_same_alleles_window(parent2, offspring, i, window_size);
			if (p1match && !p2match) {
				origins[i] = p1num;
				previous = p1num;
			} else if (p2match && !p1match) {
				origins[i] = p2num;
				previous = p2num;
			} else {
				if (certain) {
					origins[i] = 0;
				} else {
					origins[i] = previous;
				}
			}
		}
		for (; i < d->map.chr_ends[chr]; ++i) {
			origins[i] = 0;
		}
	}
	return origins;
}

/* static inline int has_same_alleles(char* p1, char* p2, int i) {
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
} */

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
 * @param d pointer to the SimData struct containing the genotypes and map under consideration.
 * @param input_file string containing the name of the file with the pairs of parents
 * and offsprings of which to calculate recombinations
 * @param output_file string containing the filename to which to save the results.
 * @param window_len an odd integer representing the number of markers to check for known parentage
 * around each marker
 * @param certain TRUE to fill locations where parentage is unknown with 0, FALSE
 * to fill locations where parentage is unknown with the most recent known parent
 * @returns 0 on success.
 */
int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file,
		int window_len, int certain) {
	struct TableSize t = get_file_dimensions(input_file, '\t');
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
	for (int j = 0; j < d->n_markers; ++j) {
		fprintf(fpo, "\t%s", d->markers[j]);
	}

	int combin_i[3];
	char* combin_genes[3];
	char buffer[3][50];
	int* r;
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s %s \n", buffer[0], buffer[1], buffer[2]);
		combin_i[0] = get_index_of_name(d->m, buffer[0]);
		combin_i[1] = get_index_of_name(d->m, buffer[1]);
		combin_i[2] = get_index_of_name(d->m, buffer[2]);
		combin_genes[0] = get_genes_of_index(d->m, combin_i[0]);
		combin_genes[1] = get_genes_of_index(d->m, combin_i[1]);
		combin_genes[2] = get_genes_of_index(d->m, combin_i[2]);

		if (window_len == 1) {
			r = calculate_min_recombinations_fw1(d, combin_genes[1],
					get_id_of_index(d->m, combin_i[1]), combin_genes[2],
					get_id_of_index(d->m, combin_i[2]), combin_genes[0], certain);
		} else {
			r = calculate_min_recombinations_fwn(d, combin_genes[1],
					get_id_of_index(d->m, combin_i[1]), combin_genes[2],
					get_id_of_index(d->m, combin_i[2]), combin_genes[0], window_len, certain);
		}
		fprintf(fpo, "\n%s", buffer[0]);
		for (int j = 0; j < d->n_markers; ++j) {
			fprintf(fpo, "\t%d", r[j]);
		}
		free(r);
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
 * from a Poisson distribution with parameter corresponding to the length of
 * the chromosome (in Morgans.).
 *
 * It generates the positions of those crossover events from a uniform distribution.
 *
 * It picks a random one of the gametes generated by picking randomly which column
 * of the parent's alleles to start with. When crossover events occur it starts reading
 * the other column.
 *
 * @param d pointer to the SimData object containing map positions for the markers
 * that make up the rows of `parent_table`. sort_markers() and locate_chromosomes()
 * should have been called previously.
 * @param parent_genome the char* containing the parent's genome as a character string
 * made up of sequential pairs of alleles for each marker in d->markers.
 * @param output the char* to which to save the gamete. It saves the alleles every second
 * character, starting at 0, so that calling generate_gamete(..., offspring_genome) &
 * generate_gamete(..., offspring_genome + 1) can be used to generate both halves of its genome.
*/
void generate_gamete(SimData* d, char* parent_genome, char* output) {
	// assumes rand is already seeded
	if (parent_genome == NULL) {
		warning( "Could not generate this gamete\n");
		return;
	}

	int num_crossovers, up_to_crossover, which;
	float crossover_where[100];
	float* p_crossover_where;

	// treat each chromosome separately.
	for (int chr = 1; chr <= d->map.n_chr; ++chr) {
		// use Poisson distribution to choose the number of crossovers in this chromosome
		num_crossovers = Rf_rpois(d->map.chr_lengths[chr - 1] / 100);

		// in the rare case where it could be >100, get enough space
		// to be able to store the crossover positions we're about to create
		if (num_crossovers <= 100) {
			p_crossover_where = crossover_where; // point at the start of array
		} else {
			p_crossover_where = get_malloc(sizeof(float) * num_crossovers);
		}

		// TASK 3: choose points where those crossovers occur
		// by randomly generating a point along the length of the chromosome
		for (int i = 0; i < num_crossovers; ++i) {
			p_crossover_where[i] = unif_rand()
				* d->map.chr_lengths[chr - 1]
				+ d->map.positions[d->map.chr_ends[chr - 1]].position;
		}

		// sort the crossover points
		if (num_crossovers > 1) {
			qsort(p_crossover_where, num_crossovers, sizeof(float),
					_ascending_float_comparer);
		}

		// pick a parent genome half at random
		which = (unif_rand() > 0.5); // if this is 0, we start with the left.

		// TASK 4: Figure out the gamete that those numbers produce.
		up_to_crossover = 0; // which crossovers we've dealt with
		for (int i = d->map.chr_ends[chr - 1]; i < d->map.chr_ends[chr]; ++i) {
			// loop through every marker for this chromosome
			if (up_to_crossover < num_crossovers &&
					d->map.positions[i].position > p_crossover_where[up_to_crossover]) {
				// if we're here then between last loop and this one we crossed over.
				// invert which and update up_to_crossover;
				which = 1 - which;
				up_to_crossover += 1;
			}
			output[2*i] = parent_genome[2*i + which];
		}

		if (num_crossovers > 100) {
			free(p_crossover_where);
		}

	}
	return;
}

/** Get the alleles of the outcome of crossing two genotypes
 *
 * Gametes are generated at the same time but are independent.
 *
 * @see generate_gamete(), on which this is based.
 *
 * @param d pointer to the SimData object that includes genetic map data
 * needed to simulate meiosis and the value of n_markers
 * @param parent1_genome a 2x(n_markers) array of characters containing the
 * alleles of the first parent
 * @param parent2_genome a 2x(n_markers) array of characters containing the
 * alleles of the second parent.
 * @param output a 2x(n_marker) array of chars which will be overwritten
 * with the offspring genome.
*/
void generate_cross(SimData* d, char* parent1_genome, char* parent2_genome, char* output) {
	// assumes rand is already seeded
	if (parent1_genome == NULL || parent2_genome == NULL) {
		warning( "Could not generate this cross\n");
		return;
	}

	int num_crossovers[2], up_to_crossover[2], which[2];
	float* p_crossover_where[2];
	float crossover_where[2][50];

	// treat each chromosome separately.
	for (int chr = 1; chr <= d->map.n_chr; ++chr) {
		// use Poisson distribution to choose the number of crossovers in this chromosome
		num_crossovers[0] = Rf_rpois(d->map.chr_lengths[chr - 1] / 100);
		num_crossovers[1] = Rf_rpois(d->map.chr_lengths[chr - 1] / 100);

		// in the rare case where it could be >100, get enough space
		// to be able to store the crossover positions we're about to create
		if (num_crossovers[0] <= 50) {
			p_crossover_where[0] = crossover_where[0]; // point to start of array
		} else {
			p_crossover_where[0] = get_malloc(sizeof(float) * num_crossovers[0]);
		}
		if (num_crossovers[1] <= 50) {
			p_crossover_where[1] = crossover_where[0]; // point to start of array
		} else {
			p_crossover_where[1] = get_malloc(sizeof(float) * num_crossovers[1]);
		}

		// TASK 3: choose points where those crossovers occur
		// by randomly generating a point along the length of the chromosome
		for (int i = 0; i < num_crossovers[0]; ++i) {
			p_crossover_where[0][i] = unif_rand()
				* d->map.chr_lengths[chr - 1]
				+ d->map.positions[d->map.chr_ends[chr - 1]].position;
		}
		for (int i = 0; i < num_crossovers[1]; ++i) {
			p_crossover_where[1][i] = unif_rand()
				* d->map.chr_lengths[chr - 1]
				+ d->map.positions[d->map.chr_ends[chr - 1]].position;
		}

		// sort the crossover points
		if (num_crossovers[0] > 1) {
			qsort(p_crossover_where[0], num_crossovers[0], sizeof(float),
					_ascending_float_comparer);
		}
		if (num_crossovers[1] > 1) {
			qsort(p_crossover_where[1], num_crossovers[1], sizeof(float),
					_ascending_float_comparer);
		}

		// pick a parent genome half at random
		which[0] = (unif_rand() > 0.5); which[1] = (unif_rand() > 0.5); // if this is 0, we start with the left.

		// TASK 4: Figure out the gamete that those numbers produce.
		up_to_crossover[0] = 0; up_to_crossover[1] = 0; // which crossovers we've dealt with
		for (int i = d->map.chr_ends[chr - 1]; i < d->map.chr_ends[chr]; ++i) {
			// loop through every marker for this chromosome
			if (up_to_crossover[0] < num_crossovers[0] &&
					d->map.positions[i].position > p_crossover_where[0][up_to_crossover[0]]) {
				// between last loop and this one we crossed over.
				// invert which and update up_to_crossover;
				which[0] = 1 - which[0];
				up_to_crossover[0] += 1;
			}
			if (up_to_crossover[1] < num_crossovers[1] &&
					d->map.positions[i].position > p_crossover_where[1][up_to_crossover[1]]) {
				which[1] = 1 - which[1];
				up_to_crossover[1] += 1;
			}
			output[2*i] = parent1_genome[2*i + which[0]];
			output[2*i + 1] = parent2_genome[2*i + which[1]];
		}

		if (num_crossovers[0] > 50) {
			free(p_crossover_where[0]);
		}
		if (num_crossovers[1] > 50) {
			free(p_crossover_where[1]);
		}

	}
	return;
}

/** Get the alleles of the outcome of producing a doubled haploid from
 * a gamete from a given parent.
 *
 * One gamete is generated, then doubled. The output will be perfectly
 * homozygous.
 *
 * @see generate_gamete(), on which this is based.
 *
 * @param d pointer to the SimData object that includes genetic map data
 * needed to simulate meiosis and the value of n_markers
 * @param parent_genome a 2x(n_markers) array of characters containing the
 * alleles of the first parent
 * @param output a 2x(n_marker) array of chars which will be overwritten
 * with the offspring genome.
*/
void generate_doubled_haploid(SimData* d, char* parent_genome, char* output) {
	// assumes rand is already seeded
	if (parent_genome == NULL) {
		warning( "Could not make this doubled haploid\n");
		return;
	}

	int num_crossovers, up_to_crossover, which;
	float crossover_where[100];
	float* p_crossover_where;

	// treat each chromosome separately.
	for (int chr = 1; chr <= d->map.n_chr; ++chr) {
		// use Poisson distribution to choose the number of crossovers in this chromosome
		num_crossovers = Rf_rpois(d->map.chr_lengths[chr - 1] / 100);

		// in the rare case where it could be >100, get enough space
		// to be able to store the crossover positions we're about to create
		if (num_crossovers <= 100) {
			p_crossover_where = crossover_where; // point at the start of array
		} else {
			p_crossover_where = get_malloc(sizeof(float) * num_crossovers);
		}

		// TASK 3: choose points where those crossovers occur
		// by randomly generating a point along the length of the chromosome
		for (int i = 0; i < num_crossovers; ++i) {
			p_crossover_where[i] = unif_rand()
				* d->map.chr_lengths[chr - 1]
				+ d->map.positions[d->map.chr_ends[chr - 1]].position;
		}

		// sort the crossover points
		if (num_crossovers > 1) {
			qsort(p_crossover_where, num_crossovers, sizeof(float),
					_ascending_float_comparer);
		}

		// pick a parent genome half at random
		which = (unif_rand() > 0.5); // if this is 0, we start with the left.

		// TASK 4: Figure out the gamete that those numbers produce.
		up_to_crossover = 0; // which crossovers we've dealt with
		for (int i = d->map.chr_ends[chr - 1]; i < d->map.chr_ends[chr]; ++i) {
			// loop through every marker for this chromosome
			if (up_to_crossover < num_crossovers &&
					d->map.positions[i].position > p_crossover_where[up_to_crossover]) {
				// if we're here then between last loop and this one we crossed over.
				// invert which and update up_to_crossover;
				which = 1 - which;
				up_to_crossover += 1;
			}
			output[2*i] = parent_genome[2*i + which];
			output[2*i + 1] = parent_genome[2*i + which];
		}

		if (num_crossovers > 100) {
			free(p_crossover_where);
		}

	}
	return;
}

/** Performs random crosses among members of a group. If the group does not
 * have at least two members, the simulation exits. Selfing/crossing an individual
 * with itself is not permitted. The resulting genotypes are allocated to a new group.
 *
 * Preferences in GenOptions are applied to this cross. The family_size parameter
 * in GenOptions allows you to repeat each particular randomly-chosen cross a
 * certain number of times.
 *
 * Parents are drawn uniformly from the group when picking which crosses to make.
 *
 * @param d pointer to the SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param from_group group number from which to draw the parents.
 * @param n_crosses number of random pairs of parents to cross.
 * @param g options for the genotypes created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
int cross_random_individuals(SimData* d, int from_group, int n_crosses, GenOptions g) {
	int g_size = get_group_size( d, from_group);
	if (g_size < 2) {
		if (g_size == 1) {
			warning("Group must contain multiple individuals to be able to perform random crossing\n");
		} else {
			warning("Group %d does not exist.\n", from_group);
		}
		return 0;
	}
	char** group_genes = get_group_genes( d, from_group, g_size);

	// create the buffer we'll use to save the output crosses before they're printed.
	AlleleMatrix* crosses;
	int n_combinations = n_crosses * g.family_size;
	int n_to_go = n_combinations;
	if (n_to_go < CONTIG_WIDTH) {
		crosses = create_empty_allelematrix(d->n_markers, n_to_go);
		n_to_go = 0;
	} else {
		crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
		n_to_go -= CONTIG_WIDTH;
	}
	int fullness = 0;
	int parent1;
	int parent2;

	// set up pedigree/id allocation, if applicable
	unsigned int cid = 0;
	unsigned int* cross_current_id;
	if (g.will_allocate_ids) {
		cross_current_id = &(d->current_id);
	} else {
		cross_current_id = &cid;
	}
	unsigned int* group_ids = NULL;
	if (g.will_track_pedigree) {
		group_ids = get_group_ids( d, from_group, g_size);
	}
	AlleleMatrix* last = NULL;
	int output_group = 0;
	if (g.will_save_to_simdata) {
		last = d->m; // for saving to simdata
		while (last->next != NULL) {
			last = last->next;
		}
		output_group = get_new_group_num( d);
	}

	// open the output files, if applicable
	char fname[NAME_LENGTH];
	FILE* fp = NULL, * fe = NULL, * fg = NULL;
	DecimalMatrix eff;
	if (g.will_save_pedigree_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-pedigree.txt");
		fp = fopen(fname, "w");
	}
	if (g.will_save_bvs_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-bv.txt");
		fe = fopen(fname, "w");
	}
	if (g.will_save_alleles_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-genotype.txt");
		fg = fopen(fname, "w");
	}

	GetRNGstate();
	// loop through each combination
	for (int i = 0; i < n_crosses; ++i) {
		// get parents, randomly.
		parent1 = round(unif_rand() * (g_size - 1));
		do {
			parent2 = round(unif_rand() * (g_size - 1));
		} while (parent1 == parent2);

		for (int f = 0; f < g.family_size; ++f, ++fullness) {
			R_CheckUserInterrupt();

			// when cross buffer is full, save these outcomes to the file.
			if (fullness >= CONTIG_WIDTH) {
				crosses->n_genotypes = CONTIG_WIDTH;
				// give the offspring their ids and names
				if (g.will_name_offspring) {
					set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
				}
				for (int j = 0; j < CONTIG_WIDTH; ++j) {
					++ *cross_current_id;
					crosses->ids[j] = *cross_current_id;
				}

				// save the offspring to files if appropriate
				if (g.will_save_pedigree_to_file) {
					save_AM_pedigree( fp, crosses, d);
				}
				if (g.will_save_bvs_to_file) {
					eff = calculate_bvs( crosses, &(d->e));
					save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
					delete_dmatrix( &eff);
				}
				if (g.will_save_alleles_to_file) {
					save_allele_matrix( fg, crosses, d->markers);
				}

				if (g.will_save_to_simdata) {
					last->next = crosses;
					last = last->next;
					if (n_to_go < CONTIG_WIDTH) {
						crosses = create_empty_allelematrix(d->n_markers, n_to_go);
						n_to_go = 0;
					} else {
						crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
						n_to_go -= CONTIG_WIDTH;
					}
				}

				fullness = 0; //reset the count and start refilling the matrix
			}
			// do the cross.
			generate_cross( d, group_genes[parent1] , group_genes[parent2] , crosses->alleles[fullness] );
			crosses->groups[fullness] = output_group;
			if (g.will_track_pedigree) {
				crosses->pedigrees[0][fullness] = group_ids[parent1];
				crosses->pedigrees[1][fullness] = group_ids[parent2];
			}
		}
	}
	PutRNGstate();

	// save the rest of the crosses to the file.
	free(group_genes);
	// give the offsprings their ids and names
	if (g.will_name_offspring) {
		set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
	}
	for (int j = 0; j < fullness; ++j) {
		++ *cross_current_id;
		crosses->ids[j] = *cross_current_id;
	}
	if (g.will_track_pedigree) {
		free(group_ids);
	}

	// save the offsprings to files if appropriate
	if (g.will_save_pedigree_to_file) {
		save_AM_pedigree( fp, crosses, d);
		fclose(fp);
	}
	if (g.will_save_bvs_to_file) {
		eff = calculate_bvs( crosses, &(d->e));
		save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
		delete_dmatrix( &eff);
		fclose(fe);
	}
	if (g.will_save_alleles_to_file) {
		save_allele_matrix( fg, crosses, d->markers);
		fclose(fg);
	}
	if (g.will_save_to_simdata) {
		last->next = crosses;
		condense_allele_matrix( d );
		return output_group;

	} else {
		delete_allele_matrix( crosses );
		return 0;
	}
}

/** Performs random crosses where the first parent comes from one group and the second from
 *  another. If each group does not have at least one member, the simulation exits.
 *  The resulting genotypes are allocated to a new group. If the user only wants one
 *  parent to be randomly chosen, the flag set_parent_gp1 or set_parent_gp2 can be
 *  set to a nonzero/truthy value, in which case the corresponding group id will instead
 *  be interpreted as an individual genotype index to which members of the other
 *  group will be crossed.
 *
 * Preferences in GenOptions are applied to this cross. The family_size parameter
 * in GenOptions allows you to repeat each particular randomly-chosen cross a
 * certain number of times.
 *
 * Parents are drawn uniformly from the group when picking which crosses to make.
 *
 * @param d pointer to the SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param group1 group number from which to draw the first parent, unless set_parent_gp1 is
 * truthy, in which case it is the index of the set/guaranteed first parent.
 * @param group2 group number from which to draw the second parent, unless set_parent_gp2 is
 * truthy, in which case it is the index of the set/guaranteed second parent.
 * @param n_crosses number of random pairs of parents to cross.
 * @param set_parent_gp1 If falsy/0, random members of group1 will be crossed, and if
 * truthy, the particular individual of index `group1` will always be the first parent of the cross.
 * @param set_parent_gp2 If falsy/0, random members of group2 will be crossed, and if
 * truthy, the particular individual of index `group2` will always be the first parent of the cross.
 * @param g options for the genotypes created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
int cross_randomly_between(SimData*d, int group1, int group2, int n_crosses, int set_parent_gp1, int set_parent_gp2, GenOptions g) {
    char* parent1_genes; char* parent2_genes;
    char** group1_genes; char** group2_genes;
    int group1_size; int group2_size;
    int parent1; int parent2;

    if (set_parent_gp1) {
        parent1_genes = get_genes_of_index( d->m, group1 );
    } else {
        group1_size = get_group_size( d, group1 );
        if (group1_size < 1) {
            warning("Group %d does not exist.\n", group1);
            return 0;
        }
        group1_genes = get_group_genes( d, group1, group1_size );
    }
    if (set_parent_gp2) {
        parent2_genes = get_genes_of_index( d->m, group2 );
    } else {
        group2_size = get_group_size( d, group2 );
        if (group2_size < 1) {
            warning("Group %d does not exist.\n", group2);
            return 0;
        }
        group2_genes = get_group_genes( d, group2, group2_size );
    }

    // create the buffer we'll use to save the output crosses before they're printed.
    AlleleMatrix* crosses;
    int n_combinations = n_crosses * g.family_size;
    int n_to_go = n_combinations;
    if (n_to_go < CONTIG_WIDTH) {
        crosses = create_empty_allelematrix(d->n_markers, n_to_go);
        n_to_go = 0;
    } else {
        crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
        n_to_go -= CONTIG_WIDTH;
    }
    int fullness = 0;

    // set up pedigree/id allocation, if applicable
    unsigned int cid = 0;
    unsigned int* cross_current_id;
    if (g.will_allocate_ids) {
        cross_current_id = &(d->current_id);
    } else {
        cross_current_id = &cid;
    }
    unsigned int parent1_id; unsigned int parent2_id;
    unsigned int* group1_ids = NULL; unsigned int* group2_ids = NULL;
    if (g.will_track_pedigree) {
        if (set_parent_gp1) {
            parent1_id = get_id_of_index( d->m, group1 );
        } else {
            group1_ids = get_group_ids( d, group1, group1_size );
        }
        if (set_parent_gp2) {
            parent2_id = get_id_of_index( d->m, group2 );
        } else {
            group2_ids = get_group_ids( d, group2, group2_size );
        }
    }
    AlleleMatrix* last = NULL;
    int output_group = 0;
    if (g.will_save_to_simdata) {
        last = d->m; // for saving to simdata
        while (last->next != NULL) {
            last = last->next;
        }
        output_group = get_new_group_num( d);
    }

    // open the output files, if applicable
    char fname[NAME_LENGTH];
    FILE* fp = NULL, * fe = NULL, * fg = NULL;
    DecimalMatrix eff;
    if (g.will_save_pedigree_to_file) {
        strcpy(fname, g.filename_prefix);
        strcat(fname, "-pedigree.txt");
        fp = fopen(fname, "w");
    }
    if (g.will_save_bvs_to_file) {
        strcpy(fname, g.filename_prefix);
        strcat(fname, "-bv.txt");
        fe = fopen(fname, "w");
    }
    if (g.will_save_alleles_to_file) {
        strcpy(fname, g.filename_prefix);
        strcat(fname, "-genotype.txt");
        fg = fopen(fname, "w");
    }

    GetRNGstate();
    // loop through each combination
    for (int i = 0; i < n_crosses; ++i) {
        // get parents, randomly.
        if (!set_parent_gp1) {
            parent1 = round(unif_rand() * (group1_size - 1));
            parent1_genes = group1_genes[parent1];
            parent1_id = group1_ids[parent1];
        }
        if (!set_parent_gp2) {
            parent2 = round(unif_rand() * (group2_size - 1));
            parent2_genes = group2_genes[parent2];
            parent2_id = group2_ids[parent2];
        }

        for (int f = 0; f < g.family_size; ++f, ++fullness) {
            R_CheckUserInterrupt();

            // when cross buffer is full, save these outcomes to the file.
            if (fullness >= CONTIG_WIDTH) {
                crosses->n_genotypes = CONTIG_WIDTH;
                // give the offspring their ids and names
                if (g.will_name_offspring) {
                    set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
                }
                for (int j = 0; j < CONTIG_WIDTH; ++j) {
                    ++ *cross_current_id;
                    crosses->ids[j] = *cross_current_id;
                }

                // save the offspring to files if appropriate
                if (g.will_save_pedigree_to_file) {
                    save_AM_pedigree( fp, crosses, d);
                }
                if (g.will_save_bvs_to_file) {
                    eff = calculate_bvs( crosses, &(d->e));
                    save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
                    delete_dmatrix( &eff);
                }
                if (g.will_save_alleles_to_file) {
                    save_allele_matrix( fg, crosses, d->markers);
                }

                if (g.will_save_to_simdata) {
                    last->next = crosses;
                    last = last->next;
                    if (n_to_go < CONTIG_WIDTH) {
                        crosses = create_empty_allelematrix(d->n_markers, n_to_go);
                        n_to_go = 0;
                    } else {
                        crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
                        n_to_go -= CONTIG_WIDTH;
                    }
                }

                fullness = 0; //reset the count and start refilling the matrix
            }
            // do the cross.
            generate_cross( d, parent1_genes , parent2_genes , crosses->alleles[fullness] );
            crosses->groups[fullness] = output_group;
            if (g.will_track_pedigree) {
                crosses->pedigrees[0][fullness] = parent1_id;
                crosses->pedigrees[1][fullness] = parent2_id;
            }
        }
    }
    PutRNGstate();

    // save the rest of the crosses to the file.
    if (!set_parent_gp1) { free(group1_genes); }
    if (!set_parent_gp2) { free(group2_genes); }
    // give the offsprings their ids and names
    if (g.will_name_offspring) {
        set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
    }
    for (int j = 0; j < fullness; ++j) {
        ++ *cross_current_id;
        crosses->ids[j] = *cross_current_id;
    }
    if (g.will_track_pedigree) {
        if (!set_parent_gp1) { free(group1_ids); }
        if (!set_parent_gp2) { free(group2_ids); }
    }

    // save the offsprings to files if appropriate
    if (g.will_save_pedigree_to_file) {
        save_AM_pedigree( fp, crosses, d);
        fclose(fp);
    }
    if (g.will_save_bvs_to_file) {
        eff = calculate_bvs( crosses, &(d->e));
        save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
        delete_dmatrix( &eff);
        fclose(fe);
    }
    if (g.will_save_alleles_to_file) {
        save_allele_matrix( fg, crosses, d->markers);
        fclose(fg);
    }
    if (g.will_save_to_simdata) {
        last->next = crosses;
        condense_allele_matrix( d );
        return output_group;

    } else {
        delete_allele_matrix( crosses );
        return 0;
    }
}

/** Performs the crosses of pairs of parents whose ids are provided in an array. The
 * resulting genotypes are allocated to a new group.
 *
 * Preferences in GenOptions are applied to this cross. The family_size parameter
 * in GenOptions allows you to repeat each particular cross a
 * certain number of times.
 *
 * @param d pointer to the SimData object that includes genetic map data and
 * allele data needed to simulate crossing.
 * @param combinations a 2D array of IDs, with the first dimension being parent 1 or 2,
 * and the second being the IDs of those parents for each combination to cross.
 * @param n_combinations the number of pairs of ids to cross/the length of `combinations`
 * @param g options for the genotypes created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
int cross_these_combinations(SimData* d, int n_combinations, int combinations[2][n_combinations],  GenOptions g) {
	if (n_combinations < 1) {
		return 0;
	}

	// create the buffer we'll use to save the output crosses before they're printed.
	AlleleMatrix* crosses;
	int n_to_go = n_combinations * g.family_size;
	if (n_to_go < CONTIG_WIDTH) {
		crosses = create_empty_allelematrix(d->n_markers, n_to_go);
		n_to_go = 0;
	} else {
		crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
		n_to_go -= CONTIG_WIDTH;
	}
	int fullness = 0, parent1id, parent2id;
	char* parent1genes, * parent2genes;

	// set up pedigree/id allocation, if applicable
	unsigned int cid = 0;
	unsigned int* cross_current_id;
	if (g.will_allocate_ids) {
		cross_current_id = &(d->current_id);
	} else {
		cross_current_id = &cid;
	}
	AlleleMatrix* last = NULL;
	int output_group = 0;
	if (g.will_save_to_simdata) {
		last = d->m; // for saving to simdata
		while (last->next != NULL) {
			last = last->next;
		}
		output_group = get_new_group_num( d);
	}

	// open the output files, if applicable
	char fname[NAME_LENGTH];
	FILE* fp = NULL, * fe = NULL, * fg = NULL;
	DecimalMatrix eff;
	if (g.will_save_pedigree_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-pedigree.txt");
		fp = fopen(fname, "w");
	}
	if (g.will_save_bvs_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-bv.txt");
		fe = fopen(fname, "w");
	}
	if (g.will_save_alleles_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-genotype.txt");
		fg = fopen(fname, "w");
	}

	GetRNGstate();
	// loop through each combination
	for (int i = 0; i < n_combinations; ++i) {
		// find the parents & do the cross. First ensure parents are valid
		if (combinations[0][i] >= 0 && combinations[1][i] >= 0) {
			parent1genes = get_genes_of_index(d->m, combinations[0][i]);
			parent2genes = get_genes_of_index(d->m, combinations[1][i]);
			if (g.will_track_pedigree) {
				parent1id = get_id_of_index(d->m, combinations[0][i]);
				parent2id = get_id_of_index(d->m, combinations[1][i]);
			}

			for (int f = 0; f < g.family_size; ++f, ++fullness) {
				R_CheckUserInterrupt();

				// when cross buffer is full, save these outcomes to the file.
				if (fullness >= CONTIG_WIDTH) {
					// give the offsprings their ids and names
					if (g.will_name_offspring) {
						set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
					}
					for (int j = 0; j < CONTIG_WIDTH; ++j) {
						++ *cross_current_id;
						crosses->ids[j] = *cross_current_id;
					}

					// save the offsprings to files if appropriate
					if (g.will_save_pedigree_to_file) {
						save_AM_pedigree( fp, crosses, d);
					}
					if (g.will_save_bvs_to_file) {
						eff = calculate_bvs( crosses, &(d->e));
						save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
						delete_dmatrix( &eff);
					}
					if (g.will_save_alleles_to_file) {
						save_allele_matrix( fg, crosses, d->markers);
					}

					if (g.will_save_to_simdata) {
						last->next = crosses;
						last = last->next;
						// get the new crosses matrix, of the right size.
						if (n_to_go < CONTIG_WIDTH) {
							crosses = create_empty_allelematrix(d->n_markers, n_to_go);
							n_to_go = 0;
						} else {
							crosses = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
							n_to_go -= CONTIG_WIDTH;
						}
					}
					fullness = 0; //reset the count and start refilling the matrix
				}

				generate_cross(d, parent1genes, parent2genes, crosses->alleles[fullness]);
				crosses->groups[fullness] = output_group;
				if (g.will_track_pedigree) {
					crosses->pedigrees[0][fullness] = parent1id;
					crosses->pedigrees[1][fullness] = parent2id;
				}
			}
		}
	}
	PutRNGstate();

	// save the rest of the crosses to the file.
	// give the offsprings their ids and names
	if (g.will_name_offspring) {
		set_names(crosses, g.offspring_name_prefix, *cross_current_id, 0);
	}
	for (int j = 0; j < fullness; ++j) {
		++ *cross_current_id;
		crosses->ids[j] = *cross_current_id;
	}

	// save the offsprings to files if appropriate
	if (g.will_save_pedigree_to_file) {
		save_AM_pedigree( fp, crosses, d);
		fclose(fp);
	}
	if (g.will_save_bvs_to_file) {
		eff = calculate_bvs( crosses, &(d->e));
		save_manual_bvs( fe, &eff, crosses->ids, crosses->names);
		delete_dmatrix( &eff);
		fclose(fe);
	}
	if (g.will_save_alleles_to_file) {
		save_allele_matrix( fg, crosses, d->markers);
		fclose(fg);
	}
	if (g.will_save_to_simdata) {
		last->next = crosses;
		condense_allele_matrix( d );
		return output_group;

	} else {
		delete_allele_matrix( crosses );
		return 0;
	}
}

/** Selfs each member of a group for a certain number of generations.
 * The resulting genotypes are allocated to a new group.
 *
 * Only the genotype after all n generations is saved. Intermediate steps will be lost.
 *
 * Preferences in GenOptions are applied to this operation. The family_size parameter
 * in GenOptions allows you to generate multiple selfed offspring from each member
 * of the group. These multiple selfed offspring all originate from independent
 * processes of selfing the parent with itself then selfing that intermediate offspring
 * with itself for the remaining n-1 generations.
 *
 * @param d pointer to the SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param n number of generations of selfing simulation to carry out.
 * @param group the genotypes on which to perform these n generations of selfing.
 * @param g options for the genotypes created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
int self_n_times(SimData* d, int n, int group, GenOptions g) {
	int group_size = get_group_size( d, group);
	if (group_size < 1) {
		warning("Group %d does not exist.\n", group);
		return 0;
	}
	if (n < 1) {
		error("That number of generations cannot be produced.\n");
	}

	AlleleMatrix* outcome;
	int n_to_go = group_size * g.family_size;
	if (n_to_go < CONTIG_WIDTH) {
		outcome = create_empty_allelematrix(d->n_markers, n_to_go);
		n_to_go = 0;
	} else {
		outcome = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
		n_to_go -= CONTIG_WIDTH;
	}
	char** group_genes = get_group_genes( d, group, group_size);
	int i, j, f, fullness = 0;

	// set up pedigree/id allocation, if applicable
	unsigned int cid = 0;
	unsigned int* cross_current_id;
	if (g.will_allocate_ids) {
		cross_current_id = &(d->current_id);
	} else {
		cross_current_id = &cid;
	}
	unsigned int* group_ids = NULL;
	if (g.will_track_pedigree) {
		group_ids = get_group_ids( d, group, group_size);
	}
	AlleleMatrix* last = NULL;
	int output_group = 0;
	if (g.will_save_to_simdata) {
		last = d->m; // for saving to simdata
		while (last->next != NULL) {
			last = last->next;
		}
		output_group = get_new_group_num( d);
	}

	// open the output files, if applicable
	char fname[NAME_LENGTH];
	FILE* fp = NULL, * fe = NULL, * fg = NULL;
	DecimalMatrix eff;
	if (g.will_save_pedigree_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-pedigree.txt");
		fp = fopen(fname, "w");
	}
	if (g.will_save_bvs_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-bv.txt");
		fe = fopen(fname, "w");
	}
	if (g.will_save_alleles_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-genotype.txt");
		fg = fopen(fname, "w");
	}

	GetRNGstate();
	for (i = 0; i < group_size; ++i) {
		// do n rounds of selfing (j-indexed loops) g.family_size times per individual (f-indexed loops)
		if (n == 1)  {
			//find the parent genes, save a shallow copy to set
			char* genes = group_genes[i];
			//int id = group_ids[i];
			for (f = 0; f < g.family_size; ++f, ++fullness) {
				R_CheckUserInterrupt();

				// when cross buffer is full, save these outcomes to the file.
				if (fullness >= CONTIG_WIDTH) {
					outcome->n_genotypes = CONTIG_WIDTH;
					// give the offspring their ids and names
					if (g.will_name_offspring) {
						set_names(outcome, g.offspring_name_prefix, *cross_current_id, 0);
					}
					for (int j = 0; j < CONTIG_WIDTH; ++j) {
						++ *cross_current_id;
						outcome->ids[j] = *cross_current_id;
					}

					// save the offspring to files if appropriate
					if (g.will_save_pedigree_to_file) {
						save_AM_pedigree( fp, outcome, d);
					}
					if (g.will_save_bvs_to_file) {
						eff = calculate_bvs( outcome, &(d->e));
						save_manual_bvs( fe, &eff, outcome->ids, outcome->names);
						delete_dmatrix( &eff);
					}
					if (g.will_save_alleles_to_file) {
						save_allele_matrix( fg, outcome, d->markers);
					}

					if (g.will_save_to_simdata) {
						last->next = outcome;
						last = last->next;
						if (n_to_go < CONTIG_WIDTH) {
							outcome = create_empty_allelematrix(d->n_markers, n_to_go);
							n_to_go = 0;
						} else {
							outcome = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
							n_to_go -= CONTIG_WIDTH;
						}
					}
					fullness = 0; //reset the count and start refilling the matrix
				}

				generate_cross( d, genes, genes, outcome->alleles[fullness] );
				outcome->groups[fullness] = output_group;
				/*outcome->groups[fullness] = output_group;
				if (g.will_track_pedigree) {
					outcome->pedigrees[0][fullness] = id;
					outcome->pedigrees[1][fullness] = id;
				}*/
			}

		} else {
			//find the parent genes, save a deep copy to set
			char* genes = get_malloc(sizeof(char) * (d->n_markers<<1));
			for (j = 0; j < d->n_markers; ++j) {
				genes[2*j] = group_genes[i][2*j];
				genes[2*j + 1] = group_genes[i][2*j + 1];
			}

			for (f = 0; f < g.family_size; ++f, ++fullness) {
				R_CheckUserInterrupt();

				// when cross buffer is full, save these outcomes to the file.
				if (fullness >= CONTIG_WIDTH) {
					outcome->n_genotypes = CONTIG_WIDTH;
					// give the offspring their ids and names
					if (g.will_name_offspring) {
						set_names(outcome, g.offspring_name_prefix, *cross_current_id, 0);
					}
					for (int j = 0; j < CONTIG_WIDTH; ++j) {
						++ *cross_current_id;
						outcome->ids[j] = *cross_current_id;
					}

					// save the offspring to files if appropriate
					if (g.will_save_pedigree_to_file) {
						save_AM_pedigree( fp, outcome, d);
					}
					if (g.will_save_bvs_to_file) {
						eff = calculate_bvs( outcome, &(d->e));
						save_manual_bvs( fe, &eff, outcome->ids, outcome->names);
						delete_dmatrix( &eff);
					}
					if (g.will_save_alleles_to_file) {
						save_allele_matrix( fg, outcome, d->markers);
					}

					if (g.will_save_to_simdata) {
						last->next = outcome;
						last = last->next;
						if (n_to_go < CONTIG_WIDTH) {
							outcome = create_empty_allelematrix(d->n_markers, n_to_go);
							n_to_go = 0;
						} else {
							outcome = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
							n_to_go -= CONTIG_WIDTH;
						}
					}
					fullness = 0; //reset the count and start refilling the matrix
				}

				for (j = 0; j < n; ++j) {
					R_CheckUserInterrupt();
					if (j % 2) {
						generate_cross( d, outcome->alleles[fullness], outcome->alleles[fullness], genes);
					} else {
						generate_cross( d, genes, genes, outcome->alleles[fullness] );
					}
				}
				if (n % 2) {
					free(outcome->alleles[fullness]);
					outcome->alleles[fullness] = genes;

					// make a new deep copy if we're still going.
					if (f + 1 < g.family_size) {
						genes = get_malloc(sizeof(char) * (d->n_markers<<1));
						for (j = 0; j < d->n_markers; ++j) {
							genes[2*j] = group_genes[i][2*j];
							genes[2*j + 1] = group_genes[i][2*j + 1];
						}
					}
				} else {
					free(genes);
				}

				outcome->groups[fullness] = output_group;

			}
		}

		if (g.will_track_pedigree) {
			fullness -= g.family_size;
			int id = group_ids[i];
			for (f = 0; f < g.family_size; ++f, ++fullness) {
				outcome->pedigrees[0][fullness] = id;
				outcome->pedigrees[1][fullness] = id;
			}
		}
	}
	PutRNGstate();

	free(group_genes);
	// save the rest of the crosses to the file.
	outcome->n_genotypes = fullness;
	// give the offspring their ids and names
	if (g.will_name_offspring) {
		set_names(outcome, g.offspring_name_prefix, *cross_current_id, 0);
	}
	for (int j = 0; j < fullness; ++j) {
		++ *cross_current_id;
		outcome->ids[j] = *cross_current_id;
	}
	if (g.will_track_pedigree) {
		free(group_ids);
	}

	// save the offspring to files if appropriate
	if (g.will_save_pedigree_to_file) {
		save_AM_pedigree( fp, outcome, d);
		fclose(fp);
	}
	if (g.will_save_bvs_to_file) {
		eff = calculate_bvs( outcome, &(d->e));
		save_manual_bvs( fe, &eff, outcome->ids, outcome->names);
		delete_dmatrix( &eff);
		fclose(fe);
	}
	if (g.will_save_alleles_to_file) {
		save_allele_matrix( fg, outcome, d->markers);
		fclose(fg);
	}
	if (g.will_save_to_simdata) {
		last->next = outcome;
		condense_allele_matrix( d );
		return output_group;
	} else {
		delete_allele_matrix( outcome );
		return 0;
	}
}

/** Creates a doubled haploid from each member of a group.
 * The resulting genotypes are allocated to a new group.
 *
 * Preferences in GenOptions are applied to this operation. The family_size parameter
 * in GenOptions allows you to generate multiple doubled haploid offspring from each member
 * of the group. These multiple offspring all originate from independent
 * processes of generating a gamete then doubling its alleles.
 *
 * @param d pointer to the SimData object that contains the genetic map and
 * genotypes of the parent group.
 * @param group the genotypes on which to perform the operation.
 * @param g options for the genotypes created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
*/
int make_doubled_haploids(SimData* d, int group, GenOptions g) {
	int group_size = get_group_size( d, group);
	if (group_size < 1) {
		warning("Group %d does not exist.\n", group);
		return 0;
	}

	AlleleMatrix* outcome;
	int n_to_go = group_size * g.family_size;
	if (n_to_go < CONTIG_WIDTH) {
		outcome = create_empty_allelematrix(d->n_markers, n_to_go);
		n_to_go = 0;
	} else {
		outcome = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
		n_to_go -= CONTIG_WIDTH;
	}
	char** group_genes = get_group_genes( d, group, group_size);
	int i, f, fullness = 0;

	// set up pedigree/id allocation, if applicable
	unsigned int cid = 0;
	unsigned int* cross_current_id;
	if (g.will_allocate_ids) {
		cross_current_id = &(d->current_id);
	} else {
		cross_current_id = &cid;
	}
	unsigned int* group_ids = NULL;
	if (g.will_track_pedigree) {
		group_ids = get_group_ids( d, group, group_size);
	}
	AlleleMatrix* last = NULL;
	int output_group = 0;
	if (g.will_save_to_simdata) {
		last = d->m; // for saving to simdata
		while (last->next != NULL) {
			last = last->next;
		}
		output_group = get_new_group_num( d);
	}

	// open the output files, if applicable
	char fname[NAME_LENGTH];
	FILE* fp = NULL, * fe = NULL, * fg = NULL;
	DecimalMatrix eff;
	if (g.will_save_pedigree_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-pedigree.txt");
		fp = fopen(fname, "w");
	}
	if (g.will_save_bvs_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-bv.txt");
		fe = fopen(fname, "w");
	}
	if (g.will_save_alleles_to_file) {
		strcpy(fname, g.filename_prefix);
		strcat(fname, "-genotype.txt");
		fg = fopen(fname, "w");
	}

	GetRNGstate();
	for (i = 0; i < group_size; ++i) {
		// do n rounds of selfing (j-indexed loops) g.family_size times per individual (f-indexed loops)
		//find the parent genes, save a shallow copy to set
		char* genes = group_genes[i];
		int id = group_ids[i];
		for (f = 0; f < g.family_size; ++f, ++fullness) {
			R_CheckUserInterrupt();

			// when cross buffer is full, save these outcomes to the file.
			if (fullness >= CONTIG_WIDTH) {
				outcome->n_genotypes = CONTIG_WIDTH;
				// give the offspring their ids and names
				if (g.will_name_offspring) {
					set_names(outcome, g.offspring_name_prefix, *cross_current_id, 0);
				}
				for (int j = 0; j < CONTIG_WIDTH; ++j) {
					++ *cross_current_id;
					outcome->ids[j] = *cross_current_id;
				}

				// save the offspring to files if appropriate
				if (g.will_save_pedigree_to_file) {
					save_AM_pedigree( fp, outcome, d);
				}
				if (g.will_save_bvs_to_file) {
					eff = calculate_bvs( outcome, &(d->e));
					save_manual_bvs( fe, &eff, outcome->ids, outcome->names);
					delete_dmatrix( &eff);
				}
				if (g.will_save_alleles_to_file) {
					save_allele_matrix( fg, outcome, d->markers);
				}

				if (g.will_save_to_simdata) {
					last->next = outcome;
					last = last->next;
					if (n_to_go < CONTIG_WIDTH) {
						outcome = create_empty_allelematrix(d->n_markers, n_to_go);
						n_to_go = 0;
					} else {
						outcome = create_empty_allelematrix(d->n_markers, CONTIG_WIDTH);
						n_to_go -= CONTIG_WIDTH;
					}
				}
				fullness = 0; //reset the count and start refilling the matrix
			}

			generate_doubled_haploid( d, genes, outcome->alleles[fullness] );
			outcome->groups[fullness] = output_group;
			if (g.will_track_pedigree) {
				outcome->pedigrees[0][fullness] = id;
				outcome->pedigrees[1][fullness] = id;
			}
		}
	}
	PutRNGstate();

	free(group_genes);
	// save the rest of the crosses to the file.
	outcome->n_genotypes = fullness;
	// give the offspring their ids and names
	if (g.will_name_offspring) {
		set_names(outcome, g.offspring_name_prefix, *cross_current_id, 0);
	}
	for (int j = 0; j < fullness; ++j) {
		++ *cross_current_id;
		outcome->ids[j] = *cross_current_id;
	}
	if (g.will_track_pedigree) {
		free(group_ids);
	}

	// save the offspring to files if appropriate
	if (g.will_save_pedigree_to_file) {
		save_AM_pedigree( fp, outcome, d);
		fclose(fp);
	}
	if (g.will_save_bvs_to_file) {
		eff = calculate_bvs( outcome, &(d->e));
		save_manual_bvs( fe, &eff, outcome->ids, outcome->names);
		delete_dmatrix( &eff);
		fclose(fe);
	}
	if (g.will_save_alleles_to_file) {
		save_allele_matrix( fg, outcome, d->markers);
		fclose(fg);
	}
	if (g.will_save_to_simdata) {
		last->next = outcome;
		condense_allele_matrix( d );
		return output_group;
	} else {
		delete_allele_matrix( outcome );
		return 0;
	}
}

/** Perform crosses between all pairs of parents
 * in the group `from_group` and allocates the resulting offspring to a new group.
 *
 * If the group has n members, there will be n * (n-1) / 2 offspring produced.
 *
 * Preferences in GenOptions are applied to this cross.
 *
 * @param d pointer to the SimData object containing or markers and parent alleles
 * @param from_group group number from which to do all these crosses.
 * @param g options for the AlleleMatrix created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
int make_all_unidirectional_crosses(SimData* d, int from_group, GenOptions g) {
	int group_size = get_group_size( d, from_group );
	if (group_size < 2) {
		if (group_size == 1) {
			warning("Group %d does not have enough members to perform crosses\n", from_group);
		} else {
			warning("Group %d does not exist.\n", from_group);
		}
		return 0;
	}
	//unsigned int* group_ids = get_group_ids( d, from_group, group_size);
    int* group_indexes = get_group_indexes( d, from_group, group_size);

	// number of crosses = number of entries in upper triangle of matrix
	//    = half of (n entries in matrix - length of diagonal)
	// 	  = half of (lmatrix * lmatrix - lmatrix);
	int n_crosses = group_size * (group_size - 1) / 2; //* g.family_size;

	int combinations[2][n_crosses];
	int cross_index = 0;
	for (int i = 0; i < group_size; ++i) {
		for (int j = i + 1; j < group_size; ++j) {
			combinations[0][cross_index] = group_indexes[i];
			combinations[1][cross_index] = group_indexes[j];

			++cross_index;
		}
	}

	return cross_these_combinations(d, n_crosses, combinations, g);

}

/** Find the top m percent of a group and perform random crosses between those
 * top individuals. The resulting offspring are allocated to a new group.
 *
 * Preferences in GenOptions are applied to this cross.
 *
 * @param d pointer to the SimData object containing or markers and parent alleles
 * @param n number of random crosses to make.
 * @param m percentage (as a decimal, i.e. 0.2 for 20percent). Take the best [this portion]
 * of the group for performing the random crosses.
 * @param group group number from which to identify top parents and do crosses
 * @param g options for the AlleleMatrix created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
int make_n_crosses_from_top_m_percent(SimData* d, int n, int m, int group, GenOptions g) {
	// move the top m% to a new group
	int group_size = get_group_size(d, group);
	int n_top_group = group_size * m / 100; //@integer division?
	Rprintf("There are %d lines in the top %d%%\n", n_top_group, m);

	int topgroup = split_by_bv(d, group, n_top_group, FALSE);

	// do the random crosses
	int gp = cross_random_individuals(d, topgroup, n, g);

	// unconvert from a group
	int to_combine[] = {group, topgroup};
	combine_groups(d, 2, to_combine);

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
 * Preferences in GenOptions are applied to this cross.
 *
 * @param d pointer to the SimData object containing or markers and parent alleles
 * @param input_file file instructing which crosses to perform
 * @param g options for the AlleleMatrix created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
int make_crosses_from_file(SimData* d, const char* input_file, GenOptions g) {
	struct TableSize t = get_file_dimensions(input_file, '\t');
	if (t.num_rows < 1) {
		warning( "No crosses exist in that file\n");
		return 0;
	}

	//open file
	FILE* fp;
	if ((fp = fopen(input_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", input_file);
	}

	int combinations[2][t.num_rows];
	char buffer[2][50];
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s \n", buffer[0], buffer[1]);
		combinations[0][i] = get_index_of_name(d->m, buffer[0]);
		combinations[1][i] = get_index_of_name(d->m, buffer[1]);
	}

	fclose(fp);
	return cross_these_combinations(d, t.num_rows, combinations, g);
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
 * Preferences in GenOptions are applied to this cross.
 *
 * @param d pointer to the SimData object containing or markers and parent alleles
 * @param input_file file instructing which crosses to perform
 * @param g options for the AlleleMatrix created. @see GenOptions
 * @returns the group number of the group to which the produced offspring were allocated.
 */
int make_double_crosses_from_file(SimData* d, const char* input_file, GenOptions g) {
	struct TableSize t = get_file_dimensions(input_file, '\t');
	if (t.num_rows < 1) {
		warning( "No crosses exist in that file\n");
		return 0;
	}

	//open file
	FILE* fp;
	if ((fp = fopen(input_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", input_file);
	}

	int combinations[2][t.num_rows];
	char buffer[4][50];
	char* to_buffer[] = {buffer[0], buffer[1], buffer[2], buffer[3]};
	unsigned int g0_id[4];
	int f1_i[2];
	// for each row in file
	for (int i = 0; i < t.num_rows; ++i) {
		// load the four grandparents
		fscanf(fp, "%s %s %s %s \n", buffer[0], buffer[1], buffer[2], buffer[3]);
		get_ids_of_names(d->m, 4, to_buffer, g0_id);
		if (g0_id[0] < 0 || g0_id[1] < 0 || g0_id[2] < 0 || g0_id[3] < 0) {
			warning( "Could not go ahead with the line %d cross - g0 names not in records\n", i);
			combinations[0][i] = -1;
			combinations[1][i] = -1;
			continue;
		}

		// identify two parents
		f1_i[0] = get_index_of_child(d->m, g0_id[0], g0_id[1]);
		f1_i[1] = get_index_of_child(d->m, g0_id[2], g0_id[3]);
		if (f1_i[0] < 0 || f1_i[1] < 0) {
			// try different permutations of the four grandparents.
			f1_i[0] = get_id_of_child(d->m, g0_id[0], g0_id[2]);
			f1_i[1] = get_id_of_child(d->m, g0_id[1], g0_id[3]);
			if (f1_i[0] < 0 || f1_i[1] < 0) {
				f1_i[0] = get_index_of_child(d->m, g0_id[0], g0_id[3]);
				f1_i[1] = get_index_of_child(d->m, g0_id[1], g0_id[2]);
				if (f1_i[0] < 0 || f1_i[1] < 0) {
					warning( "Could not go ahead with the line %d cross - f1 children do not exist for this quartet\n", i);
					combinations[0][i] = -1;
					combinations[1][i] = -1;
					continue;
				}
			}
		}

		//add them to a combinations list
		combinations[0][i] = f1_i[0];
		combinations[1][i] = f1_i[1];

	}

	fclose(fp);
	return cross_these_combinations(d, t.num_rows, combinations, g);

}


/*--------------------------------Fitness------------------------------------*/

/** Takes the `top_n` individuals in the group with the best breeding values/fitnesses
 * and puts them in a new group. The new group number is returned.
 *
 * @param d pointer to the SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param group group number from which to split the top individuals.
 * @param top_n The number of individuals to put in the new group.
 * @param lowIsBest boolean, if TRUE the `top_n` with the lowest breeding value
 * will be selected, if false the `top_n` with the highest breeding value are.
 * @returns the group number of the newly-created split-off group
 */
int split_by_bv(SimData* d, int group, int top_n, int lowIsBest) {	
	unsigned int group_size = get_group_size( d, group );
    int* group_contents = get_group_indexes( d, group, group_size);
	
	if (group_size <= top_n) {
		// well we'll just have to move em all
		int migration = split_from_group(d, group_size, group_contents);
		free(group_contents);
		return migration;
	}
	
	// This should be ordered the same as the indexes
	DecimalMatrix fits = calculate_group_bvs( d, group ); // 1 by group_size matrix
	
	// get an array of pointers to those fitnesses
	double* p_fits[fits.cols];
	for (int i = 0; i < fits.cols; i++) {
		p_fits[i] = &(fits.matrix[0][i]);
	}

	// sort descending
	if (lowIsBest) {
		qsort(p_fits, fits.cols, sizeof(double*), _ascending_double_comparer);
	} else {
		qsort(p_fits, fits.cols, sizeof(double*), _descending_double_comparer);
	}

	// save the indexes of the best n
	int top_individuals[top_n];
	for (int i = 0; i < top_n; i++) {
		top_individuals[i] = group_contents[p_fits[i] - fits.matrix[0]];
	}
	delete_dmatrix(&fits);
	free(group_contents);

	// send those n to a new group
	return split_from_group(d, top_n, top_individuals);
}

/** Calculates the fitness metric/breeding value for each genotype in the AlleleMatrix
* in a certain group, and returns the results in a DecimalMatrix.
*
* The breeding value is calculated for each genotype by taking the sum of (number of
* copies of this allele at this marker times this allele's effect at this marker)
* for each marker for each different allele. To do
* this, the vector containing each allele's effect values is multiplied by a matrix
* containing the counts of that allele at each marker.
*
* The function exits with error code 1 if no marker effect file is loaded.
*
* @param d pointer to the SimData object to which the groups and individuals belong.
* It must have a marker effect file loaded to successfully run this function.
* @param group calculate breeding values for each genotype in the group with this group number.
* @returns A DecimalMatrix containing the score for each individual in the group.
*/
DecimalMatrix calculate_group_bvs(SimData* d, unsigned int group) {
	// check that both of the items to be multiplied exist.
    if (d->e.effects.rows < 1 || d->m == NULL) {
		error( "Either effect matrix or allele matrix does not exist\n");
	}

	int group_size = get_group_size( d, group );
	DecimalMatrix sum = generate_zero_dmatrix(1, group_size);
	DecimalMatrix counts = generate_zero_dmatrix(group_size, d->n_markers);
	DecimalMatrix counts2 = generate_zero_dmatrix(group_size, d->n_markers);

    int i = 0; // highest allele index

	for (; i < d->e.effects.rows - 1; i += 2) {
		// get the allele counts in counts
		calculate_group_doublecount_matrix_of_allele(d, group, d->e.effect_names[i], &counts, d->e.effect_names[i+1], &counts2);

		// multiply counts with effects and add to bv sum
		add_doublematrixvector_product_to_dmatrix(&sum, &counts, d->e.effects.matrix[i], &counts2, d->e.effects.matrix[i+1]);
	}
	if (i < d->e.effects.rows) { // deal with the last odd-numbered allele
		calculate_group_count_matrix_of_allele(d, group, d->e.effect_names[i], &counts);
		add_matrixvector_product_to_dmatrix(&sum, &counts, d->e.effects.matrix[i]);
	}

	delete_dmatrix(&counts);
	delete_dmatrix(&counts2);

	return sum;
}

/** Calculates the fitness metric/breeding value for each genotype in the AlleleMatrix,
 * and returns the results in a DecimalMatrix struct.
 *
 * The breeding value is calculated for each genotype by taking the sum of (number of
 * copies of this allele at this marker times this allele's effect at this marker)
 * for each marker for each different allele. To do
 * this, the vector containing each allele's effect values is multiplied by a matrix
 * containing the counts of that allele at each marker.
 *
 * The function exits with error code 1 if no marker effect file is loaded.
 *
 * @param m pointer to the AlleleMatrix object to which the genotypes belong.
 * @param e pointer to the EffectMatrix that effect values have been loaded into.
 * @returns A DecimalMatrix containing the score for each individual in the group.
 */
DecimalMatrix calculate_bvs( AlleleMatrix* m, EffectMatrix* e) {
	// check that both of the items to be multiplied exist.
    if (e->effects.rows < 1 || m == NULL) {
		error( "Either effect matrix or allele matrix does not exist\n");
	}

	DecimalMatrix sum = generate_zero_dmatrix(1, m->n_genotypes);
	DecimalMatrix counts = generate_zero_dmatrix(m->n_genotypes, m->n_markers);
	DecimalMatrix counts2 = generate_zero_dmatrix(m->n_genotypes, m->n_markers);

    int i = 0; // highest allele index

	for (; i < e->effects.rows - 1; i += 2) {
		// get the allele counts in counts
		calculate_doublecount_matrix_of_allele(m, e->effect_names[i], &counts, e->effect_names[i+1], &counts2);

		// multiply counts with effects and add to bv sum
		add_doublematrixvector_product_to_dmatrix(&sum, &counts, e->effects.matrix[i], &counts2, e->effects.matrix[i+1]);
	}
	if (i < e->effects.rows) { // deal with the last odd-numbered allele
		calculate_count_matrix_of_allele(m, e->effect_names[i], &counts);
		add_matrixvector_product_to_dmatrix(&sum, &counts, e->effects.matrix[i]);
	}

	delete_dmatrix(&counts);
	delete_dmatrix(&counts2);

	return sum;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in a group.
 * Returns the result in a pre-created DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @param d pointer to the SimData.
 * @param group group number for which to calculate counts of the allele
 * @param allele the single-character allele to be counting.
 * @param counts pointer to the DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int calculate_group_count_matrix_of_allele( SimData* d, unsigned int group, char allele, DecimalMatrix* counts) {
	//DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, n_ids);
	int groupSize = get_group_size(d, group);
	if (counts->rows < groupSize || counts->cols < d->n_markers) {
        fprintf(stderr, "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                d->n_markers, groupSize, counts->rows, counts->cols);
        return 1;
	}

	char** genes = get_group_genes(d, group, groupSize);

	for (int i = 0; i < groupSize; ++i) {
		R_CheckUserInterrupt();
		if (genes[i] == NULL) {
			continue;
		}

		for (int j = 0; j < d->n_markers; ++j) {
			int cell_sum = 0;
			if (genes[i][2*j] == allele)     cell_sum += 1;
			if (genes[i][2*j + 1] == allele) cell_sum += 1;
			counts->matrix[i][j] = cell_sum;
		}
	}
	return 0;
}

/** Calculates the number of times at each marker that two particular alleles appear
 * for each genotype in a group.
 * Returns the result in two pre-created DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * This function is the same as calculate_group_count_matrix_of_allele(), but for
 * two alleles at a time, for the purpose of loop unrolling in breeding value
 * calculations.
 *
 * @param d pointer to the SimData.
 * @param group group number for which to calculate counts of the allele
 * @param allele the first single-character allele to be counting.
 * @param counts pointer to the DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @param allele2 the second single-character allele to be counting.
 * @param counts2 pointer to the DecimalMatrix into which to put the number
 * of `allele2` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int calculate_group_doublecount_matrix_of_allele( SimData* d, unsigned int group, char allele, DecimalMatrix* counts, char allele2, DecimalMatrix* counts2) {
	//DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, n_ids);
	int groupSize = get_group_size(d, group);
	if (counts->rows < groupSize || counts->cols < d->n_markers) {
        fprintf(stderr, "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                d->n_markers, groupSize, counts->rows, counts->cols);
        return 1;
	}
    if (counts2->rows < groupSize || counts2->cols < d->n_markers) {
        fprintf(stderr, "`counts2` is the wrong size to be filled: needs %d by %d but is %d by %d\n",
                d->n_markers, groupSize, counts2->rows, counts2->cols);
        return 1;
	}

	char** genes = get_group_genes(d, group, groupSize);

	for (int i = 0; i < groupSize; ++i) {
		R_CheckUserInterrupt();
		if (genes[i] == NULL) {
			continue;
		}

		for (int j = 0; j < d->n_markers; ++j) {
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
	return 0;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in a particular AlleleMatrix (does not go through the linked list).
 * Returns the result in a pre-created DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @param m pointer to the AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the single-character allele to be counting.
 * @param counts pointer to the DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the AlleleMatrix.
 * @returns 0 on success, nonzero on failure.
 * */
int calculate_count_matrix_of_allele( AlleleMatrix* m, char allele, DecimalMatrix* counts) {
	//DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, m->n_genotypes);
	if (counts->rows < m->n_genotypes || counts->cols < m->n_markers) {
        warning( "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n", m->n_markers, m->n_genotypes, counts->rows, counts->cols);
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
 * for each genotype in an AlleleMatrix (does not go through the linked list).
 * Returns the result in two pre-created DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * This function is the same as calculate_count_matrix_of_allele(), but for
 * two alleles at a time, for the purpose of loop unrolling in breeding value
 * calculations.
 *
 * @param m pointer to the AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the first single-character allele to be counting.
 * @param counts pointer to the DecimalMatrix into which to put the number
 * of `allele` occurences at each row/marker for each column/genotype in the group.
 * @param allele2 the second single-character allele to be counting.
 * @param counts2 pointer to the DecimalMatrix into which to put the number
 * of `allele2` occurences at each row/marker for each column/genotype in the group.
 * @returns 0 on success, nonzero on failure.
 * */
int calculate_doublecount_matrix_of_allele( AlleleMatrix* m, char allele, DecimalMatrix* counts, char allele2, DecimalMatrix* counts2) {
	if (counts->rows < m->n_genotypes || counts->cols < m->n_markers) {
        warning( "`counts` is the wrong size to be filled: needs %d by %d but is %d by %d\n", m->n_markers, m->n_genotypes, counts->rows, counts->cols);
        return 1;
	}
    if (counts2->rows < m->n_genotypes || counts2->cols < m->n_markers) {
        warning( "`counts2` is the wrong size to be filled: needs %d by %d but is %d by %d\n", m->n_markers, m->n_genotypes, counts2->rows, counts2->cols);
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
 * for each genotype in the linked list starting at the given AlleleMatrix.
 * Returns the result as a DecimalMatrix. Useful for multiplying to
 * effect matrix to calculate breeding values.
 *
 * @param m pointer to the AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the single-character allele to be counting.
 * @returns DecimalMatrix containing counts for every genotype in the AlleleMatrix
 * linked list.
 * */
DecimalMatrix calculate_full_count_matrix_of_allele( AlleleMatrix* m, char allele) {
	AlleleMatrix* currentm = m;
	int size = 0;
	do {
		size += currentm->n_genotypes;
	} while ((currentm = currentm->next) != NULL);

	DecimalMatrix counts = generate_zero_dmatrix(size, m->n_markers);

	currentm = m;
	int currenti = 0;
	for (int i = 0; i < size; ++i, ++currenti) {
		R_CheckUserInterrupt();
		if (currenti >= currentm->n_genotypes) {
			currenti = 0;
			currentm = currentm->next;
		}

		if (currentm->alleles[i] == NULL) {
			continue;
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
 * Remember to call the MarkerBlocks destructor delete_markerblocks() on the returned
 * struct.
 *
 * @param d pointer to the SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param n number of blocks into which to split each chromosome.
 * @returns a struct containing the markers identified as belonging to each block.
 */
MarkerBlocks create_n_blocks_by_chr(SimData* d, int n) {
	MarkerBlocks blocks;

	// count the number of chromosomes where we have markers to be allocated to blocks
	int chrs_with_contents = 0;
	for (int chr = 0; chr < d->map.n_chr; ++chr) {
		if (d->map.chr_lengths[chr] > 0) {
			++chrs_with_contents;
		}
	}

	blocks.num_blocks = n * chrs_with_contents;
	blocks.num_markers_in_block = get_malloc(sizeof(int) * blocks.num_blocks);
	blocks.markers_in_block = get_malloc(sizeof(int*) * blocks.num_blocks);
	for (int i = 0; i < blocks.num_blocks; ++i) {
		blocks.num_markers_in_block[i] = 0;
		blocks.markers_in_block[i] = NULL;
	}

	int b = 0; //index of the current block in the struct

	for (int chr = 0; chr < d->map.n_chr; ++chr) {
		if (d->map.chr_lengths[chr] <= 0) {
			// chromosome has invalid length, so we have no markers from here
			continue;
		}

		int blocks_this_chr = 1; //counter of how many blocks we have in this chr so far
		float blen = d->map.chr_lengths[chr] / n; //length of each of the n blocks
		float bend = d->map.positions[d->map.chr_ends[chr]].position + blen; // end position of the first block
		int mfirst = d->map.chr_ends[chr]; //index of first marker in the block

		// loop through each marker in this chromosome
		for (int i = d->map.chr_ends[chr]; i < d->map.chr_ends[chr + 1]; ++i) {

			// are we up to the next block yet?
			if (blocks_this_chr < n && d->map.positions[i].position > bend) {
				// save the previous block now.
				blocks.markers_in_block[b] = get_malloc(sizeof(int) * blocks.num_markers_in_block[b]);
				for (int m = 0; m < blocks.num_markers_in_block[b]; ++m) {
					blocks.markers_in_block[b][m] = mfirst + m;
				}

				// start new block
				++blocks_this_chr;
				++b;
				bend += blen;
				// check if there's any empty blocks in between previous one and this one
				while (blocks_this_chr < n && d->map.positions[i].position > bend) {
					++blocks_this_chr;
					++b;
					bend += blen;
				}

				mfirst = i;
			}

			// save marker to block
			++blocks.num_markers_in_block[b];

		}

		// save the last block of this chr
		blocks.markers_in_block[b] = get_malloc(sizeof(int) * blocks.num_markers_in_block[b]);
		for (int m = 0; m < blocks.num_markers_in_block[b]; ++m) {
			blocks.markers_in_block[b][m] = mfirst + m;
		}
		++b;

	}

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
 * @param d pointer to the SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param block_file string containing filename of the file with blocks
 * @returns a struct containing the markers identified as belonging to each block
 * according to their definitions in the file.
 */
MarkerBlocks read_block_file(SimData* d, const char* block_file) {
	struct TableSize ts = get_file_dimensions(block_file, '\t');

	MarkerBlocks blocks;
	blocks.num_blocks = ts.num_rows - 1;
	blocks.num_markers_in_block = get_malloc(sizeof(int) * blocks.num_blocks);
	blocks.markers_in_block = get_malloc(sizeof(int*) * blocks.num_blocks);

	FILE* infile;
	if ((infile = fopen(block_file, "r")) == NULL) {
		error( "Failed to open file %s.\n", block_file);
		//return blocks;
	}

	int bufferlen = d->n_markers;
	char markername[bufferlen];
	int markerbuffer[bufferlen];
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
				int markerindex = get_from_unordered_str_list(markername, d->markers, d->n_markers);
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
		blocks.markers_in_block[bi] = get_malloc(sizeof(int) * mi);
		for (int i = 0; i < mi; ++i) {
			blocks.markers_in_block[bi][i] = markerbuffer[i];
		}

		++bi;
	}

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
 * @param d pointer to the SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param b struct containing the blocks to use
 * @param output_file string containing the filename of the file to which output
 * block effects/local BVs will be saved.
 * @param group group number from which to split the top individuals.
 */
void calculate_group_local_bvs(SimData* d, MarkerBlocks b, const char* output_file, int group) {

	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}

	int bufferlen = 100;
	char buffer[bufferlen];

	int gsize = get_group_size(d, group);
	char** ggenos = get_group_genes(d, group, gsize);
	char** gnames = get_group_names(d, group, gsize);

	double beffect;

	// for each group member
	for (int i = 0; i < gsize; ++i) {
		// for each block
		sprintf(buffer, "%s_1", gnames[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);

		// for each block
		for (int j = 0; j < b.num_blocks; ++j) {
			beffect = 0;

			// calculate the local BV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
				for (int q = 0; q < d->e.effects.rows; ++q) {
					if (ggenos[i][2 * b.markers_in_block[j][k]] == d->e.effect_names[q]) {
						beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
					}
				}
			}

			// print the local BV
			fprintf(outfile, " %lf", beffect);
			fflush(outfile);
		}

		sprintf(buffer, "\n%s_2", gnames[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);

		// for each block for the second haplotype
		for (int j = 0; j < b.num_blocks; ++j) {
			beffect = 0;
			// calculate the local BV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
				for (int q = 0; q < d->e.effects.rows; ++q) {
					if (ggenos[i][2 * b.markers_in_block[j][k] + 1] == d->e.effect_names[q]) {
						beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
					}
				}
			}

			// print the local BV
			fprintf(outfile, " %lf", beffect);
			fflush(outfile);
		}
		fwrite("\n", sizeof(char), 1, outfile);
	}

	free(ggenos);
	free(gnames);

	fflush(outfile);
	fclose(outfile);
	return;
}

/** Given a set of blocks of markers in a file, for each genotype saved,
 * calculate the local BV for the first allele at each marker in the block, and
 * the local BV for the second allele at each marker in the block, then save
 * the result to a file. This gives block effects for each haplotype of each
 * individual currently saved to the SimData.
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
 * @param d pointer to the SimData object to which individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param b struct containing the blocks to use
 * @param output_file string containing the filename of the file to which output
 * block effects/local BV will be saved.
 */
void calculate_local_bvs(SimData* d, MarkerBlocks b, const char* output_file) {
	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}

	int bufferlen = 100;
	char buffer[bufferlen];

	int gsize = 0;
	AlleleMatrix* m = d->m;
	do {
		gsize += m->n_genotypes;
	} while ((m = m->next) != NULL);

	double beffect;

	// for each group member
	m = d->m;
	int total_i = 0;
	do {
		for (int i = 0; i < m->n_genotypes; ++i, ++total_i) {
			// for each block
			sprintf(buffer, "%s_1", m->names[i]);
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);

			// for each block
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;

				// calculate the local BV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
					for (int q = 0; q < d->e.effects.rows; ++q) {
						if (m->alleles[i][2 * b.markers_in_block[j][k]] == d->e.effect_names[q]) {
							beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
						}
					}
				}

				// print the local BV
				fprintf(outfile, " %lf", beffect);
				fflush(outfile);
			}

			sprintf(buffer, "\n%s_2", m->names[i]);
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);

			// for each block for the second haplotype
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;
				// calculate the local BV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {
					for (int q = 0; q < d->e.effects.rows; ++q) {
						if (m->alleles[i][2 * b.markers_in_block[j][k] + 1] == d->e.effect_names[q]) {
							beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
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

	fflush(outfile);
	fclose(outfile);
	return;
}

/** Takes a look at the currently-loaded effect values and creates a string
 * containing the allele with the highest effect value for each marker as ordered
 * in the SimData. Returns this as a null-terminated heap array of characters of
 * length d->n_markers + one null byte.
 *
 * The return value should be freed when usage is finished!
 *
 * The SimData must be initialised with marker effects for this function to succeed.
 *
 * @param d pointer to the SimData containing markers and marker effects.
 * @returns a heap array filled with a null-terminated string containing the highest
 * scoring allele at each marker.
 */
char* calculate_optimal_alleles(SimData* d) {
	if (d->e.effects.matrix == NULL || d->e.effects.rows < 1 || d->e.effect_names == NULL) {
		warning( "No effect values are loaded\n");
		return NULL;
	}

	char* optimal = get_malloc(sizeof(char)* (d->n_markers + 1));
	char best_allele;
	double best_score;

	for (int i = 0; i < d->n_markers; ++i) {
		best_allele = d->e.effect_names[0];
		best_score = d->e.effects.matrix[0][i];
		for (int a = 1; a < d->e.effects.rows; ++a) {
			if (d->e.effects.matrix[a][i] > best_score) {
				best_score = d->e.effects.matrix[a][i];
				best_allele = d->e.effect_names[a];
			}
		}
		optimal[i] = best_allele;
	}
	optimal[d->n_markers] = '\0';
	return optimal;
}

/** Takes a look at the currently-loaded effect values and returns the highest possible
 * breeding value any (diploid) genotype could have using those effect values.
 *
 * The SimData must be initialised with marker effects for this function to succeed.
 *
 * @param d pointer to the SimData containing markers and marker effects.
 * @returns the fitness metric/breeding value of the best/ideal genotype.
 */
double calculate_optimum_bv(SimData* d) {
	double best_gebv = 0;
	double best_score;

	for (int i = 0; i < d->n_markers; ++i) {
		// Find the allele with the highest effect
		best_score = d->e.effects.matrix[0][i];
		for (int a = 1; a < d->e.effects.rows; ++a) {
			if (d->e.effects.matrix[a][i] > best_score) {
				best_score = d->e.effects.matrix[a][i];
			}
		}

		// add that highest allele to the score twice over
		best_gebv += (2*best_score);
	}

	return best_gebv;
}

/** Takes a look at the currently-loaded effect values and returns the lowest possible
 * breeding value any (diploid) genotype could score using those effect values.
 *
 * The SimData must be initialised with marker effects for this function to succeed.
 *
 * @param d pointer to the SimData containing markers and marker effects.
 * @returns the fitness metric/breeding value of the worst genotype.
 */
double calculate_minimum_bv(SimData* d) {
	double worst_gebv = 0;
	double worst_score;

	for (int i = 0; i < d->n_markers; ++i) {
		// Find the allele with the highest effect
		worst_score = d->e.effects.matrix[0][i];
		for (int a = 1; a < d->e.effects.rows; ++a) {
			if (d->e.effects.matrix[a][i] < worst_score) {
				worst_score = d->e.effects.matrix[a][i];
			}
		}

		// add that highest allele to the score twice over
		worst_gebv += (2*worst_score);
	}

	return worst_gebv;
}

/*--------------------------------Printing-----------------------------------*/

/** Prints the setup data (everything except the actual genotypes) stored inside
 * a SimData to a file. Column separators are tabs.
 *
 * The printing format is:
 *
 * name	chr	pos [allele to which these effects apply, 0+ columns]
 *
 * [marker name]	[chr number]	[chr pos]	[effects, 0+ columns]
 *
 * [marker name]	[chr number]	[chr pos]	[effects, 0+ columns]
 *
 * ...
 *
 * If m->effects is NULL, m->ref_alleles is NULL, or m->genetic_map is NULL,
 * then the relevant columns are omitted.
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the SimData whose data we print
*/
void save_simdata(FILE* f, SimData* m) {
	/* Print the header. */
	//fprintf(f, "name\t");
	fwrite("name\t", sizeof(char), 5, f);
	if (m->map.positions != NULL) {
		//fprintf(f, "chr\tpos");
		fwrite("chr\tpos", sizeof(char), 7, f);
	}
	if (m->e.effect_names != NULL) {
		for (int i = 0; i < m->e.effects.rows; i++) {
			//fprintf(f, "\t%c", m->e.effect_names[i]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(m->e.effect_names + i, sizeof(char), 1, f);
		}
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

	/* Print the body. */
	for (int j = 0; j < m->n_markers; j++) {
		//fprintf(f, "%s\t", m->markers[j]);
		fwrite(m->markers[j], sizeof(char) * strlen(m->markers[j]), 1, f);
		fwrite("\t", sizeof(char), 1, f);

		if (m->map.positions != NULL) {
			//fprintf(f, "%d\t%f", m->map.positions[j].chromosome, m->map.positions[j].position);
			//fwrite(&(m->map.positions[j].chromosome), sizeof(m->map.positions[j].chromosome), 1, f);
			fprintf(f, "%d", m->map.positions[j].chromosome);
			fwrite("\t", sizeof(char), 1, f);
			//fwrite(&(m->map.positions[j].position), sizeof(m->map.positions[j].position), 1, f);
			fprintf(f, "%f", m->map.positions[j].position);
		}

		if (m->e.effects.matrix != NULL) {
			for (int i = 0; i < m->e.effects.rows; i++) {
				//fprintf(f, "\t%lf", m->e.effects.matrix[i][j]);
				fwrite("\t", sizeof(char), 1, f);
				//fwrite(m->e.effects.matrix[i] + j, sizeof(m->e.effects.matrix[i][j]), 1, f);
				fprintf(f, "%f", m->e.effects.matrix[i][j]);
			}
		}
		//fprintf(f, "\n");
		fwrite("\n", sizeof(char), 1, f);
	}
	fflush(f);
}

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
 * m7 and m9 are the names of the markers in the second block.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData whose data we print
 * @param b MarkerBlocks struct containing the groupings of markers to print.
*/
void save_marker_blocks(FILE* f, SimData* d, MarkerBlocks b) {
	const char header[] = "Chrom\tPos\tName\tClass\tMarkers\n";
	fwrite(header, sizeof(char)*strlen(header), 1, f);

	// for the moment we do not name or give locations of different blocks
	const char unspecified[] = "0\t0\tb0\tb\t";
	const int unspeci_length = strlen(unspecified);

	for (int i = 0; i < b.num_blocks; ++i) {
		fwrite(unspecified, sizeof(char)*unspeci_length, 1, f);

		for (int j = 0; j < b.num_markers_in_block[i]; ++j) {
			int k = b.markers_in_block[i][j];

			fwrite(d->markers[k], sizeof(char)*strlen(d->markers[k]), 1, f);
		}

		fwrite("\n", sizeof(char), 1, f);
	}

	fflush(f);
	return;

}

/** Prints all the genotype data saved in the linked list of AlleleMatrices
 * starting with `m` to a file. Uses the following format:
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
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the AlleleMatrix whose data we print
 * @param markers array of strings that correspond to names of the markers.
 * if this is null, the header row will be empty.
*/
void save_allele_matrix(FILE* f, AlleleMatrix* m, char** markers) {
	/* Print header */
	for (int i = 0; i < m->n_markers; ++i) {
		if (markers != NULL) { // assume all-or-nothing with marker names
			//fprintf(f, "\t%s", markers[i]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(markers[i], sizeof(char), strlen(markers[i]), f);
		}
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

	do {
		/* Print the body */
		for (int i = 0; i < m->n_genotypes; ++i) {
			// print the name or ID of the individual.
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			} else {
				//fwrite(group_contents + i, sizeof(int), 1, f);
				fprintf(f, "%d", m->ids[i]);
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
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the AlleleMatrix whose data we print
 * @param markers array of strings that correspond to names of the markers.
*/
void save_transposed_allele_matrix(FILE* f, AlleleMatrix* m, char** markers) {
	// Count number of genotypes in the AM
	AlleleMatrix* currentm = m; // current matrix
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
			fprintf(f, "%d", currentm->ids[currenti]);
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
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the genotypes of the group and
 * the marker names.
 * @param group_id group number of the group of individuals whose genotypes to print.
*/
void save_group_alleles(FILE* f, SimData* d, int group_id) {
	/* Get the stuff we'll be printing. */
	int group_size = get_group_size( d, group_id);
	char** alleles = get_group_genes( d, group_id, group_size);
	char** names = get_group_names( d, group_id, group_size);
	unsigned int* ids = get_group_ids( d, group_id, group_size);

	/* Print header */
	//fwrite(&group_id, sizeof(int), 1, f);
	fprintf(f, "%d", group_id);
	if (d->markers != NULL) {
		for (int i = 0; i < d->n_markers; ++i) {
			// assume all-or-nothing with marker names
			//fprintf(f, "\t%s", markers[i]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(d->markers[i], sizeof(char), strlen(d->markers[i]), f);
		}
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

	/* Print the body */
	for (int i = 0; i < group_size; ++i) {
		// print the name or ID of the individual.
		if (names[i] != NULL) {
			fwrite(names[i], sizeof(char), strlen(names[i]), f);
		} else {
			//fwrite(group_contents + i, sizeof(int), 1, f);
			fprintf(f, "%d", ids[i]);
		}

		for (int j = 0; j < d->n_markers; ++j) {
			//fprintf(f, "\t%c%c", m->alleles[j][2*i], m->alleles[j][2*i + 1]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(alleles[i] + 2*j, sizeof(char), 1, f);
			fwrite(alleles[i] + 2*j + 1, sizeof(char), 1, f);
		}
		///fprintf(f, "\n");
		fwrite("\n", sizeof(char), 1, f);
	}
	free(alleles);
	free(names);
	free(ids);
	fflush(f);

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
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the genotypes of the group and
 * the marker names.
 * @param group_id group number of the group of individuals whose genotypes to print.
*/
void save_transposed_group_alleles(FILE* f, SimData* d, int group_id) {
	/* Get the stuff we'll be printing. */
	int group_size = get_group_size( d, group_id);
	char** alleles = get_group_genes( d, group_id, group_size);
	char** names = get_group_names( d, group_id, group_size);
	unsigned int* ids = get_group_ids( d, group_id, group_size);

	/* Print header */
	fprintf(f, "%d", group_id);
	for (int i = 0; i < group_size; ++i) {
		fwrite("\t", sizeof(char), 1, f);
		fwrite(names[i], sizeof(char), strlen(names[i]), f);
	}
	//fprintf(f, "\n");
	fwrite("\n", sizeof(char), 1, f);

	/* Print the body */
	for (int i = 0; i < d->n_markers; ++i) {
		// print the name or ID of the individual.
		if (d->markers != NULL && d->markers[i] != NULL) {
			fwrite(d->markers[i], sizeof(char), strlen(d->markers[i]), f);
		}

		for (int j = 0; j < group_size; ++j) {
			fwrite("\t", sizeof(char), 1, f);
			fwrite(alleles[j] + 2*i, sizeof(char), 1, f);
			fwrite(alleles[j] + 2*i + 1, sizeof(char), 1, f);
		}
		///fprintf(f, "\n");
		fwrite("\n", sizeof(char), 1, f);
	}
	free(alleles);
	free(names);
	free(ids);
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
 * pedigrees field in the AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0 (which means the parent is unknown) no
 * name or id is printed for that parent.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * immediate parents of.
 */
void save_group_one_step_pedigree(FILE* f, SimData* d, int group) {
	int group_size = get_group_size( d, group);
	unsigned int* group_contents = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	unsigned int pedigree[2];
	char* name;

	for (int i = 0; i < group_size; i++) {
		/*Group member name*/
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		} else {
			//fwrite(group_contents + i, sizeof(int), 1, f);
			fprintf(f, "%d", group_contents[i]);
		}
		fwrite("\t", sizeof(char), 1, f);

		if (get_parents_of_id( d->m, group_contents[i], pedigree) == 0) {
			// Prints both parents, even if they're the same one.
			/* Parent 1 */
			name = get_name_of_id( d->m, pedigree[0]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (pedigree[0] > 0) {
				//fwrite(pedigree, sizeof(int), 1, f);
				fprintf(f, "%d", pedigree[0]);
			}
			fwrite("\t", sizeof(char), 1, f);

			/* Parent 2 */
			name = get_name_of_id( d->m, pedigree[1]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (pedigree[1] > 0) {
				fprintf(f, "%d", pedigree[1]);
				//fwrite(pedigree + 1, sizeof(int), 1, f);
			}
		}
		fwrite("\n", sizeof(char), 1, f);
	}
	free(group_names);
	free(group_contents);
	fflush(f);
}

/** Print the parents of each genotype in the SimData to a file. The following
 * tab-separated format is used:
 *
 * [genotype name]	[parent 1 name]	[parent 2 name]
 *
 * [genotype name]	[parent 1 name]	[parent 2 name]
 *
 * ...
 *
 * The parents are identified by the two ids saved in the genotype's
 * pedigrees field in the AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id is 0 (which means the parent is unknown) no
 * name or id is printed for that parent.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the genotypes and their pedigrees
 */
void save_one_step_pedigree(FILE* f, SimData* d) {
	unsigned int pedigree[2];
	char* name;
	AlleleMatrix* m = d->m;

	do {
		for (int i = 0; i < m->n_genotypes; ++i) {
			/*Group member name*/
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			} else {
				fprintf(f, "%d", m->ids[i]);
			}
			fwrite("\t", sizeof(char), 1, f);

			if (get_parents_of_id( d->m, m->ids[i], pedigree) == 0) {
				// Even if both parents are the same, print them both.
				/* Parent 1 */
				name = get_name_of_id( d->m, pedigree[0]);
				if (name != NULL) {
					fwrite(name, sizeof(char), strlen(name), f);
				} else if (pedigree[0] > 0) {
					//fwrite(pedigree, sizeof(int), 1, f);
					fprintf(f, "%d", pedigree[0]);
				}
				fwrite("\t", sizeof(char), 1, f);

				/* Parent 2 */
				name = get_name_of_id( d->m, pedigree[1]);
				if (name != NULL) {
					fwrite(name, sizeof(char), strlen(name), f);
				} else if (pedigree[1] > 0) {
					fprintf(f, "%d", pedigree[1]);
					//fwrite(pedigree + 1, sizeof(int), 1, f);
				}

			}
			fwrite("\n", sizeof(char), 1, f);
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
 * genotype's pedigrees field in the AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * pedigree of.
 */
void save_group_full_pedigree(FILE* f, SimData* d, int group) {
	int group_size = get_group_size( d, group);
	unsigned int* group_contents = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	const char newline[] = "\n";
	unsigned int pedigree[2];

	for (int i = 0; i < group_size; i++) {
		/*Group member name*/
		fprintf(f, "%d\t", group_contents[i]);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}

		if (get_parents_of_id(d->m, group_contents[i], pedigree) == 0) {
			save_parents_of(f, d->m, pedigree[0], pedigree[1]);
		}
		fwrite(newline, sizeof(char), 1, f);
	}
	free(group_names);
	free(group_contents);
	fflush(f);
}

/** Print the full known pedigree of each genotype in the SimData
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
 * genotype's pedigrees field in the AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing all genotypes to print.
 */
void save_full_pedigree(FILE* f, SimData* d) {
	const char newline[] = "\n";

	AlleleMatrix* m = d->m;

	do {
		for (int i = 0; i < m->n_genotypes; ++i) {
			/*Group member name*/
			fprintf(f, "%d\t", m->ids[i]);
			if (m->names[i] != NULL) {
				fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
			}

			if (m->pedigrees[0][i] != 0 || m->pedigrees[1][i] != 0) {
				save_parents_of(f, d->m, m->pedigrees[0][i], m->pedigrees[1][i]);
			}
			fwrite(newline, sizeof(char), 1, f);
		}
	} while ((m = m->next) != NULL);
	fflush(f);
}

/** Print the full known pedigree of each genotype in a single AlleleMatrix
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
 * genotype's pedigrees field in the AlleleMatrix struct.
 *
 * If a group member or parent has no name, the name will be
 * replaced in the output file with its session-unique id. If the parent
 * id of an individual is 0, the individual is printed without brackets or
 * parent pedigrees and recursion stops here.
 *
 * Note this does not follow through the linked list of AlleleMatrix.
 *
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the AlleleMatrix containing the genotypes to print
 * @param parents pointer to an AlleleMatrix that heads the linked list
 * containing the parents and other ancestry of the given genotypes.
 */
void save_AM_pedigree(FILE* f, AlleleMatrix* m, SimData* parents) {
	const char newline[] = "\n";

	for (int i = 0; i < m->n_genotypes; ++i) {
		/*Group member name*/
		fprintf(f, "%d\t", m->ids[i]);
		if (m->names[i] != NULL) {
			fwrite(m->names[i], sizeof(char), strlen(m->names[i]), f);
		}

        if (m->pedigrees[0][i] != 0 || m->pedigrees[1][i] != 0) {
            save_parents_of(f, parents->m, m->pedigrees[0][i], m->pedigrees[1][i]);
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
 * @param m pointer to an AlleleMatrix that heads the linked list
 * containing the parents and other ancestry of the given id.
 * @param p1 the session-unique id of the first parent to be saved.
 * @param p2 the session-unique id of the second parent to be saved.
 */
void save_parents_of(FILE* f, AlleleMatrix* m, unsigned int p1, unsigned int p2) {
	unsigned int pedigree[2];

	// open brackets
	fwrite("=(", sizeof(char), 2, f);
	char* name;

	// enables us to print only the known parent if one is unknown
	if (p1 == 0 || p2 == 0) {
		p1 = (p1 >= p2) ? p1 : p2; //max of the two
		p2 = p1;
	}

	if (p1 == p2) {
		if (p1 > 0) { //print nothing if both are unknown.
			// Selfed parent
			name = get_name_of_id( m, p1);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (p1 > 0) {
				fprintf(f, "%d", p1);
				//fwrite(pedigree, sizeof(int), 1, f);
			}

			if (get_parents_of_id(m, p1, pedigree) == 0) {
				save_parents_of(f, m, pedigree[0], pedigree[1]);
			}
		}
	} else {
		// Parent 1
		name = get_name_of_id( m, p1);
		if (name != NULL) {
			fwrite(name, sizeof(char), strlen(name), f);
		} else if (p1 > 0) {
			fprintf(f, "%d", p1);
			//fwrite(pedigree, sizeof(int), 1, f);
		}
		if (get_parents_of_id(m, p1, pedigree) == 0) {
			save_parents_of(f, m, pedigree[0], pedigree[1]);
		}

		// separator
		fwrite(",", sizeof(char), 1, f);

		// Parent 2
		name = get_name_of_id( m, p2);
		if (name != NULL) {
			fwrite(name, sizeof(char), strlen(name), f);
		} else if (p2 > 0) {
			fprintf(f, "%d", p2);
			//fwrite(pedigree + 1, sizeof(int), 1, f);
		}

		if (get_parents_of_id(m, p2, pedigree) == 0) {
			save_parents_of(f, m, pedigree[0], pedigree[1]);
		}

	}

	// close brackets
	fwrite(")", sizeof(char), 1, f);
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
 * The SimData must have loaded marker effects for this function to succeed.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * breeding values of.
 */
void save_group_bvs(FILE* f, SimData* d, int group) {
	int group_size = get_group_size( d, group);
	unsigned int* group_contents = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	DecimalMatrix effects = calculate_group_bvs(d, group);
	const char newline[] = "\n";
	const char tab[] = "\t";

	for (int i = 0; i < group_size; ++i) {
		/*Group member name*/
		//fwrite(group_contents + i, sizeof(int), 1, f);
		fprintf(f, "%d", group_contents[i]);
		fwrite(tab, sizeof(char), 1, f);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}
		fwrite(tab, sizeof(char), 1, f);
		//fwrite(effects.matrix[0], sizeof(float), 1, f);
		fprintf(f, "%f", effects.matrix[0][i]);
		fwrite(newline, sizeof(char), 1, f);
	}

	delete_dmatrix(&effects);
	free(group_names);
	free(group_contents);
	fflush(f);
}

/** Print the breeding value of each genotype in the SimData to a file. The following
 * tab-separated format is used:
 *
 * [id]	[name]	[BV]
 *
 * [id]	[name]	[BV]
 *
 * ...
 *
 * The SimData must have loaded marker effects for this function to succeed.
 *
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 */
void save_bvs(FILE* f, SimData* d) {
	AlleleMatrix* am = d->m;
	const char newline[] = "\n";
	const char tab[] = "\t";
	DecimalMatrix effects;

	do {
		effects = calculate_bvs(am, &(d->e));
		for (int i = 0; i < effects.cols; ++i) {
			/*Group member name*/
			//fwrite(group_contents + i, sizeof(int), 1, f);
			fprintf(f, "%d", am->ids[i]);
			fwrite(tab, sizeof(char), 1, f);
			if (am->names[i] != NULL) {
				fwrite(am->names[i], sizeof(char), strlen(am->names[i]), f);
			}
			fwrite(tab, sizeof(char), 1, f);
			//fwrite(effects.matrix[0], sizeof(float), 1, f);
			fprintf(f, "%f", effects.matrix[0][i]);
			fwrite(newline, sizeof(char), 1, f);
		}
		delete_dmatrix(&effects);
	} while ((am = am->next) != NULL);
	fflush(f);
}

/** Print the provided breeding values of each provided name and id to a file,
 * with the same format as a regular call to `save_bvs` or `save_group_bvs`.
 * The following tab-separated format is used:
 *
 * [id]	[name]	[BV]
 *
 * [id]	[name]	[BV]
 *
 * ...
 *
 * The function assumes the array of ids, array of names, and columns of the
 * DecimalMatrix are ordered the same, so a single index produces the corresponding
 * value from each. It is used by as-you-go savers within crossing functions. For
 * general use consider `save_bvs` or `save_group_bvs` instead.
 *
 * @param f file pointer opened for writing to put the output
 * @param e pointer to the DecimalMatrix containing the BVs in the first row.
 * @param ids array of ids to print alongside the BVs.
 * @param names array of names to print alongside the BVs.
 */
void save_manual_bvs(FILE* f, DecimalMatrix* e, unsigned int* ids, char** names) {
	char sep[] = "\t";
	char newline[] = "\n";

	for (int i = 0; i < e->cols; ++i) {
		//fwrite(ids + i, sizeof(int), 1, f);
		fprintf(f, "%d", ids[i]);
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

/** Print the number of copies of a particular allele at each marker of each genotype
 * in the SimData to a file. The following tab-separated format is used:
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
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 * @param allele the allele character to count
 */
void save_count_matrix(FILE* f, SimData* d, char allele) {
	DecimalMatrix counts = calculate_full_count_matrix_of_allele(d->m, allele);

	AlleleMatrix* currentm = d->m;
	// print the header
	for (int i = 0, currenti = 0; i < counts.rows; ++i, ++currenti) {
		if (currenti >= currentm->n_genotypes) {
			currenti = 0;
			currentm = currentm->next;
		}
		fwrite("\t", sizeof(char), 1, f);
		if (currentm->names[currenti] != NULL) {
			fwrite(currentm->names[currenti], sizeof(char), strlen(currentm->names[currenti]), f);
		}
	}

	fwrite("\n", sizeof(char), 1, f);

	// Print the body
	for (int i = 0; i < d->n_markers; ++i) { // loop through markers
		if (d->markers != NULL && d->markers[i] != NULL) {
			fwrite(d->markers[i], sizeof(char), strlen(d->markers[i]), f);
		}
		fwrite("\t", sizeof(char), 1, f);

		for (int j = 0; j < counts.rows; ++j) { // loop through genotypes
			// print the matrix entries
			fprintf(f, "%d ", (int) counts.matrix[j][i]);
		}
		//print the newline
		if (i + 1 < counts.rows) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}

	delete_dmatrix(&counts);
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
 * @param f file pointer opened for writing to put the output
 * @param d pointer to the SimData containing the group members.
 * @param group group number of the group of individuals to print the
 * allele count of.
 * @param allele the allele character to count
 */
void save_count_matrix_of_group(FILE* f, SimData* d, char allele, int group) {
	unsigned int group_size = get_group_size( d, group);
	char** group_names = get_group_names( d, group, group_size);
	DecimalMatrix counts = generate_zero_dmatrix(group_size,d->n_markers);
	calculate_group_count_matrix_of_allele(d,group,allele,&counts);

	fprintf(f, "%d", group);
	// print the header
	for (int i = 0; i < group_size; ++i) {
		fwrite("\t", sizeof(char), 1, f);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}
	}

	fwrite("\n", sizeof(char), 1, f);

	// Print the body
	for (int i = 0; i < d->n_markers; ++i) { // loop through markers
		if (d->markers != NULL && d->markers[i] != NULL) {
			fwrite(d->markers[i], sizeof(char), strlen(d->markers[i]), f);
		}
		fwrite("\t", sizeof(char), 1, f);

		for (int j = 0; j < group_size; ++j) { // loop through genotypes
			// print the matrix entries
			fprintf(f, "%d ", (int) counts.matrix[j][i]);
		}
		//print the newline
		if (i + 1 < group_size) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}

	delete_dmatrix(&counts);
	free(group_names);
	fflush(f);
}
