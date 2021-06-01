#include "utils.h"

/** Options parameter to run SimData functions in their bare-bones form.*/
const GenOptions BASIC_OPT = {
	.will_name_subjects = FALSE,
	.subject_prefix = NULL,
	.family_size = 1,
	.will_track_pedigree = FALSE,
	.will_allocate_ids = TRUE,
	.filename_prefix = NULL,
	.will_save_pedigree_to_file = FALSE,
	.will_save_effects_to_file = FALSE,
	.will_save_genes_to_file = FALSE,
	.will_save_to_simdata = TRUE
};

/** Replace calls to malloc direct with this function, which errors and exits
 * with status 2 if memory allocation fails. 
 *
 * @param size size of space to be allocated. Treat like the parameter of a 
 * regular malloc call.
 * @returns pointer to the allocated space.
 */
void* get_malloc(size_t size) {
	void* v = malloc(size);
	if (v == NULL) { 
		fprintf(stderr, "Memory allocation of size %lu failed.\n", (unsigned long)size);
		exit(3);
	}
	return v;
}


/** Creator for an empty AlleleMatrix object of a given size. Includes memory
 * allocation for `n_subjects` worth of `.alleles`.
 *
 * @param n_markers number of rows/markers to create
 * @param n_subjects number of individuals to create. This includes filling the first
 * n_subjects entries of .alleles with heap char* of length n_markers, so that the
 * alleles for these can be added without further memory allocation.
 * @returns pointer to the empty created AlleleMatrix
 */
AlleleMatrix* create_empty_allelematrix(int n_markers, int n_subjects) {
	AlleleMatrix* m = get_malloc(sizeof(AlleleMatrix));
	
	m->n_subjects = n_subjects;
	m->n_markers = n_markers;
	//m->alleles = get_malloc(sizeof(char*) * 1000);
	for (int i = 0; i < n_subjects; ++i) {
		m->alleles[i] = get_malloc(sizeof(char) * (n_markers<<1));
		memset(m->alleles[i], 0, sizeof(char) * (n_markers<<1));
		//m->ids[i] = 0;
	}

	memset(m->ids, 0, sizeof(unsigned int) * 1000);
	memset(m->pedigrees[0], 0, sizeof(unsigned int) * 1000);
	memset(m->pedigrees[1], 0, sizeof(unsigned int) * 1000);
	memset(m->groups, 0, sizeof(unsigned int) * 1000);
	memset(m->subject_names, 0, sizeof(char*) * 1000); // setting the pointers to NULL
	memset(m->alleles + n_subjects, 0, sizeof(char*) * (1000 - n_subjects + 1)); // setting the pointers to NULL
	
	/*for (int i = 0; i < 1000; ++i) {
		m->subject_names[i] = NULL;
		m->pedigrees[0][i] = 0;
		m->pedigrees[1][i] = 0;
	}*/
	
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
void set_subject_ids(SimData* d, int from_index, int to_index) {
	if (to_index < from_index) {
		error( "Bad range for setting ids\n");
	}
	
	AlleleMatrix* m = d->m;
	int i, total_i = 0;
	// find the AlleleMatrix from_index belongs to
	while (total_i + m->n_subjects <= from_index) {
		if (m->next == NULL) {
			error( "That index does not exist.");
		} else {
			total_i += m->n_subjects;
			m = m->next;
		}
	}
	total_i = from_index;
	
	// check for overflow:
	if (d->current_id > UINT_MAX - to_index + from_index) {
		warning("IDs will overflow in the process of this calculation. Past this point id-based functions (anything to do with pedigrees) have undefined behaviour.\n");
	}
	
	// allocate the ids
	while (1) {
		for (i = 0; i < m->n_subjects; ++i, ++total_i) {
			++ d->current_id;
			m->ids[i] = d->current_id;
			if (total_i >= to_index) {
				return;
			}
		}		
		
		if (m->next == NULL) {
			warning("Could not set all ids.");
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
				error("Bad columns on row %d\n", details.num_rows + 1);
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
		warning("Search failed: are you sure your list is sorted?\n");
	}
	return index;
}

/** Returns the first located index in an array of strings where the string
 * is the same as the string `target`. Returns -1 if no match was found.
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

/** Fills the designated section of the `.subject_names` array in an 
 * AlleleMatrix with the pattern `prefix`index.
 *
 * The index is padded to either 3 or 5 zeros depending on the size of 
 * `a->n_subjects`.
 *
 * WARNING: after this is run, n_subjects must be updated by n_new_subjects or
 * some of these strings will not be freed on delete.
 *
 * @param pointer to the AlleleMatrix whose `.subject_names` to modify
 * @param prefix the prefix to add to the suffix to make the new subject name
 * @param suffix suffixes start at this value and increment for each additional subject_name
 * @param from_index the new names are added to this index and all those following it in this
 * AlleleMatrix.
*/
void set_subject_names(AlleleMatrix* a, const char* prefix, int suffix, int from_index) {
	char sname[30];
	char format[30];
	if (prefix == NULL) {
		// make it an empty string instead, so it is not displayed as (null)
		prefix = ""; 
	}
	// use sname to save the number of digits to pad by:
	sprintf(sname, "%%0%dd", get_integer_digits(a->n_subjects - from_index));  // Creates: %0[n]d
	sprintf(format, "%s%s", prefix, sname);
	
	++suffix;
	for (int i = from_index; i < a->n_subjects; ++i) {
		// clear subject name if it's pre-existing
		if (a->subject_names[i] != NULL) {
			free(a->subject_names[i]);
		}
		
		// save new name
		sprintf(sname, format, suffix);
		a->subject_names[i] = get_malloc(sizeof(char) * (strlen(sname) + 1));
		strcpy(a->subject_names[i], sname);
		
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
 * @see get_top_n_fitness_subject_indexes()
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
 * @see get_top_n_fitness_subject_indexes()
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
 * to have no gaps. That is, if there are more than 1000 genotypes, the
 * first AM should be full/contain 1000 genotypes, and so forth, and the 
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
		if (filler_m->n_subjects < 1000) {
			for (filler = 0; filler < 1000; ++filler) {
				// an individual is considered to not exist if it has no genome.
				if (filler_m->alleles[filler] == NULL) {
					break; // escape for loop
				}
			}
			// assume we've found one, since n_subjects < 1000.
			break; // escape while loop
		}
		
		if (filler_m->next == NULL) {
			// No empty spaces
			return;
		} else {
			filler_m = filler_m->next;
		}
	}
	
	// loop through all subjects with checker, shifting them back when we find them.
	while (1) {
		for (checker = 0; checker < 1000; ++checker) {
			if (checker_m->alleles[checker] == NULL) {
				// check our filler has a substitute
				while (filler_m->alleles[filler] == NULL) {
					++filler;
					if (filler >= 1000) {
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
				checker_m->subject_names[checker] = filler_m->subject_names[filler];
				filler_m->subject_names[filler] = NULL;
				checker_m->ids[checker] = filler_m->ids[filler];
				filler_m->ids[filler] = -1;
				checker_m->pedigrees[0][checker] = filler_m->pedigrees[0][filler];
				checker_m->pedigrees[1][checker] = filler_m->pedigrees[1][filler];
				filler_m->pedigrees[0][filler] = 0;
				filler_m->pedigrees[1][filler] = 0;
				checker_m->groups[checker] = filler_m->groups[filler];
				filler_m->groups[filler] = 0;
				if (checker_m != filler_m) {
					++ checker_m->n_subjects;
					-- filler_m->n_subjects;
					
					if (filler_m->n_subjects == 0) {
						// find the previous matrix in the linked list
						AlleleMatrix* previous = checker_m;
						while (previous->next != filler_m) {
							previous = previous->next;
							if (previous == NULL) {
								fprintf(stderr, "Error: filler got somewhere checker can't go.\n");
								exit(1);
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
			error("During allele matrix condensing, checker got ahead of filler somehow\n");
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
 * This function assumes that ids are never reshuffled in the SimData. As of 
 * 2021.04.12 this condition is always true, but this function will need to 
 * be revisited if this is no longer the case.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose name is sought
 * @returns the name of the genotype that has id `id`, as a copy of the pointer
 * to the heap memory where the name is saved (so *don't* free the pointer returned
 * from this function)
 */
char* get_name_of_id( AlleleMatrix* start, unsigned int id) {
	if (id <= 0) {
		warning("Invalid negative ID %d\n", id);
		return NULL;	
	}
	if (start == NULL) {
		error("Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;
	
	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_subjects != 0 && id >= m->ids[0] && id <= m->ids[m->n_subjects - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_subjects);
			
			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					error("Could not find the ID %d\n", id);
				} else {
					m = m->next;
				}
			}
			
			return m->subject_names[index];

		}
		
		if (m->next == NULL) {
			error("Could not find the ID %d\n", id);
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
 * This function assumes that ids are never reshuffled in the SimData. As of 
 * 2021.04.12 this condition is always true, but this function will need to 
 * be revisited if this is no longer the case.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose alleles are sought
 * @returns the alleles of the genotype that has id `id`, as a copy of the pointer
 * to the heap memory where the genotype is saved (so *don't* free the pointer returned
 * from this function). It points to a sequence of characters, ordered according to
 * the markers in the SimData to which the AlleleMatrix belongs.
 */
char* get_genes_of_id ( AlleleMatrix* start, unsigned int id) {
	if (id <= 0) {
		warning("Invalid negative ID %d\n", id);
		return NULL;
	}
	if (start == NULL) {
		error("Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;
	
	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_subjects != 0 && id >= m->ids[0] && id <= m->ids[m->n_subjects - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_subjects);
			
			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					error("Could not find the ID %d\n", id);
				} else {
					m = m->next;
					continue;
				}
			}
			
			return m->alleles[index];

		}
		
		if (m->next == NULL) {
			error("Could not find the ID %d\n", id);
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
 * This function assumes that ids are never reshuffled in the SimData. As of 
 * 2021.04.12 this condition is always true, but this function will need to 
 * be revisited if this is no longer the case.
 *
 * @param start Pointer to the first of a linked list of AlleleMatrixes in which
 * the genotype with the provided id is assumed to be found.
 * @param id the id of the genotype whose parents are sought
 * @param output An array which the calling function can access where this function
 * will put its results.
 * @returns 0 when the id is successfully identified and at least one parent's
 * id is known, 1 otherwise. The ids of both parents if at least one parent is 
 * known/nonzero are saved to the array `output`.
 */
int get_parents_of_id( AlleleMatrix* start, unsigned int id, unsigned int output[2]) {
	if (id <= 0) {
		return 1;	
	}
	if (start == NULL) {
		error("Invalid nonexistent allelematrix\n");
	}
	AlleleMatrix* m = start;
	while (1) {
		// try to find our id. Does this AM have the right range for it?
		if (m->n_subjects != 0 && id >= m->ids[0] && id <= m->ids[m->n_subjects - 1]) {
			// perform binary search to find the exact index.
			int index = get_from_ordered_uint_list(id, m->ids, m->n_subjects);
			
			if (index < 0) {
				// search failed
				if (m->next == NULL) {
					error("Could not find the ID %d\n", id);
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
			error("Could not find the ID %d\n", id);
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
	AlleleMatrix* m = start;
	int i, j;
	
	for (i = 0; i < n_names; ++i) {
		found = FALSE;
		output[i] = -1;
		while (1) {
			// try to identify the name in this AM
			for (j = 0; j < m->n_subjects; ++j) {
				if (strcmp(m->subject_names[j], names[i]) == 0) {
					found = TRUE;
					output[i] = m->ids[j];
					break;
				}
			}
			
			if (found) {
				break;
			}
			if ((m = m->next) == NULL) {
				error("Didn't find the name %s\n", names[i]);
				break;
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
		for (j = 0; j < m->n_subjects; ++j) {
			if ((parent1id == m->pedigrees[0][j] && parent2id == m->pedigrees[1][j]) || 
					(parent1id == m->pedigrees[1][j] && parent2id == m->pedigrees[0][j])) {
				return m->ids[j];
			}
		}

		if ((m = m->next) == NULL) {
			error("Didn't find the child of %d & %d\n", parent1id, parent2id);
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
		for (j = 0; j < m->n_subjects; ++j, ++total_j) {
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
 * @param name a string to match to the subject_name of the target
 * @returns the index (0-based, starting at the start of `start`) of the first sequentially 
 * located genotype whose name is the same as the provided name.
 */
int get_index_of_name( AlleleMatrix* start, char* name) {
	AlleleMatrix* m = start;
	int j, total_j = 0;
	
	while (1) {
		// try to identify the child in this AM
		for (j = 0; j < m->n_subjects; ++j, ++total_j) {
			if (strcmp(m->subject_names[j], name) == 0) {
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
		} else if (total_j < index && total_j + m->n_subjects > index) {
			return m->ids[index - total_j];
		}
		total_j += m->n_subjects;

		if ((m = m->next) == NULL) {
			error( "Didn't find the index %d\n", index);
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
		warning("Invalid negative index %d\n", index);
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
		} else if (total_j < index && total_j + m->n_subjects > index) {
			return m->alleles[index - total_j];
		}
		total_j += m->n_subjects;

		if ((m = m->next) == NULL) {
			error( "Didn't find the index %d\n", index);
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

/*int get_number_of_subjects( AlleleMatrix* m) {
	AlleleMatrix* m = m;
	int total_subjects = m->n_subjects;
	while (m->next != NULL) {
		m = m->next;
		total_subjects += m->n_subjects;
	}
	return total_subjects;
}*/

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
		for (i = 0; i < m->n_subjects; ++i) {
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
		for (i = 0; i < m->n_subjects; ++i) {
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

/** Function to identify the next sequential integer that does not 
 * identify a group that currently has member(s). 
 *
 * This calls get_existing_groups() every time, so functions that 
 * do repeated group creation, like split_into_individuals, are 
 * recommended to use their own implementation rather than calling 
 * this function repeatedly for best speed/memory usage.
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

/** Function to count the number of genotypes that currently belong to
 * the specified group.
 * @see get_group_genes()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
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

/** Gets a shallow copy of the genes/alleles of each member of the group.
 * @see get_group_size()
 * @see get_group_names()
 * @see get_group_ids()
 * @see get_group_indexes()
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group we want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing pointers to the genes of each member of the group.
 * the vector itself is on the heap and should be freed, but its contents are only
 * shallow copies so should not be freed.
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
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing pointers to the names of each member of the group.
 * the vector itself is on the heap and should be freed, but its contents are only
 * shallow copies so should not be freed.
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
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the ids of each memmber of the group.
 * the vector itself is on the heap and should be freed.
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
 *
 * @param d the SimData struct on which to perform the operation
 * @param group_id the group number of the group you want data from
 * @param group_size if group_size has already been calculated, pass it too, otherwise
 * put in -1. This enables fewer calls to get_group_size when using multiple
 * group-data-getter functions.
 * @returns a vector containing the indexes of each member of the group.
 * the vector itself is on the heap and should be freed.
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

/** Get a DecimalMatrix composed only of row `row_index` from `m`
 *
 * Checks that `row_index` is a valid index.
 * it will throw an error if the index is not valid
 *
 * @param m pointer to the DecimalMatrix whose subset we are taking.
 * @param row_index index in m of the row from m to take as a subset
 * @returns the row of `m` at `row_index` as a new independent DecimalMatrix.
*/
DecimalMatrix subset_dmatrix_row(DecimalMatrix* m, int row_index) {
	if (row_index < 0 || row_index >= m->rows) {
		error( "Invalid index for subsetting: %d\n", row_index);
	}
	DecimalMatrix subset;
	subset.cols = m->cols;
	subset.rows = 1;
	
	subset.matrix = malloc(sizeof(double*) * subset.rows);
	subset.matrix[0] = malloc(sizeof(double) * subset.cols);
	for (int j = 0; j < subset.cols; j++) {
		subset.matrix[0][j] = m->matrix[row_index][j];
	}
	return subset;
}

/** Perform matrix multiplication on two DecimalMatrix. Result is returned as 
 * a new matrix.
 *
 * If the two matrices do not have dimensions that allow them to be multiplied,
 * it will throw an error.
 *
 * The order of parameters a and b matters.
 *
 * @param a pointer to the first DecimalMatrix.
 * @param b pointer to the second DecimalMatrix.
 * @returns the result of the matrix multiplication (*a) x (*b), as a new 
 * DecimalMatrix.
 */
DecimalMatrix multiply_dmatrices(DecimalMatrix* a, DecimalMatrix* b) {
	if (a->cols != b->rows) {
		// multiplication was not possible for these dimensions
		error("Dimensions invalid for matrix multiplication.");
	}
	
	DecimalMatrix product = generate_zero_dmatrix(a->rows, b->cols);
	
	for (int i = 0; i < product.rows; i++) {
		R_CheckUserInterrupt();
		for (int j = 0; j < product.cols; j++) {
			// we loop through every cell in the matrix.
			for (int k = 0; k < a->cols; k++) {
				// for each cell, we loop through each of the pairs adjacent to it.
				product.matrix[i][j] += (a->matrix[i][k]) * (b->matrix[k][j]);
			}
		}
	}
	
	return product;
}

/** Add the two DecimalMatrix. The result is saved in the first DecimalMatrix.
 *
 * If the two matrices do not have the same dimensions, the matrix returned has 
 * it will throw an error.
 *
 * The ordering of the parameters does not matter.
 *
 * @param a pointer to the first DecimalMatrix. Add this to b.
 * @param b pointer to the second DecimalMatrix. Add this to a.
 * @returns a new DecimalMatrix containing the result of *a + *b.
 */
void add_to_dmatrix(DecimalMatrix* a, DecimalMatrix* b) {
	if (a->cols != b->cols || a->rows != b->rows) {
		error( "Invalid dimensions for addition\n");
	}
	
	for (int i = 0; i < a->rows; i++) {
		R_CheckUserInterrupt();
		for (int j = 0; j < a->cols; j++) {
			a->matrix[i][j] += b->matrix[i][j];
		}
	}
}

/** Add the two DecimalMatrix. The result is returned as a new DecimalMatrix.
 *
 * If the two matrices do not have the same dimensions, the matrix returned has 
 * it will throw an error.
 *
 * The ordering of the parameters does not matter.
 *
 * @param a pointer to the first DecimalMatrix. Add this to b.
 * @param b pointer to the second DecimalMatrix. Add this to a.
 * @returns a new DecimalMatrix containing the result of *a + *b.
 */
DecimalMatrix add_dmatrices(DecimalMatrix* a, DecimalMatrix* b) {
	if (a->cols != b->cols || a->rows != b->rows) {
		error( "Invalid dimensions for addition\n");
	}
	
	DecimalMatrix sum = generate_zero_dmatrix(a->rows, a->cols);
	
	for (int i = 0; i < sum.rows; i++) {
		R_CheckUserInterrupt();
		for (int j = 0; j < sum.cols; j++) {
			sum.matrix[i][j] = a->matrix[i][j] + b->matrix[i][j];
		}
	}
	
	return sum;
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
		
		for (i = 0, deleted = 0; i < m->n_subjects; ++i) {
			if (m->groups[i] == group_id) {
				// delete data
				if (m->subject_names[i] != NULL) {
					free(m->subject_names[i]);
					m->subject_names[i] = NULL;
				}
				if (m->alleles[i] != NULL) {
					free(m->alleles[i]);
					m->alleles[i] = NULL;
				}
				++deleted;
			}
		}
		m->n_subjects -= deleted;
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
		if (m->alleles != NULL) {
			for (int i = 0; i < m->n_subjects; i++) {
				if (m->alleles[i] != NULL) {
					free(m->alleles[i]);
				}
				
			}
		}
		
		// free subject names
		if (m->subject_names != NULL) {
			for (int i = 0; i < m->n_subjects; i++) {
				if (m->subject_names[i] != NULL) {
					free(m->subject_names[i]);
				}
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

void SXP_delete_simdata(SEXP sd) {
	//Rprintf("Garbage collecting SimData...\n");
	SimData* d = (SimData*) R_ExternalPtrAddr(sd);
	delete_simdata(d);
	R_ClearExternalPtr(sd);
	return;
}

SEXP SXP_delete_group(SEXP exd, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	delete_group(d, group_id);
	
	return ScalarInteger(0);
}

SEXP SXP_combine_groups(SEXP exd, SEXP len, SEXP groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int n = asInteger(len);
	int *gps = INTEGER(groups); 
	if (n == NA_INTEGER) { 
		error("`len` parameter is of invalid type or `groups` vector is invalid.\n");
	}
	
	return ScalarInteger(combine_groups(d, n, gps));
}

SEXP SXP_split_individuals(SEXP exd, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	// track groups before we affect them
	int old_n_groups = 0;
	int* old_groups = get_existing_groups(d, &old_n_groups);
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	// do the actual split
	split_into_individuals(d, group_id);
	
	//track groups after the split
	int new_n_groups = 0;
	int* new_groups = get_existing_groups(d, &new_n_groups);
	
	// Get an R vector of the same length as the number of new size 1 groups created
	SEXP out = PROTECT(allocVector(INTSXP, new_n_groups - old_n_groups + 1));
	int* outc = INTEGER(out);
	int outi = 0;
	// Find the size 1 groups and save them to this vector
	// They're identified as the group nums that didn't exist before the split
	int is_new;
	for (int i = 0; i < new_n_groups; ++i) {
		is_new = TRUE;
		for (int j = 0; j < old_n_groups; ++j) {
			if (old_groups[j] == new_groups[i]) {
				is_new = FALSE;
				break;
			}
		}
		if (is_new) {
			outc[outi] = new_groups[i];
			++outi;
		}
	}
	free(old_groups);
	free(new_groups);
	UNPROTECT(1);
	return out;
}

SEXP SXP_split_familywise(SEXP exd, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	// track groups before we affect them
	int old_n_groups = 0;
	int* old_groups = get_existing_groups(d, &old_n_groups);
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	// do the actual split
	split_into_families(d, group_id);
	
	//track groups after the split
	int new_n_groups = 0;
	int* new_groups = get_existing_groups(d, &new_n_groups);
	
	// Get an R vector of the same length as the number of new size 1 groups created
	SEXP out = PROTECT(allocVector(INTSXP, new_n_groups - old_n_groups + 1));
	int* outc = INTEGER(out);
	int outi = 0;
	// Find the size 1 groups and save them to this vector
	// They're identified as the group nums that didn't exist before the split
	int is_new;
	for (int i = 0; i < new_n_groups; ++i) {
		is_new = TRUE;
		for (int j = 0; j < old_n_groups; ++j) {
			if (old_groups[j] == new_groups[i]) {
				is_new = FALSE;
				break;
			}
		}
		if (is_new) {
			outc[outi] = new_groups[i];
			++outi;
		}
	}
	free(old_groups);
	free(new_groups);
	UNPROTECT(1);
	return out;	
}

SEXP SXP_split_out(SEXP exd, SEXP len, SEXP indexes) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int n = asInteger(len);
	int *ns = INTEGER(indexes); 
	if (n == NA_INTEGER) { 
		error("`len` parameter is of invalid type or `indexes` vector is invalid.\n");
	}
	for (int i = 0; i < n; ++i) {
		if (ns[i] == NA_INTEGER || ns[i] < 0) {
			error("The `indexes` vector contains at least one invalid index.\n");
		}
	}
	
	return ScalarInteger(split_from_group(d, n, ns));
}

// get the active groups currently for a summary
SEXP SXP_get_groups(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int n_groups = 0;
	int** existing = get_existing_group_counts(d, &n_groups);
	
	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SEXP groups = PROTECT(allocVector(INTSXP, n_groups));
	int* cgroups = INTEGER(groups);
	SEXP counts = PROTECT(allocVector(INTSXP, n_groups));
	int* ccounts = INTEGER(counts);
	
	for (int i = 0; i < n_groups; ++i) {
		cgroups[i] = existing[0][i];
		ccounts[i] = existing[1][i];
	}
	SET_VECTOR_ELT(out, 0, groups);
	SET_VECTOR_ELT(out, 1, counts);
	free(existing[0]);
	free(existing[1]);
	free(existing);
	
	UNPROTECT(3);
	return out;
	
}

SEXP SXP_get_group_data(SEXP exd, SEXP group, SEXP whatData) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	char c = CHAR(asChar(whatData))[0];
	int group_size = get_group_size(d, group_id);
	if (c == 'N' || c == 'n' || c == 'G' || c == 'g') {
		char** data;
		if (c == 'N' || c == 'n') {
			data = get_group_names(d, group_id, group_size);
		} else {
			data = get_group_genes(d, group_id, group_size);
		}
		
		// save to an R vector
		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			SET_STRING_ELT(out, i, mkChar(data[i]));
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'D' || c == 'd' || c == 'I' || c == 'i') {
		unsigned int *data;
		if (c == 'D' || c == 'n') {
			data = get_group_ids(d, group_id, group_size);
		} else {
			data = get_group_indexes(d, group_id, group_size);
		}
		
		// save to an R vector
		SEXP out = PROTECT(allocVector(INTSXP, group_size));
		int* outc = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else {
		error("`data.type` parameter has unknown value.");
	}
}
