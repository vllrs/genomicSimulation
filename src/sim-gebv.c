#include "sim-gebv.h"

SEXP SXP_group_eval(SEXP exd, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	int group_size = get_group_size(d, group_id);
	unsigned int* inds = get_group_indexes(d, group_id, group_size);
	DecimalMatrix gebvs = calculate_fitness_metric_of_group(d, group_id);
	
	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SEXP index = PROTECT(allocVector(INTSXP, group_size));
	int* cindex = INTEGER(index);
	SEXP evals = PROTECT(allocVector(REALSXP, group_size));
	double* cevals = REAL(evals);
	for (int i = 0; i < group_size; ++i) {
		cindex[i] = inds[i];
		cevals[i] = gebvs.matrix[0][i];
	}
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, evals);
	free(inds);
	delete_dmatrix(&gebvs);
	
	UNPROTECT(3);
	return out;
}

SEXP SXP_simple_selection(SEXP exd, SEXP glen, SEXP groups, SEXP number, SEXP bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` vector is invalid.\n");
	}
	
	int num_to_select = asInteger(number);
	if (num_to_select == NA_INTEGER || num_to_select < 0) { 
		error("`number` parameter is of invalid type.\n");
	}
	
	int want_low = asLogical(bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type.\n"); }
	
	if (len == 1) {
		return ScalarInteger(split_group_by_fitness(d, gps[0], num_to_select, want_low));
	} else {
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = split_group_by_fitness(d, gps[i], num_to_select, want_low);
		}
		UNPROTECT(1);
		return out;
	}
	
}

SEXP SXP_simple_selection_bypercent(SEXP exd, SEXP glen, SEXP groups, SEXP percent, SEXP bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` vector is invalid.\n");
	}
	
	double pc_to_select = asReal(percent);
	if (ISNA(pc_to_select) || pc_to_select < 0) { 
		error("`percentage` parameter is of invalid type.\n");
	}
	
	int want_low = asLogical(bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type.\n"); }
	
	if (len == 1) {
		int group_size = get_group_size(d, gps[0]);
		int num_to_select = group_size * pc_to_select / 100; // integer division, so take the floor
		
		return ScalarInteger(split_group_by_fitness(d, gps[0], num_to_select, want_low));
	} else {
		int num_to_select;
		
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			num_to_select = get_group_size(d, gps[i]) * pc_to_select / 100;
			outc[i] = split_group_by_fitness(d, gps[i], num_to_select, want_low);
		}
		UNPROTECT(1);
		return out;
	}
	
}


/*--------------------------------Fitness------------------------------------*/

/** Takes the top `top_n` fitness/GEBV individuals in the group and puts them in a new group.
 * The new group number is returned.
 *
 * @param d pointer to the SimData object to which the groups and individuals belong.
 * It must have a marker effect file loaded to successfully run this function.
 * @param group group number from which to split the top individuals.
 * @param top_n The number of individuals to put in the new group.
 * @param lowIsBest boolean, if TRUE the `top_n` with the lowest fitness/GEBV score
 * will be selected, if false the `top_n` with the highest fitness/GEBV score are.
 * @returns the group number of the newly-created split-off group
 */
int split_group_by_fitness(SimData* d, int group, int top_n, int lowIsBest) {
	// get fitnesses
	unsigned int* group_contents = get_group_indexes( d, group, -1); // should be ordered same as next line would get
	DecimalMatrix fits = calculate_fitness_metric_of_group( d, group );
	
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
	int top_subjects[top_n];
	for (int i = 0; i < top_n; i++) {
		top_subjects[i] = group_contents[p_fits[i] - fits.matrix[0]];
	}
	delete_dmatrix(&fits);
	free(group_contents);
	
	// send those n to a new group
	return split_from_group(d, top_n, top_subjects);
}

/** Calculates the fitness metric/GEBV for each genotype in the AlleleMatrix
* in a certain group, and returns the results in a DecimalMatrix struct.
*
* The GEBV is calculated for each genotype by taking the sum of (number of 
* copies of this allele at this marker times this allele's effect at this marker) 
* for each marker for each different allele. To do
* this, the vector containing each allele's effect values is multiplied by a matrix 
* containing the counts of that allele at each marker. 
*
* The function exits with error code 1 if no marker effect file is loaded.
*
* @param d pointer to the SimData object to which the groups and individuals belong.
* It must have a marker effect file loaded to successfully run this function.
* @param group calculate GEBVs for each genotype in the group with this group number.
* @returns A DecimalMatrix containing the score for each individual in the group.
*/
DecimalMatrix calculate_fitness_metric_of_group(SimData* d, int group) {
	// check that both of the items to be multiplied exist.
	if (d->e.effects.rows < 1 || d->m->alleles == NULL) {
		error( "Either effect matrix or allele matrix does not exist\n");
	}
	
	int group_size = get_group_size( d, group );
	unsigned int* group_members = get_group_ids( d, group, group_size);
	DecimalMatrix sum, counts, effect_row, product;
	// sum is the sum of the product across all alleles. 
	// Its dimensions are the same as those from multiplying a row of d->e.effects by
	// a counts matrix the size of d->m.alleles
	sum = generate_zero_dmatrix(1, group_size); 
	
	for (int i = 0; i < d->e.effects.rows; i++) {
		// get the product for this allele
		counts = calculate_count_matrix_of_allele_for_ids( d->m, group_members, group_size, d->e.effect_names[i]);
		effect_row = subset_dmatrix_row(&(d->e.effects), i);
		product = multiply_dmatrices(&effect_row, &counts);
		delete_dmatrix(&counts);
		delete_dmatrix(&effect_row);
		
		// add it to the total effect GEBV sum
		add_to_dmatrix(&sum, &product);
		///counts = sum; // use counts as a temp
		//sum = add_dmatrices(&counts, &product);
		//delete_dmatrix(&counts);
		delete_dmatrix(&product);
	}
	
	free(group_members);
	return sum;
}

/** Calculates the fitness metric/GEBV for each genotype in the AlleleMatrix,
 * and returns the results in a DecimalMatrix struct.
 *
 * The GEBV is calculated for each genotype by taking the sum of (number of 
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
DecimalMatrix calculate_fitness_metric( AlleleMatrix* m, EffectMatrix* e) {
	// check that both of the items to be multiplied exist.
	if (e->effects.rows < 1 || m->alleles == NULL) {
		error("Either effect matrix or allele matrix does not exist\n");
	}
	
	DecimalMatrix sum, counts, effect_row, product;
	// sum is the sum of the product across all alleles. 
	// Its dimensions are the same as those from multiplying a row of d->e.effects by
	// a counts matrix the size of d->m.alleles
	sum = generate_zero_dmatrix(1, m->n_subjects); 
	
	for (int i = 0; i < e->effects.rows; i++) {
		// get the product for this allele
		counts = calculate_full_count_matrix_of_allele( m, e->effect_names[i]);
		effect_row = subset_dmatrix_row(&(e->effects), i);
		product = multiply_dmatrices(&effect_row, &counts);
		delete_dmatrix(&counts);
		delete_dmatrix(&effect_row);
		
		// add it to the total effect GEBV sum
		counts = sum; // use counts as a temp
		sum = add_dmatrices(&counts, &product);
		delete_dmatrix(&counts);
		delete_dmatrix(&product);
	}
	
	return sum;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype is a set of given genotype.
 * Returns the result as a DecimalMatrix. Useful for multiplying to effect matrix
 * to calculate GEBVs.
 *
 * @param m pointer to the start of a linked list that contains the alleles 
 * of each genotype in `for_ids`.
 * @param for_ids poinder to an array of ids of different genotypes. 
 * @param n_ids number of genotypes to calculate the counts of/length of 
 * `for_ids`
 * @param allele the single-character allele to be counting.
 * @returns A DecimalMatrix countaining the number of `allele` occurences at 
 * each row/marker for each column/genotype in `for_ids`.
 * */
DecimalMatrix calculate_count_matrix_of_allele_for_ids( AlleleMatrix* m, unsigned int* for_ids, unsigned int n_ids, char allele) {
	DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, n_ids);
	double cell_sum;
	char* genes;

	for (int i = 0; i < n_ids; i++) {
		R_CheckUserInterrupt();
		genes = get_genes_of_id( m, for_ids[i] );
		
		cell_sum = 0;
		if (genes == NULL) {
			continue;
		}
		
		for (int j = 0; j < m->n_markers; j++) {
			cell_sum = 0;
			if (genes[2*j] == allele) {
				cell_sum += 1;
			}
			if (genes[2*j + 1] == allele) {
				cell_sum += 1;
			}
			counts.matrix[j][i] = cell_sum;
		}
	}
	return counts;
}

/** Calculates the number of times at each marker that a particular allele appears
 * for each genotype in a particular AlleleMatrix.
 * Returns the result as a DecimalMatrix. Useful for multiplying to effect matrix
 * to calculate GEBVs.
 *
 * @param m pointer to the AlleleMatrix that contains the genotypes to count alleles.
 * @param allele the single-character allele to be counting.
 * @returns A DecimalMatrix countaining the number of `allele` occurences at 
 * each row/marker for each column/genotype in the AlleleMatrix.
 * */
DecimalMatrix calculate_full_count_matrix_of_allele( AlleleMatrix* m, char allele) {
	DecimalMatrix counts = generate_zero_dmatrix(m->n_markers, m->n_subjects);
	double cell_sum;
	char* genes;

	for (int i = 0; i < m->n_subjects; i++) {
		R_CheckUserInterrupt();
		genes = m->alleles[i];
		
		cell_sum = 0;
		if (genes == NULL) {
			continue;
		}
		
		for (int j = 0; j < m->n_markers; j++) {
			cell_sum = 0;
			if (genes[2*j] == allele) {
				cell_sum += 1;
			}
			if (genes[2*j + 1] == allele) {
				cell_sum += 1;
			}
			counts.matrix[j][i] = cell_sum;
		}
	}
	return counts;
}

void calculate_all_block_effects(SimData* d, const char* block_file, const char* output_file) {
	struct TableSize ts = get_file_dimensions(block_file, '\t');
	int n_blocks = ts.num_rows - 2;
  
	FILE* infile, * outfile;
	if ((infile = fopen(block_file, "r")) == NULL) {
		error("Failed to open file %s.\n", block_file);
	}
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error("Failed to open file %s.\n", output_file);
	}
  
	// get the group members by index
	int gsize = 0;
	AlleleMatrix* m = d->m;
	do {
		gsize += m->n_subjects;
	} while ((m = m->next) != NULL);
  
	int bufferlen = 100;
	char buffer[bufferlen];
	//char blockname[bufferlen];
	char markername[bufferlen];
  
	double beffects[gsize*2][n_blocks];
  
	// Ignore the first line
	fscanf(infile, "%*[^\n]\n");
  
	// Loop through rows of the file (each row corresponds to a block)
	//while (fscanf(infile, "%*d %*f %s %*s ", blockname) != EOF) {
	for (int bi = 0; bi < n_blocks; ++bi) { 
		// clear effect values
		for (int i = 0; i < gsize * 2; ++i) {
			beffects[i][bi] = 0;
		}
    
		fscanf(infile, "%*d %*f %s %*s ", buffer);
	
		int c, k = 0;
		while ((c = fgetc(infile)) != EOF && c !='\n') {
			if (c == ';') {
				markername[k] = '\0';
        
				// add the effects of that marker to our counts.
				int markerindex = get_from_unordered_str_list(markername, d->markers, d->n_markers);
				
				m = d->m;
				int total_i = 0;
				do {
					for (int i = 0; i < m->n_subjects; ++i, ++total_i) {
						// add the effect of that allele at that marker.
						for (int j = 0; j < d->e.effects.rows; ++j) {			
							if (m->alleles[i][2*markerindex] == d->e.effect_names[j]) {
								beffects[2*total_i + 0][bi] += d->e.effects.matrix[j][markerindex];
							}
							if (m->alleles[i][2*markerindex + 1] == d->e.effect_names[j]) {
								beffects[2*total_i + 1][bi] += d->e.effects.matrix[j][markerindex];
							}		
						}
					}
				} while ((m = m->next) != NULL);
        
				k = 0;
			} else {
				markername[k] = c;
				++k;
			}
		}
	}
	
	// Save those effects to the output file.
	m = d->m;
	int total_i = 0;
	do {
		for (int i = 0; i < m->n_subjects; ++i, ++total_i) {
			if (m->subject_names[i] != NULL) {
				sprintf(buffer, "%s", m->subject_names[i]);
			} else {
				sprintf(buffer, "G%d", total_i);
			}
			
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);
			fprintf(outfile, "_1");
			for (int j = 0; j < n_blocks; ++j) {
				fprintf(outfile, " %lf", beffects[2*total_i][j]);
			}
			fwrite("\n", sizeof(char), 1, outfile);
			
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);
			fprintf(outfile, "_2");
			for (int j = 0; j < n_blocks; ++j) {
				fprintf(outfile, " %lf", beffects[2*total_i + 1][j]);
			}
			fwrite("\n", sizeof(char), 1, outfile);
		}
	} while ((m = m->next) != NULL);
  
	fclose(infile);
	fflush(outfile);
	fclose(outfile);
	return;
}

/** Given a set of blocks of markers in a file, for each genotype in a group, 
 * calculate the local GEBV for the first allele at each marker in the block, and 
 * the local GEBV for the second allele at each marker in the block, then save
 * the result to a file. This gives block effects for each haplotype of each 
 * individual in the group.
 *
 * Note that this function is made to work on HPC, and not fixed to work on a regular
 * device. If an array as long as the number of genotypes cannot fit contiguously into
 * memory, the function will segfault. 
 * 
 * The block file is designed after the output from a call to the R SelectionTools
 * package's `st.def.hblocks` function. It should have the format (tab-separated):
 *
 * Chrom	Pos	Name	Class	Markers
 *
 * [ignored]	[ignored]	[ignored]	[ignored]	[semicolon];[separated];[list]
 * ;[of];[marker];[names];[belonging];[to];[this];[block]
 *
 * ...
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
 * @param block_file string containing filename of the file with blocks
 * @param output_file string containing the filename of the file to which output 
 * block effects/local GEBVs will be saved.
 * @param group group number from which to split the top individuals.
 */
void calculate_group_block_effects(SimData* d, const char* block_file, const char* output_file, int group) {
	struct TableSize ts = get_file_dimensions(block_file, '\t');
	int n_blocks = ts.num_rows - 2;
  
	FILE* infile, * outfile;
	if ((infile = fopen(block_file, "r")) == NULL) {
		error("Failed to open file %s.\n", block_file);
	}
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error("Failed to open file %s.\n", output_file);
	}
  
  // get the group members by index
	int gsize = get_group_size(d, group);
	char** ggenos = get_group_genes(d, group, gsize);
	char** gnames = get_group_names(d, group, gsize);
  
	int bufferlen = 100;
	char buffer[bufferlen];
	//char blockname[bufferlen];
	char markername[bufferlen];
  
	double beffects[gsize*2][n_blocks];

  
	// Ignore the first line
	fscanf(infile, "%*[^\n]\n");
  
	// Loop through rows of the file (each row corresponds to a block)
	//while (fscanf(infile, "%*d %*f %s %*s ", blockname) != EOF) {
	for (int bi = 0; bi < n_blocks; ++bi) { 
		// clear effect values
		for (int i = 0; i < gsize * 2; ++i) {
			beffects[i][bi] = 0;
		}
    
		fscanf(infile, "%*d %*f %s %*s ", buffer);
	
		int c, k = 0;
		while ((c = fgetc(infile)) != EOF && c !='\n') {
			if (c == ';') {
				markername[k] = '\0';
        
				// add the effects of that marker to our counts.
				int markerindex = get_from_unordered_str_list(markername, d->markers, d->n_markers);
				for (int i = 0; i < gsize; ++i) {
					// add the effect of that allele at that marker.
					for (int j = 0; j < d->e.effects.rows; ++j) {			
						if (ggenos[i][2*markerindex] == d->e.effect_names[j]) {
							beffects[2*i + 0][bi] += d->e.effects.matrix[j][markerindex];
						}
						if (ggenos[i][2*markerindex + 1] == d->e.effect_names[j]) {
							beffects[2*i + 1][bi] += d->e.effects.matrix[j][markerindex];
						}		
					}
				}
        
				k = 0;
			} else {
				markername[k] = c;
				++k;
			}
		}
	}
	
	// Save those effects to the output file.
	for (int i = 0; i < gsize; ++i) {
		sprintf(buffer, "%s_1", gnames[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);
		for (int j = 0; j < n_blocks; ++j) {
			fprintf(outfile, " %lf", beffects[2*i][j]);
		}

		sprintf(buffer, "\n%s_2", gnames[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);
		for (int j = 0; j < n_blocks; ++j) {
			fprintf(outfile, " %lf", beffects[2*i + 1][j]);
		}
		fwrite("\n", sizeof(char), 1, outfile);
	}
  
	free(gnames);
	free(ggenos);
	fclose(infile);
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
char* calculate_ideal_genotype(SimData* d) {
	if (d->e.effects.matrix == NULL || d->e.effects.rows < 1 || d->e.effect_names == NULL) {
		error("No effect values are loaded\n");
	}
	
char* optimal = get_malloc(sizeof(char)* (d->n_markers + 1));
	char best_allele;
	int best_score;
	
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

SEXP SXP_get_best_genotype(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	char* best_genotype = calculate_ideal_genotype(d);
	
	SEXP out = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(out, 0, mkChar(best_genotype));
	free(best_genotype);
	UNPROTECT(1);
	return out;
}