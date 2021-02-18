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

/** Takes the top `top_n` fitness individuals in the group and puts them in a new group.
 * The new group number is returned.*/
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

/** Calculates the fitness metric for each column/subject in the AlleleMatrix
* for generation `generation` in the SimData object.
*
* This is calculated as: 
* `Sum across all alleles (effect vector of this allele x count matrix of this allele)`
*
* If either the effect matrix or the allele matrix for that generation does not
* exist, the program exits with error status 2. 
*
* @param d pointer to the SimData object containing the effect vector and AlleleMatrix
* @param group integer representing the group for which to run these calculations
* @param effect_row the single-character name of the effect row to use. If it cannot 
* be found in `d->e.effect_names`, the first row is used.
* @returns A DecimalMatrix containing the score for each column of the AlleleMatrix.
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


/** Calculates the number of times in each cell that a particular allele appears.
 * Returns the result as a DecimalMatrix. Useful for multiplying to effect matrix.
 *
 * Finds the `generation`th AlleleMatrix in the linked list whose 0th element is
 * stored in `d`, and counts the number of occurences of `allele` in each
 * two-allele cell.
 *
 * @param a pointer to the AlleleMatrix whose alleles are to be counted
 * @param allele the single-character allele to be counting.
 * @returns A DecimalMatrix countaining the number of `allele` occurences in 
 * each cell of the AlleleMatrix of the correct generation.
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

/** Calculates the number of times in each cell that a particular allele appears.
 * Returns the result as a DecimalMatrix. Useful for multiplying to effect matrix.
 *
 * Finds the `generation`th AlleleMatrix in the linked list whose 0th element is
 * stored in `d`, and counts the number of occurences of `allele` in each
 * two-allele cell.
 *
 * @param a pointer to the AlleleMatrix whose alleles are to be counted
 * @param allele the single-character allele to be counting.
 * @returns A DecimalMatrix countaining the number of `allele` occurences in 
 * each cell of the AlleleMatrix of the correct generation.
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
