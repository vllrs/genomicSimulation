#include "r-wrappers.h"

/*-------------------------- Loaders -------------------------*/

SEXP SXP_load_data(SEXP alleleFile, SEXP mapFile) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile)));
	load_genmap_to_simdata(d, CHAR(asChar(mapFile)));
	
	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	
	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP SXP_load_data_weff(SEXP alleleFile, SEXP mapFile, SEXP effectFile) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile)));
	load_genmap_to_simdata(d, CHAR(asChar(mapFile)));
	load_effects_to_simdata(d, CHAR(asChar(effectFile)));

	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	
	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP SXP_load_more_genotypes(SEXP exd, SEXP alleleFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_more_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile))));
	//return ScalarInteger(load_more_transposed_genes_to_simdata(&GlobalSim, CHAR(asChar(alleleFile))));
}

SEXP SXP_load_new_effects(SEXP exd, SEXP effectFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	load_effects_to_simdata(d, CHAR(asChar(effectFile)));
	return ScalarInteger(0);
}


/*-------------------------------- Crossers ---------------------------*/
GenOptions create_genoptions(SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions go = BASIC_OPT;
	int b;
	
	b = asLogical(name);
	if (b == NA_LOGICAL) { error("`name` parameter is of invalid type.\n"); }
	go.will_name_offspring = b;
	
	go.offspring_name_prefix = CHAR(asChar(namePrefix));
	
	b = asInteger(familySize);
	if (b == NA_INTEGER) { error("`offspring` parameter is of invalid type.\n"); }
	go.family_size = b;
	
	b = asLogical(trackPedigree);
	if (b == NA_LOGICAL) { error("`track.pedigree` parameter is of invalid type.\n"); }
	go.will_track_pedigree = b;
	
	b = asLogical(giveIds);
	if (b == NA_LOGICAL) { error("`give.ids` parameter is of invalid type.\n"); }
	go.will_allocate_ids = b;
	
	go.filename_prefix = CHAR(asChar(filePrefix));
	
	b = asLogical(savePedigree);
	if (b == NA_LOGICAL) { error("`save.pedigree` parameter is of invalid type.\n"); }
	go.will_save_pedigree_to_file = b;
	b = asLogical(saveEffects);
	if (b == NA_LOGICAL) { error("`save.gebv` parameter is of invalid type.\n"); }
	go.will_save_bvs_to_file = b;
	b = asLogical(saveGenes);
	if (b == NA_LOGICAL) { error("`save.genotype` parameter is of invalid type.\n"); }
	go.will_save_alleles_to_file = b;
	
	b = asLogical(retain);
	if (b == NA_LOGICAL) { error("`retain` parameter is of invalid type.\n"); }
	go.will_save_to_simdata = b;
	
	return go;
}

SEXP SXP_cross_randomly(SEXP exd, SEXP glen, SEXP groups, SEXP crosses, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);

	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` is invalid.\n");
	}
	for (int i = 0; i < len; ++i) {
		if (gps[i] == NA_INTEGER || gps[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}
	
	int n = asInteger(crosses);
	if (n < 0 || n == NA_INTEGER) { error("`n.crosses` parameter is invalid.\n"); }
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (len == 1) {
		return ScalarInteger(cross_random_individuals(d, gps[0], n, g));
		
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = cross_random_individuals(d, gps[i], n, g);
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_cross_Rcombinations(SEXP exd, SEXP firstparents, SEXP secondparents,
		SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	if (length(firstparents) != length(secondparents)) {
		error("Parent vectors must be the same length.\n");
	}
	
	int ncrosses = length(firstparents);
	
	int combinations[2][ncrosses];
	char pname[100];
	
	if (TYPEOF(firstparents) == STRSXP) {	
		for (int i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(firstparents, i)), sizeof(char)*100);
			combinations[0][i] = get_index_of_name(d->m, pname);
		}
	} else if (TYPEOF(firstparents) == INTSXP) {
		int* indexes = INTEGER(firstparents);
		for (int i = 0; i < ncrosses; ++i) {
			combinations[0][i] = indexes[i];
		}
	} else {
		error("first.parents must be a vector of strings or integers.\n");
	}
	
	if (TYPEOF(secondparents) == STRSXP) {
		for (int i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(secondparents, i)), sizeof(char)*100);
			combinations[1][i] = get_index_of_name(d->m, pname);
		}
		
	} else if (TYPEOF(secondparents) == INTSXP) {
		int* indexes = INTEGER(secondparents);
		for (int i = 0; i < ncrosses; ++i) {
			combinations[1][i] = indexes[i];
		}
	} else {
		error("second.parents must be a vector of strings or integers.\n");
	}

	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
								 giveIds, filePrefix, savePedigree, saveEffects,
								 saveGenes, retain);

	return ScalarInteger(cross_these_combinations(d, ncrosses, combinations, g));
	
}

SEXP SXP_cross_combinations(SEXP exd, SEXP filename, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char* fname = CHAR(asChar(filename));
	
	return ScalarInteger(make_crosses_from_file(d, fname, g));
}

SEXP SXP_dcross_combinations(SEXP exd, SEXP filename, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
									 
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char* fname = CHAR(asChar(filename));
	
	return ScalarInteger(make_double_crosses_from_file(d, fname, g));
	
}

SEXP SXP_cross_unidirectional(SEXP exd, SEXP glen, SEXP groups, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` vector is invalid.\n");
	}
	for (int i = 0; i < len; ++i) {
		if (gps[i] == NA_INTEGER || gps[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	if (len == 1) {
		return ScalarInteger(make_all_unidirectional_crosses(d, gps[0], g));
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = make_all_unidirectional_crosses(d, gps[i], g);
		}
		UNPROTECT(1);
		return out;
	}
}

/*SEXP cross_top(SEXP exd, SEXP group, SEXP percent, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
	int grp = asInteger(group);
	if (grp < 0 || grp == NA_INTEGER) { error("`group` parameter is invalid.\n"); }
	
	double cpercent = asReal(percent);
	if (cpercent < 0 || cpercent > 100) { error("`threshold` parameter is invalid.\n"); }
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
	return ScalarInteger(make_crosses_from_top_m_percent(d, cpercent, grp, g)); 
}*/

SEXP SXP_selfing(SEXP exd, SEXP glen, SEXP groups, SEXP n, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` vector is invalid.\n");
	}
	for (int i = 0; i < len; ++i) {
		if (gps[i] == NA_INTEGER || gps[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}
	
	int cn = asInteger(n);
	if (cn < 0 || cn == NA_INTEGER) { error("`n` parameter is invalid.\n"); }
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	if (len == 1) {
		return ScalarInteger(self_n_times(d, cn, gps[0], g));
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = self_n_times(d, cn, gps[i], g);
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_doubled(SEXP exd, SEXP glen, SEXP groups, SEXP name, SEXP namePrefix, SEXP familySize,
		SEXP trackPedigree, SEXP giveIds, SEXP filePrefix, SEXP savePedigree,
		SEXP saveEffects, SEXP saveGenes, SEXP retain) {
	GenOptions g = create_genoptions(name, namePrefix, familySize, trackPedigree,
									 giveIds, filePrefix, savePedigree, saveEffects,
									 saveGenes, retain);
	int len = asInteger(glen);
	int *gps = INTEGER(groups); 
	if (len == NA_INTEGER) { 
		error("`groups` vector is invalid.\n");
	}
	for (int i = 0; i < len; ++i) {
		if (gps[i] == NA_INTEGER || gps[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	if (len == 1) {
		return ScalarInteger(make_doubled_haploids(d, gps[0], g));
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = make_doubled_haploids(d, gps[i], g);
		}
		UNPROTECT(1);
		return out;
	}
}


/*-----------------------------------Groups----------------------------------*/
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


/*------------------------------------Fitness------------------------------*/
SEXP SXP_group_eval(SEXP exd, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
	int group_id = asInteger(group);
	if (group_id == NA_INTEGER || group_id < 0) { 
		error("`group` parameter is of invalid type.\n");
	}
	
	int group_size = get_group_size(d, group_id);
	unsigned int* inds = get_group_indexes(d, group_id, group_size);
	DecimalMatrix gebvs = calculate_group_bvs(d, group_id);
	
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
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
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
		return ScalarInteger(split_by_bv(d, gps[0], num_to_select, want_low));
	} else {
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = split_by_bv(d, gps[i], num_to_select, want_low);
		}
		UNPROTECT(1);
		return out;
	}
	
}

SEXP SXP_simple_selection_bypercent(SEXP exd, SEXP glen, SEXP groups, SEXP percent, SEXP bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
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
		
		return ScalarInteger(split_by_bv(d, gps[0], num_to_select, want_low));
	} else {
		int num_to_select;
		
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			num_to_select = get_group_size(d, gps[i]) * pc_to_select / 100;
			outc[i] = split_by_bv(d, gps[i], num_to_select, want_low);
		}
		UNPROTECT(1);
		return out;
	}
	
}


/*----------------Data viewers---------------------*/

SEXP SXP_get_best_haplotype(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	char* best_genotype = calculate_optimal_alleles(d);
	
	SEXP out = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(out, 0, mkChar(best_genotype));
	free(best_genotype);
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_best_GEBV(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	double best_GEBV = calculate_optimum_bv(d);
	
	SEXP out = PROTECT(allocVector(REALSXP, 1));
	REAL(out)[0] = best_GEBV;
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_worst_GEBV(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	double worst_GEBV = calculate_minimum_bv(d);
	
	SEXP out = PROTECT(allocVector(REALSXP, 1));
	REAL(out)[0] = worst_GEBV;
	UNPROTECT(1);
	return out;
}

SEXP SXP_find_crossovers(SEXP exd, SEXP parentFile, SEXP outFile, SEXP windowSize, SEXP certainty) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char* load_fname = CHAR(asChar(parentFile));
	const char* save_fname = CHAR(asChar(outFile));
	
	int win = asInteger(windowSize);
	if (win < 0 || win == NA_INTEGER) { error("`window.size` parameter is invalid.\n"); }
	
	int cert = asLogical(certainty);
	if (cert == NA_LOGICAL) { error("`certainty` parameter is invalid.\n"); }
	
	return ScalarInteger(calculate_recombinations_from_file(d, load_fname, save_fname, win, cert));
}

SEXP SXP_send_map(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	MarkerPosition* mp = d->map.positions;
	
	SEXP map = PROTECT(allocVector(VECSXP, 3));
	SEXP snp = PROTECT(allocVector(STRSXP, d->n_markers));
	SEXP chr = PROTECT(allocVector(REALSXP, d->n_markers));
	SEXP pos = PROTECT(allocVector(REALSXP, d->n_markers));
	double* cchr = REAL(chr);
	double* cpos = REAL(pos);
	
	for (int i = 0; i < d->n_markers; ++i) {
		SET_STRING_ELT(snp, i, mkChar(d->markers[i]));
		cchr[i] = mp[i].chromosome;
		cpos[i] = mp[i].position;
	}
	
	SET_VECTOR_ELT(map, 0, snp);
	SET_VECTOR_ELT(map, 1, chr);
	SET_VECTOR_ELT(map, 2, pos);
	UNPROTECT(4);
	return map;
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
	char c2;
	int group_size = get_group_size(d, group_id);
	
	if (c == 'D') {
		unsigned int *data = get_group_ids(d, group_id, group_size);
		
		// save to an R vector
		SEXP out = PROTECT(allocVector(INTSXP, group_size));
		int* outc = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'X') {
		unsigned int *data = get_group_indexes(d, group_id, group_size);
		
		// save to an R vector
		SEXP out = PROTECT(allocVector(INTSXP, group_size));
		int* outc = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'B') {
		double* data = get_group_bvs(d, group_id, group_size);
		
		SEXP out = PROTECT(allocVector(REALSXP, group_size));
		double* outc = REAL(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'N') {
		char** data = get_group_names(d, group_id, group_size);
		unsigned int* backupdata = get_group_ids(d, group_id, group_size);
		char buffer[get_integer_digits(backupdata[group_size - 1]) + 1]; 
	
		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			if (data[i] != NULL) {
				SET_STRING_ELT(out, i, mkChar(data[i]));
			} else {
				sprintf(buffer, "%d", backupdata[i]);
				SET_STRING_ELT(out, i, mkChar(buffer));
			}
		}
		free(data);
		free(backupdata);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'G') {
		char** rawdata = get_group_genes(d, group_id, group_size);
		char** data;
			
		// copy over to another array, adding string terminators
		int glen = d->n_markers * 2; // genotype length
		data = get_malloc(sizeof(char*) * group_size);
		for (int i = 0; i < group_size; ++i) {
			data[i] = get_malloc(sizeof(char) * (glen + 1));
			for (int j = 0; j < glen; ++j) {
				data[i][j] = rawdata[i][j];
			}
			data[i][glen] = '\0';
		}
		
		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			SET_STRING_ELT(out, i, mkChar(data[i]));
			free(data[i]);
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'P' && (c2 = CHAR(asChar(whatData))[1]) == 'E' && CHAR(asChar(whatData))[2] == 'D') {
		char** data = get_group_pedigrees(d, group_id, group_size);
		
		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			SET_STRING_ELT(out, i, mkChar(data[i]));
			free(data[i]);
		}
		free(data);
		UNPROTECT(1);
		return out;
		
	} else if (c == 'P' && (c2 == '1' || c2 == '2')) {
		int parent = c2 == '1' ? 1 : 2;
		char** data = get_group_parent_names(d, group_id, group_size, parent);
		unsigned int* backupdata = get_group_parent_ids(d, group_id, group_size, parent);
		char buffer[get_integer_digits(backupdata[group_size - 1]) + 1]; 
		
		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			if (data[i] != NULL) {
				SET_STRING_ELT(out, i, mkChar(data[i]));
			} else {
				sprintf(buffer, "%d", backupdata[i]);
				SET_STRING_ELT(out, i, mkChar(buffer));
			}
		}
		free(data);
		free(backupdata);
		UNPROTECT(1);
		return out;
	
	} else {
		error("`data.type` parameter is not a valid option.");
	}
}

/*--------------------------------Printing-----------------------------------*/

SEXP SXP_save_simdata(SEXP exd, SEXP filename) {
	FILE* f;
	const char* fname = CHAR(asChar(filename));
	if ((f = fopen(fname, "w")) == NULL) {
		error( "Failed to open file %s.\n", fname);
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	save_simdata(f, d);
	
	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_genotypes(SEXP exd, SEXP filename, SEXP group, SEXP type) {
	FILE* f;
	const char* fname = CHAR(asChar(filename));
	if ((f = fopen(fname, "w")) == NULL) {
		error( "Failed to open file %s.\n", fname);
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char t = CHAR(asChar(type))[0];	
	if (t == 'R' || t == 'r') {
		if (isNull(group)) {
			save_allele_matrix(f, d->m, d->markers);
		} else if (asInteger(group) >= 0) {
			save_group_alleles(f, d, asInteger(group));
		} else {
			fclose(f);
			error("Supplied group number is invalid.");
		}
	} else if (t == 'T' || t == 't') {
		if (isNull(group)) {
			save_transposed_allele_matrix(f, d->m, d->markers);
		} else if (asInteger(group) >= 0) {
			save_transposed_group_alleles(f, d, asInteger(group));
		} else {
			fclose(f);
			error("Supplied group number is invalid.");
		}
	} else {
		fclose(f);
		error("Supplied printing format is invalid.");
	}
	
	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_counts(SEXP exd, SEXP filename, SEXP group, SEXP allele) {
	FILE* f;
	const char* fname = CHAR(asChar(filename));
	if ((f = fopen(fname, "w")) == NULL) {
		error( "Failed to open file %s.\n", fname);
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char t = CHAR(asChar(allele))[0];	
	if (isNull(group)) {
		save_count_matrix(f, d, t);
	} else if (asInteger(group) >= 0) {
		save_count_matrix_of_group(f, d, t, asInteger(group));
	} else {
		fclose(f);
		error("Supplied group number is invalid.");
	}
	
	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_pedigrees(SEXP exd, SEXP filename, SEXP group, SEXP type) {
	FILE* f;
	const char* fname = CHAR(asChar(filename));
	if ((f = fopen(fname, "w")) == NULL) {
		error( "Failed to open file %s.\n", fname);
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char t = CHAR(asChar(type))[0];	
	if (t == 'R' || t == 'r') { // full/recursive
		if (isNull(group)) {
			save_full_pedigree(f, d);
		} else if (asInteger(group) >= 0) {
			save_group_full_pedigree(f, d, asInteger(group));
		} else {
			error("Supplied group number is invalid.");
		}
	} else if (t == 'P' || t == 'p') { // one-step/parents
		if (isNull(group)) {
			save_one_step_pedigree(f, d);
		} else if (asInteger(group) >= 0) {
			save_group_one_step_pedigree(f, d, asInteger(group));
		} else {
			error("Supplied group number is invalid.");
		}
	} else {
		fclose(f);
		error("Supplied printing format is invalid.");
	}
	
	fclose(f);
	return ScalarInteger(0);	
}

SEXP SXP_save_GEBVs(SEXP exd, SEXP filename, SEXP group) {
	FILE* f;
	const char* fname = CHAR(asChar(filename));
	if ((f = fopen(fname, "w")) == NULL) {
		error( "Failed to open file %s.\n", fname);
	}
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
	if (isNull(group)) {
		save_bvs(f, d);
	} else if (asInteger(group) >= 0) {
		save_group_bvs(f, d, asInteger(group));
	} else {
		fclose(f);
		error("Supplied group number is invalid.");
	}
	
	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_file_block_effects(SEXP exd, SEXP filename, SEXP block_file, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
	if (isNull(group)) {
		MarkerBlocks b = read_block_file(d, CHAR(asChar(block_file)));
		calculate_local_bvs(d, b, CHAR(asChar(filename)));
		delete_markerblocks(&b);
	} else if (asInteger(group) > 0) {
		MarkerBlocks b = read_block_file(d, CHAR(asChar(block_file)));
		calculate_group_local_bvs(d, b, CHAR(asChar(filename)), asInteger(group));
		delete_markerblocks(&b);
	} else {
		error("Supplied group number is invalid.\n");
	}
	
	return ScalarInteger(0);	
}

SEXP SXP_save_chrsplit_block_effects(SEXP exd, SEXP filename, SEXP nslices, SEXP group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->e.effects.matrix == NULL) { error("Need to load effect values before running this function.\n"); } 
	
	if (isNull(group)) {
		MarkerBlocks b = create_n_blocks_by_chr(d, asInteger(nslices));
		calculate_local_bvs(d, b, CHAR(asChar(filename)));
		delete_markerblocks(&b);
	} else if (asInteger(group) > 0) {
		MarkerBlocks b = create_n_blocks_by_chr(d, asInteger(nslices));
		calculate_group_local_bvs(d, b, CHAR(asChar(filename)), asInteger(group));
		delete_markerblocks(&b);
	} else {
		error("Supplied group number is invalid.\n");
	}
	
	return ScalarInteger(0);	
}


/*--------------------------------Deletors------------------------------------*/
void SXP_delete_simdata(SEXP sd) {
	//Rprintf("Garbage collecting SimData...\n");
	SimData* d = (SimData*) R_ExternalPtrAddr(sd);
	delete_simdata(d);
	R_ClearExternalPtr(sd);
	return;
}

SEXP clear_simdata(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	delete_simdata(d);
	return ScalarInteger(0);
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
