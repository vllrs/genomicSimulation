#include "r-wrappers.h"

/*--------------------------SEXP HELPERS----------------------*/


#define GROUPNUM_IFY(n) (GroupNum){.num=n}
#define EFFECTID_IFY(n) (EffectID){.id=n}
#define LABELID_IFY(n) (LabelID){.id=n}
#define PEDIGREEID_IFY(n) (PedigreeID){.id=n}
//const enum IDENTIFIERS = {GroupNum, PedigreeID, EffectID, LabelID};

void convertVECSXP_to_GroupNum(SEXP container, GroupNum* output) {
	int len = length(container);
	int *intvec = INTEGER(container);

	for (int i = 0; i < len; ++i) {
		if (intvec[i] <= 0) {
			warning("%i is not a possible group number", intvec[i]);
			output[i] = NO_GROUP;
		} else {
			output[i] = GROUPNUM_IFY(intvec[i]);
		}
	}
}

/*
struct idVector {
	size_t len;
	unsigned int* vec;
};

struct idVector convertSEXP_to_identifier(SEXP container, int identifier) {

	switch(identifier) {
		case IDENTIFIERS.GroupNum:
			break;
		case IDENTIFIERS.PedigreeID:
			break;
		case IDENTIFIERS.EffectID:
			break;
		case IDENTIFIERS.LabelID:
			break;
		default:
			error("Unknown identifier type requested in helper function")
	}

	// Possibilities:
}*/


/*unsigned int convertSEXP_to_GroupNum(SEXP container, const char* identifier) {
	int val;
	if (isNewList(container)) {
		val = INTEGER(getListElement(container, identifier))[0];
	} else {
		val = asInteger(container);
	}

	if (val == NA_INTEGER) {
		error("group number parameter is of invalid type.\n");
	}
	return val;
}*/


/* get the list element named str, or return NULL
 *
 * From https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Handling-lists
 * Use: double g = REAL(getListElement(a, "g"))[0];
 * for a list created in R as `a <- list(f = 1, g = 2, h = 3)`
 * */
/*SEXP LISTXP_getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);


    for (int i = 0; i < length(list); i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
        }
    return elmt;
}*/


/*-------------------------- Loaders -------------------------*/

SEXP SXP_load_data(SEXP s_alleleFile, SEXP s_mapFile) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(s_alleleFile)));
	load_genmap_to_simdata(d, CHAR(asChar(s_mapFile)));

	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);

	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP SXP_load_data_weff(SEXP s_alleleFile, SEXP s_mapFile, SEXP s_effectFile) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(s_alleleFile)));
	load_genmap_to_simdata(d, CHAR(asChar(s_mapFile)));
	load_effects_to_simdata(d, CHAR(asChar(s_effectFile)));

	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);

	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP SXP_load_more_genotypes(SEXP exd, SEXP s_alleleFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_more_transposed_genes_to_simdata(d, CHAR(asChar(s_alleleFile))).num);
}

SEXP SXP_load_new_effects(SEXP exd, SEXP s_effectFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_effects_to_simdata(d, CHAR(asChar(s_effectFile))).id);
}


/*-------------------------------- Crossers ---------------------------*/
GenOptions create_genoptions(SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions go = BASIC_OPT;
	int b;

	b = asLogical(s_name);
	if (b == NA_LOGICAL) { error("`name` parameter is of invalid type.\n"); }
	go.will_name_offspring = b;

	go.offspring_name_prefix = CHAR(asChar(s_namePrefix));

	b = asInteger(s_familySize);
	if (b == NA_INTEGER) { error("`offspring` parameter is of invalid type.\n"); }
	go.family_size = b;

	b = asLogical(s_trackPedigree);
	if (b == NA_LOGICAL) { error("`track.pedigree` parameter is of invalid type.\n"); }
	go.will_track_pedigree = b;

	b = asLogical(s_giveIds);
	if (b == NA_LOGICAL) { error("`give.ids` parameter is of invalid type.\n"); }
	go.will_allocate_ids = b;

	go.filename_prefix = CHAR(asChar(s_filePrefix));

	b = asLogical(s_savePedigree);
	if (b == NA_LOGICAL) { error("`save.pedigree` parameter is of invalid type.\n"); }
	go.will_save_pedigree_to_file = b;
	b = asLogical(s_saveEffects);
	if (b == NA_LOGICAL) { error("`save.gebv` parameter is of invalid type.\n"); }
	go.will_save_bvs_to_file = EFFECTID_IFY(b);
	b = asLogical(s_saveGenes);
	if (b == NA_LOGICAL) { error("`save.genotype` parameter is of invalid type.\n"); }
	go.will_save_alleles_to_file = b;

	b = asLogical(s_retain);
	if (b == NA_LOGICAL) { error("`retain` parameter is of invalid type.\n"); }
	go.will_save_to_simdata = b;

	return go;
}

SEXP SXP_cross_randomly(SEXP exd, SEXP s_groups, SEXP s_crosses, SEXP s_cap,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	int n = asInteger(s_crosses);
	if (n == NA_INTEGER) { error("`n.crosses` parameter is invalid.\n"); }

	int cap = asInteger(s_cap);
	if (cap == NA_INTEGER) { error("`cap` parameter is invalid.\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(cross_random_individuals(d, GROUPNUM_IFY(groups[0]), n, cap, g).num);

	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (int i = 0; i < glen; ++i) {
			outc[i] = cross_random_individuals(d, GROUPNUM_IFY(groups[i]), n, cap, g).num;
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_cross_randomly_btwn(SEXP exd, SEXP s_group1, SEXP s_group2, SEXP s_cap1, SEXP s_cap2,
		SEXP s_crosses, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	int group1_c = INTEGER(s_group1)[0];
	if (group1_c == NA_INTEGER || group1_c < 0) { error("The parameter `group1` is invalid.\n"); }
	int group2_c = INTEGER(s_group2)[0];
	if (group2_c == NA_INTEGER || group2_c < 0) { error("The parameter `group2` is invalid.\n"); }

	int cap1 = asInteger(s_cap1);
	if (cap1 == NA_INTEGER) { error("The parameter `cap1` is invalid.\n"); }
	int cap2 = asInteger(s_cap2);
	if (cap2 == NA_INTEGER) { error("The parameter `cap2` is invalid.\n"); }

	int n = asInteger(s_crosses);
	if (n == NA_INTEGER) { error("`n.crosses` parameter is invalid.\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	return ScalarInteger(cross_randomly_between(d, GROUPNUM_IFY(group1_c), GROUPNUM_IFY(group2_c), n, cap1, cap2, g).num);
}

SEXP SXP_cross_Rcombinations(SEXP exd, SEXP s_firstparents, SEXP s_secondparents,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (length(s_firstparents) != length(s_secondparents)) {
		error("Parent vectors must be the same length.\n");
	}

	int ncrosses = length(s_firstparents);

	int combinations[2][ncrosses];
	char pname[NAME_LENGTH+1];

	if (TYPEOF(s_firstparents) == STRSXP) {
		for (int i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(s_firstparents, i)), sizeof(char)*NAME_LENGTH);
			combinations[0][i] = get_index_of_name(d->m, pname);
		}
	} else if (TYPEOF(s_firstparents) == INTSXP) {
		int* indexes = INTEGER(s_firstparents);
		for (int i = 0; i < ncrosses; ++i) {
			combinations[0][i] = indexes[i];
		}
	} else {
		error("first.parents must be a vector of strings or integers.\n");
	}

	if (TYPEOF(s_secondparents) == STRSXP) {
		for (int i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(s_secondparents, i)), sizeof(char)*NAME_LENGTH);
			combinations[1][i] = get_index_of_name(d->m, pname);
		}

	} else if (TYPEOF(s_secondparents) == INTSXP) {
		int* indexes = INTEGER(s_secondparents);
		for (int i = 0; i < ncrosses; ++i) {
			combinations[1][i] = indexes[i];
		}
	} else {
		error("second.parents must be a vector of strings or integers.\n");
	}

	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
								 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
								 s_saveGenes, s_retain);

	return ScalarInteger(cross_these_combinations(d, ncrosses, combinations[0], combinations[1], g).num);

}

SEXP SXP_cross_combinations(SEXP exd, SEXP s_filename, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char* filename = CHAR(asChar(s_filename));

	return ScalarInteger(make_crosses_from_file(d, filename, g).num);
}

SEXP SXP_dcross_combinations(SEXP exd, SEXP s_filename, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char* filename = CHAR(asChar(s_filename));

	return ScalarInteger(make_double_crosses_from_file(d, filename, g).num);

}

SEXP SXP_cross_unidirectional(SEXP exd, SEXP s_groups, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_all_unidirectional_crosses(d, GROUPNUM_IFY(groups[0]), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (int i = 0; i < glen; ++i) {
			outc[i] = make_all_unidirectional_crosses(d, GROUPNUM_IFY(groups[i]), g).num;
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_selfing(SEXP exd, SEXP s_groups, SEXP s_ngen, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	int cn = asInteger(s_ngen);
	if (cn < 0 || cn == NA_INTEGER) { error("`n` parameter is invalid.\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(self_n_times(d, cn, GROUPNUM_IFY(groups[0]), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (int i = 0; i < glen; ++i) {
			outc[i] = self_n_times(d, cn, GROUPNUM_IFY(groups[i]), g).num;
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_doubled(SEXP exd, SEXP s_groups, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_doubled_haploids(d, GROUPNUM_IFY(groups[0]), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (int i = 0; i < glen; ++i) {
			outc[i] = make_doubled_haploids(d, GROUPNUM_IFY(groups[i]), g).num;
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_clone(SEXP exd, SEXP s_groups, SEXP s_inherit_name, SEXP s_name,
		SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	int inherit_name = asLogical(s_inherit_name);
	if (inherit_name == NA_LOGICAL) { error("`inherit.name` parameter is of invalid type: should be logical.\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_clones(d, GROUPNUM_IFY(groups[0]), inherit_name, g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (int i = 0; i < glen; ++i) {
			outc[i] = make_clones(d, GROUPNUM_IFY(groups[i]), inherit_name, g).num;
		}
		UNPROTECT(1);
		return out;
	}
}

/*-----------------------------------Labels----------------------------------*/
SEXP SXP_create_label(SEXP exd, SEXP s_default) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int lblDefault = asInteger(s_default);
	if (lblDefault == NA_INTEGER) {
		error("`default` parameter is of invalid type: must be an integer.\n");
	}

	return ScalarInteger(create_new_label(d, lblDefault).id);
}

SEXP SXP_change_label_default(SEXP exd, SEXP s_labels, SEXP s_defaults) {

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int lblen = length(s_labels);
	int *labels = INTEGER(s_labels);

	int dlen = length(s_defaults);
	int *defaults = INTEGER(s_defaults);

	// Find the minimum length of the vectors and discard everything further along
	int matchedLen = lblen < dlen ? lblen : dlen;

	for (int i = 0; i < matchedLen; ++i) {
		if (labels[i] < 1) {
			error("entry in `label` vector is invalid: too small or large to be a label.");
		} else {
			set_label_default(d, LABELID_IFY(labels[i]), defaults[i]);
		}
	}
	Rprintf("Set the defaults of %i labels.\n", matchedLen);
	return ScalarInteger(0);
}

SEXP SXP_change_label_amount(SEXP exd, SEXP s_label, SEXP s_incr, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int label = asInteger(s_label);
	if (label == NA_INTEGER) {
		error("`label` parameter is invalid: must be an integer");
	} else if (label < 0) {
		error("`label` parameter is invalid: negative.");
	}

	int incr = asInteger(s_incr);
	if (incr == NA_INTEGER) {
		error("`amount` parameter is invalid: must be an integer");
	}

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);

	if (len == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}
		if (groups[0] < 0) {
			Rprintf("the `group` vector is invalid.");
		} else {
			increment_labels(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), incr);
		}

	} else {
		for (int i = 0; i < len; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %i in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				increment_labels(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), incr);
			}
		}
	}

	return ScalarInteger(0);
}

SEXP SXP_change_label_const(SEXP exd, SEXP s_label, SEXP s_const, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int label = asInteger(s_label);
	if (label == NA_INTEGER) {
		error("`label` parameter is invalid: must be an integer");
	} else if (label < 0) {
		error("`label` parameter is invalid: negative.");
	}

	int num = asInteger(s_const);
	if (num == NA_INTEGER) {
		error("`value` parameter is invalid: must be an integer");
	}

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);

	if (len == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			Rprintf("the `group` vector is invalid.");
		} else {
			set_labels_to_const(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), num);
		}

	} else {
		for (int i = 0; i < len; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %i in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				set_labels_to_const(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), num);
			}
		}
	}

	return ScalarInteger(0);
}

SEXP SXP_change_label_values(SEXP exd, SEXP s_label, SEXP s_values, SEXP s_group, SEXP s_start) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int label = asInteger(s_label);
	if (label == NA_INTEGER) {
		error("`label` parameter is invalid: must be an integer");
	} else if (label < 0) {
		error("`label` parameter is invalid: negative.");
	}

	int vlen = length(s_values);
	int* values = INTEGER(s_values);

	int group = asInteger(s_group);
	if (group == NA_INTEGER) {
		group = 0;
	} else if (group < 0) {
		error("`group` is invalid (wrong type or negative).\n");
	}

	int startIndex = asInteger(s_start) - 1;
	if (startIndex == NA_INTEGER || startIndex < 0) {
		error("`startIndex` is invalid (wrong type or negative).\n");
	}

	set_labels_to_values(d, GROUPNUM_IFY(group), startIndex, LABELID_IFY(label), vlen, values);

	return ScalarInteger(0);
}

/*-----------------------------------Groups----------------------------------*/
SEXP SXP_combine_groups(SEXP exd, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int len = length(s_groups);

	GroupNum groups[len];
	convertVECSXP_to_GroupNum(s_groups, groups);

	return ScalarInteger(combine_groups(d, len, groups).num);
}

SEXP SXP_split_individuals(SEXP exd, SEXP s_group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	GroupNum results[group_size];

	// do the actual split
	int nsplit = split_into_individuals(d, GROUPNUM_IFY(group_id), group_size, results);

	// Get an R vector of the same length as the number of new size 1 groups created
	SEXP out = PROTECT(allocVector(INTSXP, nsplit));
	int* outc = INTEGER(out);
	for (int i = 0; i < nsplit; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_split_familywise(SEXP exd, SEXP s_group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	GroupNum results[group_size];

	// do the actual split
	int groups_created = split_into_families(d, GROUPNUM_IFY(group_id), group_size, results);

	// Get an R vector of the same length as the number of new groups created
	SEXP out = PROTECT(allocVector(INTSXP, groups_created));
	int* outc = INTEGER(out);
	for (int i = 0; i < groups_created; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_split_halfsibwise(SEXP exd, SEXP s_group, SEXP s_parent) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int parent_num = asInteger(s_parent);
	if (parent_num == NA_INTEGER || parent_num < 0) {
		error("`parent` parameter is of invalid type.\n");
	}

	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	GroupNum results[group_size];

	// do the actual split
	int groups_created = split_into_halfsib_families(d, GROUPNUM_IFY(group_id), 
                                                  parent_num, group_size, results);

	// Get an R vector of the same length as the number of newgroups created
	SEXP out = PROTECT(allocVector(INTSXP, groups_created));
	int* outc = INTEGER(out);
	for (int i = 0; i < groups_created; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_split_randomly(SEXP exd, SEXP s_group, SEXP s_n) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int n_groups = asInteger(s_n);
	if (n_groups == NA_INTEGER || n_groups < 0) {
		error("`n` parameter is of invalid type.\n");
	}

	if (n_groups == 1) {
		return s_group;

	} else if (n_groups == 2) {
		SEXP out = PROTECT(allocVector(INTSXP, 2));
		int* outc = INTEGER(out);
		outc[0] = group_id;
		outc[1] = split_randomly_into_two(d, GROUPNUM_IFY(group_id)).num;

		UNPROTECT(1);
		return out;

	} else {
		GroupNum results[n_groups];
		n_groups = split_randomly_into_n(d, GROUPNUM_IFY(group_id), n_groups, results);

		SEXP out = PROTECT(allocVector(INTSXP, n_groups));
		int* outc = INTEGER(out);
		for (int i = 0; i < n_groups; ++i) {
			outc[i] = results[i].num;
		}

		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_split_evenly(SEXP exd, SEXP s_group, SEXP s_n) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int n_groups = asInteger(s_n);
	if (n_groups == NA_INTEGER || n_groups < 0) {
		error("`n` parameter is of invalid type.\n");
	}

	if (n_groups == 1) {
		return s_group;

	} else if (n_groups == 2) {
		SEXP out = PROTECT(allocVector(INTSXP, 2));
		int* outc = INTEGER(out);
		outc[0] = group_id;
		outc[1] = split_evenly_into_two(d, GROUPNUM_IFY(group_id)).num;

		UNPROTECT(1);
		return out;

	} else {
		GroupNum results[n_groups];
		split_evenly_into_n(d, GROUPNUM_IFY(group_id), n_groups, results);

		SEXP out = PROTECT(allocVector(INTSXP, n_groups));
		int* outc = INTEGER(out);
		for (int i = 0; i < n_groups; ++i) {
			outc[i] = results[i].num;
		}

		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_split_buckets(SEXP exd, SEXP s_group, SEXP s_buckets) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	if (TYPEOF(s_buckets) != INTSXP) { // check param name
		error("`buckets` parameter must be a vector of integers.\n");
	}

	int n_groups = length(s_buckets)+1;
	int* counts = INTEGER(s_buckets);

	GroupNum results[n_groups];
	split_by_specific_counts_into_n(d, GROUPNUM_IFY(group_id), n_groups, counts, results);

	SEXP out = PROTECT(allocVector(INTSXP, n_groups));
	int* outc = INTEGER(out);
	for (int i = 0; i < n_groups; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_split_probabilities(SEXP exd, SEXP s_group, SEXP s_probs) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	if (TYPEOF(s_probs) != REALSXP) { // check param name
		error("`probabilies` parameter must be a vector of decimals.\n");
	}

	int n_groups = length(s_probs)+1;
	double* probabilities = REAL(s_probs);

	GroupNum results[n_groups];
	memset(results, 0, sizeof(GroupNum) * n_groups);
	split_by_probabilities_into_n(d, GROUPNUM_IFY(group_id), n_groups, probabilities, results);

	SEXP out = PROTECT(allocVector(INTSXP, n_groups));
	int* outc = INTEGER(out);
	for (int i = 0; i < n_groups; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}


SEXP SXP_split_out(SEXP exd, SEXP s_indexes) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int n = length(s_indexes);
	unsigned int *ns = (unsigned int*) INTEGER(s_indexes);
	for (int i = 0; i < n; ++i) {
		if (ns[i] == NA_INTEGER || ns[i] < 0) {
			error("The `indexes` vector contains at least one invalid index.\n");
		}
	}

	return ScalarInteger(split_from_group(d, n, ns).num);
}


SEXP SXP_split_by_label_value(SEXP exd, SEXP s_label, SEXP s_value, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int label = asInteger(s_label);
	if (label == NA_INTEGER) {
		error("`label` parameter is invalid: must be an integer");
	} else if (label < 0) {
		error("`label` parameter is invalid: negative.");
	}

	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);

	int value = asInteger(s_value);
	if (value == NA_INTEGER) {
		error("`value` input must be an integer.\n");
	}

	if (glen == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			error("the `group` input is invalid (negative).");
		} else {
			return ScalarInteger(split_by_label_value(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), value).num);
		}

	} else {

		GroupNum newSubGroups[glen];
		int j = 0;
		for (int i = 0; i < glen; ++i) {
			if (groups[i] < 1) {
				Rprintf("Entry %i in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				newSubGroups[j] = split_by_label_value(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), value);
				++j;
			}
		}
		return ScalarInteger(combine_groups(d, j, newSubGroups).num);
	}
}


SEXP SXP_split_by_label_range(SEXP exd, SEXP s_label, SEXP s_lowbound, SEXP s_highbound,
		SEXP s_groups) {

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int label = asInteger(s_label);
	if (label == NA_INTEGER) {
		error("`label` parameter is invalid: must be an integer");
	} else if (label < 0) {
		error("`label` parameter is invalid: negative.");
	}

	int glen = length(s_groups);
	int *groups = INTEGER(s_groups);

	int low = asInteger(s_lowbound);
	if (low == NA_INTEGER) {
		error("`rangeLowEnd` input must be an integer.\n");
	}

	int high = asInteger(s_highbound);
	if (high == NA_INTEGER) {
		error("`rangeHighEnd` input must be an integer.\n");
	}

	if (glen == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			error("the `group` input is invalid (negative).");
		} else {
			return ScalarInteger(split_by_label_range(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), low, high).num);
		}

	} else {

		GroupNum newSubGroups[glen];
		int j = 0;
		for (int i = 0; i < glen; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %i in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				newSubGroups[j] = split_by_label_range(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), low, high);
				++j;
			}
		}
		return ScalarInteger(combine_groups(d, j, newSubGroups).num);
	}

}

/*------------------------------------Fitness------------------------------*/
SEXP SXP_group_eval(SEXP exd, SEXP s_group, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	unsigned int* inds = get_malloc(sizeof(unsigned int) * group_size);
	get_group_indexes(d, GROUPNUM_IFY(group_id), group_size, inds);
	DecimalMatrix gebvs = calculate_group_bvs(d, GROUPNUM_IFY(group_id), EFFECTID_IFY(eff_id));

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

SEXP SXP_simple_selection(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_number, SEXP s_bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	int num_to_select = asInteger(s_number);
	if (num_to_select == NA_INTEGER || num_to_select < 0) {
		error("`number` parameter is of invalid type.\n");
	}

	int want_low = asLogical(s_bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type.\n"); }

	if (len == 1) {
		return ScalarInteger(split_by_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id), num_to_select, want_low).num);
	} else {
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = split_by_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id), num_to_select, want_low).num;
		}
		UNPROTECT(1);
		return out;
	}

}

SEXP SXP_simple_selection_bypercent(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_percent, SEXP s_bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	double pc_to_select = asReal(s_percent);
	if (ISNA(pc_to_select) || pc_to_select < 0) {
		error("`percentage` parameter is of invalid type.\n");
	}

	int want_low = asLogical(s_bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type.\n"); }

	if (len == 1) {
		int group_size = get_group_size(d, GROUPNUM_IFY(groups[0]));
		int num_to_select = group_size * pc_to_select / 100; // integer division, so take the floor

		return ScalarInteger(split_by_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id), num_to_select, want_low).num);
	} else {
		int num_to_select;

		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (int i = 0; i < len; ++i) {
			num_to_select = get_group_size(d, GROUPNUM_IFY(groups[i])) * pc_to_select / 100;
			outc[i] = split_by_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id), num_to_select, want_low).num;
		}
		UNPROTECT(1);
		return out;
	}

}


/*----------------Data viewers---------------------*/

SEXP SXP_get_best_haplotype(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
		error("`effect.set` parameter is of invalid type.\n");
	}

	char* best_genotype = calculate_optimal_alleles(d, EFFECTID_IFY(eff_id));

	SEXP out = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(out, 0, mkChar(best_genotype));
	free(best_genotype);
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_best_available_haplotype(SEXP exd, SEXP s_groups, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	if (len == 1) {
		char* best_genotype = calculate_optimal_available_alleles(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id));

		SEXP out = PROTECT(allocVector(STRSXP, 1));
		SET_STRING_ELT(out, 0, mkChar(best_genotype));
		free(best_genotype);
		UNPROTECT(1);
		return out;

	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(STRSXP, len));
		for (int i = 0; i < len; ++i) {
			char* best_genotype = calculate_optimal_available_alleles(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id));
			SET_STRING_ELT(out, i, mkChar(best_genotype));
			free(best_genotype);
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_get_best_GEBV(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
		error("`effect.set` parameter is of invalid type.\n");
	}

	double best_GEBV = calculate_optimum_bv(d, EFFECTID_IFY(eff_id));

	SEXP out = PROTECT(allocVector(REALSXP, 1));
	REAL(out)[0] = best_GEBV;
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_best_available_GEBV(SEXP exd, SEXP s_groups, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	int len = length(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid.\n"); }
	}

	if (len == 1) {
		double availableBV = calculate_optimal_available_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id));

		SEXP out = PROTECT(allocVector(REALSXP, 1));
		REAL(out)[0] = availableBV;
		UNPROTECT(1);
		return out;

	} else {

		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(REALSXP, len));
		double* outc = REAL(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = calculate_optimal_available_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id));
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_get_worst_GEBV(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
		error("`effect.set` parameter is of invalid type.\n");
	}

	double worst_GEBV = calculate_minimum_bv(d, EFFECTID_IFY(eff_id));

	SEXP out = PROTECT(allocVector(REALSXP, 1));
	REAL(out)[0] = worst_GEBV;
	UNPROTECT(1);
	return out;
}

SEXP SXP_find_crossovers(SEXP exd, SEXP s_parentFile, SEXP s_outFile, SEXP s_windowSize, SEXP s_certainty) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char* load_fname = CHAR(asChar(s_parentFile));
	const char* save_fname = CHAR(asChar(s_outFile));

	int win = asInteger(s_windowSize);
	if (win < 0 || win == NA_INTEGER) { error("`window.size` parameter is invalid.\n"); }

	int cert = asLogical(s_certainty);
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

	GroupNum out_groups[d->n_groups];
	unsigned int out_counts[d->n_groups];
	int n_groups = get_existing_group_counts(d, out_groups, out_counts);

	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SEXP s_groups = PROTECT(allocVector(INTSXP, n_groups));
	int* cgroups = INTEGER(s_groups);
	SEXP counts = PROTECT(allocVector(INTSXP, n_groups));
	int* ccounts = INTEGER(counts);

	for (int i = 0; i < n_groups; ++i) {
		cgroups[i] = out_groups[i].num;
		ccounts[i] = out_counts[i];
	}
	SET_VECTOR_ELT(out, 0, s_groups);
	SET_VECTOR_ELT(out, 1, counts);

	UNPROTECT(3);
	return out;

}

SEXP SXP_get_group_data(SEXP exd, SEXP s_group, SEXP s_whatData, SEXP s_which) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	char c = CHAR(asChar(s_whatData))[0];
	char c2;
	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	if (group_size == 0) {
		error("The group is empty\n");
	}

	if (c == 'D') {
		PedigreeID data[group_size];
		memset(data, 0, sizeof(PedigreeID) * group_size);
		get_group_ids(d, GROUPNUM_IFY(group_id), group_size, data);

		// save to an R vector
		SEXP out = PROTECT(allocVector(INTSXP, group_size));
		int* outc = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i].id;
		}
		UNPROTECT(1);
		return out;

	} else if (c == 'X') {
		unsigned int data[group_size];
		get_group_indexes(d, GROUPNUM_IFY(group_id), group_size, data);

		// save to an R vector
		SEXP out = PROTECT(allocVector(INTSXP, group_size));
		int* outc = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		UNPROTECT(1);
		return out;

	} else if (c == 'B') {
		if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

		int eff_id = asInteger(s_which);
		if (eff_id == NA_INTEGER || eff_id < 1) {
		error("`which` parameter is of invalid type: needs to be effect set id.\n");
		}

		double data[group_size];
		memset(data, 0, sizeof(double) * group_size);
		get_group_bvs(d, GROUPNUM_IFY(group_id),  EFFECTID_IFY(eff_id), group_size, data);

		SEXP out = PROTECT(allocVector(REALSXP, group_size));
		double* outc = REAL(out);
		for (int i = 0; i < group_size; ++i) {
			outc[i] = data[i];
		}
		UNPROTECT(1);
		return out;

	} else if (c == 'N') {
		char* data[group_size];
		memset(data, 0, sizeof(char*) * group_size);
		get_group_names(d, GROUPNUM_IFY(group_id), group_size, data);

		PedigreeID backupdata[group_size];
		memset(backupdata, 0, sizeof(PedigreeID) * group_size);
		get_group_ids(d, GROUPNUM_IFY(group_id), group_size, backupdata);

		char buffer[get_integer_digits(backupdata[group_size - 1].id) + 1];

		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			if (data[i] != NULL) {
				SET_STRING_ELT(out, i, mkChar(data[i]));
			} else {
				sprintf(buffer, "%d", backupdata[i].id);
				SET_STRING_ELT(out, i, mkChar(buffer));
			}
		}
		UNPROTECT(1);
		return out;

	} else if (c == 'G') {
		char* rawdata[group_size];
		memset(rawdata, 0, sizeof(char*) * group_size);
		get_group_genes(d, GROUPNUM_IFY(group_id), group_size, rawdata);

		char* data[group_size];
		memset(data, 0, sizeof(char*) * group_size);
		// copy over to another array, adding string terminators
		int glen = d->n_markers * 2; // genotype length
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
		UNPROTECT(1);
		return out;

	} else if (c == 'P' && (c2 = CHAR(asChar(s_whatData))[1]) == 'E' && CHAR(asChar(s_whatData))[2] == 'D') {
		char* data[group_size];
		memset(data, 0, sizeof(char*) * group_size);
		get_group_pedigrees(d, GROUPNUM_IFY(group_id), group_size, data);

		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			SET_STRING_ELT(out, i, mkChar(data[i]));
			free(data[i]);
		}
		UNPROTECT(1);
		return out;

	} else if (c == 'P' && (c2 == '1' || c2 == '2')) {
		int parent = c2 == '1' ? 1 : 2;
		char* data[group_size];
		memset(data, 0, sizeof(char*) * group_size);
		get_group_parent_names(d, GROUPNUM_IFY(group_id), group_size, parent, data);

		PedigreeID backupdata[group_size];
		memset(backupdata, 0, sizeof(PedigreeID) * group_size);
		get_group_parent_ids(d, GROUPNUM_IFY(group_id), group_size, parent, backupdata);

		char buffer[get_integer_digits(backupdata[group_size - 1].id) + 1];

		SEXP out = PROTECT(allocVector(STRSXP, group_size));
		for (int i = 0; i < group_size; ++i) {
			if (data[i] != NULL) {
				SET_STRING_ELT(out, i, mkChar(data[i]));
			} else {
				sprintf(buffer, "%d", backupdata[i].id);
				SET_STRING_ELT(out, i, mkChar(buffer));
			}
		}
		UNPROTECT(1);
		return out;

	} else {
		error("`data.type` parameter is not a valid option.");
	}
}


SEXP SXP_get_group_gene_data(SEXP exd, SEXP s_group, SEXP s_countAllele) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type.\n");
	}

	int group_size = get_group_size(d, GROUPNUM_IFY(group_id));
	if (group_size == 0) {
		error("The group is empty\n");
	}

	if (asChar(s_countAllele) == NA_STRING) {
		// Give the pairs of alleles

		char* rawdata[group_size];
		memset(rawdata, 0, sizeof(char*) * group_size);
		get_group_genes(d, GROUPNUM_IFY(group_id), group_size, rawdata);

		SEXP out = PROTECT(allocMatrix(STRSXP, d->n_markers, group_size));
		for (int i = 0; i < group_size; ++i) {
			for (int j = 0; j < d->n_markers; ++j) {
				char alleles[] = {rawdata[i][2*j], rawdata[i][2*j + 1], '\0'};
				SET_STRING_ELT(out, j + d->n_markers*i, mkChar(alleles));
			}
		}
		UNPROTECT(1);
		return out;

	} else {
		// Give the counts of allele c
		char c = CHAR(asChar(s_countAllele))[0];

		char* rawdata[group_size];
		memset(rawdata, 0, sizeof(char*) * group_size);
		get_group_genes(d, GROUPNUM_IFY(group_id), group_size, rawdata);

		SEXP out = PROTECT(allocMatrix(INTSXP, d->n_markers, group_size));
		int* cout = INTEGER(out);
		for (int i = 0; i < group_size; ++i) {
			for (int j = 0; j < d->n_markers; ++j) {
				int count = 0;
				if (rawdata[i][2*j] == c)     { ++count; }
				if (rawdata[i][2*j + 1] == c) { ++count; }
				cout[j + d->n_markers*i] = count;
			}
		}
		UNPROTECT(1);
		return out;

	}
}


SEXP SXP_change_name_values(SEXP exd, SEXP s_values, SEXP s_group, SEXP s_start) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (TYPEOF(s_values) != STRSXP) {
		error("`values` are invalid: must be strings.");
	}
	int len = length(s_values);
	char* names[len];
	for (int i = 0; i < len; ++i) {
		names[i] = R_Calloc(sizeof(char)*(NAME_LENGTH + 1), char);
		strncpy(names[i], CHAR(STRING_ELT(s_values, i)), sizeof(char)*NAME_LENGTH);
		//names[i][NAME_LENGTH] = '\0'; // terminate 'em. Just in case they're trying
		// to spill over NAME_LENGTH
	}

	int group = asInteger(s_group);
	if (group == NA_INTEGER || group < 0) {
		error("`group` is invalid (wrong type or negative).\n");
	}

	int startIndex = asInteger(s_start);
	if (startIndex == NA_INTEGER || startIndex < 0) {
		error("`startIndex` is invalid (wrong type or negative).\n");
	}

	set_names_to_values(d, GROUPNUM_IFY(group), startIndex, len, (const char**) names);

	for (int i = 0; i < len; ++i) {
		R_Free(names[i]);
	}

	return ScalarInteger(0);
}

/*--------------------------------Printing-----------------------------------*/

SEXP SXP_save_simdata(SEXP exd, SEXP s_filename) {
	FILE* f;
	const char* filename = CHAR(asChar(s_filename));
	if ((f = fopen(filename, "w")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	save_simdata(f, d);

	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_genotypes(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_type) {
	FILE* f;
	const char* filename = CHAR(asChar(s_filename));
	if ((f = fopen(filename, "w")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char t = CHAR(asChar(s_type))[0];
	if (t == 'R' || t == 'r') {
		if (isNull(s_group)) {
			save_names_header(f, d->n_markers, (const char**) d->markers);
			save_allele_matrix(f, d->m);
		} else if (asInteger(s_group) >= 0) {
			save_group_alleles(f, d, GROUPNUM_IFY(asInteger(s_group)));
		} else {
			fclose(f);
			error("Supplied group number is invalid.");
		}
	} else if (t == 'T' || t == 't') {
		if (isNull(s_group)) {
			save_transposed_allele_matrix(f, d->m, (const char**) d->markers);
		} else if (asInteger(s_group) >= 0) {
			save_transposed_group_alleles(f, d, GROUPNUM_IFY(asInteger(s_group)));
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

SEXP SXP_save_counts(SEXP exd, SEXP s_filename, SEXP s_group, SEXP allele) {
	FILE* f;
	const char* filename = CHAR(asChar(s_filename));
	if ((f = fopen(filename, "w")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char t = CHAR(asChar(allele))[0];
	if (isNull(s_group)) {
		save_count_matrix(f, d, t);
	} else if (asInteger(s_group) >= 0) {
		save_count_matrix_of_group(f, d, t, GROUPNUM_IFY(asInteger(s_group)));
	} else {
		fclose(f);
		error("Supplied group number is invalid.");
	}

	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_pedigrees(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_type) {
	FILE* f;
	const char* filename = CHAR(asChar(s_filename));
	if ((f = fopen(filename, "w")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char t = CHAR(asChar(s_type))[0];
	if (t == 'R' || t == 'r') { // full/recursive
		if (isNull(s_group)) {
			save_full_pedigree(f, d);
		} else if (asInteger(s_group) >= 0) {
			save_group_full_pedigree(f, d, GROUPNUM_IFY(asInteger(s_group)));
		} else {
			error("Supplied group number is invalid.");
		}
	} else if (t == 'P' || t == 'p') { // one-step/parents
		if (isNull(s_group)) {
			save_one_step_pedigree(f, d);
		} else if (asInteger(s_group) >= 0) {
			save_group_one_step_pedigree(f, d, GROUPNUM_IFY(asInteger(s_group)));
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

SEXP SXP_save_GEBVs(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_eff_set) {
	FILE* f;
	const char* filename = CHAR(asChar(s_filename));
	if ((f = fopen(filename, "w")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	if (isNull(s_group)) {
		save_bvs(f, d, EFFECTID_IFY(eff_id));
	} else if (asInteger(s_group) >= 0) {
		save_group_bvs(f, d, GROUPNUM_IFY(asInteger(s_group)), EFFECTID_IFY(eff_id));
	} else {
		fclose(f);
		error("Supplied group number is invalid.");
	}

	fclose(f);
	return ScalarInteger(0);
}

SEXP SXP_save_file_block_effects(SEXP exd, SEXP s_filename, SEXP block_file, SEXP s_group, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	if (isNull(s_group)) {
		MarkerBlocks b = read_block_file(d, CHAR(asChar(block_file)));
		calculate_local_bvs(d, b, EFFECTID_IFY(eff_id), CHAR(asChar(s_filename)));
		delete_markerblocks(&b);
	} else if (asInteger(s_group) > 0) {
		MarkerBlocks b = read_block_file(d, CHAR(asChar(block_file)));
		calculate_group_local_bvs(d, b, EFFECTID_IFY(eff_id), CHAR(asChar(s_filename)), GROUPNUM_IFY(asInteger(s_group)));
		delete_markerblocks(&b);
	} else {
		error("Supplied group number is invalid.\n");
	}

	return ScalarInteger(0);
}

SEXP SXP_save_chrsplit_block_effects(SEXP exd, SEXP s_filename, SEXP s_nslices, SEXP s_group, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function.\n"); }

	int eff_id = asInteger(s_eff_set);
	if (eff_id == NA_INTEGER || eff_id < 0) {
	  error("`effect.set` parameter is of invalid type.\n");
	}

	if (isNull(s_group)) {
		MarkerBlocks b = create_n_blocks_by_chr(d, asInteger(s_nslices));
		calculate_local_bvs(d, b, EFFECTID_IFY(eff_id), CHAR(asChar(s_filename)));
		delete_markerblocks(&b);
	} else if (asInteger(s_group) > 0) {
		MarkerBlocks b = create_n_blocks_by_chr(d, asInteger(s_nslices));
		calculate_group_local_bvs(d, b, EFFECTID_IFY(eff_id), CHAR(asChar(s_filename)), GROUPNUM_IFY(asInteger(s_group)));
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

SEXP SXP_clear_simdata(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	delete_simdata(d);
	return ScalarInteger(0);
}

SEXP SXP_delete_group(SEXP exd, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int n = length(s_groups);
	int *groups = INTEGER(s_groups);

	for (int i = 0; i < n; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 1) {
			error("Entry %d in `group` parameter is of invalid type.\n", i + 1);
		}
		delete_group(d, GROUPNUM_IFY(groups[i]));
	}

	return ScalarInteger(0);
}

SEXP SXP_delete_label(SEXP exd, SEXP s_labels) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int n = length(s_labels);
	int *labels = INTEGER(s_labels);

	for (int i = 0; i < n; ++i) {
		if (labels[i] == NA_INTEGER || labels[i] < 1) {
			error("Entry %d in `labels` parameter is of invalid type.\n", i + 1);
		}
		delete_label(d, LABELID_IFY(labels[i]));
	}

	return ScalarInteger(0);
}

SEXP SXP_delete_eff_set(SEXP exd, SEXP s_eff_sets) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int n = length(s_eff_sets);
	int *eff_sets = INTEGER(s_eff_sets);

	for (int i = 0; i < n; ++i) {
		if (eff_sets[i] == NA_INTEGER || eff_sets[i] < 1) {
			error("Entry %d in `effect_sets` parameter is of invalid type.\n", i + 1);
		}
		delete_eff_set(d, EFFECTID_IFY(eff_sets[i]));
	}

	return ScalarInteger(0);
}
