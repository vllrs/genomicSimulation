#include "r-wrappers.h"

/*--------------------------SEXP HELPERS----------------------*/
#define GROUPNUM_IFY(n) (GroupNum){.num=((n == NA_INTEGER) ? 0 : n)}
#define EFFECTID_IFY(n) (EffectID){.id=n}
#define MAPID_IFY(n)  (MapID){.id=n}
#define LABELID_IFY(n) (LabelID){.id=n}
#define PEDIGREEID_IFY(n) (PedigreeID){.id=n}

// The maximum length of a buffer needed to store a base36 representation of a number of a given integer type
// Explanation: sizeof(type)*CHAR_BIT is the number of binary digits.
//     The number of binary digits is more than 5 times the number of base36 digits, but less than 6 times. 
//     So b/5 is greater than the number of base36 digits. Integer division makes this "/ 5 + 1"
//     And then add one more, to store the terminating null character of the base36 string.
#define BASE36_BUFFER_SIZE(type) (sizeof(type)*CHAR_BIT / 5 + 2)

void convertVECSXP_to_GroupNum(SEXP container, GroupNum* output) {
  R_xlen_t len = xlength(container);
	int *intvec = INTEGER(container);

	for (R_xlen_t i = 0; i < len; ++i) {
		if (intvec[i] <= 0) {
			Rprintf("NOTE! %lu is not a possible group number", intvec[i]);
			output[i] = NO_GROUP;
		} else {
			output[i] = GROUPNUM_IFY(intvec[i]);
		}
	}
}

void base36tostr(unsigned long n, char* buf) {
	// Construct the string backwards
	size_t len = 1;
	char* s = buf; 
	*s = '\0';
	do {
		++s; ++len;
		*s = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n % 36];
		n /= 36;
	} while (n > 0 && len < 15); // the 'len' condition is just for safety
	
	// turn it all around
	for (size_t i = 0; i < len/2; ++i) {
		char tmp = buf[i]; buf[i] = buf[len-1-i]; buf[len-1-i] = tmp;
	}
}

SEXP register_markerblocksEXPTR(void* b) {
	SEXP bptr = PROTECT(R_MakeExternalPtr(b, Rf_install("MarkerBlocks"), R_NilValue));
	R_RegisterCFinalizerEx(bptr, SXP_delete_markerblocks, 1);
	
	/*MarkerBlocks* b2 = (MarkerBlocks*) R_ExternalPtrAddr(bptr);
	Rprintf("Post-registration: n_blocks = %d\n", b2->num_blocks);
	Rprintf("p: %p\n", R_ExternalPtrAddr(bptr));*/
	
	UNPROTECT(1);
	return bptr;
}

void check_if_is_markerblocks(SEXP exptr) {
	SEXP btag = R_ExternalPtrTag(exptr);
	if (strcmp(CHAR(asChar(btag)), "MarkerBlocks") != 0) {
		error("The pointer is not a MarkerBlocks object");
	}
}

GSC_ID_T extract_effid(SEXP s_eff_set, SimData* d) {
	if (!d) { error("Simulation object is invalid."); }
	
	int eff_id = asInteger(s_eff_set);
	
	if (eff_id == NA_INTEGER || eff_id < 0) {
		error("`effect.set` parameter is of invalid type\n");
	} else if (eff_id == 0) {
		if (d->n_eff_sets > 0) {
			eff_id = d->eff_set_ids[0].id;
		} else {
			error("Simulation has no effect sets loaded.\n");
		}
	}
	return eff_id;
}

GSC_ID_T extract_mapid(SEXP s_mapid, SimData* d, const char* paramname) {
	if (!d) { error("Simulation object is invalid."); }
	
	int map_id = asInteger(s_mapid);
	
	if (map_id == NA_INTEGER || map_id < 0) {
		error("`%s` parameter is of invalid type\n", paramname);
	} else if (map_id == 0) {
		if (d->genome.n_maps > 0) {
			map_id = d->genome.map_ids[0].id;
		} else {
			error("Simulation has no genetic maps loaded.\n");
		}
	}
	return map_id;
}
GSC_ID_T parse_mapid(GSC_ID_T map_id, SimData* d) {
	if (!d) { error("Simulation object is invalid."); }
	
	if (map_id == 0) {
		if (d->genome.n_maps > 0) {
			map_id = d->genome.map_ids[0].id;
		} else {
			error("Simulation has no genetic maps loaded.\n");
		}
	}
	return map_id;
}

GSC_ID_T extract_labelid(SEXP s_label_id) {
	int label_id = asInteger(s_label_id);
	if (label_id == NA_INTEGER || label_id < 1) {
		error("`label` parameter is of invalid type\n");
	}
	return label_id;
}

/*-------------------------- Setup -------------------------*/
void init_format_as_matrix(FileFormatSpec* mformat) {
  if (mformat->filetype == GSC_GENOTYPEFILE_UNKNOWN) {
    mformat->filetype = GSC_GENOTYPEFILE_MATRIX;
    mformat->spec.matrix = (struct gsc_GenotypeFile_MatrixFormat) {
           .has_header = GSC_NA,
           .markers_as_rows = GSC_NA,
           .cell_style = GSC_GENOTYPECELLSTYLE_UNKNOWN
    };
    
  } else if (mformat->filetype != GSC_GENOTYPEFILE_MATRIX) {
    //warning("why are you union hacking in R.")
    error("Malformed `format` parameter");
  }
}

FileFormatSpec SXP_parse_filespec_list(SEXP s_fileSpec) {
  FileFormatSpec mformat = DETECT_FILE_FORMAT;
  
  SEXP s_listnames = getAttrib(s_fileSpec, R_NamesSymbol);
  
  for (int i = 0; i < length(s_fileSpec); i++) {
    if(strcmp(CHAR(STRING_ELT(s_listnames, i)), "has.header") == 0) {
      init_format_as_matrix(&mformat);
      
      int tmp = asLogical(VECTOR_ELT(s_fileSpec, i));
      if (tmp == NA_LOGICAL) {
        mformat.spec.matrix.has_header = GSC_NA;
      } else if (tmp) {
        mformat.spec.matrix.has_header = GSC_TRUE;
      } else {
        mformat.spec.matrix.has_header = GSC_FALSE;
      }
    }
    
    if(strcmp(CHAR(STRING_ELT(s_listnames, i)), "markers.as.rows") == 0) {
      init_format_as_matrix(&mformat);
      
      int tmp = asLogical(VECTOR_ELT(s_fileSpec, i));
      if (tmp == NA_LOGICAL) {
        mformat.spec.matrix.markers_as_rows = GSC_NA;
      } else if (tmp) {
        mformat.spec.matrix.markers_as_rows = GSC_TRUE;
      } else {
        mformat.spec.matrix.markers_as_rows = GSC_FALSE;
      }
    }
    
    if(strcmp(CHAR(STRING_ELT(s_listnames, i)), "cell.style") == 0) {
      init_format_as_matrix(&mformat);
      
      char c = CHAR(asChar(VECTOR_ELT(s_fileSpec, i)))[0];
      switch (c) {
        case 'P': mformat.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_PAIR; break;
        case '/': mformat.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_SLASHPAIR; break;
        case 'C': mformat.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_COUNT; break;
        case 'I': mformat.spec.matrix.cell_style = GSC_GENOTYPECELLSTYLE_ENCODED; break;
      }
    }
  }
  return mformat;
}

SEXP SXP_load_data(SEXP s_alleleFile, SEXP s_mapFile, SEXP s_effectFile, SEXP s_fileSpec) {
	SimData* d = create_empty_simdata();
	const char* c_alleleFile;
	if (length(s_alleleFile) == 0 || asChar(s_alleleFile) == NA_STRING) {
		c_alleleFile = NULL;
	} else {
		c_alleleFile = CHAR(asChar(s_alleleFile));
	}
	
	const char* c_mapFile;
	if (length(s_mapFile) == 0 || asChar(s_mapFile) == NA_STRING) {
		c_mapFile = NULL;
	} else {
		c_mapFile = CHAR(asChar(s_mapFile));
	}
	
	const char* c_effectFile;
	if (length(s_effectFile) == 0 || asChar(s_effectFile) == NA_STRING) {
		c_effectFile = NULL;
	} else {
		c_effectFile = CHAR(asChar(s_effectFile));
	}
	
	if (!isNewList(s_fileSpec)) {
	  error("`format` parameter is of invalid type: must be a list\n");
	}
	load_data_files(d, c_alleleFile, c_mapFile, c_effectFile, SXP_parse_filespec_list(s_fileSpec));

	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP SXP_load_genotypes(SEXP exd, SEXP s_alleleFile, SEXP s_fileSpec) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
  if (!isNewList(s_fileSpec)) {
    error("`format` parameter is of invalid type: must be a list\n");
  }
	return ScalarInteger(load_genotypefile(d, CHAR(asChar(s_alleleFile)), SXP_parse_filespec_list(s_fileSpec)).num);
}

SEXP SXP_load_map(SEXP exd, SEXP s_mapFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_mapfile(d, CHAR(asChar(s_mapFile))).id);
}

SEXP SXP_load_effects(SEXP exd, SEXP s_effectFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_effectfile(d, CHAR(asChar(s_effectFile))).id);
}

SEXP SXP_create_new_label(SEXP exd, SEXP s_default) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int lblDefault = asInteger(s_default);
	if (lblDefault == NA_INTEGER) {
		error("`default` parameter is of invalid type: must be an integer\n");
	}
	return ScalarInteger(create_new_label(d, lblDefault).id);
}

SEXP SXP_change_eff_set_centres(SEXP exd, SEXP s_values, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	GSC_ID_T eff_id = extract_effid(s_eff_set, d);
	
	// s_values: if is vector, extract as vector and call change_eff_set_centres_to_values
	//           if is list, check it has two entries of same length and right type, and call change_eff_set_centres_of_markers
	if (TYPEOF(s_values) == REALSXP) {
		double* dvalues = REAL(s_values);
		R_xlen_t n_values = xlength(s_values);
		
		gsc_change_eff_set_centres_to_values(d, EFFECTID_IFY(eff_id), n_values, dvalues);
		
	} else if (isNewList(s_values)) {
		if (xlength(s_values) != 2) {
			error("A list input for parameter `to` must have two entries, one a list of marker names, and the second a list of corresponding centring values.");
		}
		SEXP s_marker_names = VECTOR_ELT(s_values, 0);
		if (xlength(s_marker_names) != xlength(VECTOR_ELT(s_values, 1))) {
			error("A list input for parameter `to` must have two entries/two columns of equal length.");
		}
		R_xlen_t n_values = xlength(s_marker_names);
		
		char** marker_names = R_Calloc(n_values, char*);
		for (R_xlen_t i = 0; i < n_values; ++i) {
			// casting away const but promise it is still being used as readonly.
			marker_names[i] = (char*) CHAR(STRING_ELT(s_marker_names, i));
		}
		
		double* centre_values = REAL(VECTOR_ELT(s_values, 1));
		
		GSC_GENOLEN_T n_successes = gsc_change_eff_set_centre_of_markers(d, EFFECTID_IFY(eff_id), n_values, 
                                                                   (const char**) marker_names, 
                                                                   centre_values);
		
		R_Free(marker_names);
		
		if (n_successes != n_values) {
			Rprintf("%d of the %d markers in the input could not be located in the genetic map\n", n_values - n_successes, n_values);
		}
		
	} else {
		error("`to` parameter is of invalid type");
	}
	
	return ScalarInteger(0);
}

SEXP SXP_change_eff_set_centres_of_allele(SEXP exd, SEXP s_allele, SEXP s_values, SEXP s_eff_set, SEXP s_reset) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	GSC_ID_T eff_id = extract_effid(s_eff_set, d);
	
	if (asChar(s_allele) == NA_STRING) { error("`allele` parameter is of invalid type"); }
	char allele = CHAR(asChar(s_allele))[0];
	
	int reset = asLogical(s_reset);
	if (reset == NA_LOGICAL) { error("`reset.centres` parameter is of invalid type\n"); }

	// s_values: check it has two entries of same length and right type
	if (isNewList(s_values)) {
		if (xlength(s_values) != 2) {
			error("A list input for parameter `to` must have two entries, one a list of marker names, and the second a list of corresponding centring values.");
		}
		SEXP s_marker_names = VECTOR_ELT(s_values, 0);
		if (xlength(s_marker_names) != xlength(VECTOR_ELT(s_values, 1))) {
		  error("A list input for parameter `to` must have two entries/two columns of equal length.");
		}
		R_xlen_t n_values = xlength(s_marker_names);
		
		char** marker_names = R_Calloc(n_values, char*);
		for (R_xlen_t i = 0; i < n_values; ++i) {
		  // casting away const but promise it is still being used as readonly.
		  marker_names[i] = (char*) CHAR(STRING_ELT(s_marker_names, i));
		}
		
		double* centre_values = REAL(VECTOR_ELT(s_values, 1));
		
		GSC_GENOLEN_T n_successes = gsc_change_eff_set_centre_of_allele_count(d, EFFECTID_IFY(eff_id), n_values, 
                                                                        (const char**) marker_names, centre_values, 
                                                                        allele, reset);
		
		R_Free(marker_names);
		
		if (n_successes != n_values) {
			Rprintf("%d of the %d markers in the input could not be located in the genetic map\n", n_values - n_successes, n_values);
		}
		
	} else {
		error("`to` parameter is of invalid type");
	}
	
	return ScalarInteger(0);
}


SEXP SXP_create_markerblocks_df(SEXP exd, SEXP s_markernames, SEXP s_markerallocs) {
	// Assumes markernames/markerallocs are sorted according to allocation number. So blocks are always sequential 
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	// Check lengths of the two vcectors are equal
	if (xlength(s_markernames) != xlength(s_markerallocs)) {
		error("Columns passed to create.markerblocks are not the same length\n");
	} else if (xlength(s_markernames) == 0) {
		error("Empty inputs provided\n");
	}
	
	// First, convert marker names to indexes and remove missing markers
	int* c_markerallocs = INTEGER(s_markerallocs);
	
	GSC_GENOLEN_T* filt_markerix = R_Calloc(xlength(s_markernames), GSC_GENOLEN_T);
	int* filt_markerallocs = R_Calloc(xlength(s_markernames), int);
	R_xlen_t filt_ix = 0;
	
	GSC_GENOLEN_T nblocks = 1;
	int currentblock = c_markerallocs[0];
	
	for (R_xlen_t i = 0; i < xlength(s_markernames); ++i) {
		const char* markername = CHAR(STRING_ELT(s_markernames,i));
		GSC_GENOLEN_T markerix = NA_GENOLEN;
		if (get_index_of_genetic_marker(markername, d->genome, &markerix)) {
			filt_markerix[filt_ix] = markerix;
			filt_markerallocs[filt_ix] = c_markerallocs[i];
			++filt_ix;
			
			if (c_markerallocs[i] != currentblock) {
				++nblocks;
				currentblock = c_markerallocs[i];
			}
		}
	}
	
	
	// Then, save in correct MarkerBlocks struct format.
	MarkerBlocks* b = R_Calloc(1, MarkerBlocks);
	b->num_blocks=nblocks;
	b->num_markers_in_block=R_Calloc(nblocks, GSC_GENOLEN_T);
	b->markers_in_block=R_Calloc(nblocks, GSC_GENOLEN_T*);

	GSC_GENOLEN_T blockix = 0;
	GSC_GENOLEN_T startix = 0;
	while (startix < filt_ix) {
		int blockident = filt_markerallocs[startix];
		GSC_GENOLEN_T endix = startix + 1;
		while (filt_markerallocs[endix] == blockident) { ++endix; }
		GSC_GENOLEN_T nmembers = endix - startix;
		
		b->num_markers_in_block[blockix] = nmembers;
		b->markers_in_block[blockix] = R_Calloc(nmembers, GSC_GENOLEN_T);
		for (GSC_GENOLEN_T j = 0; j < nmembers; ++j) {
			b->markers_in_block[blockix][j] = filt_markerix[startix+j];
		}
		
		startix = endix;
		++blockix;
	}
	
	R_Free(filt_markerix);
	R_Free(filt_markerallocs);
	
	return register_markerblocksEXPTR(b);
}

SEXP SXP_create_markerblocks_nperchr(SEXP exd, SEXP s_nperchr, SEXP s_mapid) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	int nperchr = asInteger(s_nperchr);
	if (nperchr == NA_INTEGER || nperchr <= 0) {
		error("`n.blocks.per.chr` parameter must be a positive number\n");
	}
	
	GSC_ID_T mapid = extract_mapid(s_mapid, d, "map");
	
	MarkerBlocks b = create_evenlength_blocks_each_chr(d, MAPID_IFY(mapid), nperchr);
	MarkerBlocks* bcopy = R_Calloc(1, MarkerBlocks);
	memcpy(bcopy, &b, sizeof(b));
	
	return register_markerblocksEXPTR(bcopy);
}

/*----------------------- Calculators ---------------------*/


SEXP SXP_see_optimal_haplotype(SEXP exd, SEXP s_eff_set, SEXP s_unknown_char) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);
	
	if (asChar(s_unknown_char) == NA_STRING) { error("`unknown.allele` parameter is of invalid type"); }
	char na_char = CHAR(asChar(s_unknown_char))[0];

	char* best_genotype = (char*)R_alloc(d->genome.n_markers+1, sizeof(char));
	calculate_optimal_haplotype(d, EFFECTID_IFY(eff_id), na_char, best_genotype);

	SEXP out = PROTECT(allocVector(STRSXP, 1));
	SET_STRING_ELT(out, 0, mkChar(best_genotype));
	//free(best_genotype); // R_alloc does not need free
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_optimal_possible_haplotype(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_unknown_char) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);
	
	if (asChar(s_unknown_char) == NA_STRING) { error("`unknown.allele` parameter is of invalid type"); }
	char na_char = CHAR(asChar(s_unknown_char))[0];

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}

	if (len == 1) {
		char* best_genotype = (char*)R_alloc(d->genome.n_markers+1, sizeof(char));
		calculate_optimal_possible_haplotype(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id), na_char, best_genotype);

		SEXP out = PROTECT(allocVector(STRSXP, 1));
		SET_STRING_ELT(out, 0, mkChar(best_genotype));
		//free(best_genotype); // R_alloc does not need free
		UNPROTECT(1);
		return out;

	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(STRSXP, len));
		for (int i = 0; i < len; ++i) {
			char* best_genotype = (char*)R_alloc(d->genome.n_markers+1, sizeof(char));
			calculate_optimal_possible_haplotype(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id), na_char, best_genotype);
			SET_STRING_ELT(out, i, mkChar(best_genotype));
			//free(best_genotype); // R_alloc does not need free
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_get_optimal_GEBV(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

	double best_GEBV = calculate_optimal_bv(d, EFFECTID_IFY(eff_id));

	SEXP out = PROTECT(allocVector(REALSXP, 1));
	REAL(out)[0] = best_GEBV;
	UNPROTECT(1);
	return out;
}

SEXP SXP_get_optimal_possible_GEBV(SEXP exd, SEXP s_groups, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n",i+1); }
	}

	if (len == 1) {
		double availableBV = calculate_optimal_possible_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id));

		SEXP out = PROTECT(allocVector(REALSXP, 1));
		REAL(out)[0] = availableBV;
		UNPROTECT(1);
		return out;

	} else {

		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(REALSXP, len));
		double* outc = REAL(out);
		for (int i = 0; i < len; ++i) {
			outc[i] = calculate_optimal_possible_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id));
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_get_minimal_GEBV(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

	double worst_GEBV = calculate_minimal_bv(d, EFFECTID_IFY(eff_id));

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
	if (win < 0 || win == NA_INTEGER) { error("`window.size` parameter is invalid\n"); }

	int cert = asLogical(s_certainty);
	if (cert == NA_LOGICAL) { error("`certainty` parameter is invalid\n"); }

	return ScalarInteger(gsc_calculate_recombinations_from_file(d, load_fname, save_fname, win, cert));
}

SEXP SXP_see_map(SEXP exd, SEXP s_map) {
  SimData* d = (SimData*) R_ExternalPtrAddr(exd);
  
  GSC_ID_T mapid = extract_mapid(s_map, d, "map");
  int mapix = 0;
  if (mapid > 0) { mapix = get_index_of_map(d, MAPID_IFY(mapid)); }
  if (mapix < 0 || mapix >= d->genome.n_maps) {
	error("Invalid map identifier");
  }

	SEXP map = PROTECT(allocVector(VECSXP, 3));
	SEXP ssnp = PROTECT(allocVector(STRSXP, d->genome.n_markers));
	SEXP schr = PROTECT(allocVector(STRSXP, d->genome.n_markers));
	SEXP spos = PROTECT(allocVector(REALSXP, d->genome.n_markers));
	double* cpos = REAL(spos);
	
	char chr_name_buffer[BASE36_BUFFER_SIZE(unsigned long)];

	GSC_GENOLEN_T m_ix = 0;
	for (int chr = 0; chr < d->genome.maps[mapix].n_chr; ++chr) {
	  base36tostr(d->genome.maps[mapix].chr_names[chr], chr_name_buffer);
	  
	  if (d->genome.maps[mapix].chrs[chr].type == GSC_LINKAGEGROUP_SIMPLE) {
	    for (int i = 0; i < d->genome.maps[mapix].chrs[chr].map.simple.n_markers; ++i) {
	      SET_STRING_ELT(ssnp, m_ix, mkChar(d->genome.marker_names[i + d->genome.maps[mapix].chrs[chr].map.simple.first_marker_index]));
		  SET_STRING_ELT(schr, m_ix, mkChar(chr_name_buffer));
	      cpos[m_ix] = d->genome.maps[mapix].chrs[chr].map.simple.dists[i] * 
	        d->genome.maps[mapix].chrs[chr].map.simple.expected_n_crossovers * 100;
		  ++m_ix;
	    }
	  } else if (d->genome.maps[mapix].chrs[chr].type == GSC_LINKAGEGROUP_REORDER) {
	    for (int i = 0; i < d->genome.maps[mapix].chrs[chr].map.reorder.n_markers; ++i) {
	      SET_STRING_ELT(ssnp, m_ix, mkChar(d->genome.marker_names[d->genome.maps[mapix].chrs[chr].map.reorder.marker_indexes[i]]));
	      SET_STRING_ELT(schr, m_ix, mkChar(chr_name_buffer));
	      cpos[m_ix] = d->genome.maps[mapix].chrs[chr].map.reorder.dists[i] * 
	        d->genome.maps[mapix].chrs[chr].map.reorder.expected_n_crossovers * 100;
		  ++m_ix;
	    }
	  }
	}

	SET_VECTOR_ELT(map, 0, ssnp);
	SET_VECTOR_ELT(map, 1, schr);
	SET_VECTOR_ELT(map, 2, spos);
	
	UNPROTECT(4);
	return map;
}


SEXP SXP_see_effects(SEXP exd, SEXP s_eff_set) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T effID = extract_effid(s_eff_set, d);
	int effix = get_index_of_eff_set(d, EFFECTID_IFY(effID));
	if (effix < 0 || effix >= d->n_eff_sets) {
		error("Invalid marker effect set identifier");
	}
	MarkerEffects e = d->e[effix];
	
	// 1. create the main dataframe (marker/allele/eff)
	SEXP effset = PROTECT(allocVector(VECSXP, 3));
	
	GSC_GENOLEN_T len = e.cumn_alleles[e.n_markers - 1];
	if (e.centre != NULL) { len += e.n_markers; }
	SEXP marker = PROTECT(allocVector(STRSXP, len));
	SEXP allele = PROTECT(allocVector(STRSXP, len));
	SEXP eff = PROTECT(allocVector(REALSXP, len));
	double* ceff = REAL(eff);
	
	char an_allele[2] = {'?', '\0'};
	GSC_GENOLEN_T e_ix = 0;
	for (GSC_GENOLEN_T m_ix = 0; m_ix < e.n_markers; ++m_ix) {
		while (e_ix < e.cumn_alleles[m_ix]) {
			SET_STRING_ELT(marker, e_ix, mkChar(d->genome.marker_names[m_ix]));
			
			an_allele[0] = e.allele[e_ix];
			SET_STRING_ELT(allele, e_ix, mkChar(an_allele));
			
			ceff[e_ix] = e.eff[e_ix];
			
			++e_ix;
		}
	}
	
	if (e.centre != NULL) {
		GSC_GENOLEN_T last = e.cumn_alleles[e.n_markers - 1];
		for (GSC_GENOLEN_T m_ix = 0; m_ix < e.n_markers; ++m_ix) {
			SET_STRING_ELT(marker, last + m_ix, mkChar(d->genome.marker_names[m_ix]));
			SET_STRING_ELT(allele, last + m_ix, mkChar("(centre)"));
			ceff[last + m_ix] = e.centre[m_ix];
		}
	}
	
	SET_VECTOR_ELT(effset, 0, marker);
	SET_VECTOR_ELT(effset, 1, allele);
	SET_VECTOR_ELT(effset, 2, eff);
	
	UNPROTECT(4);
	return effset;	
}

/*----------------------- Data Access ---------------------*/

SEXP SXP_see_existing_groups(SEXP exd) {
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

enum GroupDataAccessOptions {
  ACCESS_PEDIGREEID,
  ACCESS_INDEX,
  ACCESS_NAME,
  ACCESS_BV,
  ACCESS_LABEL,
  ACCESS_GENOTYPE,
  ACCESS_PEDIGREESTRING,
  ACCESS_P1NAME,
  ACCESS_P2NAME,
  ACCESS_P1PEDIGREEID,
  ACCESS_P2PEDIGREEID,
  ACCESS_UNKNOWNTYPE
};

size_t check_group_sizes(SimData* d, R_xlen_t glen, int* groups, size_t* gsizes) {
  size_t cumulativesize = 0;
  for (R_xlen_t i = 0; i < glen; ++i) {
    gsizes[i] = 0;
    if (groups[i] == NA_INTEGER || groups[i] < 1) {
      if (glen == 1) {
        Rprintf("NOTE! `group` parameter is invalid (0 or negative)");
      } else {
        Rprintf("NOTE! %i (entry %lu in the `group` vector) is an invalid group number (0 or negative)", groups[i], i + 1);
      }
      
    } else {
      gsizes[i] = get_group_size(d, GROUPNUM_IFY(groups[i]));
      if (gsizes[i] == 0) {
        if (glen == 1) {
          Rprintf("NOTE! `group` is an empty group\n");
        } else {
          Rprintf("NOTE! Group %i (entry %lu in the `group` vector) is an empty group\n", groups[i], i + 1);
        }
      } else {
        cumulativesize += gsizes[i];
      }
    }
  }
  return cumulativesize;
}


SEXP SXP_see_group_data(SEXP exd, SEXP s_groups, SEXP s_whatData, SEXP s_eff_set_id, SEXP s_label_id) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
  
  enum GroupDataAccessOptions datatype = ACCESS_UNKNOWNTYPE;
  R_xlen_t typeselectorlen = xlength(asChar(s_whatData));
  char c = CHAR(asChar(s_whatData))[0];
  char c2, c3;
  switch (c) {
  case 'D':
    datatype = ACCESS_PEDIGREEID; break;
  case 'X':
    datatype = ACCESS_INDEX; break;
  case 'N':
    datatype = ACCESS_NAME; break;
  case 'B':
    datatype = ACCESS_BV; break;
  case 'L':
    datatype = ACCESS_LABEL; break;
  case 'G':
    datatype = ACCESS_GENOTYPE; break;
  case 'P':
    if (typeselectorlen == 1) {
      break;
    }
    c2 = CHAR(asChar(s_whatData))[1];
    c3 = typeselectorlen > 2 ? CHAR(asChar(s_whatData))[2] : '\0';
    switch (c3) {
    case 'D':
      switch (c2) {
      case 'E':
        datatype = ACCESS_PEDIGREESTRING; break;
      case '1':
        datatype = ACCESS_P1PEDIGREEID; break;
      case '2':
        datatype = ACCESS_P2PEDIGREEID; break;
      default:
        break;
      }
      break;
    default:
      switch (c2) {
      case '1':
        datatype = ACCESS_P1NAME; break;
      case '2':
        datatype = ACCESS_P2NAME; break;
      default:
        break;
      }
      break;
    }
  default: break;
  }

	if (datatype == ACCESS_UNKNOWNTYPE) {
	  error("`data.type` parameter is not a valid option");
	}
	
	// Collect group sizes
	R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	size_t* gsizes = (size_t*)R_alloc(glen, sizeof(size_t));
	size_t cumulativesize = check_group_sizes(d, glen, groups, gsizes);
	
	// Create output vector and fill.
	// Got some real mixity of methods of 'how we fill the output vector' here. Sorry.
	SEXP out; size_t outi; int* outci; double* outcr; 
	switch (datatype) {
	case ACCESS_BV: // calls a get_group function per group
	  if (d->n_eff_sets <= 0) { 
	    error("Need to load at least one set of marker effects before requesting breeding values\n"); 
	  }
	  
	  GSC_ID_T eff_id = extract_effid(s_eff_set_id, d);
	  
	  out = PROTECT(allocVector(REALSXP, cumulativesize));
	  outcr = REAL(out);
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      get_group_bvs(d, GROUPNUM_IFY(groups[i]),  EFFECTID_IFY(eff_id), gsizes[i], outcr + outi);
	      outi += gsizes[i];
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    outcr[outi] = NA_REAL;
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_GENOTYPE: // uses bidirectional iterator through each group
	  out = PROTECT(allocVector(STRSXP, cumulativesize));
	  outi = 0;
	  char* buf = R_alloc(2*d->genome.n_markers + 1, sizeof(char)); // R_alloc, not R_calloc. R will clean up the memory itself.
	  buf[d->genome.n_markers << 1] = '\0';
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        memcpy(buf,get_alleles(loc),sizeof(*buf)*(d->genome.n_markers << 1));
	        SET_STRING_ELT(out, outi, mkChar(buf));
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    SET_STRING_ELT(out, outi, NA_STRING);
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_INDEX: // uses manual iteration through each group
	  out = PROTECT(allocVector(INTSXP, cumulativesize));
	  outci = INTEGER(out);
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      AlleleMatrix* am = d->m;
	      size_t globalx = 0;
	      do {
	        for (size_t localx = 0; localx < am->n_genotypes; ++localx, ++globalx) {
	          if (am->groups[localx].num == groups[i]) {
	            outci[outi] = globalx;
	            ++outi;
	          }
	        }
	      } while ((am = am->next) != NULL);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    outci[outi] = NA_INTEGER;
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_LABEL: // uses a bidirectional iterator through each group
	  if (d->n_labels <= 0) { error("Need to create at least one custom label before requesting custom label values\n"); }
	  
	  GSC_ID_T label_id = extract_labelid(s_label_id);
	  int label_index = get_index_of_label(d, LABELID_IFY(label_id));
	  if (label_index == GSC_NA) {
	    error("`label` parameter does not match a current existing custom label id\n");
	  }
	  
	  out = PROTECT(allocVector(INTSXP, cumulativesize));
	  outci = INTEGER(out);
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        outci[outi] = get_label_value(loc, label_index);
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    outci[outi] = NA_INTEGER;
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_NAME: // uses bidirectional iterator through each group
	  out = PROTECT(allocVector(STRSXP, cumulativesize));
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        if (get_name(loc) == NULL) {
			      SET_STRING_ELT(out, outi, NA_STRING);
	          
	        } else {
	          SET_STRING_ELT(out, outi, mkChar(get_name(loc)));
	        }
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
  	  SET_STRING_ELT(out, outi, NA_STRING);
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_PEDIGREEID: // uses bidirectional iterator through each group
	  out = PROTECT(allocVector(INTSXP, cumulativesize));
	  outci = INTEGER(out);
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        outci[outi] = get_id(loc).id;
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    outci[outi] = NA_INTEGER;
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_P1PEDIGREEID: // bidirectional iterator through group
	case ACCESS_P2PEDIGREEID:
	  out = PROTECT(allocVector(INTSXP, cumulativesize));
	  outci = INTEGER(out);
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        PedigreeID parentid = c2 == '1' ? get_first_parent(loc) : get_second_parent(loc);
	        outci[outi] = parentid.id;
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    outci[outi] = NA_INTEGER;
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_P1NAME: // bidirectional iterator through group
	case ACCESS_P2NAME:
	  out = PROTECT(allocVector(STRSXP, cumulativesize));
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        PedigreeID parentid = c2 == '1' ? get_first_parent(loc) : get_second_parent(loc);
	        char* n = NULL;
	        if (parentid.id > 0) { n = gsc_get_name_of_id(d->m, parentid); }
	        if (n == NULL) {
	          SET_STRING_ELT(out, outi, NA_STRING);
	          
	        } else {
	          SET_STRING_ELT(out, outi, mkChar(n));
	        }
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    SET_STRING_ELT(out, outi, NA_STRING);
	  }
	  UNPROTECT(1);
	  return out;
	  
	case ACCESS_PEDIGREESTRING: // calls a get_group function per group
	  out = PROTECT(allocVector(STRSXP, cumulativesize));
	  outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      char* buf[gsizes[i]];
	      get_group_pedigrees(d, GROUPNUM_IFY(groups[i]), gsizes[i], buf);
	      for (int j = 0; j < gsizes[i]; ++j, ++outi) {
	        SET_STRING_ELT(out, outi, mkChar(buf[j]));
	        free(buf[j]);
	        
	      }
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    SET_STRING_ELT(out,outi,NA_STRING);
	  }
	  UNPROTECT(1);
	  return out;
	  
	default:
	  error("`data.type` parameter is not a valid option");
	}
}

SEXP SXP_see_marker_names(SEXP exd) {
  SimData* d = (SimData*) R_ExternalPtrAddr(exd);
  
  SEXP out = PROTECT(allocVector(STRSXP, d->genome.n_markers));
  for (size_t i = 0; i < d->genome.n_markers; ++i) {
    if (d->genome.marker_names == NULL) {
      SET_STRING_ELT(out,i,NA_STRING);
    } else {
      SET_STRING_ELT(out,i,mkChar(d->genome.marker_names[i]));
    }
  }
  UNPROTECT(1);
  return out;
}

SEXP SXP_see_markerblocks(SEXP exd, SEXP exb) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	check_if_is_markerblocks(exb);
	MarkerBlocks* b = (MarkerBlocks*) R_ExternalPtrAddr(exb);
	if (b->num_blocks == NA_GENOLEN) {
		error("Invalid number of blocks");
	}
	
	GSC_GENOLEN_T nrows = 0;
	for (GSC_GENOLEN_T i = 0; i < b->num_blocks; ++i) {
		nrows += b->num_markers_in_block[i];
	}
	
	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SEXP col1 = PROTECT(allocVector(STRSXP, nrows));
	SEXP col2 = PROTECT(allocVector(INTSXP, nrows));
	int* c_col2 = INTEGER(col2);
	GSC_GENOLEN_T i_overall = 0;
	for (GSC_GENOLEN_T i = 0; i < b->num_blocks; ++i) {
		for (GSC_GENOLEN_T j = 0; j < b->num_markers_in_block[i]; ++j) {
			GSC_GENOLEN_T markerix = b->markers_in_block[i][j];
			if (markerix < d->genome.n_markers && markerix >= 0) {
				SET_STRING_ELT(col1, i_overall, mkChar(d->genome.marker_names[markerix]));
			} else {
				SET_STRING_ELT(col1, i_overall, mkChar(""));
			}
			c_col2[i_overall] = i+1;
			++i_overall;
		}
	}

	SET_VECTOR_ELT(out, 0, col1);
	SET_VECTOR_ELT(out, 1, col2);
	UNPROTECT(3);
	return out;
}

SEXP SXP_see_local_gebvs(SEXP exd, SEXP exblocks, SEXP s_groups, SEXP s_eff_set_id) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	check_if_is_markerblocks(exblocks);
	MarkerBlocks* b = (MarkerBlocks*) R_ExternalPtrAddr(exblocks);
	
	if (d->n_eff_sets <= 0) { error("Need to load at least one set of marker effects before requesting breeding values\n"); }
	GSC_ID_T eff_id = extract_effid(s_eff_set_id, d);
	
	R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	size_t* gsizes = (size_t*)R_alloc(glen, sizeof(size_t));
	size_t cumulativesize = check_group_sizes(d, glen, groups, gsizes);

	SEXP out = PROTECT(allocMatrix(REALSXP, 2*cumulativesize, b->num_blocks));
	double* outc = REAL(out);
	R_xlen_t outi = 0;
	for (R_xlen_t i = 0; i < glen; ++i) {
	  if (gsizes[i] > 0) {
		gsc_DecimalMatrix m = gsc_calculate_local_bvs(d, GROUPNUM_IFY(groups[i]), *b, EFFECTID_IFY(eff_id));
		for (size_t d1 = 0; d1 < m.dim1; ++d1,++outi) {
			for (size_t d2 = 0; d2 < m.dim2; ++d2) {
				outc[outi + 2*cumulativesize*d2] = m.matrix[d1][d2];
			}
		}
		gsc_delete_dmatrix(&m);
	  }
	}
	
	UNPROTECT(1);
	return(out);
}


SEXP SXP_see_group_gene_data(SEXP exd, SEXP s_groups, SEXP s_countAllele, SEXP s_unknownAllele) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  // Collect group sizes
  R_xlen_t glen = xlength(s_groups);
  int *groups = INTEGER(s_groups);
  size_t* gsizes = (size_t*)R_alloc(glen, sizeof(size_t));
  size_t cumulativesize = check_group_sizes(d, glen, groups, gsizes);
  
  // bidirectional iterator through group
	if (asChar(s_countAllele) == NA_STRING) {
	  char printUnknown;
	  if (asChar(s_unknownAllele) == NA_STRING || !isprint(CHAR(asChar(s_unknownAllele))[0])) {
	    printUnknown = '-';
	    Rprintf("NOTE! Defaulting to printing '\\0' alleles as a dash '-'");
	  } else {
	    printUnknown = CHAR(asChar(s_unknownAllele))[0];
	  }
	  
		// Give the pairs of alleles
		SEXP out = PROTECT(allocMatrix(STRSXP, d->genome.n_markers, cumulativesize));
		size_t outi = 0;
		for (R_xlen_t i = 0; i < glen; ++i) {
		  if (gsizes[i] > 0) {
		    BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
		    GenoLocation loc = next_forwards(&it);
		    while (IS_VALID_LOCATION(loc)) {
		      char* geno = get_alleles(loc);
		      for (size_t j = 0; j < d->genome.n_markers; ++j) {
		        char alleles[] = {geno[2*j], geno[2*j + 1], '\0'};
		        if (alleles[0] == '\0') { alleles[0] = printUnknown; }
		        if (alleles[1] == '\0') { alleles[1] = printUnknown; }
		        SET_STRING_ELT(out, j + d->genome.n_markers*outi, mkChar(alleles));
		      }
		      ++outi;
		      loc = next_forwards(&it);
		    }
		    delete_bidirectional_iter(&it);
		  }
		}
		for (; outi < cumulativesize; ++outi) {
		  for (size_t j = 0; j < d->genome.n_markers; ++j) {
		    SET_STRING_ELT(out, j + d->genome.n_markers*outi, NA_STRING);
		  }
		}
		UNPROTECT(1);
		return out;

	} else {
		// Give the counts of allele c
		char c = CHAR(asChar(s_countAllele))[0];
	  
	  SEXP out = PROTECT(allocMatrix(INTSXP, d->genome.n_markers, cumulativesize));
	  int* outc = INTEGER(out);
	  size_t outi = 0;
	  for (R_xlen_t i = 0; i < glen; ++i) {
	    if (gsizes[i] > 0) {
	      BidirectionalIterator it = create_bidirectional_iter(d,GROUPNUM_IFY(groups[i]));
	      GenoLocation loc = next_forwards(&it);
	      while (IS_VALID_LOCATION(loc)) {
	        char* geno = get_alleles(loc);
	        for (size_t j = 0; j < d->genome.n_markers; ++j) {
	          int count = 0;
	          if (geno[2*j] == c)     { ++count; }
	          if (geno[2*j + 1] == c) { ++count; }
	          outc[j + d->genome.n_markers*outi] = count;
	        }
	        ++outi;
	        loc = next_forwards(&it);
	      }
	      delete_bidirectional_iter(&it);
	    }
	  }
	  for (; outi < cumulativesize; ++outi) {
	    for (size_t j = 0; j < d->genome.n_markers; ++j) {
	      outc[j + d->genome.n_markers*outi] = NA_INTEGER;
	    }
	  }
	  UNPROTECT(1);
	  return out;
	}
}


SEXP SXP_change_name_to_values(SEXP exd, SEXP s_values, SEXP s_group, SEXP s_start) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (TYPEOF(s_values) != STRSXP) {
		error("`values` are invalid: must be strings");
	}
	R_xlen_t len = xlength(s_values);
	char* names[len];
	for (R_xlen_t i = 0; i < len; ++i) {
		names[i] = R_Calloc(sizeof(char)*(NAME_LENGTH + 1), char);
		strncpy(names[i], CHAR(STRING_ELT(s_values, i)), sizeof(char)*NAME_LENGTH);
		//names[i][NAME_LENGTH] = '\0'; // terminate 'em. Just in case they're trying
		// to spill over NAME_LENGTH
	}

	int group = asInteger(s_group);
	if (group == NA_INTEGER || group < 0) {
		error("`group` is invalid (wrong type or negative)\n");
	}

	int startIndex = asInteger(s_start);
	if (startIndex == NA_INTEGER || startIndex < 0) {
		error("`startIndex` is invalid (wrong type or negative)\n");
	}

	change_names_to_values(d, GROUPNUM_IFY(group), startIndex, len, (const char**) names);

	for (int i = 0; i < len; ++i) {
		R_Free(names[i]);
	}

	return ScalarInteger(0);
}

SEXP SXP_change_allele_symbol(SEXP exd, SEXP s_markername, SEXP s_from, SEXP s_to) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char* c_markername;
	if (xlength(s_markername) == 0 || asChar(s_markername) == NA_STRING) {
		c_markername = NULL;
	} else {
		c_markername = CHAR(asChar(s_markername));
	}
	
	char from = CHAR(asChar(s_from))[0];
	char to = CHAR(asChar(s_to))[0];
	
	change_allele_symbol(d, c_markername, from, to);
	
	return ScalarInteger(0);
}


/*------------------------ Deletors ------------------------*/
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

void SXP_delete_markerblocks(SEXP sd) {
	//Rprintf("Garbage collecting MarkerBlocks...\n");
	//check_if_is_markerblocks(sd);
	MarkerBlocks* b = (MarkerBlocks*) R_ExternalPtrAddr(sd);
	delete_markerblocks(b);
	R_ClearExternalPtr(sd);
	return;
}

SEXP SXP_delete_group(SEXP exd, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t n = xlength(s_groups);
	int *groups = INTEGER(s_groups);

	for (R_xlen_t i = 0; i < n; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 1) {
			error("Entry %d in `group` parameter is of invalid type\n", i + 1);
		}
		delete_group(d, GROUPNUM_IFY(groups[i]));
	}

	return ScalarInteger(0);
}

SEXP SXP_delete_label(SEXP exd, SEXP s_labels) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t n = xlength(s_labels);
	int *labels = INTEGER(s_labels);

	for (R_xlen_t i = 0; i < n; ++i) {
		if (labels[i] == NA_INTEGER || labels[i] < 1) {
			error("Entry %d in `labels` parameter is of invalid type\n", i + 1);
		}
		delete_label(d, LABELID_IFY(labels[i]));
	}

	return ScalarInteger(0);
}

SEXP SXP_delete_eff_set(SEXP exd, SEXP s_eff_sets) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t n = xlength(s_eff_sets);
	int *eff_sets = INTEGER(s_eff_sets);

	for (R_xlen_t i = 0; i < n; ++i) {
		if (eff_sets[i] == NA_INTEGER || eff_sets[i] < 1) {
			error("Entry %d in `effect_sets` parameter is of invalid type\n", i + 1);
		}
		delete_eff_set(d, EFFECTID_IFY(eff_sets[i]));
	}

	return ScalarInteger(0);
}

SEXP SXP_delete_recombination_map(SEXP exd, SEXP s_maps) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t n = xlength(s_maps);
	int *maps = INTEGER(s_maps);

	for (R_xlen_t i = 0; i < n; ++i) {
		if (maps[i] == NA_INTEGER || maps[i] < 1) {
			error("Entry %d in `maps` parameter is of invalid type\n", i + 1);
		}
		delete_recombination_map(d, MAPID_IFY(maps[i]));
	}

	return ScalarInteger(0);
}


/*------------------------- Groups -------------------------*/
SEXP SXP_combine_groups(SEXP exd, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t len = xlength(s_groups);

	GroupNum groups[len];
	convertVECSXP_to_GroupNum(s_groups, groups);

	return ScalarInteger(combine_groups(d, len, groups).num);
}


SEXP SXP_make_group_from(SEXP exd, SEXP s_indexes) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t n = xlength(s_indexes);
	unsigned int *ns = (unsigned int*) INTEGER(s_indexes);
	/*for (R_xlen_t i = 0; i < n; ++i) {
		if (ns[i] == NA_INTEGER || ns[i] < 0) {
			error("The `indexes` vector contains at least one invalid index\n");
		}
	}*/ // <- check is carried out in R level these days

	return ScalarInteger(make_group_from(d, n, ns).num);
}


SEXP SXP_break_group_into_individuals(SEXP exd, SEXP s_group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
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

SEXP SXP_break_group_into_families(SEXP exd, SEXP s_group) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
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

SEXP SXP_break_group_into_halfsib_families(SEXP exd, SEXP s_group, SEXP s_parent) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
	}

	int parent_num = asInteger(s_parent);
	if (parent_num == NA_INTEGER || parent_num < 0) {
		error("`parent` parameter is of invalid type\n");
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

SEXP SXP_break_group_randomly(SEXP exd, SEXP s_group, SEXP s_n) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
	}

	int n_groups = asInteger(s_n);
	if (n_groups == NA_INTEGER || n_groups < 0) {
		error("`n` parameter is of invalid type\n");
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

SEXP SXP_break_group_evenly(SEXP exd, SEXP s_group, SEXP s_n) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
	}

	int n_groups = asInteger(s_n);
	if (n_groups == NA_INTEGER || n_groups < 0) {
		error("`n` parameter is of invalid type\n");
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

SEXP SXP_break_group_into_buckets(SEXP exd, SEXP s_group, SEXP s_buckets) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
	}

	if (TYPEOF(s_buckets) != INTSXP) { // check param name
		error("`buckets` parameter must be a vector of integers\n");
	}

	R_xlen_t n_groups = xlength(s_buckets)+1;
	int* counts = INTEGER(s_buckets);

	GroupNum results[n_groups];
	split_into_buckets(d, GROUPNUM_IFY(group_id), n_groups, (const unsigned int*) counts, results);

	SEXP out = PROTECT(allocVector(INTSXP, n_groups));
	int* outc = INTEGER(out);
	for (R_xlen_t i = 0; i < n_groups; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_break_group_into_probabilities(SEXP exd, SEXP s_group, SEXP s_probs) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	int group_id = asInteger(s_group);
	if (group_id == NA_INTEGER || group_id < 0) {
		error("`group` parameter is of invalid type\n");
	}

	if (TYPEOF(s_probs) != REALSXP) { // check param name
		error("`probabilies` parameter must be a vector of decimals\n");
	}

	R_xlen_t n_groups = xlength(s_probs)+1;
	double* probabilities = REAL(s_probs);

	GroupNum results[n_groups];
	memset(results, 0, sizeof(GroupNum) * n_groups);
	split_by_probabilities(d, GROUPNUM_IFY(group_id), n_groups, probabilities, results);

	SEXP out = PROTECT(allocVector(INTSXP, n_groups));
	int* outc = INTEGER(out);
	for (R_xlen_t i = 0; i < n_groups; ++i) {
		outc[i] = results[i].num;
	}

	UNPROTECT(1);
	return out;
}

SEXP SXP_break_group_by_label_value(SEXP exd, SEXP s_label, SEXP s_value, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T label = extract_labelid(s_label);

	R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);

	int value = asInteger(s_value);
	if (value == NA_INTEGER) {
		error("`value` input must be an integer\n");
	}

	if (glen == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			error("`group` parameter is invalid: negative");
		} else {
			return ScalarInteger(split_by_label_value(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), value).num);
		}

	} else {

		GroupNum newSubGroups[glen];
		int j = 0;
		for (R_xlen_t i = 0; i < glen; ++i) {
			if (groups[i] < 1) {
				Rprintf("Entry %lu in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				newSubGroups[j] = split_by_label_value(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), value);
				++j;
			}
		}
		return ScalarInteger(combine_groups(d, j, newSubGroups).num);
	}
}


SEXP SXP_break_group_by_label_range(SEXP exd, SEXP s_label, SEXP s_lowbound, SEXP s_highbound,
		SEXP s_groups) {

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T label = extract_labelid(s_label);

	R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);

	int low = asInteger(s_lowbound);
	if (low == NA_INTEGER) {
		error("`rangeLowEnd` input must be an integer\n");
	}

	int high = asInteger(s_highbound);
	if (high == NA_INTEGER) {
		error("`rangeHighEnd` input must be an integer\n");
	}

	if (glen == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			error("`group` parameter is invalid: negative");
		} else {
			return ScalarInteger(split_by_label_range(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), low, high).num);
		}

	} else {

		GroupNum newSubGroups[glen];
		int j = 0;
		for (R_xlen_t i = 0; i < glen; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %lu in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				newSubGroups[j] = split_by_label_range(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), low, high);
				++j;
			}
		}
		return ScalarInteger(combine_groups(d, j, newSubGroups).num);
	}
}


SEXP SXP_change_label_default(SEXP exd, SEXP s_labels, SEXP s_defaults) {

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

  R_xlen_t lblen = xlength(s_labels);
	int *labels = INTEGER(s_labels);

	R_xlen_t dlen = xlength(s_defaults);
	int *defaults = INTEGER(s_defaults);

	// Find the minimum length of the vectors and discard everything further along
	R_xlen_t matchedLen = lblen < dlen ? lblen : dlen;

	for (R_xlen_t i = 0; i < matchedLen; ++i) {
		if (labels[i] == NA_INTEGER || labels[i] == 0) {
			error("entry in `label` vector is invalid");
		} else {
			change_label_default(d, LABELID_IFY(labels[i]), defaults[i]);
		}
	}
	Rprintf("Set the defaults of %lu labels.\n", matchedLen);
	return ScalarInteger(0);
}

SEXP SXP_change_label_to_values(SEXP exd, SEXP s_label, SEXP s_values, SEXP s_group, SEXP s_skip) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T label = extract_labelid(s_label);
	
	R_xlen_t vlen = xlength(s_values);
	int* values = INTEGER(s_values);

	int group = asInteger(s_group);
	if (group == NA_INTEGER) {
		group = 0;
	} else if (group < 0) {
		error("`group` is invalid (wrong type or negative)\n");
	}

	int startIndex = asInteger(s_skip);
	if (startIndex == NA_INTEGER || startIndex < 0) {
		error("`skip` is invalid (wrong type or negative)\n");
	}

	change_label_to_values(d, GROUPNUM_IFY(group), (GLOBALX_T) startIndex, LABELID_IFY(label), vlen, values);

	return ScalarInteger(0);
}

SEXP SXP_change_label_by_amount(SEXP exd, SEXP s_label, SEXP s_incr, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T label = extract_labelid(s_label);

	int incr = asInteger(s_incr);
	if (incr == NA_INTEGER) {
		error("`amount` parameter is invalid: must be an integer");
	}

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);

	if (len == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}
		if (groups[0] < 0) {
			Rprintf("`group` parameter is invalid: negative");
		} else {
			change_label_by_amount(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), incr);
		}

	} else {
		for (R_xlen_t i = 0; i < len; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %lu in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				change_label_by_amount(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), incr);
			}
		}
	}

	return ScalarInteger(0);
}

SEXP SXP_change_label_to_this(SEXP exd, SEXP s_label, SEXP s_const, SEXP s_groups) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	GSC_ID_T label = extract_labelid(s_label);

	int num = asInteger(s_const);
	if (num == NA_INTEGER) {
		error("`value` parameter is invalid: must be an integer");
	}

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);

	if (len == 1) {
		if (groups[0] == NA_INTEGER) {
			groups[0] = 0;
		}

		if (groups[0] < 0) {
			Rprintf("`group` parameter is invalid: negative");
		} else {
			change_label_to(d, GROUPNUM_IFY(groups[0]), LABELID_IFY(label), num);
		}

	} else {
		for (R_xlen_t i = 0; i < len; ++i) {
			if (groups[i] < 1) {
				Rprintf("entry %lu in the `group` vector is an invalid group number (0 or negative)", i + 1);
			} else {
				change_label_to(d, GROUPNUM_IFY(groups[i]), LABELID_IFY(label), num);
			}
		}
	}

	return ScalarInteger(0);
}

SEXP SXP_break_group_by_GEBV_num(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_number, SEXP s_bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: negative at index %lu\n", i+1); }
	}

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

	int num_to_select = asInteger(s_number);
	if (num_to_select == NA_INTEGER || num_to_select < 0) {
		error("`number` parameter is of invalid type\n");
	}

	int want_low = asLogical(s_bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type\n"); }
	
	if (len == 1) {
		return ScalarInteger(split_by_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id), num_to_select, want_low).num);
	} else {
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (R_xlen_t i = 0; i < len; ++i) {
			outc[i] = split_by_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id), num_to_select, want_low).num;
		}
		UNPROTECT(1);
		return out;
	}

}

SEXP SXP_break_group_by_GEBV_percent(SEXP exd, SEXP s_groups, SEXP s_eff_set, SEXP s_percent, SEXP s_bestIsLow) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	if (d->n_eff_sets <= 0) { error("Need to load effect values before running this function\n"); }

	R_xlen_t len = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < len; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: negative at index %lu\n", i+1); }
	}

	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

	double pc_to_select = asReal(s_percent);
	if (ISNA(pc_to_select) || pc_to_select < 0) {
		error("`percentage` parameter is of invalid type\n");
	}

	int want_low = asLogical(s_bestIsLow);
	if (want_low == NA_LOGICAL) { error("`low.score.best` parameter is of invalid type\n"); }

	if (len == 1) {
		int group_size = get_group_size(d, GROUPNUM_IFY(groups[0]));
		int num_to_select = group_size * pc_to_select / 100; // integer division, so take the floor

		return ScalarInteger(split_by_bv(d, GROUPNUM_IFY(groups[0]), EFFECTID_IFY(eff_id), num_to_select, want_low).num);
	} else {
		int num_to_select;

		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, len));
		int* outc = INTEGER(out);
		for (R_xlen_t i = 0; i < len; ++i) {
			num_to_select = get_group_size(d, GROUPNUM_IFY(groups[i])) * pc_to_select / 100;
			outc[i] = split_by_bv(d, GROUPNUM_IFY(groups[i]), EFFECTID_IFY(eff_id), num_to_select, want_low).num;
		}
		UNPROTECT(1);
		return out;
	}

}

/*---------------------- Progression ----------------------*/
GenOptions SXP_create_genoptions(SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions go = BASIC_OPT;
	int b;

	b = asLogical(s_name);
	if (b == NA_LOGICAL) { error("`name` parameter is of invalid type\n"); }
	go.will_name_offspring = b;

	if (isNull(s_namePrefix)) {
		go.offspring_name_prefix = NULL;
	} else {
		go.offspring_name_prefix = CHAR(asChar(s_namePrefix));
	}

	b = asInteger(s_familySize);
	if (b == NA_INTEGER) { error("`offspring` parameter is of invalid type\n"); }
	go.family_size = b;

	b = asLogical(s_trackPedigree);
	if (b == NA_LOGICAL) { error("`track.pedigree` parameter is of invalid type\n"); }
	go.will_track_pedigree = b;

	b = asLogical(s_giveIds);
	if (b == NA_LOGICAL) { error("`give.ids` parameter is of invalid type\n"); }
	go.will_allocate_ids = b;

	if (isNull(s_filePrefix)) {
		go.filename_prefix = NULL;
	} else {
		go.filename_prefix = CHAR(asChar(s_filePrefix));
	}

	b = asLogical(s_savePedigree);
	if (b == NA_LOGICAL) { error("`save.pedigree` parameter is of invalid type\n"); }
	go.will_save_pedigree_to_file = b;
	b = asInteger(s_saveEffects);
	if (b == NA_INTEGER) { error("`save.gebv` parameter is of invalid type\n"); }
	go.will_save_bvs_to_file = EFFECTID_IFY(b);
	b = asLogical(s_saveGenes);
	if (b == NA_LOGICAL) { error("`save.genotype` parameter is of invalid type\n"); }
	go.will_save_alleles_to_file = b;

	b = asLogical(s_retain);
	if (b == NA_LOGICAL) { error("`retain` parameter is of invalid type\n"); }
	go.will_save_to_simdata = b;

	return go;
}

SEXP SXP_make_random_crosses(SEXP exd, SEXP s_groups, SEXP s_crosses, SEXP s_cap, SEXP s_map,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

  R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}
	
	R_xlen_t maplen = xlength(s_map);
	if (maplen != 1 && maplen != glen) {
		error("Cannot match up provided recombination maps to groups: lists differ in length");
	}
	int *map = INTEGER(s_map);

	int n = asInteger(s_crosses);
	if (n == NA_INTEGER) { error("`n.crosses` parameter is invalid\n"); }

	int cap = asInteger(s_cap);
	if (cap == NA_INTEGER) { error("`cap` parameter is invalid\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_random_crosses(d, GROUPNUM_IFY(groups[0]), n, cap, MAPID_IFY(parse_mapid(map[0],d)), g).num);

	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		if (maplen == 1) {
			gsc_MapID gsmap = MAPID_IFY(parse_mapid(map[0],d));
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_random_crosses(d, GROUPNUM_IFY(groups[i]), n, cap, gsmap, g).num;
			}
		} else {
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_random_crosses(d, GROUPNUM_IFY(groups[i]), n, cap, MAPID_IFY(parse_mapid(map[i],d)), g).num;
			}
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_make_random_crosses_between(SEXP exd, SEXP s_group1, SEXP s_group2, SEXP s_cap1, SEXP s_cap2, SEXP s_map1, SEXP s_map2,
		SEXP s_crosses, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	int group1_c = asInteger(s_group1);
	if (group1_c == NA_INTEGER || group1_c < 0) { error("The parameter `group1` is invalid\n"); }
	int group2_c = asInteger(s_group2);
	if (group2_c == NA_INTEGER || group2_c < 0) { error("The parameter `group2` is invalid\n"); }

	int cap1 = asInteger(s_cap1);
	if (cap1 == NA_INTEGER) { error("The parameter `cap1` is invalid\n"); }
	int cap2 = asInteger(s_cap2);
	if (cap2 == NA_INTEGER) { error("The parameter `cap2` is invalid\n"); }
	
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	if (xlength(s_map1) > 1 || xlength(s_map2) > 1) {
		Rprintf("NOTE! More than one recombination map per group provided. Only the first recombination map will be used");
	}
	GSC_ID_T map1 = extract_mapid(s_map1,d,"map1");
	GSC_ID_T map2 = extract_mapid(s_map2,d,"map2");

	int n = asInteger(s_crosses);
	if (n == NA_INTEGER) { error("`n.crosses` parameter is invalid\n"); }

	return ScalarInteger(make_random_crosses_between(d, GROUPNUM_IFY(group1_c), GROUPNUM_IFY(group2_c), n, cap1, cap2, MAPID_IFY(map1), MAPID_IFY(map2), g).num);
}

SEXP SXP_make_targeted_crosses(SEXP exd, SEXP s_firstparents, SEXP s_secondparents, SEXP s_map1, SEXP s_map2,
		SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (xlength(s_firstparents) != xlength(s_secondparents)) {
		error("Parent vectors must be the same length\n");
	}

	R_xlen_t ncrosses = xlength(s_firstparents);

	unsigned int combinations[2][ncrosses];
	char pname[NAME_LENGTH+1];

	if (TYPEOF(s_firstparents) == STRSXP) {
		for (R_xlen_t i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(s_firstparents, i)), sizeof(char)*NAME_LENGTH);
			combinations[0][i] = gsc_get_index_of_name(d->m, pname);
		}
	} else if (TYPEOF(s_firstparents) == INTSXP) {
		int* indexes = INTEGER(s_firstparents);
		for (R_xlen_t i = 0; i < ncrosses; ++i) {
			combinations[0][i] = indexes[i];
		}
	} else {
		error("first.parents must be a vector of strings or integers\n");
	}

	if (TYPEOF(s_secondparents) == STRSXP) {
		for (R_xlen_t i = 0; i < ncrosses; ++i) {
			strncpy(pname, CHAR(STRING_ELT(s_secondparents, i)), sizeof(char)*NAME_LENGTH);
			combinations[1][i] = gsc_get_index_of_name(d->m, pname);
		}
	} else if (TYPEOF(s_secondparents) == INTSXP) {
		int* indexes = INTEGER(s_secondparents);
		for (R_xlen_t i = 0; i < ncrosses; ++i) {
			combinations[1][i] = indexes[i];
		}
	} else {
		error("second.parents must be a vector of strings or integers\n");
	}
	
	R_xlen_t filleri = 0;
	for (R_xlen_t checkeri = 0; checkeri < ncrosses; ++checkeri) {
		if (combinations[0][checkeri] < 0 || combinations[1][checkeri] < 0) {
			Rprintf("NOTE! Names or indexes at row %lu of the crossing plan are invalid\n", checkeri);
		} else {
			if (checkeri != filleri) { // copy the information at "checkeri" to position "filleri"
				combinations[0][filleri] = combinations[0][checkeri];
				combinations[1][filleri] = combinations[1][checkeri];
			}
			filleri++;
		}
	}
	
	if (xlength(s_map1) > 1 || xlength(s_map2) > 1) {
		Rprintf("NOTE! More than one recombination map per group provided. Only the first recombination map each will be used");
	}
	GSC_ID_T map1 = extract_mapid(s_map1,d,"map1");
	GSC_ID_T map2 = extract_mapid(s_map2,d,"map2");

	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
								 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
								 s_saveGenes, s_retain);

	return ScalarInteger(make_targeted_crosses(d, ncrosses, combinations[0], combinations[1], MAPID_IFY(map1), MAPID_IFY(map2), g).num);

}

SEXP SXP_make_crosses_from_file(SEXP exd, SEXP s_filename, SEXP s_map1, SEXP s_map2, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char* filename = CHAR(asChar(s_filename));
	
	if (xlength(s_map1) > 1 || xlength(s_map2) > 1) {
		Rprintf("NOTE! More than one recombination map per parent provided. Only the first recombination map each will be used");
	}
	GSC_ID_T map1 = extract_mapid(s_map1,d,"map1");
	GSC_ID_T map2 = extract_mapid(s_map2,d,"map2");

	return ScalarInteger(make_crosses_from_file(d, filename, MAPID_IFY(map1),MAPID_IFY(map2), g).num);
}

SEXP SXP_make_double_crosses_from_file(SEXP exd, SEXP s_filename, SEXP s_map1, SEXP s_map2, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	const char* filename = CHAR(asChar(s_filename));
	
	if (xlength(s_map1) > 1 || xlength(s_map2) > 1) {
		Rprintf("NOTE! More than one recombination map per parent provided. Only the first recombination map each will be used");
	}
	GSC_ID_T map1 = extract_mapid(s_map1,d,"map1");
	GSC_ID_T map2 = extract_mapid(s_map2,d,"map2");

	return ScalarInteger(make_double_crosses_from_file(d, filename, MAPID_IFY(map1), MAPID_IFY(map2), g).num);

}

SEXP SXP_make_all_unidirectional_crosses(SEXP exd, SEXP s_groups, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
  R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}
	
	R_xlen_t maplen = xlength(s_map);
	if (maplen != 1 && maplen != glen) {
		error("Cannot match up provided recombination maps to groups: lists differ in length");
	}
	int *map = INTEGER(s_map);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_all_unidirectional_crosses(d, GROUPNUM_IFY(groups[0]), MAPID_IFY(parse_mapid(map[0],d)), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		if (maplen == 1) {
			gsc_MapID gsmap = MAPID_IFY(parse_mapid(map[0],d));
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_all_unidirectional_crosses(d, GROUPNUM_IFY(groups[i]), gsmap, g).num;
			}
		} else {
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_all_unidirectional_crosses(d, GROUPNUM_IFY(groups[i]), MAPID_IFY(parse_mapid(map[i],d)), g).num;
			}
	
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_self_n_times(SEXP exd, SEXP s_groups, SEXP s_ngen, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
  R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}
	
	R_xlen_t maplen = xlength(s_map);
	if (maplen != 1 && maplen != glen) {
		error("Cannot match up provided recombination maps to groups: lists differ in length");
	}
	int *map = INTEGER(s_map);

	int cn = asInteger(s_ngen);
	if (cn < 0 || cn == NA_INTEGER) { error("`n` parameter is invalid\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(self_n_times(d, cn, GROUPNUM_IFY(groups[0]), MAPID_IFY(parse_mapid(map[0],d)), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		if (maplen == 1) {
			gsc_MapID gsmap = MAPID_IFY(parse_mapid(map[0],d));
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = self_n_times(d, cn, GROUPNUM_IFY(groups[i]), gsmap, g).num;
			}
		} else {
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = self_n_times(d, cn, GROUPNUM_IFY(groups[i]), MAPID_IFY(parse_mapid(map[i],d)), g).num;
			}			
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_make_doubled_haploids(SEXP exd, SEXP s_groups, SEXP s_map, SEXP s_name, SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
  R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (R_xlen_t i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}
	
	R_xlen_t maplen = xlength(s_map);
	if (maplen != 1 && maplen != glen) {
		error("Cannot match up provided recombination maps to groups: lists differ in length");
	}
	int *map = INTEGER(s_map);

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_doubled_haploids(d, GROUPNUM_IFY(groups[0]), MAPID_IFY(parse_mapid(map[0],d)), g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		if (maplen == 1) {
			gsc_MapID gsmap = MAPID_IFY(parse_mapid(map[0],d));
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_doubled_haploids(d, GROUPNUM_IFY(groups[i]), gsmap, g).num;
			}
		} else {
			for (R_xlen_t i = 0; i < glen; ++i) {
				outc[i] = make_doubled_haploids(d, GROUPNUM_IFY(groups[i]), MAPID_IFY(parse_mapid(map[i],d)), g).num;
			}			
		}
		UNPROTECT(1);
		return out;
	}
}

SEXP SXP_make_clones(SEXP exd, SEXP s_groups, SEXP s_inherit_name, SEXP s_name,
		SEXP s_namePrefix, SEXP s_familySize,
		SEXP s_trackPedigree, SEXP s_giveIds, SEXP s_filePrefix, SEXP s_savePedigree,
		SEXP s_saveEffects, SEXP s_saveGenes, SEXP s_retain) {
	GenOptions g = SXP_create_genoptions(s_name, s_namePrefix, s_familySize, s_trackPedigree,
									 s_giveIds, s_filePrefix, s_savePedigree, s_saveEffects,
									 s_saveGenes, s_retain);
  R_xlen_t glen = xlength(s_groups);
	int *groups = INTEGER(s_groups);
	for (int i = 0; i < glen; ++i) {
		if (groups[i] == NA_INTEGER || groups[i] < 0) { error("The contents of `groups` is invalid: at index %lu\n", i+1); }
	}

	int inherit_name = asLogical(s_inherit_name);
	if (inherit_name == NA_LOGICAL) { error("`inherit.name` parameter is of invalid type: should be logical\n"); }

	SimData* d = (SimData*) R_ExternalPtrAddr(exd);

	if (glen == 1) {
		return ScalarInteger(make_clones(d, GROUPNUM_IFY(groups[0]), inherit_name, g).num);
	} else {
		// Get an R vector of the same length as the number of new size 1 groups created
		SEXP out = PROTECT(allocVector(INTSXP, glen));
		int* outc = INTEGER(out);
		for (R_xlen_t i = 0; i < glen; ++i) {
			outc[i] = make_clones(d, GROUPNUM_IFY(groups[i]), inherit_name, g).num;
		}
		UNPROTECT(1);
		return out;
	}
}


/*----------------------- Printers ----------------------*/
SEXP SXP_save_genotypes(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_markersasrows) {
	const char* filename = CHAR(asChar(s_filename));
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	GroupNum gr = GROUPNUM_IFY(asInteger(s_group));
	
	int markers_as_rows;
	if (isLogical(s_markersasrows) && asLogical(s_markersasrows) != NA_LOGICAL) {
	  markers_as_rows = asLogical(s_markersasrows);
	} else {
	  error("`markers.as.rows` parameter is invalid: must be TRUE or FALSE");
	}
	
	save_genotypes(filename, d, gr, markers_as_rows);
	return ScalarInteger(0);
}

SEXP SXP_save_allele_counts(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_allele, SEXP s_markersasrows) {
	const char* filename = CHAR(asChar(s_filename));
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	GroupNum gr = GROUPNUM_IFY(asInteger(s_group));
	
	const char allele = CHAR(asChar(s_allele))[0];
	int markers_as_rows;
	if (isLogical(s_markersasrows) && asLogical(s_markersasrows) != NA_LOGICAL) {
	  markers_as_rows = asLogical(s_markersasrows);
	} else {
	  error("`markers.as.rows` parameter is invalid: must be TRUE or FALSE");
	}
	
	save_allele_counts(filename, d, gr, allele, markers_as_rows);
	return ScalarInteger(0);
}

SEXP SXP_save_pedigrees(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_fullpedigree) {
	const char* filename = CHAR(asChar(s_filename));
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	GroupNum gr = GROUPNUM_IFY(asInteger(s_group));
	
	int full_pedigree;
	if (isLogical(s_fullpedigree) && asLogical(s_fullpedigree) != NA_LOGICAL) {
	  full_pedigree = asLogical(s_fullpedigree);
	} else {
	  error("`full.pedigree` parameter is invalid: must be TRUE or FALSE");
	}

	save_pedigrees(filename, d, gr, full_pedigree);
	return ScalarInteger(0);
}

SEXP SXP_save_GEBVs(SEXP exd, SEXP s_filename, SEXP s_group, SEXP s_eff_set) {
	const char* filename = CHAR(asChar(s_filename));
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	GroupNum gr = GROUPNUM_IFY(asInteger(s_group));
	
	if (d->n_eff_sets <= 0) { error("Need to load at least one set of marker effects before requesting breeding values\n"); }
	GSC_ID_T eff_id = extract_effid(s_eff_set, d);

  save_bvs(filename, d, gr, EFFECTID_IFY(eff_id));
	return ScalarInteger(0);
}

SEXP SXP_save_local_GEBVs(SEXP exd, SEXP s_filename, SEXP s_blocks, SEXP s_group, SEXP s_eff_set) {
	const char* filename = CHAR(asChar(s_filename));
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	GroupNum gr = GROUPNUM_IFY(asInteger(s_group));
	
	if (d->n_eff_sets <= 0) { error("Need to load at least one set of marker effects before requesting breeding values\n"); }
	GSC_ID_T eff_id = extract_effid(s_eff_set, d);
	
	check_if_is_markerblocks(s_blocks);
	MarkerBlocks* b = (MarkerBlocks*) R_ExternalPtrAddr(s_blocks);
	
	save_local_bvs(filename, d, gr, *b, EFFECTID_IFY(eff_id), 1);
	return ScalarInteger(0);
}


/*--------------------------------Deletors------------------------------------*/
