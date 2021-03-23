#include "sim-loaders.h"

SEXP clear_simdata(SEXP exd) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	delete_simdata(d);
	return ScalarInteger(0);
}

SEXP load_data(SEXP alleleFile, SEXP mapFile, SEXP groupName) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile)), CHAR(asChar(groupName)));
	load_genmap_to_simdata(d, CHAR(asChar(mapFile)));
	
	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	
	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP load_data_weff(SEXP alleleFile, SEXP mapFile, SEXP effectFile, SEXP groupName) {
	SimData* d = create_empty_simdata();
	//d->current_id = 0; // reset ID counts
	load_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile)), CHAR(asChar(groupName)));
	load_genmap_to_simdata(d, CHAR(asChar(mapFile)));
	load_effects_to_simdata(d, CHAR(asChar(effectFile)));

	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	
	SEXP sdptr = PROTECT(R_MakeExternalPtr((void*) d, R_NilValue, R_NilValue));
	R_RegisterCFinalizerEx(sdptr, SXP_delete_simdata, 1);
	UNPROTECT(1);
	return sdptr;
}

SEXP load_more_genotypes(SEXP exd, SEXP alleleFile, SEXP groupName) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	return ScalarInteger(load_more_transposed_genes_to_simdata(d, CHAR(asChar(alleleFile)), 
																CHAR(asChar(groupName))));
	//return ScalarInteger(load_more_transposed_genes_to_simdata(&GlobalSim, CHAR(asChar(alleleFile))));
}

SEXP load_new_effects(SEXP exd, SEXP effectFile) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	load_effects_to_simdata(d, CHAR(asChar(effectFile)));
	return ScalarInteger(0);
}

/*-------------------------------SimData loaders-----------------------------*/



/** Populates a SimData combination with marker allele data.
 *
 * Assumes it is starting from a clean/empty SimData.
 *
 * Since the lines are as columns, this function cannot load files of more than 
 * 1000 lines. If you want to load that, please transpose the file and then
 * call @see load_genes_to_simdata()
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
 * @param d pointer to SimData to be populated
 * @param filename string containing name/path of file containing SNP marker 
 * allele data.
*/
int load_transposed_genes_to_simdata(SimData* d, const char* filename, const char* groupName) {	
	struct TableSize t = get_file_dimensions(filename, '\t');
	FILE* fp;
	const int gp = 1;
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.
	
	// discard the column title that is not a line
	char word[30];

	fscanf(fp, "%s", word);
	
	// now we want to read the header columns. 
	// There are num_columns-1 of these because of the 'name' entry	
	// this will also create our unique ids
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < 1000) {
		current_am = create_empty_allelematrix(t.num_rows - 1, n_to_go);
		d->m = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(t.num_rows - 1, 1000);
		d->m = current_am;
		n_to_go -= 1000;
		while (n_to_go) {
			if (n_to_go < 1000) {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, 1000);
				current_am = current_am->next;
				n_to_go -= 1000;
			}
		}
		current_am = d->m;
	}
	
	// load in the subject names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, "\t%s", word);
		
		if (i_am >= current_am->n_subjects) {
			i_am = 0;
			current_am = current_am->next;
		}
		
		current_am->subject_names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->subject_names[i_am], word);
		
	}
	
	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");
	
	// set the ids for the subjects we loaded
	set_subject_ids(d, 0, t.num_columns - 2);
	
	// get space to put marker names and data we gathered
	d->n_markers = t.num_rows - 1;
	d->markers = get_malloc(sizeof(char*) * (t.num_rows-1));
	//memset(d->markers, '\0', sizeof(char*) * (t.num_rows-1));

	// now read the rest of the table.
	char word2[30];
	int badRows = 0;
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.
		
		// get the row name, store in markers
		fscanf(fp, "%s", word);
		d->markers[j] = get_malloc(sizeof(char) * strlen(word) + 1);
		strncpy(d->markers[j], word, strlen(word) + 1);
		
		current_am = d->m;
		//d->m->alleles[j] = get_malloc(sizeof(char) * d->m[0].n_subjects * 2);
		for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
			// looping through the remaining columns in this row.
			fscanf(fp, "\t%s", word2);
			
			// save the two alleles.
			if (strlen(word2) != 2) {
				++badRows;
				//warning("This file is invalid, but nothing will be done about it.\n");
			}
			
			if (i_am >= current_am->n_subjects) {
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

int load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename, const char* groupName) {
	
	struct TableSize t = get_file_dimensions(filename, '\t');
	FILE* fp;
	const int gp = 1;
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.
	
	// discard the column title that is not a line
	char word[30];

	fscanf(fp, "%s", word);
	
	// now we want to read the header columns. 
	// There are num_columns-1 of these because of the 'name' entry	
	// this will also create our unique ids
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < 1000) {
		current_am = create_empty_allelematrix(t.num_rows - 1, n_to_go);
		d->m = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(t.num_rows - 1, 1000);
		d->m = current_am;
		n_to_go -= 1000;
		while (n_to_go) {
			if (n_to_go < 1000) {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(t.num_rows - 1, 1000);
				current_am = current_am->next;
				n_to_go -= 1000;
			}
		}
		current_am = d->m;
	}
	
	
	// load in the subject names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, "\t%s", word);
		
		if (i_am >= current_am->n_subjects) {
			i_am = 0;
			current_am = current_am->next;
		}
		
		current_am->subject_names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->subject_names[i_am], word);
	}
	
	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");
	
	// set the ids for the subjects we loaded
	set_subject_ids(d, 0, t.num_columns - 2);
	
	// get space to put marker names and data we gathered
	d->n_markers = t.num_rows - 1;
	d->markers = get_malloc(sizeof(char*) * (t.num_rows-1));
	//memset(d->markers, '\0', sizeof(char*) * (t.num_rows-1));

	// now read the rest of the table.
	GetRNGstate();
	char c, decoded[2];
	int r;
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.
		
		// get the row name, store in markers
		fscanf(fp, "%s", word);
		d->markers[j] = get_malloc(sizeof(char) * strlen(word) + 1);
		strncpy(d->markers[j], word, strlen(word) + 1);
		
		current_am = d->m;
		//d->m->alleles[j] = get_malloc(sizeof(char) * d->m[0].n_subjects * 2);
		for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
			// looping through the remaining columns in this row.
			fscanf(fp, "\t%c", &c);
			
			if (i_am >= current_am->n_subjects) {
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

/*Can't add to the number of markers available. */
int load_more_transposed_genes_to_simdata(SimData* d, const char* filename, const char* groupName) {
	struct TableSize t = get_file_dimensions(filename, '\t');
	FILE* fp;
	int gp = get_new_group_num(d);
	if ((fp = fopen(filename, "r")) == NULL) {
		error( "Failed to open file %s.\n", filename);
	}
	// we have successfully opened the file.
	
	// discard the column title that is not a line
	char word[30];

	fscanf(fp, "%s", word);
	
	// now we want to read the header columns. 
	// There are num_columns-1 of these because of the 'name' entry	
	// this will also create our unique ids
	
	// find the end of the AM chain so far
	AlleleMatrix* last_am = d->m; 
	int last_n_subjects = last_am->n_subjects;
	while (last_am->next != NULL) {
		last_am = last_am->next;
		last_n_subjects += last_am->n_subjects;
	}	
	// Create new AMs that will be populated from the file.
	AlleleMatrix* current_am;
	int n_to_go = t.num_columns - 1;
	if (n_to_go < 1000) {
		current_am = create_empty_allelematrix(d->n_markers, n_to_go);
		last_am->next = current_am;
		n_to_go = 0;
	} else {
		current_am = create_empty_allelematrix(d->n_markers, 1000);
		last_am->next = current_am;
		n_to_go -= 1000;
		while (n_to_go) {
			if (n_to_go <= 1000) {
				current_am->next = create_empty_allelematrix(d->n_markers, n_to_go);
				n_to_go = 0;
			} else {
				current_am->next = create_empty_allelematrix(d->n_markers, 1000);
				current_am = current_am->next;
				n_to_go -= 1000;
			}
		}
		current_am = last_am->next;
	}
	
	// set the ids for the subjects we loaded
	set_subject_ids(d, last_n_subjects, last_n_subjects + t.num_columns - 2);
	
	// load in the subject names from the header
	for (int i = 0, i_am = 0; i < (t.num_columns-1); ++i, ++i_am) {
		R_CheckUserInterrupt();
		fscanf(fp, "\t%s", word);
		
		if (i_am >= current_am->n_subjects) {
			i_am = 0;
			current_am = current_am->next;
		}
		
		current_am->subject_names[i_am] = get_malloc(sizeof(char) * strlen(word) + 1);
		strcpy(current_am->subject_names[i_am], word);
	}
	
	// get the rest of the line, to be clean
	fscanf(fp, "%*[^\n]\n");
	
	// get space to put marker names and data we gathered
	// now read the rest of the table.
	char word2[30];
	int markeri;
	current_am = last_am->next;
	for (int j = 0; j < (t.num_rows - 1); ++j) {
		R_CheckUserInterrupt();
		// looping through rows in the table.
		
		// get the row name, store in markers
		fscanf(fp, "%s", word);
		markeri = get_from_unordered_str_list(word, d->markers, d->n_markers);
		
		if (markeri >= 0) {
			for (int i = 0, i_am = 0; i < (t.num_columns - 1); ++i, ++i_am) {
				// looping through the remaining columns in this row.
				fscanf(fp, "\t%s", word2);
				
				// save the two alleles.
				if (strlen(word2) != 2) {
					warning("This file is invalid, but nothing will be done about it.\n");
				}
				
				if (i_am >= current_am->n_subjects) {
					i_am = 0;
					current_am = current_am->next;
				}
				
				//strncpy(d->m->alleles[i] + (2*j), word2, 2);
				current_am->alleles[i_am][2*markeri] = word2[0];
				current_am->alleles[i_am][2*markeri + 1] = word2[1];
				current_am->groups[i_am] = gp;
			}
		} else {
			error( "Could not find the marker %s\n", word);
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
	if (d->m != NULL && d->m->alleles != NULL) {
		//temp = get_malloc(sizeof(char) * ((actual_n_markers * 2)));
		AlleleMatrix* am = d->m;
		
		do {
			R_CheckUserInterrupt();
			/*for (int i = 0; i < actual_n_markers; ++i) {
				location_in_old = sortable[i] - d->map.positions; 
				
				for (int j = 0; j < am->n_subjects; ++j) {
					strncpy(temp, am->alleles[j], am->n_subjects * sizeof(char));
					am->alleles[j][2*i] = temp[2*location_in_old];
					am->alleles[j][2*i + 1] = temp[2*location_in_old + 1];
				}
			}*/
			
			for (int i = 0; i < am->n_subjects; ++i) {
				//strncpy(temp, am->alleles[i], sizeof(char) * ((am->n_markers * 2)));
				temp = get_malloc(sizeof(char) * ((actual_n_markers * 2)));
				
				for (int j = 0; j < actual_n_markers; ++j) {
					location_in_old = sortable[j] - d->map.positions; 
					temp[2*j] = am->alleles[i][2*location_in_old];
					temp[2*j + 1] = am->alleles[i][2*location_in_old + 1];
					
					//am->alleles[i][2*j] = temp[2*location_in_old];
					//am->alleles[i][2*j + 1] = temp[2*location_in_old + 1];
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
	
	MarkerPosition* new_map = get_malloc(sizeof(MarkerPosition) * actual_n_markers);
	for (int i = 0; i < actual_n_markers; ++i) {
		R_CheckUserInterrupt();
		location_in_old = sortable[i] - d->map.positions; 
		new_map[i].chromosome = d->map.positions[location_in_old].chromosome;
		new_map[i].position = d->map.positions[location_in_old].position;
	}
	delete_genmap(&(d->map));
	d->map.positions = new_map;
		
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
		R_CheckUserInterrupt();
		d->map.chr_lengths[i] = d->map.positions[d->map.chr_ends[i+1] - 1].position 
				- d->map.positions[d->map.chr_ends[i]].position;
	}
}

/** Populates a SimData combination with effect values. The SimData must already 
 * have its allele data and map data loaded (so that it has an ordered `markers`
 * list and no markers that will not be used for simulation.
 *
 * It loads in the file as rows of effects for each allele that appears in the
 * allele data out of 'A', 'C', 'G', 'T', encoded in the file as 1, 2, 3, 4
 * respectively.
 *
 * The file should have format:
 *
 * [n]
 *
 * [marker] [allele] [effect]
 *
 * [marker] [allele] [effect]
 *
 * ...
 *
 * The value of [n] in the file is ignored.
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
		error( "Failed to open file %s.\n", filename);
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
			int symbol_index; 
			char* symbol_location = strchr(alleles_loaded, allele);
			if (symbol_location == NULL) {
				symbol_index = n_allele;
				++n_allele;
				alleles_loaded[symbol_index] = allele;
			} else {
				symbol_index = symbol_location - alleles_loaded; // difference between the pointers
			}
			
			// the marker is in our list and the allele value is valid
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
	int gp = load_transposed_genes_to_simdata(d, data_file, NULL);
	
	load_genmap_to_simdata(d, map_file);
	load_effects_to_simdata(d, effect_file);
	
	get_sorted_markers(d, d->n_markers);
	get_chromosome_locations(d);
	return gp;
}




/*--------------------------Recombination counts-----------------------------*/

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

// window should be odd
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


SEXP find_crossovers(SEXP exd, SEXP parentFile, SEXP outFile, SEXP windowSize, SEXP certainty) {
	SimData* d = (SimData*) R_ExternalPtrAddr(exd);
	
	const char* load_fname = CHAR(asChar(parentFile));
	const char* save_fname = CHAR(asChar(outFile));
	
	int win = asInteger(windowSize);
	if (win < 0 || win == NA_INTEGER) { error("`window.size` parameter is invalid.\n"); }
	
	int cert = asLogical(certainty);
	if (cert == NA_LOGICAL) { error("`certainty` parameter is invalid.\n"); }
	
	return ScalarInteger(calculate_recombinations_from_file(d, load_fname, save_fname, win, cert));
}

SEXP send_map(SEXP exd) {
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
