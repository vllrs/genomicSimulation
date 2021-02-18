#include "sim-printers.h"

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
			save_full_pedigree(f, d->m, d->m);
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
		save_all_fitness(f, d);
	} else if (asInteger(group) >= 0) {
		save_group_fitness(f, d, asInteger(group));
	} else {
		fclose(f);
		error("Supplied group number is invalid.");
	}
	
	fclose(f);
	return ScalarInteger(0);
}

/*--------------------------------Printing-----------------------------------*/


/** Prints the data inside a SimData to a file. Column separators are '\t'. Note:
 * does not print gene data saved in AlleleMatrix of the SimData
 *
 * The printing format is:
 *
 * name	chr	pos [effect row names]	
 *
 * [marker name]	[chr number]	[chr pos]	[effects]
 *
 * [marker name]	[chr number]	[chr pos]	[effects]
 *
 * ...
 *
 * If m->effects is NULL, m->ref_alleles is NULL, or m->genetic_map is NULL, 
 * then the relevant columns are omitted.
 *
 * @see save_allele_matrix() for details on allele table formats
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
		fwrite(m->markers[j], strlen(m->markers[j]), 1, f);
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

/** Prints all the gene data saved in the linked list starting with `m` to the 
 * file. Uses the following format:
 *
 * 		[marker name]	[marker name]
 *
 * [subject id]:[subject name]	[allele pairs for each marker]
 *
 * [subject id]:[subject name] 	[allele pairs for each marker]
 *
 * ...
 * 
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the AlleleMatrix whose data we print
 * @param markers array of strings that correspond to names of the markers.
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
		for (int i = 0; i < m->n_subjects; ++i) {
			// print the name or ID of the individual.
			if (m->subject_names[i] != NULL) {
				fwrite(m->subject_names[i], sizeof(char), strlen(m->subject_names[i]), f);
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
 * file. Has markers as columns and individual genomes as rows, in the following format:
 *
 * 		[marker name]	[marker name]
 *
 * [subject id]:[subject name]	[allele pairs for each marker]
 *
 * [subject id]:[subject name] 	[allele pairs for each marker]
 *
 * ...
 * 
 * @param f file pointer opened for writing to put the output
 * @param m pointer to the AlleleMatrix whose data we print
 * @param markers array of strings that correspond to names of the markers.
*/
void save_transposed_allele_matrix(FILE* f, AlleleMatrix* m, char** markers) {
	// Count number of genotypes in the AM
	AlleleMatrix* currentm = m; // current matrix
	int tn_subjects = 0;
	do {
		tn_subjects += currentm->n_subjects;
	} while ((currentm = currentm->next) != NULL);
		
	currentm = m;
	
	/* Print header */
	for (int i = 0, currenti = 0; i < tn_subjects; ++i, ++currenti) {
		if (currenti >= currentm->n_subjects) {
			currenti = 0;
			currentm = currentm->next;
		}
		if (currentm->subject_names[currenti] != NULL) { // assume all-or-nothing with marker names
			//fprintf(f, "\t%s", markers[i]);
			fwrite("\t", sizeof(char), 1, f);
			fwrite(currentm->subject_names[currenti], sizeof(char), strlen(currentm->subject_names[currenti]), f);
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
		
		for (int i = 0, currenti = 0; i < tn_subjects; ++i, ++currenti) {
			if (currenti >= currentm->n_subjects) {
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
			} else if (pedigree[0] > 0) {
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

void save_one_step_pedigree(FILE* f, SimData* d) {
	unsigned int pedigree[2];
	char* name;
	AlleleMatrix* m = d->m;
	
	do {
		for (int i = 0; i < m->n_subjects; ++i) {
			/*Group member name*/
			if (m->subject_names[i] != NULL) {
				fwrite(m->subject_names[i], sizeof(char), strlen(m->subject_names[i]), f);
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
				} else if (pedigree[0] > 0) {
					fprintf(f, "%d", pedigree[1]);
					//fwrite(pedigree + 1, sizeof(int), 1, f);
				}
				
			}
			fwrite("\n", sizeof(char), 1, f);
		}
	} while ((m = m->next) != NULL);
	fflush(f);
}


void save_group_full_pedigree(FILE* f, SimData* d, int group) {
	int group_size = get_group_size( d, group);
	unsigned int* group_contents = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	const char newline[] = "\n";
	
	for (int i = 0; i < group_size; i++) {
		/*Group member name*/
		fprintf(f, "%d\t", group_contents[i]);
		if (group_names[i] != NULL) {
			fwrite(group_names[i], sizeof(char), strlen(group_names[i]), f);
		}
		
		save_parents_of(f, d->m, group_contents[i]);
		fwrite(newline, sizeof(char), 1, f);
	}
	free(group_names);
	free(group_contents);
	fflush(f);	
}

/*Saves pedigree of all in m.*/
void save_full_pedigree(FILE* f, AlleleMatrix* m, AlleleMatrix* parents) {
	const char newline[] = "\n";
	const char tab[] = "\t";
	
	for (int i = 0; i < m->n_subjects; i++) {
		/*Group member name*/
		//fwrite(m->ids + i, sizeof(int), 1, f);
		fprintf(f, "%d", m->ids[i]);
		fwrite(tab, sizeof(char), 1, f);
		if (m->subject_names[i] != NULL) {
			fwrite(m->subject_names[i], sizeof(char), strlen(m->subject_names[i]), f);
		}
		
		// save our parents
		// open brackets
		fwrite("=(", sizeof(char), 2, f);
		char* name;
		if (m->pedigrees[0][i] == m->pedigrees[1][i]) {
			// Selfed parent
			name = get_name_of_id( parents, m->pedigrees[0][i]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (m->pedigrees[0][i] > 0) {
				fprintf(f, "%d", m->pedigrees[0][i]);
				//fwrite(m->pedigrees[0] + i, sizeof(int), 1, f);
			}	
			save_parents_of(f, parents, m->pedigrees[0][i]);
		} else {
			// Parent 1
			name = get_name_of_id( parents, m->pedigrees[0][i]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (m->pedigrees[0][i] > 0) {
				fprintf(f, "%d", m->pedigrees[0][i]);
				//fwrite(m->pedigrees[0] + i, sizeof(int), 1, f);
			}
			save_parents_of(f, parents, m->pedigrees[0][i]);
			// separator
			fwrite(",", sizeof(char), 1, f);
			// Parent 2
			name = get_name_of_id( parents, m->pedigrees[1][i]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (m->pedigrees[1][i] > 0) {
				fprintf(f, "%d", m->pedigrees[1][i]);
				//fwrite(m->pedigrees[0] + i + 1, sizeof(int), 1, f);
			}
			save_parents_of(f, parents, m->pedigrees[1][i]);
		}
		// close brackets
		fwrite(")", sizeof(char), 1, f);
		
		fwrite(newline, sizeof(char), 1, f);
	}
	fflush(f);	
}

/*Recursive, saves all parents and their parents.*/
void save_parents_of(FILE* f, AlleleMatrix* m, unsigned int id) {
	unsigned int pedigree[2];
	
	if (get_parents_of_id(m, id, pedigree) == 0) {
		// open brackets
		fwrite("=(", sizeof(char), 2, f);
		char* name;
		
		if (pedigree[0] == pedigree[1]) {
			// Selfed parent
			name = get_name_of_id( m, pedigree[0]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (pedigree[0] > 0) {
				fprintf(f, "%d", pedigree[0]);
				//fwrite(pedigree, sizeof(int), 1, f);
			}
			save_parents_of(f, m, pedigree[0]);
					
		} else {
			// Parent 1
			name = get_name_of_id( m, pedigree[0]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (pedigree[0] > 0) {
				fprintf(f, "%d", pedigree[0]);
				//fwrite(pedigree, sizeof(int), 1, f);
			}
			save_parents_of(f, m, pedigree[0]);
			
			// separator
			fwrite(",", sizeof(char), 1, f);
			
			// Parent 2
			name = get_name_of_id( m, pedigree[1]);
			if (name != NULL) {
				fwrite(name, sizeof(char), strlen(name), f);
			} else if (pedigree[1] > 0) {
				fprintf(f, "%d", pedigree[1]);
				//fwrite(pedigree + 1, sizeof(int), 1, f);
			}
			save_parents_of(f, m, pedigree[1]);
			
		}

		// close brackets
		fwrite(")", sizeof(char), 1, f);
	}
}


void save_group_fitness(FILE* f, SimData* d, int group) {
	int group_size = get_group_size( d, group);
	unsigned int* group_contents = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	DecimalMatrix effects = calculate_fitness_metric_of_group(d, group);
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

void save_all_fitness(FILE* f, SimData* d) {
	AlleleMatrix* am = d->m;
	const char newline[] = "\n";
	const char tab[] = "\t";
	DecimalMatrix effects;
	
	do {
		effects = calculate_fitness_metric(am, &(d->e));
		for (int i = 0; i < effects.cols; ++i) {
			/*Group member name*/
			//fwrite(group_contents + i, sizeof(int), 1, f);
			fprintf(f, "%d", am->ids[i]);
			fwrite(tab, sizeof(char), 1, f);
			if (am->subject_names[i] != NULL) {
				fwrite(am->subject_names[i], sizeof(char), strlen(am->subject_names[i]), f);
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

void save_fitness(FILE* f, DecimalMatrix* e, unsigned int* ids, char** names) {
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

void save_count_matrix(FILE* f, SimData* d, char allele) {
	DecimalMatrix counts = calculate_full_count_matrix_of_allele(d->m, allele);
	
	AlleleMatrix* currentm = d->m;
	// print the header
	for (int i = 0, currenti = 0; i < counts.cols; ++i, ++currenti) {
		if (currenti >= currentm->n_subjects) {
			currenti = 0;
			currentm = currentm->next;
		}
		fwrite("\t", sizeof(char), 1, f);
		if (currentm->subject_names[currenti] != NULL) {
			fwrite(currentm->subject_names[currenti], sizeof(char), strlen(currentm->subject_names[currenti]), f);
		}
	}
	
	fwrite("\n", sizeof(char), 1, f);
	
	// Print the body
	for (int i = 0; i < d->n_markers; ++i) { // loop through markers
		if (d->markers != NULL && d->markers[i] != NULL) {
			fwrite(d->markers[i], sizeof(char), strlen(d->markers[i]), f);
		}
		fwrite("\t", sizeof(char), 1, f);
	
		for (int j = 0; j < counts.cols; ++j) { // loop through subjects
			// print the matrix entries
			fprintf(f, "%f ", counts.matrix[i][j]);
		}
		//print the newline
		if (i + 1 < counts.cols) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}
	
	delete_dmatrix(&counts);
	fflush(f);
}

void save_count_matrix_of_group(FILE* f, SimData* d, char allele, int group) {
	unsigned int group_size = get_group_size( d, group);
	unsigned int* group_ids = get_group_ids( d, group, group_size);
	char** group_names = get_group_names( d, group, group_size);
	DecimalMatrix counts = calculate_count_matrix_of_allele_for_ids(d->m, group_ids, group_size, allele);
	
	fprintf(f, "%d", group);
	// print the header
	for (int i = 0; i < counts.cols; ++i) {
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
	
		for (int j = 0; j < counts.cols; ++j) { // loop through subjects
			// print the matrix entries
			fprintf(f, "%f ", counts.matrix[i][j]);
		}
		//print the newline
		if (i + 1 < counts.cols) {
			fwrite("\n", sizeof(char), 1, f);
		}
	}
	
	delete_dmatrix(&counts);
	free(group_ids);
	free(group_names);
	fflush(f);
}
