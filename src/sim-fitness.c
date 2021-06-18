#include "sim-fitness.h"

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
		error("Failed to open file %s.\n", block_file);
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
 * calculate the local GEBV for the first allele at each marker in the block, and 
 * the local GEBV for the second allele at each marker in the block, then save
 * the result to a file. This gives block effects for each haplotype of each 
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
 * block effects/local GEBVs will be saved.
 * @param group group number from which to split the top individuals.
 */
void calculate_group_block_effects(SimData* d, MarkerBlocks b, const char* output_file, int group) {
	
	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error("Failed to open file %s.\n", output_file);
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
			
			// calculate the local GEBV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {	
				for (int q = 0; q < d->e.effects.rows; ++q) {			
					if (ggenos[i][2 * b.markers_in_block[j][k]] == d->e.effect_names[q]) {
						beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
					}
				}	
			}
			
			// print the local GEBV
			fprintf(outfile, " %lf", beffect);
			fflush(outfile);
		}
		
		sprintf(buffer, "\n%s_2", gnames[i]);
		fwrite(buffer, sizeof(char), strlen(buffer), outfile);
		
		// for each block for the second haplotype
		for (int j = 0; j < b.num_blocks; ++j) {
			beffect = 0;
			// calculate the local GEBV
			for (int k = 0; k < b.num_markers_in_block[j]; ++k) {	
				for (int q = 0; q < d->e.effects.rows; ++q) {			
					if (ggenos[i][2 * b.markers_in_block[j][k] + 1] == d->e.effect_names[q]) {
						beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
					}		
				}	
			}
			
			// print the local GEBV
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
 * calculate the local GEBV for the first allele at each marker in the block, and 
 * the local GEBV for the second allele at each marker in the block, then save
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
 * block effects/local GEBVs will be saved.
 */
void calculate_all_block_effects(SimData* d, MarkerBlocks b, const char* output_file) {
	FILE* outfile;
	if ((outfile = fopen(output_file, "w")) == NULL) {
		error( "Failed to open file %s.\n", output_file);
	}
	
	int bufferlen = 100;
	char buffer[bufferlen];
	
	int gsize = 0;
	AlleleMatrix* m = d->m;
	do {
		gsize += m->n_subjects;
	} while ((m = m->next) != NULL);
	
	double beffect;
	
	// for each group member
	m = d->m;
	int total_i = 0;
	do {
		for (int i = 0; i < m->n_subjects; ++i, ++total_i) {
			// for each block
			sprintf(buffer, "%s_1", m->subject_names[i]);
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);
			
			// for each block
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;
				
				// calculate the local GEBV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {	
					for (int q = 0; q < d->e.effects.rows; ++q) {			
						if (m->alleles[i][2 * b.markers_in_block[j][k]] == d->e.effect_names[q]) {
							beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
						}
					}	
				}
				
				// print the local GEBV
				fprintf(outfile, " %lf", beffect);
				fflush(outfile);
			}
			
			sprintf(buffer, "\n%s_2", m->subject_names[i]);
			fwrite(buffer, sizeof(char), strlen(buffer), outfile);
			
			// for each block for the second haplotype
			for (int j = 0; j < b.num_blocks; ++j) {
				beffect = 0;
				// calculate the local GEBV
				for (int k = 0; k < b.num_markers_in_block[j]; ++k) {	
					for (int q = 0; q < d->e.effects.rows; ++q) {			
						if (m->alleles[i][2 * b.markers_in_block[j][k] + 1] == d->e.effect_names[q]) {
							beffect += d->e.effects.matrix[q][b.markers_in_block[j][k]];
						}		
					}	
				}
				
				// print the local GEBV
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
char* calculate_ideal_genotype(SimData* d) {
	if (d->e.effects.matrix == NULL || d->e.effects.rows < 1 || d->e.effect_names == NULL) {
		fprintf(stderr, "No effect values are loaded\n");
		return NULL;
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

/** Takes a look at the currently-loaded effect values returns the highest possible
 * GEBV any genotype could score using those effect values.
 *
 * The SimData must be initialised with marker effects for this function to succeed.
 * 
 * @param d pointer to the SimData containing markers and marker effects.
 * @returns the GEBV of the best/ideal genotype.
 */
double calculate_optimal_gebv(SimData* d) {
	char* best_alleles = calculate_ideal_genotype(d);
	double best_gebv = 0;
	
	DecimalMatrix counts = generate_zero_dmatrix(d->n_markers, 1);
	DecimalMatrix effect_row, product;
	
	for (int i = 0; i < d->e.effects.rows; ++i) {
		// fill the count matrix for this allele.
		for (int j = 0; j < d->n_markers; ++j) {
			if (best_alleles[j] == d->e.effect_names[i]) {
				counts.matrix[j][0] = 2;
			} else {
				counts.matrix[j][0] = 0;
			}
		}
		
		// calculate the GEBV contribution from this allele
		effect_row = subset_dmatrix_row(&(d->e.effects), i);
		product = multiply_dmatrices(&effect_row, &counts);
		
		best_gebv += product.matrix[0][0];
		delete_dmatrix(&effect_row);
		delete_dmatrix(&product);
	}
	
	delete_dmatrix(&counts);
	free(best_alleles);
	
	return best_gebv;
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