#ifndef SIM_OPERATIONS_H
#define SIM_OPERATIONS_H

#include <string.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

//#define PI 3.1415926535897932384626433832795028841971693993751
#define TRUE 1
#define FALSE 0


/** \page Guides Guides to simulating particular situations
 *
 * To be added
 *
 */
 
 /** \page Methods Simulation Methodology
 *
 * For the moment, please refer to the R package vignette for 
 * descriptions of simulation methodologies and assumptions. 
 *
 * [genomicSimulation (R package)](https://kiravill.github.io/genomicSimulation/)
 *
 */
 
 /* This section contains settings for the simulation that users can modify if they have the need.
 * To apply the modified settings the simulation tool must be re-compiled.
 * To modify a setting, replace the number (eg 1000) with the new value without modifying the name
 */
 
/** The largest contiguous block of memory that could be requested in the
 * process of simulation is CONTIG_WIDTH integers. This setting's default value 
 * is 1000.
 *
 * This could be decreased to help long simulations or lower-end machines. 
 * Increasing this may provide some speed gain. 
 */
#define CONTIG_WIDTH 1000
 
 /** The maximum number of characters allowed in a name field.
 * These include names of SNPs, names of genotypes loaded from files, 
 * names of generated genotypes, and save-as-you-go filenames. Default is 30.
 *
 * Increase this if there is a risk some names may be longer than 
 * this value.
 */
#define NAME_LENGTH 30
 

/** @defgroup structs Data Structures
 *
 * How the simulation stores data.
 *
 * genomicSimulation is a state-based package. These data structures
 * are used to store the library's data/state. Many use dynamically 
 * allocated memory so call the relevant delete_ function if one exists
 * when you are finished with them.
 *
 * SimData is the central state/data
 * storage struct, and so is a required parameter to most user-facing functions. 
 * It contains pointers to a GeneticMap (storing the loaded genome map), an
 * EffectMatrix (storing the loaded allele effects), and an AlleleMatrix
 * (storing metadata and genotypes of founders and simulated offspring).
 *
 * Other structs in this group (TableSize, MarkerBlocks, GenOptions) represented
 * specific types of data and are used as parameters and return value of certain
 * functions.
 *
 * @{
 */

/** A struct representing a single marker location. the attribute
 * `chromosome` represents the chromosome number, and `position` the
 * position on the chromosome in centiMorgans.
 */
typedef struct {
	int chromosome; /**< The chromosome number */
	float position; /**< The distance in centiMorgans along the chromosome */
} MarkerPosition;

/** A simple struct used for returning the dimensions of a matrix or table.*/
struct TableSize {
	int num_columns;
	int num_rows;
};

/** A struct used to store a set of blocks of markers.
 */
typedef struct {
    /** The number of blocks whose details are stored here. */
	int num_blocks; 
    
    /** Pointer to a heap array of length num_blocks
      * containing the number of markers that make up each block */
	int* num_markers_in_block; 
    
    /** Pointer to a heap array of length num_blocks, each
     * entry in which is a pointer to a heap array with length corresponding to
     * the value of the corresponding entry in num_markers_in_block whose values
     * are the indexes in the SimData of the markers that make up that block. */
	int** markers_in_block; 
} MarkerBlocks;


/** A row-major heap matrix that contains floating point numbers. `dmatrix` functions
 * are designed to act on this matrix.
 *
 * Rows make the first index of the matrix and columns the second.
 */
typedef struct {
	double** matrix; /**< The actual matrix and contents */
	int rows;        /**< Number of rows in the matrix */
	int cols;        /**< number of columns in the matrix */
} DecimalMatrix;


/** A type that contains choices of settings for SimData functions that create a
 * new AlleleMatrix/generation.
 *
 * The family_size parameter will affect how many offspring are produced.
 *
 * The will_name_offspring, will_track_pedigree, and will_allocate_ids parameters
 * affect how much extra detail about the offspring is generated/saved.
 *
 * The will_save_to_simdata toggle allows you the option of generating offspring without
 * saving them in memory. This may be useful in combination with save-as-you-go toggles
 * will_save_pedigree_to_file, will_save_bvs_to_file, and will_save_alleles_to_file,
 * to generate a larger number of offspring than will fit in memory.
*/
typedef struct {
	int will_name_offspring; /**< A boolean: whether generated offspring should be given names. */
	char* offspring_name_prefix; /**< If `will_name_offspring` is true, generated 
                           * offspring are named [offspring_name_prefix][index]. */

	int family_size; /**< The number of offspring to produce from each cross.*/

	int will_track_pedigree; /**< A boolean: whether to track parentage of generated offspring.*/
	int will_allocate_ids; /**< A boolean: whether to allocate generated offspring session-
                            * unique IDs. IDs are used for pedigree tracking. The 
                            * offspring of an anonymous individual (one without an ID)
                            * cannot identify that individual as their parent. */

	char* filename_prefix; /**< A string used in save-as-you-go file names. */
	int will_save_pedigree_to_file; /**< A boolean. If true, the full/recursive 
                            * pedigrees of every offspring generated in the cross
                            * are saved to "[filename_prefix}-pedigree.txt", even 
                            * if `will_save_to_simdata` is false. 
                            * Pedigrees are saved in the format of save_full_pedigree()*/
	int will_save_bvs_to_file; /**< A boolean. If true, the breeding values
                            * of every offspring generated in the cross
                            * are saved to "[filename_prefix}-bv.txt", even 
                            * if `will_save_to_simdata` is false. 
                            * BVs are saved in the format of save_bvs() */
	int will_save_alleles_to_file; /**< A boolean. If true, the set of alleles
                            * of every offspring generated in the cross
                            * are saved to "[filename_prefix}-genotype.txt", even 
                            * if `will_save_to_simdata` is false. 
                            * Genotypes are saved in the format of save_group_alleles()*/
	int will_save_to_simdata; /**< A boolean. If true, the generated offspring exist 
                            * in the SimData struct after the function executes. 
                            * If false, they are discarded after creation. */
} GenOptions;


/** A type that stores the genetic map for a set of markers.
 *
 * To get all markers belonging to a particular chromosome, use the following rule:
 * Chr n includes all markers in `positions` starting at index chr_ends[n-1] up
 * but not including the marker at index chr_ends[n]
 *
 * Chromosomes must be numbered. All chromosomes from 1 up to and including the 
 * highest chromosome number found in the loaded map are represented in these
 * arrays.
 *
*/
typedef struct {
	int n_chr; /**< The number of chromosomes represented in the map. This 
                * corresponds to the highest numbered chromosome with a tracked
                * marker (some chromosomes in between may be empty) */
	int* chr_ends; /**< An array of ints. The entry at index i is the index in
                * `positions` of the first marker that belongs to Chr(i + 1). 
                * The array is n_chr + 1 integers long.*/
	float* chr_lengths; /**< An array of floats. The entry at index i is the length
                * of Chr(i + 1), calculated by 
                * `position of last marker - position of first marker`.
                * The array is n_chr entries long. */

	MarkerPosition* positions; /**< An array of MarkerPositions, ordered from lowest to highest. */
} GeneticMap;

/** A linked list entry that stores a matrix of alleles for a set of SNP markers
 * and genotypes.
*/
typedef struct AlleleMatrix AlleleMatrix;
struct AlleleMatrix {
    /** A matrix of SNP markers by lines/genotypes containing pairs of alleles
     * eg TT, TA. Use `alleles[line index][marker index * 2]` to get the 
     * first allele and `alleles[lines index][marker index * 2 + 1]` to 
     * get the second. If CONTIG_WIDTH lines are saved here, another 
     * AlleleMatrix is added to the linked list when there's a need to save more.*/
	char* alleles[CONTIG_WIDTH];
    
	int n_genotypes; /**< Number of genotypes currently loaded in this matrix.*/
	int n_markers; /**< Number of markers across which genotypes are tracked. This has
                    * redundancy with SimData and other members of its linked list
                    * but it's good to know how big your own `alleles` array is.*/

    char* names[CONTIG_WIDTH]; /**< Array of dynamically allocated strings 
                    * containing the names of the lines/genotypes in this matrix. 
                    * Guaranteed to be NULL if they do not have names. */
	unsigned int ids[CONTIG_WIDTH]; /**< Unique ID for each genotype. */
	unsigned int pedigrees[2][CONTIG_WIDTH]; /**< Two lists of integer IDs of the 
                    * parents of this genotype (if tracked), or 0 if we don't know/care.*/
	unsigned int groups[CONTIG_WIDTH]; /**< Group allocation of each genotype. */
    
	AlleleMatrix* next; /**< Pointer to the next AlleleMatrix in the linked list, 
                         * or NULL if this entry is the last. */
};

/** A type that stores a matrix of effect values and their names.
 */
typedef struct {
	DecimalMatrix effects; /**< Effect on breeding value of alleles at markers. 
        * Rows correspond to `effect_names`/alleles, columns to markers. */
	char* effect_names; /**< Character array containing allele characters ordered 
        * to match rows of `effects`. */
} EffectMatrix;

/** Composite type that is used to run crossing simulations.
 *
 * The core of this type is a list of markers. These are used to index the rows
 * of the allele matrix and the position map, and the columns of the effect matrix.
 */
typedef struct {
	int n_markers;  /**< The number of markers/length of `markers`. */
	char** markers; /**< Array of strings containing the names of markers. */

	GeneticMap map; /**< A GeneticMap. If this is set, then `markers`
                     * will be ordered and all markers have a known position.*/
	AlleleMatrix* m; /**< Pointer to an AlleleMatrix, which stores data and
                      * metadata of founders and simulated offspring. The 
                      * AlleleMatrix is start of a linked list if there are
                      * many genotypes. */
	EffectMatrix e; /**< An EffectMatrix, optional for the use of the simulation.
                     * Used for calculating breeding values from which alleles
                     * a genotype has at each marker.*/

	unsigned int current_id; /**< Highest SimData-unique ID that has been generated
                              * so far. Used to track which IDs have already been 
                              * given out.*/
} SimData;

//const SimData EMPTY_SIMDATA;
const GenOptions BASIC_OPT;
/** @} */

/** @defgroup maths Mathematical functions
 *
 * For mathematical and statistical operations as required by the package. 
 *
 * Includes matrix operations defined on a DecimalMatrix struct, and 
 * draws from certain random distributions.
 *
 * @{
 */
DecimalMatrix generate_zero_dmatrix(int r, int c);
DecimalMatrix subset_dmatrix_row(DecimalMatrix* m, int row_index);
DecimalMatrix add_dmatrices(DecimalMatrix* a, DecimalMatrix* b);
void add_to_dmatrix(DecimalMatrix* a, DecimalMatrix* b);
DecimalMatrix multiply_dmatrices(DecimalMatrix* a, DecimalMatrix* b);
/** @} */

/** @defgroup supporters Utils/Supporting Functions
 *
 * @{
 */
struct TableSize get_file_dimensions(const char* filename, char sep);
int get_from_ordered_uint_list(unsigned int target, unsigned int* list, unsigned int list_len);
int get_from_unordered_str_list(char* target, char** list, int list_len) ;

void get_sorted_markers(SimData* d, int actual_n_markers);
void get_chromosome_locations(SimData *d);

void set_names(AlleleMatrix* a, char* prefix, int suffix, int from_index);
void set_ids(SimData* d, int from_index, int to_index);
int get_integer_digits(int i);
int get_new_group_num( SimData* d);
void set_group_list( SimData* d, int by_n, int new_group);
void condense_allele_matrix( SimData* d);
void* get_malloc(size_t size);

int _simdata_pos_compare(const void *pp0, const void *pp1);
int _descending_double_comparer(const void* pp0, const void* pp1);
int _ascending_double_comparer(const void* pp0, const void* pp1);
int _ascending_float_comparer(const void* p0, const void* p1);
int _ascending_int_comparer(const void* p0, const void* p1);
int _ascending_int_dcomparer(const void* pp0, const void* pp1);
/**@}*/

/** @defgroup getters Data Access and Search Functions
 *
 * For non-persistent access to simulation results and contents.
 *
 * @{
 */
char* get_name_of_id( AlleleMatrix* start, unsigned int id);
char* get_genes_of_id ( AlleleMatrix* start, unsigned int id);
int get_parents_of_id( AlleleMatrix* start, unsigned int id, unsigned int output[2]);
void get_ids_of_names( AlleleMatrix* start, int n_names, char* names[n_names], unsigned int* output);
unsigned int get_id_of_child( AlleleMatrix* start, unsigned int parent1id, unsigned int parent2id);
int get_index_of_child( AlleleMatrix* start, unsigned int parent1id, unsigned int parent2id);
int get_index_of_name( AlleleMatrix* start, char* name);
unsigned int get_id_of_index( AlleleMatrix* start, int index);
char* get_genes_of_index( AlleleMatrix* start, int index);

int get_group_size( SimData* d, int group_id);
char** get_group_genes( SimData* d, int group_id, int group_size);
char** get_group_names( SimData* d, int group_id, int group_size);
unsigned int* get_group_ids( SimData* d, int group_id, int group_size);
unsigned int* get_group_indexes(SimData* d, int group_id, int group_size);
double* get_group_bvs( SimData* d, int group_id, int group_size);
unsigned int* get_group_parent_ids( SimData* d, int group_id, int group_size, int parent);
char** get_group_parent_names( SimData* d, int group_id, int group_size, int parent);
char** get_group_pedigrees( SimData* d, int group_id, int group_size);

int* get_existing_groups( SimData* d, int* n_groups);
int** get_existing_group_counts( SimData* d, int* n_groups);
/**@}*/

/** @defgroup groupmod Seletion/Group Modification Functions
 *
 * For simulation of selection or structure in breeding programs.
 *
 * @{
 */
int combine_groups( SimData* d, int list_len, int group_ids[list_len]);
void split_into_individuals( SimData* d, int group_id);
void split_into_families(SimData* d, int group_id);
int split_from_group( SimData* d, int n, int indexes_to_split[n]);
/**@}*/

/** @defgroup deletors Deletor Functions
 *
 * For deleting and free associated memory of data structures.
 *
 * @ingroup structs
 * @{
 */
void delete_group(SimData* d, int group_id);
void delete_genmap(GeneticMap* m);
void delete_allele_matrix(AlleleMatrix* m);
void delete_effect_matrix(EffectMatrix* m);
void delete_simdata(SimData* m);
void delete_markerblocks(MarkerBlocks* b);
void delete_dmatrix(DecimalMatrix* m);
/**@}*/

/** @defgroup loaders Setup Functions
 *
 * For setup of the simulation (loading founders, genetic maps, and optionally allele effects).
 *
 * @{
 */
AlleleMatrix* create_empty_allelematrix(int n_markers, int n_genotypes);
SimData* create_empty_simdata();
 
int load_transposed_genes_to_simdata(SimData* d, const char* filename);
int load_more_transposed_genes_to_simdata(SimData* d, const char* filename);
int load_genes_to_simdata(SimData* d, const char* filename); //@ add
int load_more_genes_to_simdata(SimData* d, const char* filename); //@ add
int load_transposed_encoded_genes_to_simdata(SimData* d, const char* filename);
void load_genmap_to_simdata(SimData* d, const char* filename);
void load_effects_to_simdata(SimData* d, const char* filename);
int load_all_simdata(SimData* d, const char* data_file, const char* map_file, const char* effect_file);
/** @} */

/** @defgroup recomb Recombination Calculators 
 *
 * Experimental functions for retroactively calculating number of recombinations. 
 *
 * This functionality is for interest only. It is not clear, or tidy,
 * or checked against real data.
 *
 * @{
 */
int* calculate_min_recombinations_fw1(SimData* d, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int certain); // forward filling, window size 1
int* calculate_min_recombinations_fwn(SimData* d, char* parent1, unsigned int p1num, char* parent2,
		unsigned int p2num, char* offspring, int window_size, int certain); // forward filling, window size n


// int has_same_alleles(char* p1, char* p2, int i);
// int has_same_alleles_window(char* g1, char* g2, int start, int w);
/** Simple operator to determine if at marker i, two genotypes share at least
 * one allele. Checks only 3 of four possible permutations because assumes
 * there cannot be more than two alleles at a given marker.
 *
 * @param p1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param p2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param i index of the marker at which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles(char* p1, char* p2, int i) {
	return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}
// w is window length, i is start value
/** Simple operator to determine if at markers with indexes i to i+w inclusive, two genotypes
 * share at least one allele. Checks only 3 of four possible permutations at each marker
 * because assumes there cannot be more than two alleles at a given marker. For the return value
 * to be true, there must be at least one match at every one of the markers in the window.
 *
 * @param g1 pointer to a character array genotype of the type stored in an AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param g2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param start index of the first marker in the window over which to perform the check
 * @param w length of the window over which to perform the check
 * @returns boolean result of the check
 */
static inline int has_same_alleles_window(char* g1, char* g2, int start, int w) {
	int same = TRUE;
	int i;
	for (int j = 0; j < w; ++j) {
		i = start + j;
		same = same && (g1[i<<1] == g2[i<<1] || g1[(i<<1) + 1] == g2[i] || g1[i] == g2[(i<<1) + 1]);
	}
	return same;
}

int calculate_recombinations_from_file(SimData* d, const char* input_file, const char* output_file,
		int window_len, int certain);
/**@}*/

/** @defgroup crossers Crossing and Progression Simulation Functions
 *
 * For simulation of progression steps in breeding programs.
 * 
 * @{
 */
 
    /** @defgroup meiosis Meiosis Simulation Functions
     * 
     * For simulation of meiosis.
     *
     * @{
     */
void generate_gamete(SimData* d, char* parent_genome, char* output);
void generate_cross(SimData* d, char* parent1_genome, char* parent2_genome, char* output);
void generate_doubled_haploid(SimData* d, char* parent_genome, char* output);
    /**@}*/

int cross_random_individuals(SimData* d, int from_group, int n_crosses, GenOptions g);
int cross_these_combinations(SimData* d, int n_combinations, int combinations[2][n_combinations],  GenOptions g);
int self_n_times(SimData* d, int n, int group, GenOptions g);
int make_doubled_haploids(SimData* d, int group, GenOptions g);

int make_all_unidirectional_crosses(SimData* d, int from_group, GenOptions g);
int make_n_crosses_from_top_m_percent(SimData* d, int n, int m, int group, GenOptions g);
int make_crosses_from_file(SimData* d, const char* input_file, GenOptions g);
int make_double_crosses_from_file(SimData* d, const char* input_file, GenOptions g);
/**@}*/

/** @defgroup calculators Breeding Value and Allele Count Calculators 
 *
 * For calculations related to the loaded allele effects and internal additive breeding value model.
 *
 * @{
 */
int split_by_bv(SimData* d, int group, int top_n, int lowIsBest);
DecimalMatrix calculate_group_bvs(SimData* d, int group);
DecimalMatrix calculate_bvs( AlleleMatrix* m, EffectMatrix* e);
DecimalMatrix calculate_count_matrix_of_allele_for_ids( AlleleMatrix* m, unsigned int* for_ids, unsigned int n_ids, char allele);
DecimalMatrix calculate_full_count_matrix_of_allele( AlleleMatrix* m, char allele);

MarkerBlocks create_n_blocks_by_chr(SimData* d, int n);
MarkerBlocks read_block_file(SimData* d, const char* block_file);
void calculate_group_local_bvs(SimData* d, MarkerBlocks b, const char* output_file, int group);
void calculate_local_bvs(SimData* d, MarkerBlocks b, const char* output_file);

char* calculate_optimal_alleles(SimData* d);
double calculate_optimum_bv(SimData* d);
double calculate_minimum_bv(SimData* d);
/**@}*/

/** @defgroup savers Saving Functions
 * 
 * For saving persistent simulation results.
 * 
 * @{
 */
void save_simdata(FILE* f, SimData* m);

void save_marker_blocks(FILE* f, SimData* d, MarkerBlocks b);

void save_allele_matrix(FILE* f, AlleleMatrix* m, char** markers);
void save_transposed_allele_matrix(FILE* f, AlleleMatrix* m, char** markers);
void save_group_alleles(FILE* f, SimData* d, int group_id);
void save_transposed_group_alleles(FILE* f, SimData* d, int group_id);

void save_group_one_step_pedigree(FILE* f, SimData* d, int group);
void save_one_step_pedigree(FILE* f, SimData* d);
void save_group_full_pedigree(FILE* f, SimData* d, int group);
void save_full_pedigree(FILE* f, SimData* d);
void save_AM_pedigree(FILE* f, AlleleMatrix* m, SimData* parents);
void save_parents_of(FILE* f, AlleleMatrix* m, unsigned int p1, unsigned int p2);

void save_group_bvs(FILE* f, SimData* d, int group);
void save_manual_bvs(FILE* f, DecimalMatrix* e, unsigned int* ids, char** names);
void save_bvs(FILE* f, SimData* d);

void save_count_matrix(FILE* f, SimData* d, char allele);
void save_count_matrix_of_group(FILE* f, SimData* d, char allele, int group);

/**@}*/
#endif
