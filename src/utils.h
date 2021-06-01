#ifndef SIM_UTILS_H
#define SIM_UTILS_H

#include <string.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#define TRUE 1
#define FALSE 0

/** A struct representing a single marker location. the attribute 
 * `chromosome` represents the chromosome number, and `position` the 
 * position on the chromosome in centiMorgans. 
 * This should probably be used in an array indexed by marker.
 */
typedef struct {
	int chromosome;
	float position;
} MarkerPosition;

/** A struct used to store a set of blocks of markers. 
 *
 * @param num_blocks the number of blocks whose details are stored here.
 * @param num_markers_in_block pointer to a heap array of length num_blocks
 * containing the number of markers that make up each block
 * @param markers_in_block pointer to a heap array of length num_blocks, each
 * entry in which is a pointer to a heap array with length corresponding to 
 * the value of the corresponding entry in num_markers_in_block whose values
 * are the indexes in the SimData of the markers that make up that block.
 */
typedef struct {
	int num_blocks;
	int* num_markers_in_block;
	int** markers_in_block;
} MarkerBlocks;

/** A simple struct used for returning the dimensions of a matrix or table.*/
struct TableSize {
	int num_columns;
	int num_rows;
};

/** A heap matrix that contains floating point numbers. `dmatrix` functions
 * are designed to act on this matrix.
 *
 * Rows make the first index of the matrix and columns the second.
 */
typedef struct {
	double** matrix; // the actual matrix with its contents
	int rows; // number of rows
	int cols; // number of columns
} DecimalMatrix;


/** A type that contains choices of settings for SimData functions that create a 
 * new AlleleMatrix/generation.
 *
 * The family_size parameter will affect how many offspring are produced.
 *
 * The will_name_subjects, will_track_pedigree, and will_allocate_ids parameters
 * affect how much extra detail about the offspring is generated/saved.
 *
 * The will_save_to_simdata toggle allows you the option of generating offspring without
 * saving them in memory. This may be useful in combination with save-as-you-go toggles 
 * will_save_pedigree_to_file, will_save_effects_to_file, and will_save_genes_to_file,
 * to generate a larger number of offspring than will fit in memory.
 *
 * @param will_name_subjects a boolean representing if subject_names will be 
 * filled or not.
 * @param subject_prefix If `will_name_subjects` is true, subjects are named 
 * [subject_prefix][index].
 * @param family_size the number of offspring to produce from each cross.
 * @param will_track_pedigree a boolean representing whether to bother to track
 * parentage of individual outcomes of the cross or not.
 * @param will_allocate_ids a boolean representing whether to allocate individuals
 * session-unique ids used for pedigree tracking.
 * @param filename_prefix a string
 * @param will_save_pedigree_to_file a boolean. If true, the full/recursive pedigrees
 * of every genotype generated in the cross are saved to "[filename_prefix]-pedigree",
 * even if the genotypes are not later saved to SimData.
 * @param will_save_effects_to_file a boolean. If true, the GEBVs
 * of every genotype generated in the cross are saved to "[filename_prefix]-eff",
 * even if the genotypes are not later saved to SimData.
 * @param will_save_genes_to_file a boolean. If true, the set of alleles
 * of every genotype generated in the cross are saved to "[filename_prefix]-genome",
 * even if the genotypes are not later saved to SimData.
 * @param will_save_to_simdata a boolean. If true, the offspring are retained in the 
 * SimData as a new group. If false, they are discarded after creation.
*/
typedef struct {
	int will_name_subjects;
	const char* subject_prefix;
	
	int family_size;
	
	int will_track_pedigree;
	int will_allocate_ids;
	
	const char* filename_prefix;
	int will_save_pedigree_to_file;
	int will_save_effects_to_file;
	int will_save_genes_to_file;
	int will_save_to_simdata;
} GenOptions;


/** A type that stores the genetic map for a set of markers.
 *
 * To get all markers belonging to a particular chromosome, use the following rule:
 * Chr n includes all markers in `positions` starting at index chr_ends[n-1] up
 * but not including the marker at index chr_ends[n]
 *
 * @param n_chr the number of chromosomes represented in the map.
 * @param chr_ends An array of ints. The entry at index i is the index in 
 * `positions` of the first marker that belongs to Chr(i + 1). The array is
 * n_chr + 1 integers long.
 * @param chr_lengths An array of floats. The entry at index i is the length
 * of Chr(i + 1), calculated by `position of last marker - position of first marker`.
 * The array is n_chr entries long.
 * @param positions An array of MarkerPositions, ordered from lowest to highest.
*/
typedef struct {
	int n_chr;
	int* chr_ends;
	float* chr_lengths;
	
	MarkerPosition* positions;
} GeneticMap;

/** A type used for containing the basic genotype and parentage of a new cross.
 *
 * @param gametes the alleles inherited from each parent as character sequences 
 * whose order matches the order of the markers in the main SimData object that they
 * are alleles for. There are two of these sequences, one from each parent.
 * @param p_index the indexes/column numbers of each parent. The generation AlleleMatrix
 * will be visible in the call to the function that returns a NewCross
 */
typedef struct {
	char* gametes[2];
	int p_index[2];
} NewCross;

/** A linked list entry that stores a matrix of alleles for a set of SNP markers
 * and subjects.
 *
 * @param alleles a matrix of SNP markers by subjects containing pairs of alleles
 * eg TT, TA. Use `alleles[subject index][marker index * 2]` to get the first allele
 * and `alleles[subject index][marker index * 2 + 1]` to get the second. In usual use
 * the char** alleles will point to an 1000-ptr-long piece of heap memory, and if you
 * need to add more than 100 subjects you will add another AlleleMatrix to the linked list
 * @param subject_names array of strings containing the names of the lines/subjects
 * whose data is stored in `alleles`. NULL if they do not have names.
 * @param ids unique ids for each subject in the matrix. If you use a creator function,
 * these will be generated for you.
 * @param n_subjects the number of subjects currently loaded in `alleles`
 * @param n_markers the number of markers in `alleles`
 * @param pedigrees two lists of integer IDs of the parents of this subject if tracked,
 * or 0 if we don't know/care about this pedigree.
 * @param next pointer to the next AlleleMatrix in the linked list, or NULL
 * if this entry is the last.
*/
typedef struct AlleleMatrix AlleleMatrix;
struct AlleleMatrix {
	
	char* alleles[1000];
	char* subject_names[1000];
	unsigned int ids[1000];
	int n_subjects;
	int n_markers; // slight redundancy but allows this to stand alone
	
	unsigned int pedigrees[2][1000]; 
	unsigned int groups[1000];
	AlleleMatrix* next;
};

/** A type that stores a matrix of effect values and their names.
 *
 * @param effects the effect of each marker. rows correspond to `effect_names`
 * columns correspond to markers.
 * @param effect_names the one-character name of each row in `effects`.
 * Eg 'A' for allele A, 'r' for range between homozygous ref and alternate.
 */
typedef struct {
	DecimalMatrix effects;
	char* effect_names;
} EffectMatrix;

/** Composite type that is used to run crossing simulations.
 *
 * The core of this type is a list of markers. These are used to index the rows
 * of the allele matrix and the position map, and the cols of the effect matrix.
 * 
 * @param n_markers the number of markers/length of `markers`
 * @param markers an array of strings containing the names of markers
 * @param map GeneticMap with positions of markers on chromosomes. If this is 
  * set, then `markers` will be ordered and all markers have a known position.
 * @param m pointer to an AlleleMatrix containing data for the parent generation.
 * Successive generations are added to the end of the linked list that starts here.
 * @param e EffectMatrix containing the effects at all markers.
 * @param current_id integer denoting the highest id that has been allocated to a 
 * subject. Used to track where we are in generating unique ids.
 */
typedef struct {
	int n_markers;
	char** markers; // marker names
	
	GeneticMap map;
	AlleleMatrix* m;
	EffectMatrix e;
	
	unsigned int current_id;
} SimData; 

typedef struct {
	unsigned int gid;
	char* name; // optional field for a group name
} LineGroup;

const GenOptions BASIC_OPT;

AlleleMatrix* create_empty_allelematrix(int n_markers, int n_subjects);
SimData* create_empty_simdata();

/* Supporters */
struct TableSize get_file_dimensions(const char* filename, char sep);
int get_from_ordered_uint_list(unsigned int target, unsigned int* list, unsigned int list_len);
int get_from_unordered_str_list(char* target, char** list, int list_len) ;

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

void set_subject_names(AlleleMatrix* a, const char* prefix, int suffix, int from_index);
void set_subject_ids(SimData* d, int from_index, int to_index);

int get_integer_digits(int i);
int get_new_group_num( SimData* d);
LineGroup get_new_group(SimData* d, const char* groupName);

void set_group_list( SimData* d, int by_n, int new_group); 
void condense_allele_matrix( SimData* d); 
void* get_malloc(size_t size);

int _simdata_pos_compare(const void *pp0, const void *pp1);
int _descending_double_comparer(const void* pp0, const void* pp1);
int _ascending_double_comparer(const void* pp0, const void* pp1);
int _ascending_float_comparer(const void* p0, const void* p1);
int _ascending_int_comparer(const void* p0, const void* p1);
int _ascending_int_dcomparer(const void* pp0, const void* pp1);


/* Group Modification */
int combine_groups( SimData* d, int list_len, int group_ids[list_len]);
void split_into_individuals( SimData* d, int group_id);
void split_into_families(SimData* d, int group_id);
int split_from_group( SimData* d, int n, int indexes_to_split[n]);
int* get_existing_groups( SimData* d, int* n_groups);
int** get_existing_group_counts( SimData* d, int* n_groups);

/* Deletors */
void delete_group(SimData* d, int group_id);
void delete_genmap(GeneticMap* m);
void delete_allele_matrix(AlleleMatrix* m);
void delete_effect_matrix(EffectMatrix* m);
void delete_simdata(SimData* m);
void delete_markerblocks(MarkerBlocks* b);
void SXP_delete_simdata(SEXP sd);

SEXP SXP_delete_group(SEXP exd, SEXP group);
SEXP SXP_combine_groups(SEXP exd, SEXP len, SEXP groups);
SEXP SXP_split_individuals(SEXP exd, SEXP group);
SEXP SXP_split_familywise(SEXP exd, SEXP group);
SEXP SXP_split_out(SEXP exd, SEXP len, SEXP indexes);
SEXP SXP_get_groups(SEXP exd);
SEXP SXP_get_group_data(SEXP exd, SEXP group, SEXP whatData);

DecimalMatrix generate_zero_dmatrix(int r, int c);
DecimalMatrix subset_dmatrix_row(DecimalMatrix* m, int row_index);
DecimalMatrix add_dmatrices(DecimalMatrix* a, DecimalMatrix* b);
void add_to_dmatrix(DecimalMatrix* a, DecimalMatrix* b);
DecimalMatrix multiply_dmatrices(DecimalMatrix* a, DecimalMatrix* b);
void delete_dmatrix(DecimalMatrix* m);


#endif

