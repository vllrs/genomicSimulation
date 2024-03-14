#ifndef SIM_OPERATIONS_H
#define SIM_OPERATIONS_H
/* 
genomicSimulationC v0.2.4.005 

    Last edit: 13 Mar 2024 
	License: MIT License

Copyright (c) 2021 Kira Villiers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

	Project repository: https://github.com/vllrs/genomicSimulationC

*/

/* Compile-time options:

	#define GSC_NO_SHORT_NAMES
		
		If simple names (eg `make_group_from`) might cause namespace conflicts, 
		then add the above definition before the import, and only module-prefixed 
		names (eg `gsc_make_group_from`) will be exposed. 
		
	#define GSC_MALLOC(ptr,size)  better_malloc
	#define GSC_FREE(ptr)         better_free
	
		By default, genomicSimulation uses stdlib malloc() and free(). These can
		be replaced with other equivalent memory management functions by redefining
		both GSC_MALLOC and GSC_FREE. Do not redefine only one of the pair.
	
	#define CONTIG_WIDTH 1000	
	
		The largest contiguous block of memory that could be requested in the 
		process of simulation is CONTIG_WIDTH integers. This setting's default
		value is 1000 if not redefined.
		
		This could be decreased to help long simulations or lower-end machines.
		Increasing this may provide some speed gain (likely very minor. I haven't
		tested).
		
	#define NAME_LENGTH 45
	
		The maximum number of characters allowed in a name field.
		These include names of SNPs, names of genotypes loaded from files,
		names of generated genotypes, and save-as-you-go filenames.
		
		@TODO render this setting meaningless. 
	
*/

#ifdef SIM_OPERATIONS
    #define RND_IMPLEMENTATION
#endif

#include <string.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

#define GSC_TRUE 1
#define GSC_FALSE 0
#define GSC_UNINIT -1  // uninitialised

#ifndef CONTIG_WIDTH
	#define CONTIG_WIDTH 1000
#endif
#ifndef NAME_LENGTH
	#define NAME_LENGTH 45
#endif

#if defined(GSC_MALLOC) && !defined(GSC_FREE) || !defined(GSC_MALLOC) && defined(GSC_FREE)
#error "You must define both GSC_MALLOC and GSC_FREE, or neither."
#endif
#if !defined(GSC_MALLOC) && !defined(GSC_FREE)
	#define GSC_MALLOC(size) malloc(size)
	#define GSC_FREE(ptr) free(ptr)
#endif

/** @defgroup shortnames Functions with Short Names Available
 *
 * List of functions and data structures that have short names available.
 *
 * Everything in this library is prefixed with `gsc_`. It might however
 * be more appealing to use the short names (without `gsc_` prefix). Short 
 * names are provided (listed below) for all data structures and user-friendly
 * functions. 
 *
 * Each function's documentation page will list its short name, if it has one.
 * Some functions are exposed but are not given short names: they may be useful 
 * for expert users trying to do specific things, but otherwise you can leave 
 * them alone. 
 *
 * If, due to namespace collision issues, you do not wish to have the short 
 * names available, simply `#define GSC_NO_SHORT_NAMES` before importing 
 * genomicSimulationC. 
 *
 * @{
 */

#ifndef GSC_NO_SHORT_NAMES
#define PedigreeID                 gsc_PedigreeID
#define NO_PEDIGREE                gsc_NO_PEDIGREE
#define GroupNum                   gsc_GroupNum
#define NO_GROUP                   gsc_NO_GROUP
#define EffectID                   gsc_EffectID
#define NOT_AN_EFFECT_SET          gsc_NOT_AN_EFFECT_SET
#define LabelID                    gsc_LabelID
#define NOT_A_LABEL                gsc_NOT_A_LABEL
#define GroupAndEffectSet          gsc_GroupAndEffectSet

#define TableSize                  gsc_TableSize
#define DecimalMatrix              gsc_DecimalMatrix
#define GenOptions                 gsc_GenOptions
#define BASIC_OPT                  gsc_BASIC_OPT
#define MarkerPosition             gsc_MarkerPosition
#define MarkerBlocks               gsc_MarkerBlocks
#define GeneticMap                 gsc_GeneticMap
#define AlleleMatrix               gsc_AlleleMatrix
#define EffectMatrix               gsc_EffectMatrix
#define SimData                    gsc_SimData
#define GenoLocation               gsc_GenoLocation
#define INVALID_GENO_LOCATION      gsc_INVALID_GENO_LOCATION
#define IS_VALID_LOCATION          gsc_IS_VALID_LOCATION
#define BidirectionalIterator      gsc_BidirectionalIterator
#define GappyIterator              gsc_GappyIterator
#define RandomAccessIterator       gsc_RandomAccessIterator

#define create_new_label           gsc_create_new_label
#define change_label_default       gsc_change_label_default
#define change_label_to            gsc_change_label_to
#define change_label_by_amount     gsc_change_label_by_amount
#define change_label_to_values     gsc_change_label_to_values
#define change_names_to_values     gsc_change_names_to_values

#define create_empty_simdata       gsc_create_empty_simdata
#define clear_simdata              gsc_clear_simdata
#define load_all_data              gsc_load_all_data
#define load_genmap                gsc_load_genmap
#define load_effects               gsc_load_effects
// #define load_genotypes             gsc_load_genotypes
#define load_genotypes_transposed  gsc_load_genotypes_transposed
#define load_genotypes_encoded_and_transposed  gsc_load_genotypes_encoded_and_transposed

#define create_bidirectional_iter  gsc_create_bidirectional_iter 
#define create_randomaccess_iter   gsc_create_randomaccess_iter
#define set_bidirectional_iter_to_start gsc_set_bidirectional_iter_to_start
#define set_bidirectional_iter_to_end   gsc_set_bidirectional_iter_to_end
#define next_forwards              gsc_next_forwards
#define next_backwards             gsc_next_backwards
#define next_get_nth               gsc_next_get_nth
#define get_name                   gsc_get_name
#define set_name                   gsc_set_name
#define get_alleles                gsc_get_alleles
#define get_first_parent           gsc_get_first_parent
#define get_second_parent          gsc_get_second_parent
#define get_id                     gsc_get_id
#define get_group                  gsc_get_group
#define set_group                  gsc_set_group
#define get_label_value            gsc_get_label_value
#define get_existing_groups        gsc_get_existing_groups
#define get_existing_group_counts  gsc_get_existing_group_counts
#define get_group_size             gsc_get_group_size
#define get_group_names            gsc_get_group_names
#define get_group_ids              gsc_get_group_ids
#define get_group_indexes          gsc_get_group_indexes
#define get_group_genes            gsc_get_group_genes
#define get_group_bvs              gsc_get_group_bvs
#define get_group_parent_names     gsc_get_group_parent_names
#define get_group_parent_ids       gsc_get_group_parent_ids
#define get_group_pedigrees        gsc_get_group_pedigrees

#define combine_groups             gsc_combine_groups
#define make_group_from            gsc_make_group_from
#define split_by_label_value       gsc_split_by_label_value
#define split_by_label_range       gsc_split_by_label_range
#define split_into_individuals     gsc_split_into_individuals
#define split_into_families        gsc_split_into_families
#define split_into_halfsib_families gsc_split_into_halfsib_families
#define split_evenly_into_two      gsc_split_evenly_into_two 
#define split_evenly_into_n        gsc_split_evenly_into_n
#define split_into_buckets         gsc_split_into_buckets
#define split_randomly_into_two    gsc_split_randomly_into_two
#define split_randomly_into_n      gsc_split_randomly_into_n
#define split_by_probabilities     gsc_split_by_probabilities 

#define generate_gamete            gsc_generate_gamete
#define generate_doubled_haploid   gsc_generate_doubled_haploid
#define generate_clone             gsc_generate_clone
#define make_random_crosses        gsc_make_random_crosses
#define make_random_crosses_between gsc_make_random_crosses_between
#define make_targeted_crosses      gsc_make_targeted_crosses
#define self_n_times               gsc_self_n_times
#define make_doubled_haploids      gsc_make_doubled_haploids
#define make_clones                gsc_make_clones
#define make_all_unidirectional_crosses   gsc_make_all_unidirectional_crosses
#define make_n_crosses_from_top_m_percent gsc_make_n_crosses_from_top_m_percent
#define make_crosses_from_file            gsc_make_crosses_from_file
#define make_double_crosses_from_file     gsc_make_double_crosses_from_file

#define split_by_bv                gsc_split_by_bv
#define calculate_group_bvs        gsc_calculate_group_bvs
#define calculate_bvs              gsc_calculate_bvs
#define calculate_group_count_matrix      gsc_calculate_group_count_matrix
#define calculate_count_matrix            gsc_calculate_count_matrix
#define calculate_full_count_matrix       gsc_calculate_full_count_matrix 
#define calculate_group_count_matrix_pair gsc_calculate_group_count_matrix_pair
#define calculate_count_matrix_pair       gsc_calculate_count_matrix_pair
#define create_evenlength_blocks_each_chr gsc_create_evenlength_blocks_each_chr
#define load_blocks                       gsc_load_blocks
#define calculate_group_local_bvs         gsc_calculate_group_local_bvs
#define calculate_local_bvs               gsc_calculate_local_bvs
#define calculate_optimal_haplotype          gsc_calculate_optimal_haplotype
#define calculate_optimal_possible_haplotype gsc_calculate_optimal_possible_haplotype
#define calculate_optimal_bv                 gsc_calculate_optimal_bv
#define calculate_optimal_possible_bv        gsc_calculate_optimal_possible_bv
#define calculate_minimal_bv                 gsc_calculate_minimal_bv

#define delete_group               gsc_delete_group
#define delete_label               gsc_delete_label
#define delete_genmap              gsc_delete_genmap
#define delete_eff_set             gsc_delete_eff_set
#define delete_simdata             gsc_delete_simdata
#define delete_markerblocks        gsc_delete_markerblocks
#define delete_dmatrix             gsc_delete_dmatrix
#define delete_bidirectional_iter  gsc_delete_bidirectional_iter
#define delete_randomaccess_iter   gsc_delete_randomaccess_iter

#define save_marker_blocks         gsc_save_marker_blocks
#define save_names_header          gsc_save_names_header
#define save_allele_matrix         gsc_save_allele_matrix
#define save_transposed_allele_matrix   gsc_save_transposed_allele_matrix
#define save_group_genotypes       gsc_save_group_genotypes
#define save_transposed_group_genotypes gsc_save_transposed_group_genotypes
#define save_count_matrix          gsc_save_count_matrix
#define save_group_count_matrix    gsc_save_group_count_matrix
#define save_one_step_pedigree     gsc_save_one_step_pedigree 
#define save_group_one_step_pedigree gsc_save_group_one_step_pedigree
#define save_full_pedigree         gsc_save_full_pedigree
#define save_group_full_pedigree   gsc_save_group_full_pedigree
#define save_bvs                   gsc_save_bvs
#define save_group_bvs             gsc_save_group_bvs
#endif

/** @} */


/** @defgroup structs Data Structures
 *
 * How the simulation stores data.
 *
 * genomicSimulation is a state-based package. These data structures
 * are used to store the library's data/state. Many use dynamically
 * allocated memory so call the relevant delete_ function if one exists
 * when you are finished with them.
 *
 * gsc_SimData is the central state/data
 * storage struct, and so is a required parameter to most user-facing functions.
 * It contains pointers to a gsc_GeneticMap (storing the loaded genome map), an
 * gsc_EffectMatrix (storing the loaded allele effects), and an gsc_AlleleMatrix
 * (storing metadata and genotypes of founders and simulated offspring).
 *
 * Other structs in this group (gsc_TableSize, gsc_MarkerBlocks, gsc_GenOptions) represented
 * specific types of data and are used as parameters and return value of certain
 * functions.
 *
 * @{
 */

/** A struct representing a single marker location. the attribute
 * `chromosome` represents the chromosome number, and `position` the
 * position on the chromosome in centiMorgans.
 *
 * @shortnamed{MarkerPosition}
 */
typedef struct {
	int chromosome; /**< The chromosome number */
	float position; /**< The distance in centiMorgans along the chromosome */
} gsc_MarkerPosition;

/** A simple struct used for returning the dimensions of a matrix or table.
 *
 * @shortnamed{TableSize}
 */
struct gsc_TableSize {
	int num_columns;
	int num_rows;
};

/** A struct used to store a set of blocks of markers.
 *
 * @shortnamed{MarkerBlocks}
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
     * are the indexes in the gsc_SimData of the markers that make up that block. */
	int** markers_in_block;
} gsc_MarkerBlocks;

/** A row-major heap matrix that contains floating point numbers. `dmatrix` functions
 * are designed to act on this matrix.
 *
 * Rows make the first index of the matrix and columns the second.
  *
 * @shortnamed{DecimalMatrix}
 */
typedef struct {
	double** matrix; /**< The actual matrix and contents */
	int rows;        /**< Number of rows in the matrix */
	int cols;        /**< number of columns in the matrix */
} gsc_DecimalMatrix;


/** A type representing a program-lifetime-unique identifier for a genotype,
 *  to be used in tracking pedigree.
 *
 * @shortnamed{PedigreeID}
 */
typedef struct {
    unsigned int id;
} gsc_PedigreeID;
/** Empty/null value for pedigree fields.
 *
 * @shortnamed{NO_PEDIGREE} */
#define gsc_NO_PEDIGREE (gsc_PedigreeID){.id=0}


/** A type representing the identifier of a group of genotypes
 *
 * @shortnamed{GroupNum}
 */
typedef struct {
    unsigned int num;
} gsc_GroupNum;
/** Empty/null value for group allocations.
 *
 * @shortnamed{NO_GROUP} */
#define gsc_NO_GROUP (gsc_GroupNum){.num=0}

/** A type representing a particular loaded set of marker effects
 *
 * @shortnamed{EffectID}
 */
typedef struct {
    int id;
} gsc_EffectID;
/** Empty/null value for effect set identifiers.
 *
 * @shortnamed{NOT_AN_EFFECT_SET} */
#define gsc_NOT_AN_EFFECT_SET (gsc_EffectID){.id=0}

/** A type representing a particular integer label
 *
 * @shortnamed{LabelID}
 */
typedef struct {
    int id;
} gsc_LabelID;
/** Empty/null value for custom label identifiers.
 *
 * @shortnamed{NOT_A_LABEL} */
#define gsc_NOT_A_LABEL (gsc_LabelID){.id=0}

/** Simple crate (stores a GroupNum and an EffectID, nothing more).
 *
 * @shortnamed{GroupAndEffectSet} */
struct gsc_GroupAndEffectSet {
    gsc_GroupNum group;
    gsc_EffectID effectSet;
};

/** A type that contains choices of settings for gsc_SimData functions that create a
 * new gsc_AlleleMatrix/generation.
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
 *
 * @shortnamed{GenOptions}
*/
typedef struct {
	int will_name_offspring; /**< A boolean: whether generated offspring should be given names. */
    const char* offspring_name_prefix; /**< If `will_name_offspring` is true, generated
                           * offspring are named [offspring_name_prefix][index]. */

	int family_size; /**< The number of offspring to produce from each cross.*/

	int will_track_pedigree; /**< A boolean: whether to track parentage of generated offspring.*/
	int will_allocate_ids; /**< A boolean: whether to allocate generated offspring session-
                            * unique IDs. IDs are used for pedigree tracking. The
                            * offspring of an anonymous individual (one without an ID)
                            * cannot identify that individual as their parent. */

    const char* filename_prefix; /**< A string used in save-as-you-go file names. */
	int will_save_pedigree_to_file; /**< A boolean. If true, the full/recursive
                            * pedigrees of every offspring generated in the cross
                            * are saved to "[filename_prefix}-pedigree.txt", even
                            * if `will_save_to_simdata` is false.
                            * Pedigrees are saved in the format of gsc_save_full_pedigree()*/
    gsc_EffectID will_save_bvs_to_file; /**< If equal to NOT_AN_EFFECT_SET, no bvs are calculated or saved.
                            * Otherwise, for each offspring in the cross,
                            * the breeding values according
                            * to the marker effect set with this gsc_EffectID
                            * are saved to "[filename_prefix}-bv.txt", even
                            * if `will_save_to_simdata` is false.
                            * BVs are saved in the format of gsc_save_bvs() */
	int will_save_alleles_to_file; /**< A boolean. If true, the set of alleles
                            * of every offspring generated in the cross
                            * are saved to "[filename_prefix}-genotype.txt", even
                            * if `will_save_to_simdata` is false.
                            * Genotypes are saved in the format of gsc_save_allele_matrix()*/
	int will_save_to_simdata; /**< A boolean. If true, the generated offspring exist
                            * in the gsc_SimData struct after the function executes.
                            * If false, they are discarded after creation. */
} gsc_GenOptions;


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
 * @shortnamed{GeneticMap}
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

	gsc_MarkerPosition* positions; /**< An array of MarkerPositions, ordered from lowest to highest. */
} gsc_GeneticMap;

/** A linked list entry that stores a matrix of alleles for a set of SNP markers
 * and genotypes.
 *
 * The simulation stores its genotypes in a list of AlleleMatrix nodes. Each node can 
 * store up to CONTIG_WIDTH genotypes.
 *
 * @shortnamed{AlleleMatrix}
*/
typedef struct gsc_AlleleMatrix gsc_AlleleMatrix;
struct gsc_AlleleMatrix {
    /** A matrix of SNP markers by lines/genotypes containing pairs of alleles
     * eg TT, TA. Use `alleles[line index][marker index * 2]` to get the
     * first allele and `alleles[lines index][marker index * 2 + 1]` to
     * get the second. If CONTIG_WIDTH lines are saved here, another
     * gsc_AlleleMatrix is added to the linked list when there's a need to save more.*/
	char* alleles[CONTIG_WIDTH];

	int n_genotypes; /**< Number of genotypes currently loaded in this matrix.*/
	int n_markers; /**< Number of markers across which genotypes are tracked. This has
                    * redundancy with gsc_SimData and other members of its linked list
                    * but it's good to know how big your own `alleles` array is.*/

    char* names[CONTIG_WIDTH]; /**< Array of dynamically allocated strings
                    * containing the names of the lines/genotypes in this matrix.
                    * Guaranteed to be NULL if they do not have names. */
    gsc_PedigreeID ids[CONTIG_WIDTH]; /**< Unique ID for each genotype. */
    gsc_PedigreeID pedigrees[2][CONTIG_WIDTH]; /**< Two lists of integer IDs of the
                    * parents of this genotype (if tracked), or 0 if we don't know/care.*/
    gsc_GroupNum groups[CONTIG_WIDTH]; /**< Group allocation of each genotype. */

    int n_labels; /**< Number of custom labels currently available to this gsc_AlleleMatrix. This has
                    * redundancy with gsc_SimData and other members of its linked list
                    * but it's good to know how big your own `labels` array is.*/
    int** labels; /**< Pointer to list of labels. Size of first dimension is n_labels,
                           * of second dimension is arrays of labels of length CONTIG_WIDTH*/

	gsc_AlleleMatrix* next; /**< Pointer to the next gsc_AlleleMatrix in the linked list,
                         * or NULL if this entry is the last. */
};

/** A type that stores a matrix of effect values and their names.
 *
 * @shortnamed{EffectMatrix}
 */
typedef struct {
	gsc_DecimalMatrix effects; /**< Effect on breeding value of alleles at markers.
        * Rows correspond to `effect_names`/alleles, columns to markers. */
	char* effect_names; /**< Character array containing allele characters ordered
        * to match rows of `effects`. */
} gsc_EffectMatrix;

/** Composite type that is used to run crossing simulations.
 *
 * The core of this type is a list of markers. These are used to index the rows
 * of the allele matrix and the position map, and the columns of the effect matrix.
 *
 * @shortnamed{SimData}
 */
typedef struct {
	int n_markers;  /**< The number of markers/length of `markers`. */
	char** markers; /**< Array of strings containing the names of markers. */

    int n_labels; /**< The number of custom labels in the simulation.*/
    gsc_LabelID* label_ids; /**< The identifier number of each label in the simulation, in order
                     * of their lookup index. */
    int* label_defaults; /**< Array containing the default (birth) value of each
                          * custom label. */

	gsc_GeneticMap map; /**< A gsc_GeneticMap. If this is set, then `markers`
                     * will be ordered and all markers have a known position.*/
	gsc_AlleleMatrix* m; /**< Pointer to an gsc_AlleleMatrix, which stores data and
                      * metadata of founders and simulated offspring. The
                      * gsc_AlleleMatrix is start of a linked list if there are
                      * many genotypes. */

    int n_eff_sets; /**< The number of sets of allele effects in the simulation **/
    gsc_EffectID* eff_set_ids; /**< The identifier number of each set of allele effects in the simulation,
                     * ordered by their lookup index. */
    gsc_EffectMatrix* e; /**< Array of n_eff_sets gsc_EffectMatrix, optional for the use of the simulation.
                     * Used for calculating breeding values from which alleles
                     * a genotype has at each marker.*/

    //CRANDOMGENERATOR /**< Random number generator working memory. */
    gsc_PedigreeID current_id; /**< Highest SimData-unique ID that has been generated
                              * so far. Used to track which IDs have already been
                              * given out.*/
    unsigned int n_groups; /**< Number of groups currently existing in simulation. It is
                        * guaranteed to never be less than the number of groups in simulation
                        * even if not perfectly accurate. */
} gsc_SimData;

extern const gsc_GenOptions gsc_BASIC_OPT;
/** @} */

/** @defgroup maths Mathematical functions
 *
 * For mathematical and statistical operations as required by the package.
 *
 * Includes matrix operations defined on a gsc_DecimalMatrix struct, and
 * draws from certain random distributions.
 *
 * @{
 */
gsc_DecimalMatrix gsc_generate_zero_dmatrix(const int r, const int c);
int 		  gsc_add_matrixvector_product_to_dmatrix(gsc_DecimalMatrix* result, const gsc_DecimalMatrix* a, const double* b);
int			  gsc_add_doublematrixvector_product_to_dmatrix(gsc_DecimalMatrix* result, const gsc_DecimalMatrix* amat, const double* avec,
                                              const gsc_DecimalMatrix* bmat, const double* bvec);

/** @} */

/** @defgroup supporters Utils/Supporting Functions
 *
 * @{
 */
struct gsc_TableSize gsc_get_file_dimensions(const char* filename, const char sep);
int gsc_get_from_ordered_uint_list(const unsigned int target, const unsigned int listLen, const unsigned int list[listLen]);
int gsc_get_from_ordered_pedigree_list(const gsc_PedigreeID target, const unsigned int listLen, const gsc_PedigreeID list[listLen]);
int gsc_get_from_unordered_str_list(const char* target, const int listLen, const char* list[listLen]);
void gsc_shuffle_up_to( unsigned int* sequence, const unsigned int total_n, const unsigned int n_to_shuffle);
unsigned int gsc_randomdraw_replacementrules(gsc_SimData* d, unsigned int max, unsigned int cap, unsigned int* member_uses, unsigned int noCollision);

gsc_LabelID gsc_create_new_label(gsc_SimData* d, const int setTo);
void gsc_change_label_default(gsc_SimData* d, const gsc_LabelID whichLabel, const int newDefault);
void gsc_change_label_to(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel, const int setTo);
void gsc_change_label_by_amount(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel, const int byValue);
void gsc_change_label_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const int startIndex, const gsc_LabelID whichLabel,
                          const int n_values, const int values[n_values]);

//static void gsc_get_sorted_markers(gsc_SimData* d, int actual_n_markers);
//static void gsc_get_chromosome_locations(gsc_SimData *d);

void gsc_change_names_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const int startIndex, const int n_values, const char* values[n_values]);
//static void gsc_set_names(gsc_AlleleMatrix* a, const char* prefix, const int suffix, const int from_index);
//static void gsc_set_ids(gsc_SimData* d, const int from_index, const int to_index);
int gsc_get_integer_digits(const int i);
int gsc_get_index_of_label( const gsc_SimData* d, const gsc_LabelID label );
int gsc_get_index_of_eff_set( const gsc_SimData* d, const gsc_EffectID eff_set_id );

gsc_LabelID gsc_get_new_label_id( const gsc_SimData* d );
gsc_EffectID gsc_get_new_eff_set_id( const gsc_SimData* d );
gsc_GroupNum gsc_get_next_free_group_num( const int n_existing_groups, const gsc_GroupNum* existing_groups, int* cursor,  gsc_GroupNum previous);
gsc_GroupNum gsc_get_new_group_num( gsc_SimData* d );
void gsc_get_n_new_group_nums( gsc_SimData* d, const int n, gsc_GroupNum* result);
void gsc_condense_allele_matrix( gsc_SimData* d);
//static void* gsc_malloc_wrap(const unsigned int size);

//static int gsc_helper_simdata_pos_compare(const void *pp0, const void *pp1);
//static int gsc_helper_descending_double_comparer(const void* pp0, const void* pp1);
//static int gsc_helper_ascending_double_comparer(const void* pp0, const void* pp1);
//static int gsc_helper_ascending_float_comparer(const void* p0, const void* p1);
/**@}*/


/** @defgroup loaders Setup Functions
 *
 * For setup of the simulation (loading founders, genetic maps, and optionally allele effects).
 *
 * @{
 */
gsc_AlleleMatrix* gsc_create_empty_allelematrix(const int n_markers, const int n_labels, const int labelDefaults[n_labels], const int n_genotypes);
gsc_SimData* gsc_create_empty_simdata();
void gsc_clear_simdata(gsc_SimData* d);

gsc_GroupNum gsc_load_genotypes_transposed(gsc_SimData* d, const char* filename);
gsc_GroupNum gsc_load_more_genotypes_transposed(gsc_SimData* d, const char* filename);
gsc_GroupNum gsc_load_genotypes_encoded_and_transposed(gsc_SimData* d, const char* filename);
void gsc_load_genmap(gsc_SimData* d, const char* filename);
gsc_EffectID gsc_load_effects(gsc_SimData* d, const char* filename);
struct gsc_GroupAndEffectSet gsc_load_all_data(gsc_SimData* d, const char* data_file, const char* map_file, const char* effect_file);
/** @} */


/** @defgroup getters Data Access and Search Functions
 *
 * For non-persistent access to simulation results and contents.
 *
 * @{
 */

    /** @defgroup iterators Genotype Iterators
     *
     * For iterating through the genotypes in the simulation.
     * It is possible to iterate through one group or through
     * every genotype in the simulation.
     *
     * @{
     */

/** An gsc_AlleleMatrix/gsc_AlleleMatrix index coordinate of a particular
 *  genotype in the simulation. To be used to look up details of
 *  that genotype using the `get_` family of functions.
 *
 * @shortnamed{GenoLocation}
*/
typedef struct {
    gsc_AlleleMatrix* localAM; /**< Pointer to the gsc_AlleleMatrix in which
                            * the genotype can be found. */
    int localPos; /**< Index in the localAM where the genotype can be
                   * found (min value: 0. Max value: CONTIG_WIDTH-1). */
} gsc_GenoLocation;

/** Constant representing a nonexistent location in the simulation.
 *
 * @shortnamed{INVALID_GENO_LOCATION}
 */
#define gsc_INVALID_GENO_LOCATION (gsc_GenoLocation){.localAM=0,.localPos=-1}
/** Check if a @ref GenoLocation is @ref INVALID_GENO_LOCATION
 *
 * @shortnamed{IS_VALID_LOCATION} */
#define gsc_IS_VALID_LOCATION(g) (g.localAM != 0 && g.localPos != -1)	

/** Identify whether a gsc_GenoLocation is INVALID_GENO_LOCATION
 *
 * @see IS_VALID_LOCATION IS_VALID_LOCATION has the same function. 
 * @param g location to check.
 * @return FALSE if g has either of the attributes of
 * INVALID_GENO_LOCATION, TRUE otherwise
 */
static inline int gsc_isValidLocation(const gsc_GenoLocation g) {
    // Either entry of INVALID_GENO_LOCATION is inappropriate in a valid gsc_GenoLocation
    return (g.localAM != INVALID_GENO_LOCATION.localAM &&
            g.localPos != INVALID_GENO_LOCATION.localPos);
}

/** A structure to iterate forwards and backwards through all
 *  genotypes in a gsc_SimData or through only the members of a group.
 *
 * @shortnamed{BidirectionalIterator}
 *
 *  @see gsc_create_bidirectional_iter
 */
typedef struct {
    gsc_SimData* d; /**< Simulation data through which to iterate */
    const gsc_GroupNum group; /**< Group through which to iterate. If it is 0,
                          * then iterate through all genotypes in the simulation.
                          * Otherwise, iterate through members of the group with
                          * this as their group number. */
    unsigned int localPos; /**< Local index (index within the cachedAM) of the genotype in the linked list
                       * of gsc_AlleleMatrix beginning at `d->m` where the
                       * iterator's 'cursor' currently sits. */

    gsc_AlleleMatrix* cachedAM; /**< Pointer to the gsc_AlleleMatrix from the linked list
                              * of gsc_AlleleMatrix beginning at `d->m` where the
                              * iterator's 'cursor' currently sits. Contains
                              * the genotype at `localPos`. */
    unsigned int cachedAMIndex; /**< Index of `cachedAM` in the linked list of
                                  * gsc_AlleleMatrix beginning at `d->m`. `d->m`
                                  * is considered to be index 0. */

    char atEnd; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                  * the last genotype (genotype with the highest index in the
                  * gsc_SimData) that fulfils the `group` critera of this iterator. */
    char atStart; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                    * the first genotype (genotype with the lowest index in the
                    * gsc_SimData) that fulfils the `group` critera of this iterator. */

} gsc_BidirectionalIterator;

/** A structure to iterate forwards through all
 *  positions in the gsc_AlleleMatrix linked list in gsc_SimData. Used
 *  in @see gsc_condense_allele_matrix. Internal, not recommended for end users.
 *
 * @shortnamed{GappyIterator}
 */
struct gsc_GappyIterator {
    gsc_GenoLocation cursor;
    unsigned int cursorAMIndex;
};

/** A structure to search and cache indexes of all
 *  genotypes in a gsc_SimData or of all the members of a group.
 *  @see gsc_create_randomaccess_iter
 *
 * @shortnamed{RandomAccessIterator}
 */
typedef struct {
    gsc_SimData* d; /**< Simulation data through which to iterate */
    const gsc_GroupNum group; /**< Group through which to iterate. If it is 0,
                          * then iterate through all genotypes in the simulation.
                          * Otherwise, iterate through members of the group with
                          * this as their group number. */

    unsigned int cacheSize; /**< Length in gsc_GenoLocations of `cache` */
    gsc_GenoLocation* cache; /**< Array iteratively updated with the known
                           * genotypes in the simulation that fulfil the
                           * `group` criteria of the iterator as they
                           * are discovered during calls to next_ functions */

    int largestCached; /**< Local/group index (that is, index in `cache`) of the
                         * highest cell in `cache` that has been filled. */
    int groupSize; /**< If the number of genotypes in the simulation that fulfil
                     * the iterator's `group` criteria is known, it is saved here.
                     * This value is left uninitialised until then. */
} gsc_RandomAccessIterator;

gsc_BidirectionalIterator gsc_create_bidirectional_iter( gsc_SimData* d, const gsc_GroupNum group);
gsc_RandomAccessIterator gsc_create_randomaccess_iter( gsc_SimData* d, const gsc_GroupNum group);

gsc_AlleleMatrix* gsc_get_nth_AlleleMatrix( gsc_AlleleMatrix* listStart, const unsigned int n);

gsc_GenoLocation gsc_set_bidirectional_iter_to_start(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_set_bidirectional_iter_to_end(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_forwards(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_backwards(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_get_nth(gsc_RandomAccessIterator* it, const unsigned int n);

//static gsc_GenoLocation gsc_nextgappy_get_gap(struct gsc_GappyIterator* it);
//static gsc_GenoLocation gsc_nextgappy_get_nongap(struct gsc_GappyIterator* it);
//static gsc_GenoLocation gsc_nextgappy_valid_pos(struct gsc_GappyIterator* it);
    /**@}*/

    /** @defgroup liteget Getting data from an Iterator
     *
     * These functions take as a parameter a gsc_GenoLocation (the output
     * of iterators), and access the data for the single genotype that
     * the iterator has found.
     *
     * They are very lightweight: they have no error-checking and
     * implemented inline.
     *
     * @{
     */
/** Get the name of a genotype
 *
 * @shortnamed{get_name}
 *
 * @param loc location of the relevant genotype
 * @return shallow copy of the name of the
 * genotype at location `loc`
 */
static inline char* gsc_get_name(const gsc_GenoLocation loc) {
    return loc.localAM->names[loc.localPos];
}

/** Set the name of a genotype
 *
 * @shortnamed{set_name}
 *
 * @param loc location of the relevant genotype
 * @param name name of the genotype
 * at location `loc` after this call
 */
static inline void gsc_set_name(const gsc_GenoLocation loc, char* name) {
    char* oldname = loc.localAM->names[loc.localPos];
    if (oldname != NULL) free(oldname);
    loc.localAM->names[loc.localPos] = name;
}

/** Get the alleles of a genotype
 *
 * @shortnamed{get_alleles}
 *
 * @param loc location of the relevant genotype
 * @return shallow copy of the allele string of the
 * genotype at location `loc`. The loci are ordered
 * according to the genetic map (by chromosome, then
 * by location). The entire string is 2*n characters
 * long, where n is the number of loci. The two
 * alleles of each locus are presented side-by-side,
 * at positions 2i and 2i+1.
 */
static inline char* gsc_get_alleles(const gsc_GenoLocation loc) {
    return loc.localAM->alleles[loc.localPos];
}

/** Get the first/left parent of a genotype
 *
 * @shortnamed{get_first_parent}
 *
 * @param loc location of the relevant genotype
 * @return id of the left parent of the genotype
 * at location `loc`
 */
static inline gsc_PedigreeID gsc_get_first_parent(const gsc_GenoLocation loc) {
    return loc.localAM->pedigrees[0][loc.localPos];
}

/** Get the second/right parent of a genotype
 *
 * @shortnamed{get_second_parent}
 *
 * @param loc location of the relevant genotype
 * @return id of the right parent of the genotype
 * at location `loc`
 */
static inline gsc_PedigreeID gsc_get_second_parent(const gsc_GenoLocation loc) {
    return loc.localAM->pedigrees[1][loc.localPos];
}

/** Get the persistent id of a genotype
 *
 *  The persistent id is the number used in pedigree
 *  tracing, and is unique for the lifetime of the simulation.
 *
 * @shortnamed{get_id}
 *
 * @param loc location of the relevant genotype
 * @return id of the genotype
 * at location `loc`
 */
static inline gsc_PedigreeID gsc_get_id(const gsc_GenoLocation loc) {
    return loc.localAM->ids[loc.localPos];
}

/** Get the current group membership of a genotype
 *
 * @shortnamed{get_group}
 *
 * @param loc location of the relevant genotype
 * @return group number of the group affiliation
 * of the genotype at location `loc`
 */
static inline gsc_GroupNum gsc_get_group(const gsc_GenoLocation loc) {
    return loc.localAM->groups[loc.localPos];
}

/** Set the current group membership of a genotype
 *
 * @shortnamed{set_group}
 *
 * @param loc location of the relevant genotype
 * @param group gsc_GroupNum of the group affiliation
 * of the genotype at location `loc` after this call
 */
static inline void gsc_set_group(const gsc_GenoLocation loc, const gsc_GroupNum group) {
    loc.localAM->groups[loc.localPos] = group;
}

//static inline int get_bv(GenotypeLocation loc) {}

/** Get the value of a specific label of a genotype
 *
 * @shortnamed{get_label_value}
 *
 * @param loc location of the relevant genotype
 * @param labelIndex index of the relevant label.
 * @return value of the `labelIndex`th label
 * of the genotype at location `loc`
 */
static inline int gsc_get_label_value(const gsc_GenoLocation loc, const int labelIndex) {
    return loc.localAM->labels[labelIndex][loc.localPos];
}
    /**@}*/

    /** @defgroup search Data Searching Functions
     *
     * These functions search the simulation data for the genotype
     * that matches a particular known piece of information, eg a
     * name, id, global index, or set of parents. Depending on
     * the size of the simulation, they may not be fast.
     *
     * @{
     */
char* gsc_get_name_of_id( const gsc_AlleleMatrix* start, const gsc_PedigreeID id);
int gsc_get_parents_of_id( const gsc_AlleleMatrix* start, const gsc_PedigreeID id, gsc_PedigreeID output[2]);
void gsc_get_ids_of_names( const gsc_AlleleMatrix* start, const int n_names, const char* names[n_names], gsc_PedigreeID* output);
int gsc_get_index_of_child( const gsc_AlleleMatrix* start, const gsc_PedigreeID parent1id, const gsc_PedigreeID parent2id);
int gsc_get_index_of_name( const gsc_AlleleMatrix* start, const char* name);
gsc_PedigreeID gsc_get_id_of_index( const gsc_AlleleMatrix* start, const int index);
char* gsc_get_genes_of_index( const gsc_AlleleMatrix* start, const int index);
    /**@}*/

    /** @defgroup collgetters Collective Data Access Functions
     *
     * These functions return vector data, rather than data on
     * a single genotype or single group.
     *
     * The `get_group_` family are left here for legacy purposes: using an
     * iterator is the new and less memory intensive way to do the tasks
     * these were used for.
     *
     * @{
     */
int gsc_get_group_size( const gsc_SimData* d, const gsc_GroupNum group_id);
int gsc_get_group_genes( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output);
int gsc_get_group_names( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output);
int gsc_get_group_ids( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, gsc_PedigreeID* output);
int gsc_get_group_indexes( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, unsigned int* output);
int gsc_get_group_bvs( const gsc_SimData* d, const gsc_GroupNum group_id, const gsc_EffectID effID, int group_size, double* output);
int gsc_get_group_parent_ids( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, const int whichParent, gsc_PedigreeID* output);
int gsc_get_group_parent_names( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, const int whichParent, char** output);
int gsc_get_group_pedigrees( const gsc_SimData* d, const gsc_GroupNum group_id, int group_size, char** output);

int gsc_get_existing_groups( gsc_SimData* d, gsc_GroupNum* output);
int gsc_get_existing_group_counts( gsc_SimData* d, gsc_GroupNum* out_groups, unsigned int* out_sizes);
    /**@}*/
/**@}*/


/** @defgroup groupmod Seletion/Group Modification Functions
 *
 * For simulation of selection or structure in breeding programs.
 *
 * @{
 */
gsc_GroupNum gsc_combine_groups( gsc_SimData* d, const int list_len, const gsc_GroupNum group_ids[list_len]);
gsc_GroupNum gsc_make_group_from( gsc_SimData* d, const int n, const unsigned int genotype_indexes[n]);
gsc_GroupNum gsc_split_by_label_value( gsc_SimData* d, const gsc_GroupNum group, const gsc_LabelID whichLabel, const int valueToSplit);
gsc_GroupNum gsc_split_by_label_range( gsc_SimData* d, const gsc_GroupNum group, const gsc_LabelID whichLabel, const int valueLowBound, const int valueHighBound);

// GENERIC
unsigned int gsc_scaffold_split_by_somequality( gsc_SimData* d, const gsc_GroupNum group_id,
        void* somequality_data,
        gsc_GroupNum (*somequality_tester)(gsc_GenoLocation, void*, unsigned int, unsigned int, gsc_GroupNum*),
        unsigned int maxentries_results, gsc_GroupNum* results);
// APPLICATIONS
unsigned int gsc_split_into_individuals( gsc_SimData* d, const gsc_GroupNum group_id, unsigned int maxentries_results, gsc_GroupNum results[maxentries_results]);
    //static gsc_GroupNum gsc_helper_split_by_quality_individuate(gsc_GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum** results);
unsigned int gsc_split_into_families(gsc_SimData* d, const gsc_GroupNum group_id, unsigned int maxentries_results, gsc_GroupNum results[maxentries_results]);
    //static gsc_GroupNum gsc_helper_split_by_quality_family(gsc_GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum** results);
unsigned int gsc_split_into_halfsib_families( gsc_SimData* d, const gsc_GroupNum group_id, const int parent, unsigned int maxentries_results, gsc_GroupNum results[maxentries_results]);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsib1(gsc_GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsib2(gsc_GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsibtemplate(gsc_GenoLocation loc, void** datastore, unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results, gsc_PedigreeID (*getparent)(gsc_GenoLocation));

// GENERIC
unsigned int gsc_scaffold_split_by_someallocation( gsc_SimData* d, const gsc_GroupNum group_id, void* someallocator_data,
        gsc_GroupNum (*someallocator)(gsc_GenoLocation, gsc_SimData*, void*, unsigned int, unsigned int*, gsc_GroupNum*),
        unsigned int n_outgroups, gsc_GroupNum outgroups[n_outgroups]);
// APPLICATIONS
gsc_GroupNum gsc_split_evenly_into_two(gsc_SimData* d, const gsc_GroupNum group_id);
    //static gsc_GroupNum gsc_helper_split_by_allocator_knowncounts(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups)
unsigned int gsc_split_evenly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, const int n, gsc_GroupNum* results);
unsigned int gsc_split_into_buckets(gsc_SimData* d, const gsc_GroupNum group_id, const int n, const int* counts, gsc_GroupNum* results);
gsc_GroupNum gsc_split_randomly_into_two(gsc_SimData* d, const gsc_GroupNum group_id);
unsigned int gsc_split_randomly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, const int n, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_allocator_equalprob(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups);
unsigned int gsc_split_by_probabilities(gsc_SimData* d, const gsc_GroupNum group_id, const int n, const double* probs, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_allocator_unequalprob(gsc_GenoLocation loc, gsc_SimData* d, void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups);

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
void gsc_generate_gamete(gsc_SimData* d, const char* parent_genome, char* output);
void gsc_generate_doubled_haploid(gsc_SimData* d, const char* parent_genome, char* output);
void gsc_generate_clone(gsc_SimData* d, const char* parent_genome, char* output);
    /**@}*/

// SUPPORTER FUNCS FOR GENOPTIONS OPTIONS
// static FILE* gsc_helper_genoptions_save_pedigrees_setup(const gsc_GenOptions g);
// static FILE* gsc_helper_genoptions_save_bvs_setup(const gsc_SimData* d, const gsc_GenOptions g, int* effIndexp);
// static FILE* gsc_helper_genoptions_save_genotypes_setup(const gsc_SimData* d, const gsc_GenOptions g);
// static void gsc_helper_genoptions_save_pedigrees(FILE* fp, gsc_SimData* d, gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_save_bvs(FILE* fe, gsc_EffectMatrix* effMatrices, int effIndex, gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_save_genotypes(FILE* fg, gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_give_names_and_ids(gsc_AlleleMatrix* am, gsc_SimData* d, const gsc_GenOptions g);

// PARAMETER FUNCTIONS FOR THE FOLLOWING GENERIC
// static void gsc_helper_make_offspring_cross(gsc_SimData* d, void* datastore, gsc_GenoLocation parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_self_n_times(gsc_SimData* d, void* datastore, gsc_GenoLocation parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_doubled_haploids(gsc_SimData* d, void* datastore, gsc_GenoLocation parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_clones(gsc_SimData* d, void* datastore, gsc_GenoLocation parents[static 2], gsc_GenoLocation putHere);

// static int gsc_helper_parentchooser_cross_randomly(void* parentIterator, void* datastore, unsigned int* counter, gsc_GenoLocation parents[static 2]);
// static int gsc_helper_parentchooser_cross_randomly_between(void* parentIterator, void* datastore, unsigned int* counter, gsc_GenoLocation parents[static 2]);
// static int gsc_helper_parentchooser_cross_targeted(void* parentIterator, void* datastore, unsigned int* counter, gsc_GenoLocation parents[static 2]);
// static int gsc_helper_parentchooser_selfing(void* parentIterator, void* datastore, unsigned int* counter, gsc_GenoLocation parents[static 2]);
// static int gsc_helper_parentchooser_cloning(void* parentIterator, void* datastore, unsigned int* counter, gsc_GenoLocation parents[static 2]);

// static int gsc_helper_random_cross_checks(gsc_SimData* d, const gsc_GroupNum from_group, const int n_crosses, const int cap);

// GENERIC
gsc_GroupNum gsc_scaffold_make_new_genotypes(gsc_SimData* d, const gsc_GenOptions g,
        void* parentIterator, void* datastore,
        int (*parentChooser)(void*, void*, unsigned int*, gsc_GenoLocation[static 2]),
        void (*offspringGenerator)(gsc_SimData*, void*, gsc_GenoLocation[static 2], gsc_GenoLocation) );
// APPLICATIONS
gsc_GroupNum gsc_make_random_crosses(gsc_SimData* d, const gsc_GroupNum from_group, const int n_crosses, const int cap, const gsc_GenOptions g);
gsc_GroupNum gsc_make_random_crosses_between(gsc_SimData*d, const gsc_GroupNum group1, const gsc_GroupNum group2, const int n_crosses, const int cap1, const int cap2, const gsc_GenOptions g);
gsc_GroupNum gsc_make_targeted_crosses(gsc_SimData* d, const int n_combinations, const int* firstParents, const int* secondParents, const gsc_GenOptions g);
gsc_GroupNum gsc_self_n_times(gsc_SimData* d, const unsigned int n, const gsc_GroupNum group, const gsc_GenOptions g);
gsc_GroupNum gsc_make_doubled_haploids(gsc_SimData* d, const gsc_GroupNum group, const gsc_GenOptions g);
gsc_GroupNum gsc_make_clones(gsc_SimData* d, const gsc_GroupNum group, const int inherit_names, const gsc_GenOptions g);

gsc_GroupNum gsc_make_all_unidirectional_crosses(gsc_SimData* d, const gsc_GroupNum from_group, const gsc_GenOptions g);
gsc_GroupNum gsc_make_n_crosses_from_top_m_percent(gsc_SimData* d, const int n, const int m, const gsc_GroupNum group, const gsc_EffectID effID, const gsc_GenOptions g);
gsc_GroupNum gsc_make_crosses_from_file(gsc_SimData* d, const char* input_file, const gsc_GenOptions g);
gsc_GroupNum gsc_make_double_crosses_from_file(gsc_SimData* d, const char* input_file, const gsc_GenOptions g);
/**@}*/


/** @defgroup calculators Breeding Value and Allele Count Calculators
 *
 * For calculations related to the loaded allele effects and internal additive breeding value model.
 *
 * @{
 */
gsc_GroupNum gsc_split_by_bv(gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID, const int top_n, const int lowIsBest);
gsc_DecimalMatrix gsc_calculate_group_bvs(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID);
gsc_DecimalMatrix gsc_calculate_bvs( const gsc_AlleleMatrix* m, const gsc_EffectMatrix* e);
int gsc_calculate_group_count_matrix( const gsc_SimData* d, const gsc_GroupNum group, const char allele, gsc_DecimalMatrix* counts);
int gsc_calculate_group_count_matrix_pair( const gsc_SimData* d, const gsc_GroupNum group, const char allele, gsc_DecimalMatrix* counts, const char allele2, gsc_DecimalMatrix* counts2);
int gsc_calculate_count_matrix( const gsc_AlleleMatrix* m, const char allele, gsc_DecimalMatrix* counts);
int gsc_calculate_count_matrix_pair( const gsc_AlleleMatrix* m , const char allele, gsc_DecimalMatrix* counts, const char allele2, gsc_DecimalMatrix* counts2);
gsc_DecimalMatrix gsc_calculate_full_count_matrix( const gsc_AlleleMatrix* m, const char allele);

gsc_MarkerBlocks gsc_create_evenlength_blocks_each_chr(const gsc_SimData* d, const int n);
gsc_MarkerBlocks gsc_load_blocks(const gsc_SimData* d, const char* block_file);
void gsc_calculate_group_local_bvs(const gsc_SimData* d, const gsc_MarkerBlocks b, const gsc_EffectID effID, const char* output_file, const gsc_GroupNum group);
void gsc_calculate_local_bvs(const gsc_SimData* d, const gsc_MarkerBlocks b, const gsc_EffectID effID, const char* output_file);

char* gsc_calculate_optimal_haplotype(const gsc_SimData* d, const gsc_EffectID effID);
char* gsc_calculate_optimal_possible_haplotype(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID);
double gsc_calculate_optimal_bv(const gsc_SimData* d, const gsc_EffectID effID);
double gsc_calculate_optimal_possible_bv(const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID);
double gsc_calculate_minimal_bv(const gsc_SimData* d, const gsc_EffectID effID);
/**@}*/


/** @defgroup deletors Deletor Functions
 *
 * For deleting and free associated memory of data structures.
 *
 * @ingroup structs
 * @{
 */
void gsc_delete_group(gsc_SimData* d, const gsc_GroupNum group_id);
void gsc_delete_label(gsc_SimData* d, const gsc_LabelID whichLabel);
void gsc_delete_genmap(gsc_GeneticMap* m);
void gsc_delete_allele_matrix(gsc_AlleleMatrix* m);
void gsc_delete_effect_matrix(gsc_EffectMatrix* m);
void gsc_delete_eff_set(gsc_SimData* d, gsc_EffectID whichID);
void gsc_delete_simdata(gsc_SimData* m);
void gsc_delete_markerblocks(gsc_MarkerBlocks* b);
void gsc_delete_dmatrix(gsc_DecimalMatrix* m);
void gsc_delete_bidirectional_iter(gsc_BidirectionalIterator* it);
void gsc_delete_randomaccess_iter(gsc_RandomAccessIterator* it);

void gsc_move_genotype(gsc_GenoLocation from, gsc_GenoLocation to, int* label_defaults);

/**@}*/


/** @defgroup savers Saving Functions
 *
 * For saving persistent simulation results.
 *
 * @{
 */

void gsc_save_marker_blocks(FILE* f, const gsc_SimData* d, const gsc_MarkerBlocks b);

void gsc_save_names_header(FILE* f, unsigned int n, const char* names[n]);
void gsc_save_allele_matrix(FILE* f, const gsc_AlleleMatrix* m);
void gsc_save_transposed_allele_matrix(FILE* f, const gsc_AlleleMatrix* m, const char** markers);

void gsc_save_group_genotypes(FILE* f, gsc_SimData* d, const gsc_GroupNum group_id);
void gsc_save_transposed_group_genotypes(FILE* f, const gsc_SimData* d, const gsc_GroupNum group_id);

void gsc_save_count_matrix(FILE* f, const gsc_SimData* d, const char allele);
void gsc_save_group_count_matrix(FILE* f, const gsc_SimData* d, const char allele, const gsc_GroupNum group);

void gsc_save_one_step_pedigree(FILE* f, const gsc_SimData* d);
void gsc_save_group_one_step_pedigree(FILE* f, const gsc_SimData* d, const gsc_GroupNum group);
void gsc_save_full_pedigree(FILE* f, const gsc_SimData* d);
void gsc_save_group_full_pedigree(FILE* f, const gsc_SimData* d, const gsc_GroupNum group);
void gsc_save_allelematrix_full_pedigree(FILE* f, const gsc_AlleleMatrix* m, const gsc_SimData* parents);
void gsc_save_parents_of(FILE* f, const gsc_AlleleMatrix* m, gsc_PedigreeID p1, gsc_PedigreeID p2);

void gsc_save_bvs(FILE* f, const gsc_SimData* d, const gsc_EffectID effID);
void gsc_save_group_bvs(FILE* f, const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID);
void gsc_save_manual_bvs(FILE* f, const gsc_DecimalMatrix* e, const gsc_PedigreeID* ids, const char** names);


/**@}*/


/** @defgroup recomb Recombination Calculators
 *
 * Experimental functions for retroactively calculating number of recombinations.
 *
 * This functionality is for interest only. It is not clear, or tidy,
 * or checked against real data.
 *
 * @{
 */
int* gsc_calculate_min_recombinations_fw1(gsc_SimData* d, char* parent1, unsigned int p1num, char* parent2,
        unsigned int p2num, char* offspring, int certain); // forward filling, window size 1
int* gsc_calculate_min_recombinations_fwn(gsc_SimData* d, char* parent1, unsigned int p1num, char* parent2,
        unsigned int p2num, char* offspring, int window_size, int certain); // forward filling, window size n

/** Simple operator to determine if at marker i, two genotypes share at least
 * one allele. Checks only 3 of four possible permutations because assumes
 * there cannot be more than two alleles at a given marker.
 *
 * @param p1 pointer to a character array genotype of the type stored in an gsc_AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param p2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param i index of the marker at which to perform the check
 * @returns boolean result of the check
 */
static inline int gsc_has_same_alleles(const char* p1, const char* p2, const int i) {
    return (p1[i<<1] == p2[i<<1] || p1[(i<<1) + 1] == p2[i] || p1[i] == p2[(i<<1) + 1]);
}
// w is window length, i is start value
/** Simple operator to determine if at markers with indexes i to i+w inclusive, two genotypes
 * share at least one allele. Checks only 3 of four possible permutations at each marker
 * because assumes there cannot be more than two alleles at a given marker. For the return value
 * to be true, there must be at least one match at every one of the markers in the window.
 *
 * @param g1 pointer to a character array genotype of the type stored in an gsc_AlleleMatrix
 * (2*n_markers long, representing the two alleles at a marker consecutively) for the first
 * of the genotypes to compare.
 * @param g2 pointer to a character array genotype for the second of the genotypes to compare.
 * @param start index of the first marker in the window over which to perform the check
 * @param w length of the window over which to perform the check
 * @returns boolean result of the check
 */
static inline int gsc_has_same_alleles_window(const char* g1, const char* g2, const int start, const int w) {
    int same = GSC_TRUE;
    int i;
    for (int j = 0; j < w; ++j) {
        i = start + j;
        same = same && (g1[i<<1] == g2[i<<1] || g1[(i<<1) + 1] == g2[i] || g1[i] == g2[(i<<1) + 1]);
    }
    return same;
}

int gsc_calculate_recombinations_from_file(gsc_SimData* d, const char* input_file, const char* output_file,
        int window_len, int certain);
/**@}*/

#endif
