#ifndef SIM_OPERATIONS_H
#define SIM_OPERATIONS_H
/* 
genomicSimulationC v0.2.6.17
// Converted using Rconversion.sh v2

    Last edit: 12 June 2025
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
        
    #define GSC_ID_T unsigned int
    #define GSC_LOCALX_T unsigned int
    #define GSC_GLOBALX_T unsigned int
    #define GSC_GENOLEN_T unsigned int
        
        These type aliases are used for different types of internal identifiers.
        They can be redefined to any signed or unsigned integer type, if needed.
    
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
#include <ctype.h> // for isdigit etc.
#include <stdint.h> // for SIZE_MAX etc.

/** genomicSimulation's "logical value" type
 *
 * A type representing one of three possible states: true, false, or invalid/uninitialised.
 * This type is only used when the invalid/uninitialised value is needed: otherwise, use _Bool
 */
typedef enum {
    GSC_TRUE = 1,
    GSC_FALSE = 0,
    GSC_NA = -1
} GSC_LOGICVAL;

#ifndef GSC_ID_T
    /** genomicSimulation's "ID" type 
    *
    * A type representing the maximum number of unique session IDs.
    *
    * It is used by the various ID types (PedigreeID, MapID, LabelID, etc.)
    * as well as being used as the type to represent the index of the stored
    * map/label/effect set, because the maximum potential size of these values
    * is pinned to the size of the corresponding ID type.
    *
    * Can be redefined to another unsigned or signed integer type.
    */
    #define GSC_ID_T unsigned int
#endif
/** For unique session IDs, the INVALID/UNINITIALISED value is 0. */
#define GSC_NA_ID 0
/** When accessing the current array index of a unique session ID,
 * the "ID not found"/failure value is -1 (for signed types) or 
 * the maximum value of the type (for unsigned types).
 * i.e. if GSC_ID_T is an unsigned type, the .id field of a GroupNum,
 * PedigreeID, MapID, etc., could hold this value, but it is the 
 * failure/NA value for an index into a list of GroupNums, PedigreeIDs, etc.
 */
#define GSC_NA_IDX (GSC_ID_T)-1
// ^ Equivalent but less flexible definition: #define GSC_NA_IDX UINT_MAX

#ifndef GSC_GLOBALX_T
    /** genomicSimulation's "Candidate global index" type 
    *
    * A type representing the maximum number of candidates stored in simulation. 
    *
    * It is used to represent the maximum number of candidates that could exist
    * or be produced in the simulation, and to represent "global" candidate indexes
    * counted cumulatively from the first AlleleMatrix in the SimData linked list.
    *
    * This type should be the same size as, or larger than, GSC_LOCALX_T.
    *
    * Can be redefined to another unsigned or signed integer type.
    */
    #define GSC_GLOBALX_T unsigned int
#endif
/** For candidate global indexes, the INVALID/UNINITIALISED value is 
 * -1 (for signed types) or the maximum value of the type (for unsigned types). */
#define GSC_NA_GLOBALX (GSC_GLOBALX_T)-1

#ifndef GSC_LOCALX_T
    /** genomicSimulation's "Candidate local index" type
    *
    * A type at least large enough to store the number of candidates in a single
    * AlleleMatrix (i.e., a type whose maximum value is greater than or equal to
    * CONTIG_WIDTH).
    *
    * It is used to represent the index of a candidate within its current AlleleMatrix.
    *
    * Can be redefined to another unsigned or signed integer type. This type should have a
    * maximum value less than CONTIG_WIDTH.
    */ 
    #define GSC_LOCALX_T unsigned int
#endif
/** For candidate local indexes, the INVALID/UNINITIALISED value is 
 * -1 (for signed types) or the maximum value of the type (for unsigned types). */
#define GSC_NA_LOCALX (GSC_LOCALX_T)-1

#ifndef GSC_GENOLEN_T
    /** genomicSimulation's "Genotype length" type
    *
    * A type representing the maximum number of genetic markers tracked by the
    * simulation, or the index of a particular genetic marker inside the map of 
    * markers in the simulation.
    *
    * Can be redefined to another unsigned or signed integer type. 
    */
    #define GSC_GENOLEN_T unsigned int
#endif
/** For genetic marker indexes, the INVALID/UNINITIALISED value is
 * -1 (for signed types) or the maximum value of the type (for unsigned types). */
#define GSC_NA_GENOLEN (GSC_GENOLEN_T)-1

#ifndef CONTIG_WIDTH
    #define CONTIG_WIDTH 1000
#else
    #if CONTIG_WIDTH > SIZE_MAX
        #error "This C compiler cannot handle a CONTIG_WIDTH that large. Please redefine."
    #elif CONTIG_WIDTH > UINT_MAX
        #undef GSC_LOCALX_T
        #undef GSC_NA_LOCALX
        #define GSC_LOCALX_T unsigned int
        #define GSC_NA_LOCALX SIZE_MAX
        #undef GSC_GLOBALX_T
        #undef GSC_NA_GLOBALX
        #define GSC_GLOBAL_T unsigned int
        #define GSC_NA_GLOBALX SIZE_MAX
    #endif //else LOCALX is already defined as unsigned int
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

#define GENOLEN_T                  GSC_GENOLEN_T
#define LOCALX_T                   GSC_LOCALX_T
#define GLOBALX_T                  GSC_GLOBALX_T
#define ID_T                       GSC_ID_T
#define LOGICVAL                   GSC_LOGICVAL
#define NA_GENOLEN                 GSC_NA_GENOLEN
#define NA_LOCALX                  GSC_NA_LOCALX
#define NA_GLOBALX                 GSC_NA_GLOBALX
#define NA_IDX                     GSC_NA_IDX
#define TRUE                       GSC_TRUE 
#define FALSE                      GSC_FALSE
#define NA                         GSC_NA

#define PedigreeID                 gsc_PedigreeID
#define NO_PEDIGREE                GSC_NO_PEDIGREE
#define GroupNum                   gsc_GroupNum
#define NO_GROUP                   GSC_NO_GROUP
#define EffectID                   gsc_EffectID
#define NO_EFFECTSET               GSC_NO_EFFECTSET
#define LabelID                    gsc_LabelID
#define NO_LABEL                   GSC_NO_LABEL
#define MapID                      gsc_MapID
#define NO_MAP                     GSC_NO_MAP
#define MultiIDSet                 gsc_MultiIDSet

#define TableSize                  gsc_TableSize
#define DecimalMatrix              gsc_DecimalMatrix
#define GenOptions                 gsc_GenOptions
#define BASIC_OPT                  GSC_BASIC_OPT
#define MarkerBlocks               gsc_MarkerBlocks
#define KnownGenome                gsc_KnownGenome
#define RecombinationMap           gsc_RecombinationMap
#define AlleleMatrix               gsc_AlleleMatrix
#define MarkerEffects              gsc_MarkerEffects
#define SimData                    gsc_SimData
#define GenoLocation               gsc_GenoLocation
#define INVALID_GENO_LOCATION      GSC_INVALID_GENO_LOCATION
#define IS_VALID_LOCATION          GSC_IS_VALID_LOCATION
#define BidirectionalIterator      gsc_BidirectionalIterator
#define RandomAccessIterator       gsc_RandomAccessIterator

#define create_new_label           gsc_create_new_label
#define change_label_default       gsc_change_label_default
#define change_label_to            gsc_change_label_to
#define change_label_by_amount     gsc_change_label_by_amount
#define change_label_to_values     gsc_change_label_to_values
#define change_names_to_values     gsc_change_names_to_values
#define change_allele_symbol       gsc_change_allele_symbol
#define change_eff_set_centres_to_values        gsc_change_eff_set_centres_to_values
#define change_eff_set_centre_of_markers        gsc_change_eff_set_centre_of_markers
#define change_eff_set_centre_of_allele_count   gsc_change_eff_set_centre_of_allele_count

#define create_empty_simdata       gsc_create_empty_simdata
#define clear_simdata              gsc_clear_simdata
#define load_mapfile               gsc_load_mapfile
#define load_effectfile            gsc_load_effectfile
#define load_genotypefile          gsc_load_genotypefile
#define load_data_files            gsc_load_data_files
#define DETECT_FILE_FORMAT         GSC_DETECT_FILE_FORMAT
#define FileFormatSpec             gsc_FileFormatSpec
#define define_matrix_format_details  gsc_define_matrix_format_details

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
#define get_index_of_label         gsc_get_index_of_label
#define get_index_of_eff_set       gsc_get_index_of_eff_set
#define get_index_of_map           gsc_get_index_of_map
#define get_index_of_genetic_marker gsc_get_index_of_genetic_marker
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
#define make_crosses_from_file            gsc_make_crosses_from_file
#define make_double_crosses_from_file     gsc_make_double_crosses_from_file

#define split_by_bv                gsc_split_by_bv
#define calculate_bvs              gsc_calculate_bvs
#define calculate_allele_counts    gsc_calculate_allele_counts
#define create_evenlength_blocks_each_chr gsc_create_evenlength_blocks_each_chr
#define load_blocks                       gsc_load_blocks
#define calculate_local_bvs               gsc_calculate_local_bvs
#define calculate_optimal_haplotype          gsc_calculate_optimal_haplotype
#define calculate_optimal_possible_haplotype gsc_calculate_optimal_possible_haplotype
#define calculate_optimal_bv                 gsc_calculate_optimal_bv
#define calculate_optimal_possible_bv        gsc_calculate_optimal_possible_bv
#define calculate_minimal_bv                 gsc_calculate_minimal_bv

#define delete_group               gsc_delete_group
#define delete_label               gsc_delete_label
#define delete_recombination_map   gsc_delete_recombination_map
#define delete_eff_set             gsc_delete_eff_set
#define delete_dmatrix             gsc_delete_dmatrix
#define delete_simdata             gsc_delete_simdata
#define delete_markerblocks        gsc_delete_markerblocks
#define delete_bidirectional_iter  gsc_delete_bidirectional_iter
#define delete_randomaccess_iter   gsc_delete_randomaccess_iter

#define save_markerblocks          gsc_save_markerblocks
#define save_genotypes             gsc_save_genotypes
#define save_allele_counts         gsc_save_allele_counts
#define save_pedigrees             gsc_save_pedigrees
#define save_bvs                   gsc_save_bvs
#define save_local_bvs             gsc_save_local_bvs

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
 * It contains pointers to a gsc_KnownGenome (storing the loaded markers and any
 * recombination maps), a
 * gsc_MarkerEffects (storing the loaded allele effects), and an gsc_AlleleMatrix
 * (storing metadata and genotypes of founders and simulated offspring).
 *
 * Other structs in this group (gsc_TableSize, gsc_MarkerBlocks, gsc_GenOptions) represented
 * specific types of data and are used as parameters and return value of certain
 * functions.
 *
 * @{
 */
 
 
/** Macro to create a stretchy buffer of any type and some length
 *
 * After this macro is run, a buffer with the requested type and capacity
 * will exist in the scope under the requested name. 
 *
 * This macro will also create two helper variables. Their names will be 
 * generated based on the name of the buffer:
 * - {name}cap will contain the current capacity of the buffer.
 * - {name}stack will be a stack array of size CONTIG_WIDTH. If the buffer 
 * length is not greater than CONTIG_WIDTH, then {name}stack will point to
 * the same array as the buffer.
 *
 * Use this buffer only within one local scope. These functions won't work
 * for buffers that have escaped their scope and so left their helper variables 
 * behind.
 *
 * The buffer will be allocated on the stack if its length will not make it 
 * exceed CONTIG_WIDTH*sizeof(int) in size in bytes. Otherwise, it will be allocated
 * on the heap. For safety, you 
 * should call @a GSC_DELETE_BUFFER on the buffer once you have finished using it,
 * even if you believe it is small enough to have been allocated on the stack.
 *
 * @see GSC_DELETE_BUFFER
 * @see GSC_STRETCH_BUFFER
 *
 * @param n name for the buffer.
 * @param type type of each entry in the buffer (eg int).
 * @param length number of entries the buffer should be able to hold.
 */
#define GSC_CREATE_BUFFER(n,type,length) \
type n##stack[sizeof(int)*CONTIG_WIDTH/sizeof(type)]; unsigned int n##cap = length; \
type* n = (n##cap >= sizeof(n##stack)/sizeof(type)) ? gsc_malloc_wrap(sizeof(type)*n##cap,GSC_TRUE) : n##stack;

/** For debugging purposes.
 *
 * @see GSC_CREATE_BUFFER
 * @see GSC_STRETCH_BUFFER
 * @see GSC_DELETE_BUFFER
 */
#define GSC_BUFFER_ISHEAP(n) n##cap >= sizeof(n##stack)/sizeof(n##stack[0])

/** Macro to convert a stretchy buffer to a solid heap vector.
 *
 *  @see GSC_DELETE_BUFFER
 *
 * The buffer named {n}, and its assistant variable {n}cap and {n}stack, must exist
 * in the current scope. They would be created by @a GSC_CREATE_BUFFER
 *
 * This is an alternative to GSC_DELETE_BUFFER, if you want to keep the results.
 *
 * @param n name of the buffer.
 * @param as name of the finalised buffer.
 * @param nentries number of entries to copy, if less than buffer capacity
 */
#define GSC_FINALISE_BUFFER(n,as,nentries) do { if (n##cap >= sizeof(n##stack)/sizeof(n##stack[0])) { as = n; } else \
{ unsigned int len = nentries > n##cap ? n##cap : nentries; as = gsc_malloc_wrap(sizeof(n##stack[0])*len,GSC_TRUE); memcpy(as,n,sizeof(n##stack[0])*len); } } while (0)

/** Macro to delete a stretchy buffer
 *
 * @see GSC_CREATE_BUFFER
 * @see GSC_STRETCH_BUFFER
 * @see GSC_FINALISE_BUFFER
 *
 * The buffer named {n}, and its assistant variable {n}cap and {n}stack, must exist 
 * in the current scope. They would be created by @a GSC_CREATE_BUFFER
 *
 * @param n name of the buffer.
 */
#define GSC_DELETE_BUFFER(n) do { if (n##cap >= sizeof(n##stack)/sizeof(n##stack[0])) { GSC_FREE(n); } \
n = NULL; n##cap = 0; } while (0)

/** Macro to expand the capacity of a stretchy buffer
 *
 * @see GSC_CREATE_BUFFER
 * @see GSC_DELETE_BUFFER
 *
 * The buffer named {n}, and its assistant variables {n}cap and {n}stack, must exist 
 * in the current scope. They would be created by @a GSC_CREATE_BUFFER
 *
 * After this macro executes, the buffer named {n} will have the capacity 
 * to hold {n}cap entries. Unless memory allocation failed, {n}cap will be 
 * greater than or equal to the requested new length. Check the value of 
 * {n}cap to check that resizing succeeded. 
 *
 * @param n name of the buffer.
 * @param newlen after execution, the buffer should be able to hold this 
 * many entries, unless memory allocation failed (can be checked with
 * n{cap} >= newlen )
 */
#define GSC_STRETCH_BUFFER(n,newlen) do {    \
    if (newlen < n##cap) { } \
    else if (n##cap >= sizeof(n##stack)/sizeof(n##stack[0])) { \
        void* tmp = gsc_malloc_wrap(sizeof(n##stack[0])*newlen,GSC_FALSE); \
        if (tmp != NULL) { \
        memcpy(tmp,n,sizeof(n[0])*n##cap); \
        GSC_FREE(n); n = tmp; n##cap = newlen; }} \
    else if (newlen >= sizeof(n##stack)/sizeof(n##stack[0])) { \
        n = gsc_malloc_wrap(sizeof(n##stack[0])*newlen,GSC_FALSE); \
        if (n != NULL) { \
        memcpy(n,n##stack,sizeof(n##stack[0])*n##cap); n##cap = newlen; }} \
    else if (newlen < CONTIG_WIDTH) { n##cap = newlen; } \
} while (0)


/* A simple struct used for returning the dimensions of a matrix or table.
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
    GSC_GENOLEN_T num_blocks;

    /** Pointer to a heap array of length num_blocks
      * containing the number of markers that make up each block */
    GSC_GENOLEN_T* num_markers_in_block;

    /** Pointer to a heap array of length num_blocks, each
     * entry in which is a pointer to a heap array with length corresponding to
     * the value of the corresponding entry in num_markers_in_block whose values
     * are the indexes in the gsc_SimData of the markers that make up that block. */
    GSC_GENOLEN_T** markers_in_block;
} gsc_MarkerBlocks;

/** A row-major heap matrix that contains floating point numbers.
 *
 * Rows make the first index of the matrix and columns the second.
 *
 * @shortnamed{DecimalMatrix}
 */
typedef struct {
    double** matrix; /**< The actual matrix and contents */
    unsigned int dim1;        /**< Number of rows in the matrix */
    unsigned int dim2;        /**< number of columns in the matrix */
} gsc_DecimalMatrix;


/** A type representing a program-lifetime-unique identifier for a genotype,
 *  to be used in tracking pedigree.
 *
 * @shortnamed{PedigreeID}
 */
typedef struct {
    GSC_ID_T id;
} gsc_PedigreeID;
/** Empty/null value for pedigree fields.
 *
 * @shortnamed{NO_PEDIGREE} */
#define GSC_NO_PEDIGREE (gsc_PedigreeID){.id=GSC_NA_ID}


/** A type representing the identifier of a group of genotypes
 *
 * @shortnamed{GroupNum}
 */
typedef struct {
    GSC_ID_T num;
} gsc_GroupNum;
/** Empty/null value for group allocations.
 *
 * @shortnamed{NO_GROUP} */
#define GSC_NO_GROUP (gsc_GroupNum){.num=GSC_NA_ID}

/** A type representing a particular loaded set of marker effects
 *
 * @shortnamed{EffectID}
 */
typedef struct {
    GSC_ID_T id;
} gsc_EffectID;
/** Empty/null value for effect set identifiers.
 *
 * @shortnamed{NO_EFFECTSET} */
#define GSC_NO_EFFECTSET (gsc_EffectID){.id=GSC_NA_ID}

/** A type representing a particular custom label
 *
 * @shortnamed{LabelID}
 */
typedef struct {
    GSC_ID_T id;
} gsc_LabelID;
/** Empty/null value for custom label identifiers.
 *
 * @shortnamed{NO_LABEL} */
#define GSC_NO_LABEL (gsc_LabelID){.id=GSC_NA_ID}

/** A type representing a particular loaded recombination map
 *
 * @shortnamed{MapID}
 */
typedef struct {
    GSC_ID_T id;
} gsc_MapID;
/** Empty/null value for recombination map identifiers.
 *
 * @shortnamed{NO_MAP} */
#define GSC_NO_MAP (gsc_MapID){.id=GSC_NA_ID}

/** Simple crate that stores a GroupNum, a MapID, and an EffectID.
 *
 * @shortnamed{MultiIDSet} */
struct gsc_MultiIDSet {
    gsc_GroupNum group;
    gsc_MapID map;
    gsc_EffectID effSet;
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
    _Bool will_name_offspring; /**< A boolean: whether generated offspring should be given names. */
    const char* offspring_name_prefix; /**< If `will_name_offspring` is true, generated
                           * offspring are named with the concatenation {offspring_name_prefix}{index}. */

    GSC_GLOBALX_T family_size; /**< The number of offspring to produce from each cross.*/

    _Bool will_track_pedigree; /**< A boolean: whether to track parentage of generated offspring.*/
    _Bool will_allocate_ids; /**< A boolean: whether to allocate generated offspring session-
                            * unique IDs. IDs are used for pedigree tracking. The
                            * offspring of an anonymous individual (one without an ID)
                            * cannot identify that individual as their parent. */

    const char* filename_prefix; /**< A string used in save-as-you-go file names. */
    _Bool will_save_pedigree_to_file; /**< A boolean. If true, the full/recursive
                            * pedigrees of every offspring generated in the cross
                            * are saved to "{filename_prefix}-pedigree.txt", even
                            * if `will_save_to_simdata` is false.
                            * Pedigrees are saved in the format of gsc_save_full_pedigree()*/
    gsc_EffectID will_save_bvs_to_file; /**< If equal to NO_EFFECTSET, no bvs are calculated or saved.
                            * Otherwise, for each offspring in the cross,
                            * the breeding values according
                            * to the marker effect set with this gsc_EffectID
                            * are saved to "{filename_prefix}-bv.txt", even
                            * if `will_save_to_simdata` is false.
                            * BVs are saved in the format of gsc_save_bvs() */
    _Bool will_save_alleles_to_file; /**< A boolean. If true, the set of alleles
                            * of every offspring generated in the cross
                            * are saved to "{filename_prefix}-genotype.txt", even
                            * if `will_save_to_simdata` is false.
                            * Genotypes are saved in the format of gsc_save_allele_matrix()*/
    /* _Bool will_save_recombinations_to_file; **< A boolean. If true, will log recombination
                            * events to "{filename_prefix}-recomb.txt". This information
                            * is not accessible by another saver function - it is only accessible
                            * if logged by the function generating a new set of offspring.
                            * The columns of the log file are:
                            */
    _Bool will_save_to_simdata; /**< A boolean. If true, the generated offspring exist
                            * in the gsc_SimData struct after the function executes.
                            * If false, they are discarded after creation. */
} gsc_GenOptions;

extern const gsc_GenOptions GSC_BASIC_OPT;

/** Parameters for simulating meiosis on a linkage group whose markers are stored contiguously
 *  in the simulation.
 *
 * @shortnamed{SimpleLinkageGroup}
 */
typedef struct {
    double expected_n_crossovers; /**< Expected value of the Poisson distribution from which
                                   * the number of crossovers in this linkage group will be
                                   * drawn when simulating meiosis. Probably corresponds to
                                   * the length of the chromosome/linkage group in Morgans. */
    GSC_GENOLEN_T n_markers; /**< The number of markers in this chromosome/linkage group. Their indexes in
                       * the simulation's corresponding @a gsc_KnownGenome and @a gsc_AlleleMatrix are
                       * indexes [ @a first_marker_index ] to [ @a first_marker_index + @a n_markers - 1 ] */
    GSC_GENOLEN_T first_marker_index; /**< The index of the first marker in this chromosome/linkage group
                          * in the simulation's corresponding @a gsc_KnownGenome and @a gsc_AlleleMatrix. */
    double* dists; /**< Array with @a n_markers entries, containing at position i
               * the distance in centimorgans along the linkage group of the i-th marker,
               * divided by the length of the linkage group in centimorgans.
               * Effectively, positions of markers when the length of the linkage
               * group is normalised to 1. */
} gsc_SimpleLinkageGroup;

/** Parameters for simulating meiosis on a linkage group whose markers are re-ordered
 *  compared to the first recombination map.
 *
 *  Could represent chromosomal inversion, or different allocations of markers to linkage groups.
 *
 * @shortnamed{ReorderedLinkageGroup}
 */
typedef struct {
    double expected_n_crossovers; /**< Expected value of the Poisson distribution from which
                                   * the number of crossovers in this linkage group will be
                                   * drawn when simulating meiosis. Probably corresponds to
                                   * the length of the chromosome/linkage group in Morgans. */
    GSC_GENOLEN_T n_markers; /**< The number of markers in this chromosome/linkage group. */
    GSC_GENOLEN_T* marker_indexes; /**< Array with @a n_markers entries. Each entry is the index of a marker
                          * in the simulation's corresponding @a gsc_KnownGenome and @a gsc_AlleleMatrix.
                          * This set may represent a different linkage group or different ordering
                          * of markers to the default gsc_KnownGenome set. Positions in @a dists correspond
                          * to the markers at the same index in this vector. */
    double* dists; /**< Array with @a n_markers entries, containing at position i
               * the distance in centimorgans along the linkage group of marker i
               * divided by the length of the linkage group in centimorgans.
               * Effectively, positions of markers when the length of the linkage
               * group is normalised to 1. */
} gsc_ReorderedLinkageGroup;

/** A generic store for a linkage group, used to simulate meiosis on
 *  a certain subset of markers.
 *
 * @shortnamed{LinkageGroup}
 */
typedef struct {
    enum gsc_LinkageGroupType {
        GSC_LINKAGEGROUP_SIMPLE, /**< @see gsc_SimpleLinkageGroup */
        GSC_LINKAGEGROUP_REORDER /**< @see gsc_ReorderedLinkageGroup */
    } type;
    union {
        gsc_SimpleLinkageGroup simple;
        gsc_ReorderedLinkageGroup reorder;
    } map;
} gsc_LinkageGroup;

/** A type that stores linkage groups and crossover probabilities for simulating meiosis.
 *
 * @shortnamed{RecombinationMap}
 */
typedef struct {
    unsigned int n_chr; /**< The number of chromosomes/linkage groups represented in the map. **/
    char** chr_names; /**< An identifying code for each chromosome/linkage group in the map. Optional. */
    gsc_LinkageGroup* chrs; /**< Vector of @a n_chr recombination maps, one for each chromosome/linkage group
              * in this recombination map. */

} gsc_RecombinationMap;

/** A type that stores the genome structure used in simulation.
 *
 * @shortnamed{KnownGenome}
 */
typedef struct {
    GSC_GENOLEN_T n_markers; /**< The total number of markers.**/
    char** marker_names; /**< A vector of @a n_markers strings containing the names of markers, ordered
                          * according to their index in an AlleleMatrix. */
    char*** names_alphabetical; /**< A vector of @a n_markers pointers to names in @a marker_names, ordered
                          * in alphabetical order of the names. For speeding up functions involving loading
                          * new genetic maps, or observing genotypes at particular markers. */

    GSC_ID_T n_maps; /**< The number of recombination maps currently stored. */
    gsc_MapID* map_ids; /**< A vector of @a n_maps identifiers for each of the recombination maps
                         * currently stored. These IDs correspond to the map in the corresponding
                         * position in the vector @a maps. */
    gsc_RecombinationMap* maps; /**< A vector of @a n_maps recombination maps, to use for simulating meiosis. */
} gsc_KnownGenome;


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
     * eg TT, TA. Use `alleles[genotype index][marker index * 2]` to get the
     * first allele and `alleles[genotype index][marker index * 2 + 1]` to
     * get the second. If CONTIG_WIDTH lines are saved here, another
     * gsc_AlleleMatrix is added to the linked list when there's a need to save more.*/
    char* alleles[CONTIG_WIDTH];

    GSC_LOCALX_T n_genotypes; /**< Number of genotypes currently loaded in this matrix.*/
    GSC_GENOLEN_T n_markers; /**< Number of markers across which genotypes are tracked. This has
                    * redundancy with gsc_SimData and other members of its linked list
                    * but it's good to know how big your own `alleles` array is.*/

    char* names[CONTIG_WIDTH]; /**< Array of dynamically allocated strings
                    * containing the names of the lines/genotypes in this matrix.
                    * Guaranteed to be NULL if they do not have names. */
    gsc_PedigreeID ids[CONTIG_WIDTH]; /**< Unique ID for each genotype. */
    gsc_PedigreeID pedigrees[2][CONTIG_WIDTH]; /**< Two lists of integer IDs of the
                    * parents of this genotype (if tracked), or 0 if we don't know/care.*/
    gsc_GroupNum groups[CONTIG_WIDTH]; /**< Group allocation of each genotype. */

    GSC_ID_T n_labels; /**< Number of custom labels currently available to this gsc_AlleleMatrix. This has
                    * redundancy with gsc_SimData and other members of its linked list
                    * but it's good to know how big your own `labels` array is.*/
    int** labels; /**< Pointer to list of labels. Size of first dimension is n_labels,
                           * of second dimension is arrays of labels of length CONTIG_WIDTH*/

    gsc_AlleleMatrix* next; /**< Pointer to the next gsc_AlleleMatrix in the linked list,
                         * or NULL if this entry is the last. */
};

/** A type that stores the information needed to calculate breeding values from 
 * alleles at markers.
 *
 * The equation to calculate breeding values for a genotype g across markers m
 * using the marker effects in this struct is:
 *
 * \f$\sum_{m}\left[\sum_{a_m}\left(c_{g,m,a_m} \cdot eff_{m,a_m}\right) - centre_{m}\right]\f$
 *
 * where this struct stores effects for a set of alleles \f$a_m\f$ for each marker m, 
 * \f$c_[g,m,a_m}$\f$ is the number of copies of allele \f$a_m\f$ that genotype g has
 * at marker m, \f$eff_{m,a_m}\f$ is the effect of that allele at that marker (stored in .eff), 
 * and \f$centre_{m}\f$ is a centring value applied per marker (stored in .centre).
 *
 * The formula for calculating the breeding value of a haplotype is the same, except 
 * \f$c_[g,m,a_m}$\f$ will have possible values 0 or 1 instead of 0, 1 or 2. 
 *
 * @shortnamed{MarkerEffects}
 */
typedef struct {
    GSC_GENOLEN_T n_markers; /**< Number of markers across which genotypes are tracked. 
                              * This number has redundancy with gsc_SimData. */
    double* centre; /**< Vector of length @a n_markers, containing a value for each marker 
        * which represents the value to subtract from the breeding value
        * contribution of the alleles at that marker. 
        * Markers are ordered according to the SimData's KnownGenome.marker_names.
        * If this value is NULL, no values are to be subtracted
        * from the products of allele count and allele effect.
        */

    GSC_GENOLEN_T* cumn_alleles; /**< A vector of length @a n_markers holding the 
        * cumulative number of alleles that have effects on breeding values at each marker.
        * Markers are ordered according to the SimData's KnownGenome.marker_names.
        * At position i, contains the number of alleles with effects at the ith marker
        * plus the number of alleles at previous markers. So at index 0, it will have a nonzero value. */
    char* allele; /**< A vector holding the symbol/character representing each allele at each marker.
        * Index into with (MarkerEffects.n_markers * markerix + alleleix). */
    double* eff; /**< A vector holding the effect on breeding value of each allele at each marker.
        * Index into with (MarkerEffects.n_markers * markerix + alleleix). */
} gsc_MarkerEffects;

/** Composite type that is used to run crossing simulations.
 *
 * The core of this type is a list of markers. These are used to index the rows
 * of the allele matrix and the position map, and the columns of the effect matrix.
 *
 * @shortnamed{SimData}
 */
typedef struct {
    GSC_ID_T n_labels; /**< The number of custom labels in the simulation.*/
    gsc_LabelID* label_ids; /**< The identifier number of each label in the simulation, in order
                     * of their lookup index. */
    int* label_defaults; /**< Array containing the default (birth) value of each
                          * custom label. */

    gsc_KnownGenome genome; /**< A gsc_KnownGenome, which stores the information of
                             * known markers and linkage groups, as well as one or more
                             * recombination maps for use in simulating meiosis.*/
    gsc_AlleleMatrix* m; /**< Pointer to an gsc_AlleleMatrix, which stores data and
                      * metadata of founders and simulated offspring. The
                      * gsc_AlleleMatrix is start of a linked list if there are
                      * many genotypes. */

    GSC_ID_T n_eff_sets; /**< The number of sets of allele effects in the simulation **/
    gsc_EffectID* eff_set_ids; /**< The identifier number of each set of allele effects in the simulation,
                     * ordered by their lookup index. */
    gsc_MarkerEffects* e; /**< Array of n_eff_sets gsc_MarkerEffects, optional for the use 
                     * of the simulation. Used for calculating breeding values for genotypes. */

    //CRANDOMGENERATOR /**< Random number generator working memory. */
    gsc_PedigreeID current_id; /**< Highest SimData-unique ID that has been generated
                              * so far. Used to track which IDs have already been
                              * given out.*/
    GSC_ID_T n_groups; /**< Number of groups currently existing in simulation. It is
                        * guaranteed to never be less than the number of groups in simulation
                        * even if not perfectly accurate. */
} gsc_SimData;


/** Represent possible states of the cursor of a @a gsc_TableFileReader */
enum gsc_TableFileCurrentStatus {
    GSC_TABLEFILE_NEWLINE,
    GSC_TABLEFILE_COLUMNGAP,
    GSC_TABLEFILE_CONTENTS,
    GSC_TABLEFILE_ERROR_EOF,
    GSC_TABLEFILE_ERROR_EOBUF
};

/** Stream reader for files of some tabular format.
 *
 * Expected usage is the use of @a gsc_tablefilereader_create to
 * create an instance of this struct, repeated calls to
 * @a gsc_tablefilereader_get_next to parse the file,
 * and @a gsc_tablefilereader_close to close the file after use. The
 * internals of this struct are not intended to be manually manipulated.
 *
 * A table file consists of one or more lines (separated by \n, \r,
 * or \r\n) and one or more columns (separated by any combination of \t,
 * ' ', or ','). Line and column
 * separators do not need to be consistent across different line/column
 * gaps in the file.
 *
 * In the current implementation, the maximum length of a cell that can be
 * read by this file reader is 8192 characters.
 *
 * Multi-byte characters are not handled correctly in this current implementation.
 *
 * @see gsc_tablefilereader_create
 * @see gsc_tablefilereader_close
 * @see gsc_tablefilereader_get_next
 */
typedef struct {
    FILE* fp; /**< File being read. @a gsc_tablefilereader_create ensures this is a file successfully opened for reading. */

    char buf[8192]; /**< A window of characters from the file, loaded into memory for current processing. In current implementation,
                     * has a maximum capacity of 8192 characters. */
    int buf_fill; /**< Number of characters from the file that are currently loaded in @a buf. */
    int cursor; /**< Index in @a buf of the first character that the file reader has not yet parsed. */
} gsc_TableFileReader;

/** Represent a cell read by a @a gsc_TableFileReader */
typedef struct {
    _Bool isCellShallow; /**< is the string in 'cell' a shallow copy or deep copy? */
    char* cell; /**< deep copy of the cell contents, or NULL */
    unsigned int cell_len; /**< length of cell contents (because a shallow copy may not be null-terminated) */
    int predCol; /**< since last read, how many column gaps have there been? */
    int predNewline; /**< since last read, how many newlines have there been? */
    _Bool eof; /**< are we (this cell) at end of file */
} gsc_TableFileCell;


gsc_TableFileReader gsc_tablefilereader_create(const char* filename);
void                gsc_tablefilereader_close(gsc_TableFileReader* tbl);
void gsc_helper_tablefilereader_refill_buffer(gsc_TableFileReader* tbl);
enum gsc_TableFileCurrentStatus gsc_helper_tablefilereader_classify_char(gsc_TableFileReader* tbl);

void gsc_tablefilecell_deep_copy(gsc_TableFileCell* c);
gsc_TableFileCell gsc_tablefilereader_get_next_cell(gsc_TableFileReader* tbl);


/** Represent possible representations of alleles at a marker in a genotype file */
enum gsc_GenotypeFileCellStyle {
    GSC_GENOTYPECELLSTYLE_PAIR,
    GSC_GENOTYPECELLSTYLE_COUNT,
    GSC_GENOTYPECELLSTYLE_ENCODED,
    //GSC_GENOTYPECELLSTYLE_SPACEDPAIR, // to come later maybe
    GSC_GENOTYPECELLSTYLE_SLASHPAIR,
    GSC_GENOTYPECELLSTYLE_UNKNOWN
};


/** Unprocessed data for one marker (linkage group and position) loaded from a map file. */
struct gsc_MapfileUnit {
    char* name;
    char* chr;
    double pos;
};

/** Unprocessed data for one marker effect loaded from an effect file. */
struct gsc_EffectfileUnit {
    GSC_GENOLEN_T markerix;
    char allele;
    double centre;
    double eff;
};

/** Enumerate types of genotype files that the simulation knows how to load.
 *
 * The format of the file cannot be automatically in all cases. This type exists so that
 * users can specify the format of input files.
 *
 * @shortnamed{GenotypeFileType}
 *
 * The format of the file is decoded separately from the format of the alleles
 * of each genotype at each marker. In the templates of each format, "[alleles]"
 * can represent:
 * - a pair of ASCII characters (in which case the two characters are interpreted as
 * the two alleles, with their ordering representing their phase),
 * - a single character from the standard IUPAC nucleotide encoding (in which case the
 * character is decoded to represent the alleles observed at this marker. The phase of
 * the alleles at this marker is randomly chosen if the genotype is heterozygous at that
 * marker), or
 * - a single digit from {0,1,2} (in which case the digit represents the number of copies of the
 * major allele. The phase of the alleles at this marker if the digit is 1 is randomly chosen.)
 *
 * IUPAC nucleotide encoding: Code => Alleles key:
 * A => AA    ; C => CC    ; G => GG    ; T => TT   ;
 * R => AG    ; Y => CT    ; S => CG    ; W => AT   ; K => GT   ; M => AC
 *
 * The format of the "[alleles]" cells in a file can be automatically determined
 * out of these options.
 *
 */
enum gsc_GenotypeFileType {
    GSC_GENOTYPEFILE_UNKNOWN,
    /** Either a marker-by-line matrix, where each marker is a row, or a
     *  line-by-marker matrix, where each marker is a column.
     *
     * The other axis represents lines/organisms/founders.
     *
     * @section Template (marker-as-row form)
     *
     * [corner] [line] [line] [line] ... [line]
     *
     * [marker] [alleles] [alleles] [alleles] ... [alleles]
     *
     * [marker] [alleles] [alleles] [alleles] ... [alleles]
     *
     * ...
     *
     * @section Template (marker-as-column form)
     *
     * [corner] [marker] [marker] [marker] ... [marker]
     *
     * [line] [alleles] [alleles] [alleles] ... [alleles]
     *
     * [line] [alleles] [alleles] [alleles] ... [alleles]
     *
     * ...
     *
     * @section Format Details
     *
     * The corner cell may or may not be filled. Its value is ignored.
     *
     * Any combination of spaces or tabs between non-space/non-tab characters
     * is interpreted as a column separator. Length and order of spaces and tabs
     * do not need to be consistent between column separators in the file.
     *
     * Any one or two consecutive characters from {'\n', '\r'}, in any order,
     * will be interpreted as a single line break.
     *
     * The default, when no map is present in simulation, is to assume markers are
     * rows in this file. However, if
     * any of the column headers of a matrix file are names of markers being
     * tracked by the simulation, then that file is interpreted as having markers
     * as columns.
    */
    GSC_GENOTYPEFILE_MATRIX,
    GSC_GENOTYPEFILE_BED,
    GSC_GENOTYPEFILE_PED,
    GSC_GENOTYPEFILE_VCF,
    // GSC_GENOTYPEFILE_FLATFILE
};

/** Variants in the format of a genotype matrix file. */
struct gsc_GenotypeFile_MatrixFormat {
    GSC_LOGICVAL has_header; /** < Boolean: Is the first row of the file a header row? (Note: genotype matrix files must have
                    * row headers, so this setting only applies to the presence/absense of column headers).
                    * Set to GSC_TRUE or GSC_FALSE if this detail of the formatting is known, or a negative number like GSC_UNINIT if
                    * the detail is not known and file loaders should auto-detect headers. */
    GSC_LOGICVAL markers_as_rows; /** < Boolean: Are genetic markers the rows of the matrix (GSC_TRUE) or the columns of the matrix (GSC_FALSE)?
                    * Set to GSC_TRUE or GSC_FALSE if this detail of the formatting is known, or a negative number like GSC_UNINIT if
                    * the detail is not known and file loaders should auto-detect matrix orientation. */
    enum gsc_GenotypeFileCellStyle cell_style; /** < How are the alleles for a genotype and marker encoded, in the body of the matrix?
                    * Set to GSC_GENOTYPECELLSTYLE_UNKNOWN if this detail is not known and file loaders should auto-detect cell format. */
};

/** File format specifier for the genotype input file.
 *
 * If @a filetype is set to @a GSC_GENOTYPEFILE_UNKNOWN, then
 * no specification is defined and the union @a spec will not be
 * accessed. Otherwise, the values should correspond. If @a filetype is not @a GSC_GENOTYPEFILE_UNKNOWN, then
 * the struct in @a spec should be of the matching file type
 * (eg a struct gsc_GenotypeFile_MatrixFormat if the file type is GSC_GENOTYPEFILE_MATRIX).
**/
typedef struct {
    enum gsc_GenotypeFileType filetype;
    union {
        struct gsc_GenotypeFile_MatrixFormat matrix;
        //struct gsc_GenotypeFile_plinkFormat plink;
        //struct gsc_GenotypeFile_vcfFormat vcf;
    } spec;
} gsc_FileFormatSpec;

/** File format specifier to instruct genomicSimulation loaders to auto-detect all details of the file format. */
#define GSC_DETECT_FILE_FORMAT ((gsc_FileFormatSpec){.filetype=GSC_GENOTYPEFILE_UNKNOWN})

gsc_FileFormatSpec gsc_define_matrix_format_details(const GSC_LOGICVAL has_header, 
        const GSC_LOGICVAL markers_as_rows, const enum gsc_GenotypeFileCellStyle cell_style);

/** @} */
/** @defgroup supporters Utils/Supporting Functions
 *
 * @{
 */
struct gsc_TableSize gsc_get_file_dimensions(const char* filename, const char sep);
//int gsc_get_from_ordered_uint_list(const unsigned int target, const unsigned int listLen, const unsigned int* list);
GSC_LOCALX_T gsc_get_from_ordered_pedigree_list(const gsc_PedigreeID target, const GSC_LOCALX_T listLen, const gsc_PedigreeID* list);
unsigned int gsc_get_from_unordered_str_list(const char* target, const unsigned int listLen, const char** list);
unsigned int gsc_get_from_ordered_str_list(const char* target, const unsigned int listLen, const char** list);


void gsc_shuffle_up_to( void* sequence, const unsigned int item_size, const unsigned int total_n, const unsigned int n_to_shuffle);
GSC_GLOBALX_T gsc_randomdraw_replacementrules(gsc_SimData* d, GSC_GLOBALX_T max, GSC_GLOBALX_T cap, 
                                              GSC_GLOBALX_T* member_uses, GSC_GLOBALX_T noCollision);

gsc_LabelID gsc_create_new_label(gsc_SimData* d, const int setTo);
void gsc_change_label_default(gsc_SimData* d, const gsc_LabelID whichLabel, const int newDefault);
void gsc_change_label_to(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel, const int setTo);
void gsc_change_label_by_amount(gsc_SimData* d, const gsc_GroupNum whichGroup, const gsc_LabelID whichLabel,
                                const int byValue);
void gsc_change_label_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const GSC_GLOBALX_T startIndex, 
                                const gsc_LabelID whichLabel, const unsigned int n_values, const int* values);
void gsc_change_names_to_values(gsc_SimData* d, const gsc_GroupNum whichGroup, const GSC_GLOBALX_T startIndex, 
                                const unsigned int n_values, const char** values);
void gsc_change_allele_symbol(gsc_SimData* d, const char* which_marker, const char from, const char to);

_Bool gsc_change_eff_set_centres_to_values(gsc_SimData* d,
                                           const gsc_EffectID effset,
                                           const GSC_GENOLEN_T n_values,
                                           const double* values);
GSC_GENOLEN_T gsc_change_eff_set_centre_of_markers(gsc_SimData* d,
                                                   const gsc_EffectID effset,
                                                   const GSC_GENOLEN_T n_markers,
                                                   const char** marker_names,
                                                   const double* centres);
GSC_GENOLEN_T gsc_change_eff_set_centre_of_allele_count(gsc_SimData* d,
                                                        const gsc_EffectID effset,
                                                        const GSC_GENOLEN_T n_markers,
                                                        const char** marker_names,
                                                        const double* centres,
                                                        const char allele,
                                                        const _Bool reset_centres);

//static void gsc_set_names(gsc_AlleleMatrix* a, const char* prefix, const int suffix, const GSC_LOCALX_T from_index);
int gsc_get_integer_digits(const int i);
GSC_ID_T gsc_get_index_of_label( const gsc_SimData* d, const gsc_LabelID label );
GSC_ID_T gsc_get_index_of_eff_set( const gsc_SimData* d, const gsc_EffectID eff_set_id );
GSC_ID_T gsc_get_index_of_map( const gsc_SimData* d, const gsc_MapID map );
_Bool gsc_get_index_of_genetic_marker(const char* target, gsc_KnownGenome g, GSC_GENOLEN_T* out);

gsc_LabelID gsc_get_new_label_id( const gsc_SimData* d );
gsc_EffectID gsc_get_new_eff_set_id( const gsc_SimData* d );
gsc_MapID gsc_get_new_map_id( const gsc_SimData* d);
gsc_GroupNum gsc_get_next_free_group_num( const unsigned int n_existing_groups, const gsc_GroupNum* existing_groups, 
                                          unsigned int* cursor, gsc_GroupNum previous);
gsc_GroupNum gsc_get_new_group_num( gsc_SimData* d );
void gsc_get_n_new_group_nums( gsc_SimData* d, const unsigned int n, gsc_GroupNum* result);
void gsc_condense_allele_matrix( gsc_SimData* d);
//static void* gsc_malloc_wrap(const unsigned int size, char exitonfail);

//static int gsc_helper_simdata_pos_compare(const void *pp0, const void *pp1);
//static int gsc_helper_descending_double_comparer(const void* pp0, const void* pp1);
//static int gsc_helper_ascending_double_comparer(const void* pp0, const void* pp1);
// static int gsc_helper_alphabetical_str_comparer(const void* p0, const void* p1);
// static int gsc_helper_indirect_alphabetical_str_comparer(const void* p0, const void* p1);
// static int gsc_helper_mapfileunit_ascending_chr_comparer(const void* p0, const void* p1);
// static int gsc_helper_mapfileunit_ascending_d_comparer(const void* p0, const void* p1);
/**@}*/


/** @defgroup loaders Setup Functions
 *
 * For setup of the simulation (loading founders, genetic maps, and optionally allele effects).
 *
 * @{
 */
gsc_AlleleMatrix* gsc_create_empty_allelematrix(const GSC_GENOLEN_T n_markers, const GSC_ID_T n_labels, 
                                                const int* labelDefaults, const GSC_LOCALX_T n_genotypes);
gsc_SimData* gsc_create_empty_simdata();
void gsc_clear_simdata(gsc_SimData* d);

gsc_GroupNum gsc_load_genotypefile(SimData* d, const char* filename, const gsc_FileFormatSpec format);
// static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_orientation(const SimData* d, 
        //const gsc_TableFileCell* cellqueue, const unsigned int firstrowlen, const unsigned int queuelen, 
        //struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog);
// static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_cellstyle(
        //const gsc_TableFileCell* cellqueue, const unsigned int firstrowlen, const unsigned int queuelen, 
        // struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog);
// static enum gsc_GenotypeFileCellStyle gsc_helper_genotype_matrix_identify_cell_style(char* cell);
// static struct gsc_GenotypeFile_MatrixFormat gsc_helper_genotypefile_matrix_detect_header(
    //const gsc_TableFileCell* cellqueue, const unsigned int firstrowlen, const unsigned int queuelen, 
    // struct gsc_GenotypeFile_MatrixFormat format, const char* filenameforlog);
// static int gsc_helper_genotypefile_matrix_detect_cornercell_presence(const unsigned int ncellsfirstrow, 
    // const unsigned int ncellssecondrow, const int secondrowheaderisempty);


// static struct gsc_MultiIDSet gsc_load_genotypefile_matrix(gsc_SimData* d, const char* filename, const char* mapfile)
// static void gsc_helper_parse_genofile(char* filename, enum gsc_GenotypeFileType type, 
        //GSC_GENOLEN_T* n_markers, char*** marker_names, GSC_GLOBALX_T* n_lines, char*** line_names, 
        //struct gsc_GenotypeFileCell** cells)

// static int gsc_helper_genotypefile_matrix_check_markers_are_rows(gsc_SimData* d, int hasheader, 
        //gsc_TableFileCell* firstrow, unsigned int firstrowlen, gsc_TableFileCell secondrowcellone)
//static void gsc_helper_genotypecell_to_allelematrix(GenoLocation loc, GSC_GENOLEN_T markerix, 
        //enum gsc_GenotypeFileCellStyle style, char* cell, gsc_SimData* forrng)

gsc_MapID gsc_load_mapfile(gsc_SimData* d, const char* filename);
// static GSC_GENOLEN_T gsc_helper_sort_markerlist(GSC_GENOLEN_T n_markers, struct gsc_MapfileUnit* markerlist) 
// static GSC_GENOLEN_T gsc_helper_parse_mapfile(const char* filename, struct gsc_MultiTypeUnit** out)
// static gsc_MapID gsc_helper_insert_recombmap_into_simdata(gsc_SimData* d, gsc_RecombinationMap map)
// static GSC_GENOLEN_T gsc_helper_str_markerlist_leftjoin(gsc_KnownGenome g, GSC_GENOLEN_T n_markers_in_list, 
        //struct gsc_MapfileUnit** markerlist)
gsc_MapID gsc_create_recombmap_from_markerlist(gsc_SimData* d, GSC_GENOLEN_T n_markers, struct gsc_MapfileUnit* markerlist);
gsc_MapID gsc_create_uniformspaced_recombmap(gsc_SimData* d, GSC_GENOLEN_T n_markers, char** markernames, 
                                             double expected_n_recombinations);
gsc_MapID gsc_create_unlinked_recombmap(gsc_SimData* d, GSC_GENOLEN_T n_markers, char** markernames);

gsc_EffectID gsc_load_effectfile(gsc_SimData* d, const char* filename);

struct gsc_MultiIDSet gsc_load_data_files(gsc_SimData* d, const char* data_file, const char* map_file, 
                                          const char* effect_file, const gsc_FileFormatSpec format);
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
    GSC_LOCALX_T localPos; /**< Index in the localAM where the genotype can be
                   * found (min value: 0. Max value: CONTIG_WIDTH-1). */
} gsc_GenoLocation;

/** Constant representing a nonexistent location in the simulation.
 *
 * @shortnamed{INVALID_GENO_LOCATION}
 */
#define GSC_INVALID_GENO_LOCATION (gsc_GenoLocation){.localAM=0,.localPos=GSC_NA_LOCALX}
/** Check if a @ref GenoLocation is @ref INVALID_GENO_LOCATION
 *
 * @shortnamed{IS_VALID_LOCATION} */
#define GSC_IS_VALID_LOCATION(g) (g.localAM != 0 && g.localPos != GSC_NA_LOCALX)

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


typedef struct {
    gsc_GenoLocation loc; /**< Location in the simulation where this parent is stored. */
    GSC_ID_T mapindex; /**< Index in d->genome.maps of the recombination map to use when producing gametes 
                                from this parent. Stores an index, not an ID */
} gsc_ParentChoice;

/** A structure to iterate forwards and backwards through all
 *  genotypes in a gsc_SimData or through only the members of a group.
 *
 * @shortnamed{BidirectionalIterator}
 *
 *  @see gsc_create_bidirectional_iter
 */
typedef struct {
    gsc_AlleleMatrix* am; /**< Simulation genotypes through which to iterate */
    const gsc_GroupNum group; /**< Group through which to iterate. If it is 0,
                          * then iterate through all genotypes in the simulation.
                          * Otherwise, iterate through members of the group with
                          * this as their group number. */
    GSC_LOCALX_T localPos; /**< Local index (index within the cachedAM) of the genotype in the linked list
                       * of gsc_AlleleMatrix beginning at `d->m` where the
                       * iterator's 'cursor' currently sits. */

    gsc_AlleleMatrix* cachedAM; /**< Pointer to the gsc_AlleleMatrix from the linked list
                              * of gsc_AlleleMatrix beginning at `d->m` where the
                              * iterator's 'cursor' currently sits. Contains
                              * the genotype at `localPos`. */
    unsigned int cachedAMIndex; /**< Index of `cachedAM` in the linked list of
                                  * gsc_AlleleMatrix beginning at `d->m`. `d->m`
                                  * is considered to be index 0. */

    _Bool atEnd; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                  * the last genotype (genotype with the highest index in the
                  * gsc_SimData) that fulfils the `group` critera of this iterator. */
    _Bool atStart; /**< Boolean that is TRUE if the iterator's 'cursor' is on
                    * the first genotype (genotype with the lowest index in the
                    * gsc_SimData) that fulfils the `group` critera of this iterator. */

} gsc_BidirectionalIterator;

/** A structure to iterate forwards through all
 *  positions in the gsc_AlleleMatrix linked list in gsc_SimData. Used
 *  in @a gsc_condense_allele_matrix. Internal, not recommended for end users.
 */
struct gsc_GappyIterator {
    gsc_GenoLocation cursor;
    unsigned int cursorAMIndex;
};

/** A structure to hold an initially empty AlleleMatrix list
 * whose genotypes can be accessed sequentially by storage order. 
 * For internal use by @a gsc_load_genotypefile.
 */
struct gsc_EmptyListNavigator {
    gsc_SimData* d;
    gsc_GroupNum alloctogroup;
    gsc_PedigreeID currentid;
    gsc_AlleleMatrix* firstAM;
    gsc_AlleleMatrix* localAM;
    GSC_LOCALX_T localPos;
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

    GSC_GLOBALX_T cacheSize; /**< Length in gsc_GenoLocations of `cache` */
    gsc_GenoLocation* cache; /**< Array iteratively updated with the known
                           * genotypes in the simulation that fulfil the
                           * `group` criteria of the iterator as they
                           * are discovered during calls to next_ functions */

    GSC_GLOBALX_T largestCached; /**< Local/group index (that is, index in `cache`) of the
                         * highest cell in `cache` that has been filled. */
    GSC_GLOBALX_T groupSize; /**< If the number of genotypes in the simulation that fulfil
                     * the iterator's `group` criteria is known, it is saved here.
                     * This value is left uninitialised until then. */
} gsc_RandomAccessIterator;

gsc_BidirectionalIterator gsc_create_bidirectional_iter( gsc_SimData* d, const gsc_GroupNum group);
gsc_BidirectionalIterator gsc_create_bidirectional_iter_fromAM( gsc_AlleleMatrix* am, const gsc_GroupNum group);
gsc_RandomAccessIterator gsc_create_randomaccess_iter( gsc_SimData* d, const gsc_GroupNum group);

gsc_AlleleMatrix* gsc_get_nth_AlleleMatrix( gsc_AlleleMatrix* listStart, const unsigned int n);

gsc_GenoLocation gsc_set_bidirectional_iter_to_start(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_set_bidirectional_iter_to_end(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_forwards(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_backwards(gsc_BidirectionalIterator* it);
gsc_GenoLocation gsc_next_get_nth(gsc_RandomAccessIterator* it, const GSC_GLOBALX_T n);

//static gsc_GenoLocation gsc_nextgappy_get_gap(struct gsc_GappyIterator* it);
//static gsc_GenoLocation gsc_nextgappy_get_nongap(struct gsc_GappyIterator* it);
//static gsc_GenoLocation gsc_nextgappy_valid_pos(struct gsc_GappyIterator* it);

// static struct gsc_EmptyListNavigator gsc_create_emptylistnavigator(SimData* d)
// static gsc_GenoLocation gsc_emptylistnavigator_get_first(struct gsc_EmptyListNavigator* it)
// static gsc_GenoLocation gsc_emptylistnavigator_get_next(struct gsc_EmptyListNavigator* it)
// static void gsc_emptylistnavigator_finaliselist(struct gsc_EmptyListNavigator* it)
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
int gsc_get_parents_of_id( const gsc_AlleleMatrix* start, const gsc_PedigreeID id, 
                                    gsc_PedigreeID output[static 2]);
void gsc_get_ids_of_names( const gsc_AlleleMatrix* start, const unsigned int n_names, const char** names, 
                          gsc_PedigreeID* output);
GSC_GLOBALX_T gsc_get_index_of_child( const gsc_AlleleMatrix* start, const gsc_PedigreeID parent1id, 
                                      const gsc_PedigreeID parent2id);
GSC_GLOBALX_T gsc_get_index_of_name( const gsc_AlleleMatrix* start, const char* name);
gsc_PedigreeID gsc_get_id_of_index( const gsc_AlleleMatrix* start, const GSC_GLOBALX_T index);
char* gsc_get_genes_of_index( const gsc_AlleleMatrix* start, const GSC_GLOBALX_T index);
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
GSC_GLOBALX_T gsc_get_group_size( const gsc_SimData* d, const gsc_GroupNum group_id);
GSC_GLOBALX_T gsc_get_group_genes( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, char** output);
GSC_GLOBALX_T gsc_get_group_names( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, char** output);
GSC_GLOBALX_T gsc_get_group_ids( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, gsc_PedigreeID* output);
GSC_GLOBALX_T gsc_get_group_indexes( const gsc_SimData* d, const gsc_GroupNum group_id, 
                                    GSC_GLOBALX_T group_size, GSC_GLOBALX_T* output);
GSC_GLOBALX_T gsc_get_group_bvs( const gsc_SimData* d, const gsc_GroupNum group_id, const gsc_EffectID effID, 
                       GSC_GLOBALX_T group_size, double* output);
GSC_GLOBALX_T gsc_get_group_parent_ids( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, 
                              const int whichParent, gsc_PedigreeID* output);
GSC_GLOBALX_T gsc_get_group_parent_names( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, 
                                const int whichParent, char** output);
GSC_GLOBALX_T gsc_get_group_pedigrees( const gsc_SimData* d, const gsc_GroupNum group_id, GSC_GLOBALX_T group_size, char** output);

unsigned int gsc_get_existing_groups( gsc_SimData* d, gsc_GroupNum* output);
unsigned int gsc_get_existing_group_counts( gsc_SimData* d, gsc_GroupNum* out_groups, GSC_GLOBALX_T* out_sizes);
    /**@}*/
/**@}*/


/** @defgroup groupmod Seletion/Group Modification Functions
 *
 * For simulation of selection or structure in breeding programs.
 *
 * @{
 */
gsc_GroupNum gsc_combine_groups( gsc_SimData* d, const unsigned int list_len, const gsc_GroupNum* grouplist);
gsc_GroupNum gsc_make_group_from( gsc_SimData* d, const unsigned int index_list_len, const GSC_GLOBALX_T* genotype_indexes);
gsc_GroupNum gsc_split_by_label_value( gsc_SimData* d, const gsc_GroupNum group, 
                                       const gsc_LabelID whichLabel, const int valueToSplit);
gsc_GroupNum gsc_split_by_label_range( gsc_SimData* d, const gsc_GroupNum group, const gsc_LabelID whichLabel, 
                                       const int valueLowBound, const int valueHighBound);

// GENERIC
unsigned int gsc_scaffold_split_by_somequality( gsc_SimData* d, const gsc_GroupNum group_id,
        void* somequality_data,
        gsc_GroupNum (*somequality_tester)(gsc_GenoLocation, void*, unsigned int, unsigned int, gsc_GroupNum*),
        unsigned int maxentries_results, gsc_GroupNum* results);
// APPLICATIONS
unsigned int gsc_split_into_individuals( gsc_SimData* d, const gsc_GroupNum group_id, 
                                        unsigned int maxentries_results, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_individuate(gsc_GenoLocation loc, void** datastore, 
        //unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum** results);
unsigned int gsc_split_into_families(gsc_SimData* d, const gsc_GroupNum group_id, unsigned int maxentries_results, 
                                     gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_family(gsc_GenoLocation loc, void** datastore, 
        //unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum** results);
unsigned int gsc_split_into_halfsib_families( gsc_SimData* d, const gsc_GroupNum group_id, 
                                              const int parent, unsigned int maxentries_results, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsib1(gsc_GenoLocation loc, void** datastore, 
        // unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsib2(gsc_GenoLocation loc, void** datastore, 
        // unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_quality_halfsibtemplate(gsc_GenoLocation loc, void** datastore, 
        //unsigned int maxgroups, unsigned int groupsfound, gsc_GroupNum* results, 
        //gsc_PedigreeID (*getparent)(gsc_GenoLocation));

// GENERIC
unsigned int gsc_scaffold_split_by_someallocation(gsc_SimData* d, const gsc_GroupNum group_id, void* someallocator_data,
        gsc_GroupNum (*someallocator)(gsc_GenoLocation, gsc_SimData*, void*, unsigned int, unsigned int*, gsc_GroupNum*),
        unsigned int n_outgroups, gsc_GroupNum* outgroups);
// APPLICATIONS
gsc_GroupNum gsc_split_evenly_into_two(gsc_SimData* d, const gsc_GroupNum group_id);
    //static gsc_GroupNum gsc_helper_split_by_allocator_knowncounts(gsc_GenoLocation loc, gsc_SimData* d, 
        //void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups)
unsigned int gsc_split_evenly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, 
                                     const unsigned int n, gsc_GroupNum* results);
unsigned int gsc_split_into_buckets(gsc_SimData* d, const gsc_GroupNum group_id, const unsigned int n, 
                                    const GSC_GLOBALX_T* counts, gsc_GroupNum* results);
gsc_GroupNum gsc_split_randomly_into_two(gsc_SimData* d, const gsc_GroupNum group_id);
unsigned int gsc_split_randomly_into_n(gsc_SimData* d, const gsc_GroupNum group_id, 
                                       const unsigned int n, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_allocator_equalprob(gsc_GenoLocation loc, gsc_SimData* d, 
        //void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups);
unsigned int gsc_split_by_probabilities(gsc_SimData* d, const gsc_GroupNum group_id, 
                                        const unsigned int n, const double* probs, gsc_GroupNum* results);
    //static gsc_GroupNum gsc_helper_split_by_allocator_unequalprob(gsc_GenoLocation loc, gsc_SimData* d, 
        //void* datastore, unsigned int n_outgroups, unsigned int* subgroupsfound, gsc_GroupNum* outgroups);

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
void gsc_generate_gamete(gsc_SimData* d, const char* parent_genome, char* output, const GSC_ID_T mapindex);
void gsc_generate_doubled_haploid(gsc_SimData* d, const char* parent_genome, char* output, const GSC_ID_T mapindex);
void gsc_generate_clone(gsc_SimData* d, const char* parent_genome, char* output);
    /**@}*/

// SUPPORTER FUNCS FOR GENOPTIONS OPTIONS
// static FILE* gsc_helper_genoptions_save_pedigrees_setup(const gsc_GenOptions g);
// static FILE* gsc_helper_genoptions_save_bvs_setup(const gsc_SimData* d, const gsc_GenOptions g, int* effIndexp);
// static FILE* gsc_helper_genoptions_save_genotypes_setup(const gsc_SimData* d, const gsc_GenOptions g);
// static void gsc_helper_genoptions_save_pedigrees(FILE* fp, gsc_SimData* d, gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_save_bvs(FILE* fe, gsc_MarkerEffects* effMatrices, int effIndex, 
        //gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_save_genotypes(FILE* fg, gsc_AlleleMatrix* tosave);
// static void gsc_helper_genoptions_give_names_and_ids(gsc_AlleleMatrix* am, gsc_SimData* d, 
        //const gsc_GenOptions g);

// PARAMETER FUNCTIONS FOR THE FOLLOWING GENERIC
// static void gsc_helper_make_offspring_cross(gsc_SimData* d, union gsc_datastore_make_genotype* datastore, 
        //gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_self_n_times(gsc_SimData* d, union gsc_datastore_make_genotype* datastore, 
        //gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_doubled_haploids(gsc_SimData* d, union gsc_datastore_make_genotype* datastore, 
        //gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere);
// static void gsc_helper_make_offspring_clones(gsc_SimData* d, union gsc_datastore_make_genotype* datastore, 
        //gsc_ParentChoice parents[static 2], gsc_GenoLocation putHere);

// static int gsc_helper_parentchooser_cross_randomly(void* parentIterator, 
        //union gsc_datastore_make_genotype* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]);
// static int gsc_helper_parentchooser_cross_randomly_between(void* parentIterator, 
        //union gsc_datastore_make_genotype* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]);
// static int gsc_helper_parentchooser_cross_targeted(void* parentIterator, 
        //union gsc_datastore_make_genotype* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]);
// static int gsc_helper_parentchooser_selfing(void* parentIterator, 
        //union gsc_datastore_make_genotype* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]);
// static int gsc_helper_parentchooser_cloning(void* parentIterator, 
        //union gsc_datastore_make_genotype* datastore, unsigned int* counter, gsc_ParentChoice parents[static 2]);

// static int gsc_helper_random_cross_checks(gsc_SimData* d, const gsc_GroupNum from_group, 
        // const unsigned int n_crosses, const int cap);

// GENERIC
union gsc_datastore_make_genotypes;
gsc_GroupNum gsc_scaffold_make_new_genotypes(gsc_SimData* d, const gsc_GenOptions g,
        void* parentIterator, union gsc_datastore_make_genotypes* datastore,
        int (*parentChooser)(void*, union gsc_datastore_make_genotypes*, GSC_GLOBALX_T*, gsc_ParentChoice[static 2]),
        void (*offspringGenerator)(gsc_SimData*, union gsc_datastore_make_genotypes*, gsc_ParentChoice[static 2], 
                                   gsc_GenoLocation));
// APPLICATIONS
gsc_GroupNum gsc_make_random_crosses(gsc_SimData* d, const gsc_GroupNum from_group, const GSC_GLOBALX_T n_crosses,
                                     const GSC_GLOBALX_T cap, const gsc_MapID which_map, const gsc_GenOptions g);
gsc_GroupNum gsc_make_random_crosses_between(gsc_SimData*d, const gsc_GroupNum group1, const gsc_GroupNum group2,
                                             const GSC_GLOBALX_T n_crosses, const GSC_GLOBALX_T cap1, 
                                             const GSC_GLOBALX_T cap2,
                                             const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g);
gsc_GroupNum gsc_make_targeted_crosses(gsc_SimData* d, const unsigned int n_combinations,
                                       const GSC_GLOBALX_T* firstParents, const GSC_GLOBALX_T* secondParents,
                                       const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g);
gsc_GroupNum gsc_self_n_times(gsc_SimData* d, const unsigned int n, const gsc_GroupNum group, 
                              const gsc_MapID which_map, const gsc_GenOptions g);
gsc_GroupNum gsc_make_doubled_haploids(gsc_SimData* d, const gsc_GroupNum group, 
                                       const gsc_MapID which_map, const gsc_GenOptions g);
gsc_GroupNum gsc_make_clones(gsc_SimData* d, const gsc_GroupNum group, 
                             const _Bool inherit_names, const gsc_GenOptions g);

gsc_GroupNum gsc_make_all_unidirectional_crosses(gsc_SimData* d, const gsc_GroupNum from_group, 
                                                 const gsc_MapID mapID, const gsc_GenOptions g);
gsc_GroupNum gsc_make_crosses_from_file(gsc_SimData* d, const char* input_file, 
                                        const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g);
gsc_GroupNum gsc_make_double_crosses_from_file(gsc_SimData* d, const char* input_file, 
                                        const gsc_MapID map1, const gsc_MapID map2, const gsc_GenOptions g);
/**@}*/


/** @defgroup calculators Breeding Value and Allele Count Calculators
 *
 * For calculations related to the loaded allele effects and internal additive breeding value model.
 *
 * @{
 */
gsc_GroupNum gsc_split_by_bv(gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID,
                             const GSC_GLOBALX_T top_n, const _Bool lowIsBest);
gsc_DecimalMatrix gsc_calculate_bvs( const gsc_SimData* d, const gsc_GroupNum group, const gsc_EffectID effID);
gsc_DecimalMatrix gsc_calculate_utility_bvs(gsc_BidirectionalIterator* targets, const gsc_MarkerEffects* effset);
gsc_DecimalMatrix gsc_calculate_utility_local_bvs(gsc_BidirectionalIterator* targets, 
                                                  gsc_MarkerBlocks b,
                                                  gsc_MarkerEffects e);
gsc_DecimalMatrix gsc_calculate_allele_counts(const gsc_SimData* d, 
                                              const gsc_GroupNum group, 
                                              const char allele);
gsc_MarkerBlocks gsc_create_evenlength_blocks_each_chr(const gsc_SimData* d, const gsc_MapID mapid, const GSC_ID_T n);
gsc_MarkerBlocks gsc_load_blocks(const gsc_SimData* d, const char* block_file);
gsc_DecimalMatrix gsc_calculate_local_bvs(const gsc_SimData* d, 
                                                const gsc_GroupNum group,
                                                const gsc_MarkerBlocks b, 
                                                const gsc_EffectID effID); 
void gsc_calculate_optimal_haplotype(const gsc_SimData* d, 
                                      const gsc_EffectID effID, 
                                      const char symbol_na,
                                      char* opt_haplotype);
void gsc_calculate_optimal_possible_haplotype(const gsc_SimData* d, 
                                               const gsc_GroupNum group, 
                                               const gsc_EffectID effID,
                                               const char symbol_na,
                                               char* opt_haplotype);
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
void gsc_delete_genome(gsc_KnownGenome* g);
void gsc_delete_recombination_map(SimData* d, const gsc_MapID whichMap);
void gsc_delete_recombination_map_nointegrity(gsc_RecombinationMap* m);
void gsc_delete_allele_matrix(gsc_AlleleMatrix* m);
void gsc_delete_effects_table(gsc_MarkerEffects* m);
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

// User-facing saving functions
void gsc_save_markerblocks(const char* fname, const gsc_SimData* d, const gsc_MarkerBlocks b, 
                           const gsc_MapID labelMapID);
void gsc_save_genotypes(const char* fname, const gsc_SimData* d, const gsc_GroupNum groupID, 
                        const _Bool markers_as_rows);
void gsc_save_allele_counts(const char* fname, 
                            const gsc_SimData* d, 
                            const gsc_GroupNum groupID, 
                            const char allele, 
                            const _Bool markers_as_rows);
void gsc_save_pedigrees(const char* fname, const gsc_SimData* d, const gsc_GroupNum groupID, 
                        const _Bool full_pedigree);
void gsc_save_bvs(const char* fname, const gsc_SimData* d, const gsc_GroupNum groupID, const gsc_EffectID effID);
void gsc_save_local_bvs(const char* fname, 
						const gsc_SimData* d, 
						const gsc_GroupNum groupID, 
						const gsc_MarkerBlocks b,
						const gsc_EffectID effID, 
						const _Bool headers);

// Utility and helper saving functions
void gsc_save_utility_markerblocks(FILE* f, const gsc_MarkerBlocks b, const GSC_GENOLEN_T n_markers, 
        char** const marker_names, const gsc_RecombinationMap* map);
void gsc_save_utility_genotypes(FILE* f, gsc_BidirectionalIterator* targets, const GSC_GENOLEN_T n_markers,
        char** const marker_names, const _Bool markers_as_rows);
void gsc_save_utility_allele_counts(FILE* f, 
                                    gsc_BidirectionalIterator* targets,
                                    GSC_GENOLEN_T n_markers, 
                                    char** const marker_names, 
                                    const _Bool markers_as_rows, 
                                    const char allele);
void gsc_save_utility_pedigrees(FILE* f, gsc_BidirectionalIterator* targets,
        const _Bool full_pedigree, const gsc_AlleleMatrix* parent_pedigree_store);
void gsc_save_utility_bvs(FILE* f, gsc_BidirectionalIterator* targets, const gsc_MarkerEffects* eff);

void gsc_save_utility_dmatrix(FILE* f, DecimalMatrix* dec, char** row_headers, char** col_headers, _Bool dim1_is_columns);

// static GSC_LOGICVAL gsc_helper_is_marker_in_chr(const GSC_GENOLEN_T markerix, const gsc_LinkageGroup chr, double* pos);

// static void gsc_scaffold_save_genotype_info(FILE* f, gsc_BidirectionalIterator* targets, 
//        GSC_GENOLEN_T n_markers, char** const marker_names, const _Bool markers_as_rows,
//        void (*bodycell_printer)(FILE*, gsc_GenoLocation, unsigned int, void*), void* bodycell_printer_data);
// static void gsc_helper_output_genotypematrix_cell(FILE* f, gsc_GenoLocation loc, GSC_GENOLEN_T markerix, void* NA);
// static void gsc_helper_output_countmatrix_cell(FILE* f, gsc_GenoLocation loc, GSC_GENOLEN_T markerix, void* data);

// static void gsc_scaffold_save_ancestry_of(const gsc_AlleleMatrix* m, gsc_PedigreeID p1, gsc_PedigreeID p2,
//        void (*strprinter)(char*, unsigned int, void*), void (*intprinter)(int, void*), void* printer_data);
// static void gsc_helper_ancestry_strprinter_file(char* str, unsigned int strlen, void* data);
// static void gsc_helper_ancestry_intprinter_file(int i, void* data);

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
int* gsc_calculate_min_recombinations_fw1(gsc_SimData* d, gsc_MapID mapid, 
                                          char* parent1, unsigned int p1num, 
                                          char* parent2, unsigned int p2num, 
                                          char* offspring, int certain); // forward filling, window size 1
int* gsc_calculate_min_recombinations_fwn(gsc_SimData* d, gsc_MapID mapid, 
                                          char* parent1, unsigned int p1num, 
                                          char* parent2, unsigned int p2num, 
                                          char* offspring, int window_size, 
                                          int certain); // forward filling, window size n

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
static inline int gsc_has_same_alleles(const char* p1, const char* p2, const unsigned int i) {
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
static inline int gsc_has_same_alleles_window(const char* g1, const char* g2, const unsigned int start, const unsigned int w) {
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
