/*
 * sw_align.h
 * A header file for Smith-Waterman alignment in C
 * Implemented by Thomas Sunghoon Heo
 */
 
#ifndef _SW_ALIGN_THOMAS_H
#define _SW_ALIGN_THOMAS_H
#ifdef __cplusplus
extern "C"
{
#endif

// Declaration of Smith-Waterman alignment result
typedef struct sw_align
{
	int map_pos;
	int map_end;
	char* cigar;
	char* reference;
	char* query;
}sw_align_t;

// Declaration of Smith-Waterman alignment engine
typedef struct sw_aligner
{
	sw_align_t *alignment;
	int **score_matrix;
	int match_score;
	int mismatch_score;
	int gap_open;
	int gap_ext;
	int rows; // size of matrix
	int cols; // size of matrix
}sw_aligner_t;

// Declaration of aux linked list to store chars
typedef struct sw_aux_ll_char_node
{
	char letter;
	struct sw_aux_ll_char_node *next;
}sw_aux_char_node_t;
typedef struct sw_aux_ll
{
	int N; // number of elements
	sw_aux_char_node_t* start_node; // start of iterator
	sw_aux_char_node_t* end_node; // end of iterator
}sw_aux_char_ll_t;


typedef struct sw_aux_ll_numeric_node
{
	char letter;
	struct sw_aux_ll_numeric_node *next;
}sw_aux_numeric_node_t;
typedef struct sw_aux_ll2
{
	int N; // number of elements
	sw_aux_numeric_node_t* start_node; // start of iterator
	sw_aux_numeric_node_t* end_node; // end of iterator
}sw_aux_numeric_ll_t;

// iterators
typedef sw_aux_char_node_t* sw_char_node_it;
typedef sw_aux_numeric_node_t* sw_numric_node_it;

// Constructor and destructor
sw_aligner_t *sw_aligner_init(char *reference, char* query, int m, int mm, int go, int ge);
void sw_aligner_destroy(sw_aligner_t *aligner);

sw_aux_char_ll_t *sw_init_char_vector();
void sw_destroy_char_vector(sw_aux_char_ll_t* v);

sw_aux_numeric_ll_t *sw_init_numeric_vector();
void sw_destroy_numeric_vector(sw_aux_numeric_ll_t* v);


// Functions for three major steps of alignment
void sw_aligner_fill_matrix(sw_aligner_t *aligner);
void sw_aligner_back_trace(sw_aligner_t *aligner);
sw_align_t* sw_aligner_align(sw_aligner_t *aligner);

// Helper function
int sw_aligner_similarity_score(sw_aligner_t *aligner, char a, char b);
int sw_aligner_affine_penalty(sw_aligner_t *aligner, int k);
int sw_aligner_max_value(int L, int D, int U);
int sw_aligner_max_value_loc(int L, int D, int U);

// Helpers for linked list
sw_char_node_it sw_get_next_char_node(sw_char_node_it node);
sw_char_node_it sw_get_first_char_node(sw_char_ll_t *vector);
sw_char_node_it sw_get_end_char_node(sw_char_ll_t *vector);
void sw_push_char_vector(char a, sw_char_ll_t *v);

sw_numeric_node_it sw_get_next_numeric_node(sw_numeric_node_it node);
sw_numeric_node_it sw_get_first_numeric_node(sw_numeric_ll_t *vector);
sw_numeric_node_it sw_get_end_numeric_node(sw_numeric_ll_t *vector);
void sw_push_numeric_vector(int a, sw_numeric_ll_t *v);

#ifdef __cplusplus
}
#endif
#endif /*_SW_ALIGN_THOMAS_H*/