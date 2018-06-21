/*
 * sw_align.c
 * A source file for sw_align.h
 * Implemented by Thomas Sunghoon Heo
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ref = reference sequence allocated on memory
// query = query sequence allocated on memory
// m = match score in integer > 0
// mm = mismatch score in integer < 0
// go = gap open penalty in integer > 0
// ge = gap extension penalty in integer > 0
sw_aligner_t *sw_aligner_init(char *ref, char* query, int m, int mm, int go, int ge)
{
	// assumes reference and query are already allocated on mem
	sw_aligner_t *A = (sw_aligner_t*)malloc(sizeof(sw_aligner_t));
	A->alignment = (sw_align_t *)malloc(sizeof(sw_align_t));
	size_t rows = strlen(ref) + 1;
	size_t cols = strlen(query) + 1;
	int i;
	A->rows = rows;
	A->cols = cols;
	A->matrix = (int **)malloc(sizeof(int*) * rows);
	for(i = 0; i < rows; i++)
	{
		A->matrix[i] = (int *)malloc(sizeof(int) * cols);
		memset(A->matrix[i] , 0 , sizeof(int) * cols);
	}
	A->gap_open = go;
	A->gap_ext = ge;
	A->match_score = m;
	A->mismatch_score = mm;
	A->alignment->reference = ref;
	A->alignment->query = query;
	A->alignment->map_pos = -1;
	A->alignment->map_end = -1;
	A->alignment->cigar = NULL;
	return A;
}

void sw_aligner_destroy(sw_aligner_t *aligner)
{
	int i;
	for(i = 0; i < aligner->rows; i++)
	{
		free(aligner->matrix[i]);
	}
	free(alinger->matrix);
	free(alinger->alignment->cigar);
	free(aligner->alignment);
	free(aligner);
	*(&aligner) = NULL;
}


sw_aux_char_ll_t *sw_init_char_vector()
{
	sw_aux_char_ll_t* l = (sw_aux_char_ll_t*)malloc(sizeof(sw_aux_char_ll_t));
	l->start_node = NULL;
	l->end_node=NULL;
	l->N = 0;
	return l;
}

void sw_destroy_char_vector(sw_aux_char_ll_t* v)
{
	sw_char_node_it it2, it1 = v->start_node;
	while(it1 != NULL)
	{
		it2 = it1->next;
		free(it1);
		it1 = it2;
	}
}

sw_aux_numeric_ll_t *sw_init_numeric_vector()
{
	sw_aux_numeric_ll_t* l = (sw_aux_numeric_ll_t*)malloc(sizeof(sw_aux_numeric_ll_t));
	l->start_node = NULL;
	l->end_node=NULL;
	l->N = 0;
	return l;
}

void sw_destroy_numeric_vector(sw_aux_numeric_ll_t* v)
{
	sw_numeric_node_it it2, it1 = v->start_node;
	while(it1 != NULL)
	{
		it2 = it1->next;
		free(it1);
		it1 = it2;
	}
}

sw_char_node_it sw_get_next_char_node(sw_char_node_it node)
{
	if(node != NULL)
	{
		return node->next;
	}
	else
	{
		return NULL;
	}
}

sw_char_node_it sw_get_first_char_node(sw_char_ll_t *vector)
{
	if(vector != NULL)
	{
		return vector->start_node;
	}
	else
	{
		return NULL;
	}
}

sw_numeric_node_it sw_get_next_numeric_node(sw_numeric_node_it node)
{
	if(node != NULL)
	{
		return node->next;
	}
	else
	{
		return NULL;
	}
}

sw_numeric_node_it sw_get_first_numeric_node(sw_numeric_ll_t *vector)
{
	if(vector != NULL)
	{
		return vector->start_node;
	}
	else
	{
		return NULL;
	}
}

sw_char_node_it sw_get_end_char_node(sw_char_ll_t *vector)
{
	if(vector != NULL)
	{
		return vector->end_node;
	}
	else
	{
		return NULL;
	}
}

void sw_push_char_vector(char a, sw_char_ll_t *v)
{
	if(v == NULL)
	{
		return;
	}
	if(v->start_node == NULL)
	{
		v->start_node = (sw_aux_char_node_t *)malloc(sizeof(sw_aux_char_node_t));
		v->start_node->next = NULL;
		v->start_node->letter = a;
		v->end_node = v->start_node;
	}
	else
	{
		v->end_node->next = (sw_aux_char_node_t *)malloc(sizeof(sw_aux_char_node_t));
		v->end_node->letter= a;
		v->end_node->next->next = NULL;
		v->end_node = v->end_node->next;
	}
}

sw_numeric_node_it sw_get_end_numeric_node(sw_numeric_ll_t *vector)
{
	if(vector != NULL)
	{
		return vector->end_node;
	}
	else
	{
		return NULL;
	}
}

void sw_push_numeric_vector(numeric a, sw_numeric_ll_t *v)
{
	if(v == NULL)
	{
		return;
	}
	if(v->start_node == NULL)
	{
		v->start_node = (sw_aux_numeric_node_t *)malloc(sizeof(sw_aux_numeric_node_t));
		v->start_node->next = NULL;
		v->start_node->letter = a;
		v->end_node = v->start_node;
	}
	else
	{
		v->end_node->next = (sw_aux_numeric_node_t *)malloc(sizeof(sw_aux_numeric_node_t));
		v->end_node->letter= a;
		v->end_node->next->next = NULL;
		v->end_node = v->end_node->next;
	}
}

// helper functions 

int sw_aligner_similarity_score(sw_aligner_t * aligner, char a, char b)
{
	if(a == b)
	{
		return aligner->match_score;
	}
	else
	{
		return aligner->mismatch_score;
	}
}

int sw_aligner_affine_penalty(sw_aligner_t *aligner, int k)
{
	return aligner->gap_open * (k-1) + aligner->gap_ext;
}

int sw_aligner_max_value(int L, int D, int U)
{
	int max = -1;
	max = (max > L ? max : L);
	max = (max > D ? max : D);
	max = (max > U ? max : U);
	return max;
}

int sw_aligner_max_value_loc(int L, int D, int U)
{
	int max = sw_aligner_max_value(L,D,U);
	if(max == D)
	{
		// diagonal
		return 0;
	}
	else if(max == U)
	{
		// upper
		return 1;
	}
	else
	{
		// left
		return -1;
	}
}

void sw_aligner_fill_matrix(sw_aligner_t *aligner)
{
	char *ref = aligner->alignment->reference;
	char *query = aligner->alignment->reference;
	size_t rows = aligner->rows;
	size_t cols = aligner->cols;
	size_t i,j,k;
	int *col_tmp;
	int *row_tmp;
	for(i = 1; i < rows; i++)
	{
		for(j = 1; j < cols; j++)
		{
			
		}
	}
}
