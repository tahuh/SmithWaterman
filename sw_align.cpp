/*
 * sw_align.cpp
 * A source file for sw_align.h
 * Implemented by Thomas Sunghoon Heo
 */
 
// usage is described below
// or see sw_main.cpp at this repo
// #include <string.h> // for strdup
// #include <stdlib.h> // for free of example below
// #include <iostream>
// #include "sw_align.h"

// int main()
// {
	// const char *ref = "TGTTACGG";
	// const char *q = "GGTTGACTA";
	// char *ref_copy = strdup(ref);
	// char *query_copy = strdup(q);
	
	// // creates aligner object
	// // always important that use pointers for both
	// // aligner and alignment result
	
	// SWAligner* aligner = new SWAligner(ref_copy, query_copy);
	// SWAlignment* alignment = aligner->align();
	// std::cout << "Reference : " << alignment->get_refseq() << std::endl;
	// std::cout << "Query : " << alignment->get_query() << std::endl;
	// std::cout << "cigar string : " << alignment->get_cigar() << std::endl;
	
	// // delete from memory
	// delete alignment;
	// delete aligner;
	// return 0;
// }
#include <stdlib.h>
#include <string.h> // memcpy
#include <string>
#include <sstream>
#include <algorithm>
#include "sw_align.h"
SWAlignment::SWAlignment()
{
	reference = NULL;
	query = NULL;
	cigar_string = NULL;
	map_pos = -1;
	map_end = -1;
}
SWAlignment::~SWAlignment()
{
	if(reference != NULL)
	{
		free(reference);
	}
	if(query != NULL)
	{
		free(query);
	}
	if(cigar_string != NULL)
	{
		free(cigar_string);
	}
}

void SWAlignment::set_reference(char *s)
{
	reference = strdup(s);
}

void SWAlignment::set_query(char *s)
{
	query = strdup(s);
}

void SWAlignment::set_mapping_pos(int v)
{
	map_pos = v;
}

void SWAlignment::set_map_end(int v)
{
	map_end = v;
}
char *SWAlignment::get_refseq()
{
	return reference;
}

char *SWAlignment::get_query()
{
	return query;
}

char *SWAlignment::get_cigar()
{
	return cigar_string;
}

int SWAlignment::get_mapping_pos()
{
	return map_pos;
}

void SWAlignment::build_cigar(std::vector<char> &pset)
{
	
	std::vector<int> numbers;
	std::vector<char> cigars;
	for(auto it = pset.begin(); it != pset.end(); it++)
	{
		char letter = *it;
		if(letter == '^')
		{
			// query has deletion
			if(cigars.size() == 0)
			{
				numbers.push_back(1);
				cigars.push_back('D');
			}
			else
			{
				if(cigars.back() == 'D')
				{
					int last = numbers.back();
					int last2 = last + 1;
					numbers.pop_back();
					numbers.push_back(last2);
				}
				else
				{
					numbers.push_back(1);
					cigars.push_back('D');
				}
			}
		}
		else if(letter =='<')
		{
			// query has insertion
			if(cigars.size() == 0)
			{
				numbers.push_back(1);
				cigars.push_back('I');
			}
			else
			{
				if(cigars.back() == 'I')
				{
					int last = numbers.back();
					int last2 = last + 1;
					numbers.pop_back();
					numbers.push_back(last2);
				}
				else
				{
					numbers.push_back(1);
					cigars.push_back('I');
				}
			}
		}
		else
		{
			// matched
			if(cigars.size() == 0)
			{
				numbers.push_back(1);
				cigars.push_back('M');
			}
			else
			{
				if(cigars.back() == 'M')
				{
					int last = numbers.back();
					int last2 = last + 1;
					numbers.pop_back();
					numbers.push_back(last2);
				}
				else
				{
					numbers.push_back(1);
					cigars.push_back('M');
				}
			}
		}
	}
	// TODO: build final cigar string
	auto cigar_it = cigars.begin();
	auto num_it = numbers.begin();
	if(map_pos != 1)
	{
		num_it = numbers.insert(num_it, map_pos);
		cigar_it = cigars.insert(cigar_it, 'H');
	}
	if(map_end != strlen(query))
	{
		int diff = strlen(query) - map_end;
		numbers.push_back(diff);
		cigars.push_back('H');
	}
	
	std::ostringstream stream;
	for(int idx = 0; idx < cigars.size(); idx++)
	{
		stream << numbers[idx];
		stream << cigars[idx];
	}
	
	std::string str = stream.str();
	const char *tmp = str.c_str();
	cigar_string = strdup(tmp);
}
SWAligner::SWAligner()
{
	reference = NULL;
	query = NULL;
	matrix = NULL;
	row_size = 0;
	col_size = 0;
	init_basic();
}
SWAligner::SWAligner(char *r, char *q)
{
	reference = strdup(r);
	query = strdup(q);
	row_size = strlen(r) + 1;
	col_size = strlen(q) + 1;
	matrix = (int **)malloc(sizeof(int *) * row_size);
	for(int i = 0; i < col_size; i++)
	{
		matrix[i] = (int *)malloc(sizeof(int) * col_size);
		memset(matrix[i] , 0 , sizeof(int) * col_size);
	}
	init_basic();
}
void SWAligner::init_basic()
{
	match_score = 1;
	mismatch_score = -1;
	gap_open = 1;
	gap_ext = 1;
}
SWAligner::~SWAligner()
{
	if(matrix != NULL)
	{
		for(int i = 0; i < row_size; i++)
		{
			free(matrix[i]);
		}
		free(matrix);
	}
	
	if(reference != NULL)
	{
		free(reference);
	}
	
	if(query != NULL)
	{
		free(query);
	}
}

void SWAligner::set_reference(char *s)
{
	reference = strdup(s);
}

void SWAligner::set_query(char *s)
{
	query = strdup(s);
}

void SWAligner::set_gap_ext_penalty(int v)
{
	gap_ext = v;
}

void SWAligner::set_gap_open_penalty(int v)
{
	gap_open= v;
}

void SWAligner::set_match_score(int v)
{
	match_score = v;
}

void SWAligner::set_mismatch_score(int v)
{
	mismatch_score = v;
}
int SWAligner::affine_penalty(int k)
{
	return gap_open * (k-1) + gap_ext;
}

int SWAligner::find_max(int L, int D, int U)
{
	std::vector<int> v = {L,D,U};
	auto max_it = std::max_element(v.begin(), v.end());
	int max = *max_it;
	return max;
}

int SWAligner::find_max_loc(int L, int D, int U)
{
	int max = find_max(L,D,U);
	if(max == D)
	{
		return 0;
	}
	else if(max == L)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}
void SWAligner::forward_pass()
{
	// Assumes matrix is initialized
	int i,j,k,l,max_col, max_row,L,D,U;
	for(i = 1; i < row_size; i++)
	{
		for(j = 1; j < col_size; j++)
		{
			char ref_base = reference[i-1];
			char query_base = query[j-1];
			int score = (ref_base == query_base ? 1 : -1);
			
			std::vector<int> column_way_weights;
			std::vector<int> row_way_weights;
			
			for(k = 1; k < i; k++)
			{
				int wk = affine_penalty(k);
				int h = matrix[i-k][j];
				column_way_weights.push_back(h-wk);
			}
			
			for(l = 1; l < j; l++)
			{
				int wl = affine_penalty(k);
				int h2 = matrix[i][j-l];
				row_way_weights.push_back(h2-wl);
			}
			if(column_way_weights.size() == 0)
			{
				max_col = 0;
			}
			else
			{
				max_col = *std::max_element(column_way_weights.begin(), column_way_weights.end());
			}
			
			if(row_way_weights.size() == 0)
			{
				max_row = 0;
			}
			else
			{
				max_row = *std::max_element(row_way_weights.begin(), row_way_weights.end());
			}
			L = max_row;
			D = matrix[i-1][j-1];
			U = max_col;
			matrix[i][j] = find_max(L,D,U);
		}
	}
}
SWAlignment* SWAligner::align()
{
	// perform forward pass first
	forward_pass();
	// must use delete function to remove alignment result from mem
	SWAlignment* alignment = new SWAlignment();
	std::vector<char> trace_result;
	// before performing backtrace algorithm
	// we must find the maximum location
	int map_pos;
	int map_end;
	int max_i = -1;
	int max_j = -1;
	int max_val = -1;
	int direction, left, up, diag, current_value = -1;
	for(int i = 1; i <row_size; i++)
	{
		for(int j = 1; j < col_size; j++)
		{
			if(matrix[i][j] > max_val)
			{
				max_val = matrix[i][j];
				max_i = i;
				max_j = j;
			}
		}
	}
	map_end = max_j;
	// Here we start backtrace algorithm
	alignment->set_map_end(map_end);
	map_pos = 0;
	
	while(current_value != 0)
	{
		left =matrix[max_i][max_j-1];
		up = matrix[max_i-1][max_j];
		diag = matrix[max_i-1][max_j-1];
		direction = find_max_loc(left, diag, up);
		if(direction == -1)
		{
			// query has insertion
			max_j--;
			trace_result.push_back('<');
		}
		else if(direction == 0)
		{
			// matched
			max_i -= 1;
			max_j -= 1;
			trace_result.push_back('.');
		}
		else
		{
			// query has deletion
			max_i--;
			trace_result.push_back('^');
		}
		current_value = find_max(left, diag, up);
	}
	
	// done alignment
	map_pos = max_i + 1;
	alignment->set_mapping_pos(map_pos);
	alignment->set_reference(reference);
	alignment->set_query(query);
	alignment->build_cigar(trace_result);
	
	return alignment;
}