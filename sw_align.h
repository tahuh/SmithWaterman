/*
 * sw_align.h
 * A C++ version of Smith-Waterman alignment
 * This is header file for aligner
 * We will use affine penalty for gap
 * Implemented by Thomas Sunghoon Heo
 */
 
#ifndef SW_ALIGN_H
#define SW_ALIGN_H

#include <vector>
#include <utility>

class SWAlignment
{
private:
	char *reference;
	char *query;
	char *cigar_string;
	int map_pos;
	int map_end;
public:
	SWAlignment();
	~SWAlignment();
	void set_reference(char *s);
	void set_query(char *s);
	void set_mapping_pos(int v);
	void set_map_end(int v);
	void build_cigar(std::vector<char> &pset);
	char* get_refseq();
	char* get_query();
	char *get_cigar();
	int get_mapping_pos();
};
class SWAligner
{
private:
	char* reference;
	char* query;
	int gap_ext;
	int gap_open;
	int match_score;
	int mismatch_score;
	int **matrix;
	int row_size;
	int col_size;
public:
	SWAligner();
	SWAligner(char* ref, char* query);
	~SWAligner();
	void set_reference(char *s);
	void set_query(char *s);
	void set_gap_ext_penalty(int v);
	void set_gap_open_penalty(int v);
	void set_match_score(int v);
	void set_mismatch_score(int v);
	void init_basic();
	int affine_penalty(int k);
	int find_max(int L,int D,int U);
	int find_max_loc(int L,int D,int U);
	void forward_pass();
	SWAlignment* align();
};

#endif