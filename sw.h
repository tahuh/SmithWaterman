/*
 * sw.h
 * Smith-Waterman alignment algorithm
 * Implemented by Thomas Sunghoon Heo
 */

#ifndef _SW_H
#define _SW_H

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
class SWAlignment
{
private:
    char *ref;
    char *query;
    char *cigar_string;
    int map_pos;
    int mem_alloc;
public:
    SWAlignment();
    SWAlignment(char *ref, char *query, int pos, const char *cigar);
    ~SWAlignment()
    {
        if(mem_alloc)
        {
            free(ref); free(query); free(cigar_string);
        }
    }
    int get_mapping_pos();
    char *get_cigar_string();
    char *get_reference();
    char *get_query();
}; // end of class SWalignment
class SWAligner
{
private:
    char *ref;
    char *query;
    double **matrix;
    int rows;
    int cols;
    int gap_open;
    int gap_ext;
    int similarity;
    int mismatch;
    bool is_matrix_init;
    bool is_seq_init;
public:
    SWAligner();
    SWAligner(std::string ref, std::string query);
    SWAligner(std::string ref, std::string query, int go, int ge, int m, int M);
    SWAligner(char *ref, char* query);
    SWAligner(char *ref, char* query, int go, int ge, int m, int M);
    SWAligner(const char *ref, const char *query);
    SWAligner(const char *ref, const char *query, int go, int ge, int m, int M);
    ~SWAligner()
    {
        if(is_seq_init)
        {
            free(ref); free(query);
        }
        if(is_matrix_init)
        {
            for(int i = 0; i < rows; i++)
            {
                free(matrix[i]);
            }
            free(matrix);
        }
    }
    void init_initial();
    /* modifiers */
    void set_reference(std::string ref);
    void set_reference(char *ref);
    void set_reference(const char *ref);
    void set_query(std::string query);
    void set_query(char *q);
    void set_query(const char *query);
    void set_gap_open(int go);
    void set_gap_ext(int ge);
    void set_similarity(int s);
    void set_mismatch(int m);
    /* accessors */
    char *get_reference();
    char *get_query();
    int get_gap_open_penalty();
    int get_gap_ext_penalty();
    int get_similarity_score();
    int get_mismatch_penalty();

    // Algorithm functions
    void _build_matrix();
    void _fill_matrix();
    SWAlignment* _backtrace();
    // The main algorithm
    SWAlignment* align();

    // auxilaries
    double compute_similarity(char a, char b);
    double compute_gap_score(int k); // affine
    double find_max_value(double L, double D, double U);
    int find_max_direction(double L, double D, double U);
}; // end of class SWAligner

#endif
