/*
 * sw.cpp
 * Algorithm source code for sw.h
 * Compiler requirements : C++11 or higher support
 * Implemented by Thomas Sunghoon Heo
 */

#include <sstream> // string stream
#include <algorithm> // std::max_element
#include "sw.h"

SWAlignment::SWAlignment()
{}
SWAlignment::SWAlignment(char *r, char*q, int p, const char *c)
{
    ref = strdup(r);
    query = strdup(q);
    map_pos = p;
    cigar_string = strdup(c);
    mem_alloc = 1;
}
int SWAlignment::get_mapping_pos()
{
    return map_pos;
}
char *SWAlignment::get_cigar_string()
{
    return cigar_string;
}
char *SWAlignment::get_reference()
{
    return ref;
}

char *SWAlignment::get_query()
{
    return query;
}

SWAligner::SWAligner()
{
    is_seq_init = false;
    is_matrix_init = false;
}

SWAligner::SWAligner(std::string r, std::string q)
{
    ref = strdup(r.c_str());
    query = strdup(q.c_str());
    is_seq_init = true;
    init_initial();
}
SWAligner::SWAligner(std::string r, std::string q, int o, int e, int m, int M)
{
    ref = strdup(r.c_str()); query = strdup(q.c_str());
    gap_open = o; gap_ext = e; similarity = m; mismatch= M;
    rows = strlen(ref) + 1; cols = strlen(query) + 1;
    is_seq_init = true; is_matrix_init = false;
}
SWAligner::SWAligner(char *r, char *q)
{
    ref = strdup(r); query = strdup(q); is_seq_init = true;
    init_initial();
}
SWAligner::SWAligner(char *r, char *q, int o , int e, int m , int M)
{
    ref = strdup(r); query = strdup(q); gap_open = 0; gap_ext = e;
    similarity = m; mismatch = M;
    rows = strlen(ref) + 1; cols = strlen(query) + 1;
    is_seq_init = true; is_matrix_init = false;
}
SWAligner::SWAligner(const char *r, const char *q)
{
    ref = strdup(r); query = strdup(q); is_seq_init = true;
    init_initial();
}
SWAligner::SWAligner(const char *r, const char *q, int o , int e, int m , int M)
{
    ref = strdup(r); query = strdup(q); gap_open = 0; gap_ext = e;
    similarity = m; mismatch = M;
    rows = strlen(ref) + 1; cols = strlen(query) + 1;
    is_seq_init = true; is_matrix_init = false;
}


void SWAligner::init_initial()
{
    rows = strlen(ref) + 1;
    cols = strlen(query) + 1;
    gap_open = 2;
    gap_ext = 3;
    similarity = 1;
    mismatch = -1;
    is_matrix_init = false;
}

void SWAligner::set_reference(std::string r)
{
    if(is_seq_init)
    {
        free(ref); 
        ref = strdup(r.c_str());
    }
    else
    {
        ref = strdup(r.c_str());
    }
}

void SWAligner::set_reference(char *r)
{
    if(is_seq_init)
    {
        free(ref); ref = strdup(r); 
    }
    else
    {
        ref = strdup(r);
    }
}

void SWAligner::set_reference(const char *r)
{
     if(is_seq_init)
    {
       free(ref);
    }
    ref = strdup(r);
}
void SWAligner::set_query(std::string r)
{
    if(is_seq_init)
    {
        free(query); 
        query = strdup(r.c_str());
    }
    else
    {
        query = strdup(r.c_str());
    }
}

void SWAligner::set_query(char *r)
{
    if(is_seq_init)
    {
        free(query); query = strdup(r); 
    }
    else
    {
        query = strdup(r);
    }
}

void SWAligner::set_query(const char *r)
{
     if(is_seq_init)
    {
       free(query);
    }
    query = strdup(r);
}
void SWAligner::set_gap_open(int v)
{
    gap_open = v;
}
void SWAligner::set_gap_ext(int v)
{
    gap_ext = v;
}
void SWAligner::set_similarity(int s)
{
    similarity = s;
}
void SWAligner::set_mismatch(int v)
{
    mismatch = v;
}
char* SWAligner::get_reference()
{
    return ref;
}
char *SWAligner::get_query()
{
    return query;
}
int SWAligner::get_gap_open_penalty()
{
    return gap_open;
}
int SWAligner::get_gap_ext_penalty()
{
    return gap_ext;
}
int SWAligner::get_similarity_score()
{
    return similarity;
}
int SWAligner::get_mismatch_penalty()
{
    return mismatch;
}
void SWAligner::_build_matrix()
{
    size_t per_cols = sizeof(double) * cols;
    size_t per_rows = sizeof(double *) * rows;
    matrix = (double **)malloc(per_rows);
    for(int i = 0; i < rows; i++)
    {
        matrix[i] = (double *)malloc(per_cols);
        memset(matrix[i] , 0.0, per_cols);
    }
}
void SWAligner::_fill_matrix()
{
    double max_row, max_col, current_maxima;
    for(int i = 1; i < rows; i++)
    {
        for(int j = 1; j < cols; j++)
        {
            char ref_base = ref[i-1];
            char query_base = query[j-1];
            double s = compute_similarity(ref_base, query_base);
            std::vector<double> column_way_weights;
            std::vector<double> row_way_weights;
            for(int k = 1; k < i; k++)
            {
                double wk = compute_gap_score(k);
                double h = matrix[i-k][j];
                column_way_weights.push_back(h - wk);
            }
            for(int l = 1; l < j; l++)
            {
                double wl = compute_gap_score(l);
                double hl = matrix[i][j-l];
                row_way_weights.push_back(hl-wl);
            }
            if(column_way_weights.size() == 0)
            {
                max_col = 0.0;
            }
            else
            {
                max_col = *std::max_element(column_way_weights.begin(),
                                           column_way_weights.end());
            }
            if(row_way_weights.size() == 0)
            {
                max_row = 0.0;
            }
            else
            {
                max_row = *std::max_element(row_way_weights.begin(),
                                            row_way_weights.end());
            }
            double L = max_row;
            double U = max_col;
            double D = matrix[i-1][j-1] + s;
            current_maxima = find_max_value(L,D,U);
            matrix[i][j] = current_maxima;
        } // end of index j for loop
    } // end of index i loop
}
SWAlignment* SWAligner::_backtrace()
{
    SWAlignment *alignment;
    int max_i = -1;
    int max_j = -1;
    double max_val = -1;
    double current_value = -1;
    // find local alignment score max
    for(int i = 1; i < rows; i++)
    {
        for(int j = 1; j < cols; j++)
        {
            double mat_val = matrix[i][j];
            if(mat_val > max_val)
            {
                max_val = mat_val; max_i = i; max_j = j;
            }
        }
    } 

    int alignment_end_pos = max_j;
    int alignment_start_pos = 0;
    std::vector<char> trace_result;
    while(current_value != 0)
    {
        double left = matrix[max_i][max_j-1];
        double up = matrix[max_i-1][max_j];
        double diag = matrix[max_i-1][max_j-1];
        int max_loc = find_max_direction(left, diag, up);
        if(max_loc == -1)
        {
            max_j--;
            trace_result.push_back('<');
        }
        else if (max_loc == 0)
        {
             max_i--; max_j--;
             trace_result.push_back('.');
        }
        else
        {
             max_i--;
             trace_result.push_back('^');
        }
        current_value = find_max_value(left, up, diag);
    }
    alignment_start_pos = max_i + 1;
    std::vector<int> numbers;
    std::vector<char> alphabets;
    for(auto it = trace_result.begin(); it != trace_result.end(); it++)
    {
        char letter = *it;
        if(letter == '<')
        {
            if(alphabets.size() == 0)
            {
                alphabets.push_back('I');
                numbers.push_back(1);
            }
            else
            {
                if(alphabets.back() == 'I')
                {
                    int back = numbers.back();
                    numbers.pop_back();
                    back = back + 1;
                    numbers.push_back(back);
                }
                else
                {
                    alphabets.push_back('I');
                    numbers.push_back(1);
                }
            }
        }
        else if(letter == '.')
        {
            if(alphabets.size() == 0)
            {
                alphabets.push_back('M');
                numbers.push_back(1);
            }
            else
            {
                if(alphabets.back() == 'M')
                {
                    int back = numbers.back();
                    numbers.pop_back();
                    back = back + 1;
                    numbers.push_back(back);
                }
                else
                {
                    alphabets.push_back('M');
                    numbers.push_back(1);
                }
            }
        }
        else
        {
            if(alphabets.size() == 0)
            {
                alphabets.push_back('D');
                numbers.push_back(1);
            }
            else
            {
                if(alphabets.back() == 'D')
                {
                    int back = numbers.back();
                    numbers.pop_back();
                    back = back + 1;
                    numbers.push_back(back);
                }
                else
                {
                    alphabets.push_back('D');
                    numbers.push_back(1);
                }
            }
        }
    } // end of cigar string generation
    if(alignment_start_pos != 1)
    {
        auto nit = numbers.begin();
        numbers.insert(nit , alignment_start_pos);
        auto char_it = alphabets.begin();
        alphabets.insert(char_it, 'H');
    }
    if(alignment_end_pos != cols-1)
    {
        size_t diff = strlen(query) - alignment_end_pos;
        numbers.push_back(diff);
        alphabets.push_back('H');
    }
    std::stringstream ss;
    for(int index = 0; index < numbers.size(); index++)
    {
        ss << numbers[index]; ss << alphabets[index];
    }
    std::string cpp_cig_str = ss.str();
    const char *cigar = cpp_cig_str.c_str();
    char *my_ref = get_reference();
    char *my_query = get_query();
    alignment = new SWAlignment(my_ref, my_query,
                                alignment_start_pos, cigar);
    return alignment;
}
SWAlignment* SWAligner::align()
{
    _build_matrix();
    _fill_matrix();
    SWAlignment *ptr = _backtrace();
    return ptr;
}

double SWAligner::compute_similarity(char a, char b)
{
    if(a == b)
    {
        return similarity;
    }
    else
    {
        return mismatch;
    }
}
double SWAligner::compute_gap_score(int k)
{
    int value = gap_open * (k - 1) + gap_ext;
    return (double)value;
}
double SWAligner::find_max_value(double L, double D, double U)
{
    std::vector<double> v;
    v.push_back(L); v.push_back(D); v.push_back(U);
    auto it = std::max_element(v.begin(), v.end());
    return *it;
}
int SWAligner::find_max_direction(double L, double D, double U)
{
    double max = find_max_value(L,D,U);
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
