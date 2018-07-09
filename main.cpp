/*
 * Test code for sw.cpp
 * compilation follows
 * g++ -std=c++11 -g -c sw.cpp
 * g++ -std=c++11 -g -c main.cpp
 * g++ -std=c++11 -g -o main main.o sw.o
 */
#include <iostream>
#include "sw.h"

const char *ref = "TGTTACGG";
const char *query = "GGTTGACTA";
int gap_e = 1;
int gap_o = 1;
int match = 1;
int mismatch = -1;
int main()
{
    SWAligner *aligner = new SWAligner(ref, query, gap_o, gap_e, match, mismatch);
    SWAlignment* alignment = aligner->align();
    std::cout << "CIGAR string : " << alignment->get_cigar_string() << std::endl;
    delete alignment;
    delete aligner;
}
