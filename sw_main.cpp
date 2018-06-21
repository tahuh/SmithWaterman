#include <string.h> // for strdup
#include <stdlib.h> // for free of example below
#include <iostream>
#include "sw_align.h"

int main()
{
	const char *ref = "TGTTACGG";
	const char *q = "GGTTGACTA";
	char *ref_copy = strdup(ref);
	char *query_copy = strdup(q);
	
	// creates aligner object
	// always important that use pointers for both
	// aligner and alignment result
	
	SWAligner* aligner = new SWAligner(ref_copy, query_copy);
	SWAlignment* alignment = aligner->align();
	std::cout << "Reference : " << alignment->get_refseq() << std::endl;
	std::cout << "Query : " << alignment->get_query() << std::endl;
	std::cout << "cigar string : " << alignment->get_cigar() << std::endl;
	
	// delete from memory
	delete alignment;
	delete aligner;
	free(ref_copy);
	free(query_copy);
	return 0;
}