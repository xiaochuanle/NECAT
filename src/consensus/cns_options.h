#ifndef CNS_OPTONS_H
#define CNS_OPTONS_H

#include "../common/ontcns_defs.h"

typedef struct {
	int min_align_size;
	int min_cov;
	int max_cov;
	int min_size;
	int full_consensus;
	double error;
	double mapping_ratio;
	int num_threads;
	int rescue_long_indels;
	int use_fixed_ident_cutoff;
	int small_memory;
} CnsOptions;

void
describe_CnsOptions();

BOOL
parse_CnsOptions(int argc, char* argv[], CnsOptions* options);

void
print_CnsOptions(const CnsOptions* options);

#endif // CNS_OPTONS_H
