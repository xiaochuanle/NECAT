#include "cns_options.h"

#include <stdio.h>
#include <getopt.h>

#include "../common/ontcns_aux.h"

static const char* argn_list = "a:x:y:l:f:e:p:t:r:u:s:";

static const CnsOptions sDefaultCnsOptions = {
	400,  	// minimal align size
	4, 		// minimal coverage 
	12, 	// maximal coverage 
	500, 	// minimal corrected read size
	0,		// full consensus
	0.5, 	// sequencing error
	0.8,	// mapping ratio
	1, 		// number of threads
	0, 		// rescue long indels
	0,		// use fixed ident cutoff
	0		// small memory
};

void
print_CnsOptions(const CnsOptions* options)
{
	FILE* out = stderr;
	fprintf(out, "-a %d ", options->min_align_size);
	fprintf(out, "-x %d ", options->min_cov);
	fprintf(out, "-y %d ", options->max_cov);
	fprintf(out, "-l %d ", options->min_size);
	fprintf(out, "-f %d ", options->full_consensus);
	fprintf(out, "-e %f ", options->error);
	fprintf(out, "-p %f ", options->mapping_ratio);
	fprintf(out, "-t %d ", options->num_threads);
	fprintf(out, "-r %d ", options->rescue_long_indels);
	fprintf(out, "-u %d ", options->use_fixed_ident_cutoff);
	fprintf(out, "-s %d	", options->small_memory);
	fprintf(out, "\n");
}

BOOL
parse_CnsOptions(int argc, char* argv[], CnsOptions* options)
{
	*options = sDefaultCnsOptions;
	int c;
	while ((c = getopt(argc, argv, argn_list)) != -1) {
		switch (c) {
			case 'a':
				options->min_align_size = atoi(optarg);
				break;
			case 'x':
				options->min_cov = atoi(optarg);
				break;
			case 'y':
				options->max_cov = atoi(optarg);
				break;
			case 'l':
				options->min_size = atoi(optarg);
				break;
			case 'f':
				options->full_consensus = atoi(optarg);
				break;
			case 'e':
				options->error = atof(optarg);
				break;
			case 'p':
				options->mapping_ratio = atof(optarg);
				break;
			case 't':
				options->num_threads = atoi(optarg);
				break;
			case 'r':
				options->rescue_long_indels = atoi(optarg);
				break;
			case 'u':
				options->use_fixed_ident_cutoff = atoi(optarg);
				break;
			case 's':
				options->small_memory = atoi(optarg);
				break;
			case '?':
				fprintf(stderr, "invalid option '%c'\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			case ':':
				fprintf(stderr, "argument to option '%c' is not provided\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			default:
				break;
		}
	}
	return ARG_PARSE_SUCCESS;
}

void
describe_CnsOptions()
{
	FILE* out = stderr;
	fprintf(out, "-a <Integer>\talign length cutoff\n");
	fprintf(out, "-x <Integer>\tminimal coverage\n");
	fprintf(out, "-y <Integer>\tmaximal coverage\n");
	fprintf(out, "-l <Integer>\tminimal length of corrected reads.\n");
	fprintf(out, "-f <0 or 1>\tfull consensus or not: 1 = yes, 0 = no\n");
	fprintf(out, "-e <Real>\tsequencing error\n");
	fprintf(out, "-p <Real>\tminimal mapping ratio\n");
	fprintf(out, "-t <Integer>\tnumber of cpu threads\n");
	fprintf(out, "-r <0 or 1>\trescue long indels or not: 1 = yes, 0 = no\n");
	fprintf(out, "-u <0 or 1>\tuse dynamic or fixed ident cutoff: 1 = fixed, 0 = dynamic\n");
	fprintf(out, "-s <0 or 1>\tuse small memoty");
	
	fprintf(out, "\n");
	fprintf(out, "DEFAULT OPTIONS:\n");
	print_CnsOptions(&sDefaultCnsOptions);
}
