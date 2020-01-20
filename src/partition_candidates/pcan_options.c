#include "pcan_options.h"

#include <getopt.h>

#include "../common/ontcns_aux.h"

static const char*
argn_list = "p:f:t:";

static const PcanOptions sDefaultPcanOptions = {
	100000,
	100,
	1
};

int
parse_PcanOptions(int argc, char* argv[], PcanOptions* options)
{
	*options = sDefaultPcanOptions;
	int c;
	while ((c = getopt(argc, argv, argn_list)) != -1) {
		switch (c) {
			case 'p':
				options->batch_size = atoi(optarg);
				break;
			case 'f':
				options->num_output_files = atoi(optarg);
				break;
			case 't':
				options->num_threads = atoi(optarg);
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
print_PcanOptions(const PcanOptions* options)
{
	FILE* out = stderr;
	fprintf(out, "-p %d ", options->batch_size);
	fprintf(out, "-f %d ", options->num_output_files);
	fprintf(out, "-t %d ", options->num_threads);
	fprintf(out, "\n");
}

void
describe_PcanOptions()
{
	FILE* out = stderr;
	fprintf(out, "-p <Integer> batch size\n");
	fprintf(out, "-f <Integer> number of partition files\n");
	fprintf(out, "-t <Integer> number threads\n");
	
	fprintf(out, "\n");
	fprintf(out, "DEFAULT OPTIONS:\n");
	print_PcanOptions(&sDefaultPcanOptions);
}
