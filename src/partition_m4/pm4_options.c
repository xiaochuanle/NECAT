#include "pm4_options.h"

#include <getopt.h>

#include "../common/ontcns_aux.h"

static const char*
argn_list = "p:f:";

static const PM4Options sDefaultPM4Options = {
	100000,
	10
};

int
parse_PM4Options(int argc, char* argv[], PM4Options* options)
{
	*options = sDefaultPM4Options;
	int c;
	while ((c = getopt(argc, argv, argn_list)) != -1) {
		switch (c) {
			case 'p':
				options->batch_size = atoi(optarg);
				break;
			case 'f':
				options->num_output_files = atoi(optarg);
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
print_PM4Options(const PM4Options* options)
{
	FILE* out = stderr;
	fprintf(out, "-p %d ", options->batch_size);
	fprintf(out, "-f %d ", options->num_output_files);
	fprintf(out, "\n");
}

void
describe_PM4Options()
{
	FILE* out = stderr;
	fprintf(out, "-p <Integer> batch size\n");
	fprintf(out, "-f <Integer> number of partition files\n");
	
	fprintf(out, "\n");
	fprintf(out, "DEFAULT OPTIONS:\n");
	print_PM4Options(&sDefaultPM4Options);
}
