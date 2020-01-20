#include "pcan.h"
#include "pcan_aux.h"
#include "pcan_options.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk_dir candidates\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_PcanOptions();
}

int main(int argc, char* argv[])
{
	PcanOptions options;
	if (argc < 3 || parse_PcanOptions(argc - 2, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	const char* wrk_dir = argv[argc - 2];
	const char* can_path = argv[argc - 1];
	pcan_main(&options, wrk_dir, can_path);
	return 0;
}
