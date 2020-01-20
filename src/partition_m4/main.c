#include "pm4.h"
#include "pm4_aux.h"
#include "pm4_options.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk_dir candidates\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_PM4Options();
}

int main(int argc, char* argv[])
{
	PM4Options options;
	if (argc < 3 || parse_PM4Options(argc - 2, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	const char* wrk_dir = argv[argc - 2];
	const char* m4_path = argv[argc - 1];
	pm4_main(&options, wrk_dir, m4_path);
	return 0;
}
