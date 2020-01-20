#include "../common/map_options.h"
#include "rm_worker.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] reads reference output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions(&sDefaultReferenceMapingOptions);
}

int main(int argc, char* argv[])
{
	TIMING_START(__func__);
	
	MapOptions options = sDefaultReferenceMapingOptions;
	if (argc < 4 || (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* reads_path = argv[argc - 3];
	const char* reference_path = argv[argc - 2];
	const char* output = argv[argc - 1];
	rm_main(&options, reads_path, reference_path, output);
	
	TIMING_END(__func__);
}
