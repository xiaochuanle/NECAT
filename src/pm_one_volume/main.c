#include "pm_worker.h"
#include "../common/map_options.h"
#include "../common/makedb_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [options] wrk-dir volume-id output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions(&sDefaultPairwiseMapingOptions);
}

int
find_solution(const char* wrk_dir, const int num_threads)
{
	VolumesInfo* volumes = load_volumes_info(wrk_dir);
	int ng = (num_threads + GroupThreadSize - 1) / GroupThreadSize;
	int use_mg = 1;
	if (ng > volumes ->num_volumes) use_mg = 0;
	destroy_volumes_info(volumes);
	return use_mg;
}

int main(int argc, char* argv[])
{
	if (argc < 4) {
		print_usage(argv[0]);
		return 1;
	}
	MapOptions options;
	if (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	const char* wrk_dir = argv[argc - 3];
	int vid = atoi(argv[argc - 2]);
	const char* output = argv[argc - 1];
	if (0) {
		pm_multi_group(&options, vid, wrk_dir, output);
	} else {
		pm_main(&options, vid, wrk_dir, output);
	}
	return 0;
}
