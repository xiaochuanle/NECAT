#include "pm4_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s wrk_dir m4 error_cutoff num_threads\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		return 1;
	}
	const char* wrk_dir = argv[1];
	const char* m4_path = argv[2];
	const double min_ident_perc = 100.0 - 100.0 * atof(argv[3]);
	int num_threads = atoi(argv[4]);
	num_threads = OC_MIN(num_threads, 8);
	int num_dumpped_files = 100;
	int partition_size = 100000;
	pm4_main(wrk_dir, m4_path, num_dumpped_files, partition_size, min_ident_perc, num_threads);
	return 0;
}
