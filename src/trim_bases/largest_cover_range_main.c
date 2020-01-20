#include "../common/makedb_aux.h"
#include "pm4_aux.h"
#include "largest_cover_range.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s m4 packed_reads_dir error_cutoff min_ovlp_size min_cov min_read_size num_threads lcrv_path\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 9) {
		print_usage(argv[0]);
		return 1;
	}

	const char* m4_path = argv[1];
	const char* reads_dir = argv[2];
	const double min_ident_perc = 100.0 - 100.0 * atof(argv[3]);
	const int min_ovlp_size = atoi(argv[4]);
	const int min_cov = atoi(argv[5]);
	const int min_size = atoi(argv[6]);
	const int num_threads = atoi(argv[7]);
	const char* output = argv[8];

	int num_reads = load_num_reads(reads_dir);
	LargestCoverRange* lcrv = (LargestCoverRange*)calloc((num_reads+2), sizeof(LargestCoverRange));
	int num_partitions = load_num_partitions(m4_path);
	for (int i = 0; i < num_partitions; ++i) {
		get_largest_cover_range_for_one_partition(m4_path, 
				i, 
				lcrv,
				min_ident_perc,
				min_ovlp_size,
				min_cov,
				num_threads);
	}

	DFOPEN(out, output, "w");
	fprintf(out, "0\t0\t0\t0\n");
	for (int i = 1; i <= num_reads; ++i) {
		if (lcrv[i].size == 0) {
			invalidate_lcr(lcrv[i]);
		} else if (lcrv[i].right - lcrv[i].left < min_size) {
			invalidate_lcr(lcrv[i]);
		}
		fprintf(out, "%d\t%d\t%d\t%d\n", i, lcrv[i].left, lcrv[i].right, lcrv[i].size);
	}
	FCLOSE(out);
}
