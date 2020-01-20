#include "../common/makedb_aux.h"
#include "pm4_aux.h"
#include "largest_cover_range.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s m4 packed_reads_dir error_cutoff min_ovlp_size min_cov min_read_size output num_threads\n", prog);
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
	const char* output = argv[7];
	const int num_threads = atoi(argv[8]);
	PackedDB* reads = NULL;
	reads = merge_volumes(reads_dir);
printf("Number of reads: %lu\n", pdb_num_seqs(reads));
	DFOPEN(out, output, "w");
	int read_id = 1;
	pthread_mutex_t out_lock;
	pthread_mutex_init(&out_lock, NULL);
	
	int num_partitions = load_num_partitions(m4_path);
	for (int i = 0; i < num_partitions; ++i) {
		DFOPEN(log_out, "lcr_log.txt", "a+");
		fprintf(log_out, "processing partition %d\n", i);
		FCLOSE(log_out);
		get_largest_cover_range_for_one_partition(m4_path, 
				i, 
				reads,
				&read_id,
				out,
				&out_lock,
				min_ident_perc,
				min_ovlp_size,
				min_cov,
				num_threads);
	}
	FCLOSE(out);
	free_PackedDB(reads);
}
