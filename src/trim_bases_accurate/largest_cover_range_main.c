#include "../common/makedb_aux.h"
#include "pm4_aux.h"
#include "largest_cover_range.h"

#include <stdio.h>

static const double kErrorCutoff = 0.1;
static const int kMinOvlpSize = 1;
static const int kMinCov = 1;
static const int kMinReadSize = 1000;

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s m4 packed_reads_dir small_memoty num_threads output\n", prog);
	fprintf(out, "\n");
	fprintf(out, "If Multiple Nodes Are Used:\n");
	fprintf(out, "%s m4 packed_reads_dir small_memoty num_threads output -mn node_id num_nodes\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 6 && argc != 9) {
		print_usage(argv[0]);
		return 1;
	}
	int spid = 0;
	int n_node = 1;
	if (argc == 9) {
		oc_assert(strcmp(argv[6], "-mn") == 0);
		spid = atoi(argv[7]);
		n_node = atoi(argv[8]);
		oc_assert(spid >= 0);
		oc_assert(n_node >= 0);
	}

	const char* m4_path = argv[1];
	const char* reads_dir = argv[2];
	const int small_memoty = atoi(argv[3]);
	const double min_ident_perc = 100.0 - 100.0 * kErrorCutoff;
	const int min_ovlp_size = kMinOvlpSize;
	const int min_cov = kMinCov;
	const int min_size = kMinReadSize;
	const int num_threads = atoi(argv[4]);
	const char* output = argv[5];
	PackedDB* reads = NULL;
	if (!small_memoty) reads = merge_volumes(reads_dir);
	DFOPEN(out, output, "w");
	int read_id = 1;
	pthread_mutex_t out_lock;
	pthread_mutex_init(&out_lock, NULL);
	
	int num_partitions = load_num_partitions(m4_path);
	for (int i = spid; i < num_partitions; i += n_node) {
		DFOPEN(log_out, "lcr_log.txt", "a+");
		fprintf(log_out, "processing partition %d\n", i);
		FCLOSE(log_out);
		get_largest_cover_range_for_one_partition(m4_path, 
				i, 
				reads,
				reads_dir,
				&read_id,
				out,
				&out_lock,
				min_ident_perc,
				min_ovlp_size,
				min_cov,
				num_threads);
	}
	FCLOSE(out);
	if (reads) free_PackedDB(reads);
}
