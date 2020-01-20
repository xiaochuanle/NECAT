#include "read_id_pool.h"

#include "cns_options.h"
#include "../common/makedb_aux.h"
#include "../partition_candidates/pcan_aux.h"
#include "consensus_one_partition.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stderr;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s [options] wrk_dir candidates cns_out raw_out\n", prog);
	fprintf(out, "\n");
	fprintf(out, "If Multiple Nodes Are Used:\n");
	fprintf(out, "%s [options] wrk_dir candidates cns_out raw_out -mn node_id num_nodes\n", prog);
	fprintf(out, "\n");
	fprintf(out, "OPTIONS AND DESCRIPTIONS:\n");
	describe_CnsOptions();
}

static int parse_params(int argc, 
			char* argv[], 
			CnsOptions* options, 
			const char** wrk_dir, 
			const char** candidates_path, 
			const char** cns_out_path, 
			const char** raw_out_path, 
			int* spid, 
			int* nnode)
{
	if (argc < 5) return -1;
	*spid = 0;
	*nnode = 1;
	if (strcmp(argv[argc - 3], "-mn") == 0) {
		*spid = atoi(argv[argc-2]);
		*nnode = atoi(argv[argc-1]);
		argc -= 3;
	}	
	if (argc < 5) return -1;
	if (parse_CnsOptions(argc - 4, argv, options) != ARG_PARSE_SUCCESS) return -1;
	*wrk_dir = argv[argc - 4];
	*candidates_path = argv[argc - 3];
	*cns_out_path = argv[argc - 2];
	*raw_out_path = argv[argc - 1];
	return 0;
}

int main(int argc, char* argv[])
{
	CnsOptions options;
	const char* wrk_dir = NULL;
	const char* candidates_path = NULL;
	const char* cns_out_path = NULL;
	const char* raw_out_path = NULL;
	int spid = 0;
	int nnode = 1;
	if (parse_params(argc, argv, &options, &wrk_dir, &candidates_path, &cns_out_path, &raw_out_path, &spid, &nnode) != 0) {
		print_usage(argv[0]);
		return 1;
	}

	PackedDB* reads = NULL;
	if (!options.small_memory) reads = merge_volumes(wrk_dir);
	int num_partitions = load_num_partitions(candidates_path);
	DFOPEN(cns_out, cns_out_path, "w");
	DFOPEN(raw_out, raw_out_path, "w");
	
	for (int i = spid; i < num_partitions; i += nnode) {
		consensus_one_partition(wrk_dir, reads, candidates_path, &options, cns_out, raw_out, i);
	}
	
	FCLOSE(cns_out);
	FCLOSE(raw_out);
	if (reads) free_PackedDB(reads);
}
