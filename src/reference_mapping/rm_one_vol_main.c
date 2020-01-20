#include "../common/map_options.h"
#include "../common/makedb_aux.h"
#include "../lookup_table/lookup_table.h"
#include "../common/map_aux.h"
#include "rm_worker.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s [OPTIONS] wrk-dir reference output\n", prog);
	fprintf(stderr, "\n");
	fprintf(stderr, "If Multiple Nodes Are Used:\n");
    fprintf(stderr, "%s [OPTIONS] wrk_dir reference output -mn node_id num_nodes", prog);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
	fprintf(stderr, "OPTIONS AND DESCRIPTIONS:\n");
	describe_MapOptions(&sDefaultPairwiseMapingOptions);
}

int main(int argc, char* argv[])
{
    TIMING_START(__func__);

    int svid = 0;
    int num_nodes = 1;
    if (argc >= 7 && strcmp(argv[argc-3], "-mn") == 0) {
        svid = atoi(argv[argc-2]);
        num_nodes = atoi(argv[argc-1]);
        argc -= 3;
    }

    MapOptions options = sDefaultReferenceMapingOptions;
	if (argc < 4 || (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
    const char* reads_path = argv[argc - 3];
	const char* reference_path = argv[argc - 2];
	const char* output = argv[argc - 1];

    PackedDB* reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_NANOPORE);
	LookupTable* lktbl = build_lookup_table(reference, options.kmer_size, options.kmer_cnt_cutoff, options.num_threads);
	const int reference_start_id = 0;
    OcMutex out_lock;
	pthread_mutex_init(&out_lock, NULL);
	OcMutex chunk_lock;
	pthread_mutex_init(&chunk_lock, NULL);
	const int num_threads = options.num_threads;
	pthread_t tids[num_threads];
	int read_start_id = 0;
	int chunk_id;
	DFOPEN(out, output, "w");
	MappingThreadData* mdata[num_threads];

    VolumesInfo* vis = load_volumes_info(reads_path);
	for (int i = svid; i < vis->num_volumes; i += num_nodes) {
		PackedDB* volume = new_PackedDB();
		pdb_load(volume, vi_volume_name(vis, i), TECH_NANOPORE);
        chunk_id = 0;
        oc_assert(i < kv_size(vis->read_start_id));
        read_start_id = kv_A(vis->read_start_id, i);
		for (int t = 0; t < num_threads; ++t) {
		    mdata[t] = new_MappingThreadData(t,
										   &options,
										   volume,
										   read_start_id,
										   reference,
										   reference_start_id,
										   lktbl,
										   out,
										   &out_lock,
										   500,
										   &chunk_id,
										   &chunk_lock);
            pthread_create(tids + t, NULL, rm_search_one_volume, mdata[t]);
        }
        for (int t = 0; t < num_threads; ++t) {
            pthread_join(tids[t], NULL);
            free_MappingThreadData(mdata[t]);
        }
		free_PackedDB(volume);
	}
	
	destroy_volumes_info(vis);
    FCLOSE(out);
    destroy_lookup_table(lktbl);
    free_PackedDB(reference);

    TIMING_END(__func__);
}
