#include "asm_pm_common.c"

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s [map options] wrk_dir volume_id output\n", prog);
}

int main(int argc, char* argv[])
{
	MapOptions options = sDefaultPairwiseMapingOptions;
	options.num_candidates = options.num_output = MAXC;
	if (argc < 4 || (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
	print_MapOptions(&options);
	
	const char* wrk_dir = argv[argc - 3];
	const int vid = atoi(argv[argc - 2]);
	const char* output = argv[argc - 1];
	out = fopen(output, "w");
	alloc_m4_sink();
	pthread_t jobs[options.num_threads];
	char job_name[1024];
	
	VolumesInfo* volumes = load_volumes_info(wrk_dir);
	reference_start_id = kv_A(volumes->read_start_id, vid);
	const char* reference_path = vi_volume_name(volumes, vid);
	reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_NANOPORE);
	lktbl = build_lookup_table(reference, options.kmer_size, options.kmer_cnt_cutoff, options.num_threads);
	seed_len = options.kmer_size;
	num_extended_can = options.num_candidates;
	binary_output = options.binary_output;
	BC = options.scan_window;
	
	for (int i = vid; i < volumes->num_volumes; ++i) {
		//if (i != 14) continue;
		sprintf(job_name, "pairwise mapping v%d vs v%d", vid, i);
		TIMING_START(job_name);
		const char* reads_path = vi_volume_name(volumes, i);
		reads = new_PackedDB();
		pdb_load(reads, reads_path, TECH_NANOPORE);
		read_start_id = kv_A(volumes->read_start_id, i);
		runnumber=0;
		readcount = PDB_NUM_SEQS(reads);
		terminalnum = (readcount + PLL - 1) / PLL;
		
		for (int k = 0; k < options.num_threads; ++k) {
			pthread_create(jobs + k, NULL, pairwise_mapping, NULL);
		}
		for (int k = 0; k < options.num_threads; ++k) {
			pthread_join(jobs[k], NULL);
		}
		dump_m4_sink();
		reads = free_PackedDB(reads);
		TIMING_END(job_name);
	}
	
	free_m4_sink();
	fclose(out);
	lktbl = destroy_lookup_table(lktbl);
	reference = free_PackedDB(reference);
	volumes = destroy_volumes_info(volumes);
}
