#include "../common/makedb_aux.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static const size_t kVolSize = 2000000000;

int
pack_one_file(const char* path,
			  PackedDB* pdb,
			  const char* wrk_dir,
			  int* vid,
			  int* read_start_id,
			  FILE* vn_out)
{
	int num_reads = 0;
	idx num_bps = 0;
	idx cvs = PDB_SIZE(pdb);
	DGZ_OPEN(in, path, "r");
	kseq_t* read = kseq_init(in);
	while (kseq_read(read) >= 0) {
		pdb_add_one_seq(pdb, read, TECH_PACBIO);
		size_t s = kseq_size(*read);
		++num_reads;
		num_bps += s;
		cvs += s;
		if (cvs >= kVolSize) {
			new_kstring(vname);
			make_volume_name(wrk_dir, *vid, &vname);
			pdb_dump(pdb, kstr_str(vname));
			fprintf(vn_out, "%s\t%d\t%lu\n", kstr_str(vname), *read_start_id, PDB_NUM_SEQS(pdb));
			(*read_start_id) += PDB_NUM_SEQS(pdb);
			pdb_clear(pdb);
			cvs = 0;
			++(*vid);
			free_kstring(vname);
		}
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
	
	OC_LOG("pack %s: %d reads, %lu bps", path, num_reads, num_bps);
	return num_reads;
}

int
pack_one_list(const char* list_path,
			  PackedDB* pdb,
			  int* vid,
			  const char* wrk_dir,
			  int* read_start_id,
			  FILE* vn_out)
{
	DFOPEN(in, list_path, "r");
	char line[2048];
	int num_reads = 0;
	while (ontcns_getline(line, 2048, in)) {
		num_reads += pack_one_file(line, pdb, wrk_dir, vid, read_start_id, vn_out);
	}
	FCLOSE(in);
	return num_reads;
}

void
print_usage(const char* prog)
{
	fprintf(stderr, "USAGE:\n");
	fprintf(stderr, "%s wrk-dir file-list [file-list]\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc < 3) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* wrk_dir = argv[1];
	if (access(wrk_dir, F_OK) == -1) {
		if (mkdir(wrk_dir, 0755) == -1) {
			OC_ERROR("failed to create folder %s\n", wrk_dir);
		}
	}
	
	PackedDB pdb;
	init_packed_db(&pdb);
	pdb_enlarge_size(&pdb, kVolSize + 2000000);
	new_kstring(vnn);
	make_volume_name_name(wrk_dir, &vnn);
	DFOPEN(vn_out, kstr_str(vnn), "w");
	int vid = 0;
	int num_reads = 0;
	int read_start_id = 0;
	int c = 2;
	while (c < argc) {
		const char* file_list = argv[c];
		fprintf(stdout, "file_list: %s\n", file_list);
		num_reads += pack_one_list(file_list, &pdb, &vid, wrk_dir, &read_start_id, vn_out);
		++c;
	}
	
	if (PDB_SIZE(&pdb)) {
		new_kstring(vname);
		make_volume_name(wrk_dir, vid, &vname);
		pdb_dump(&pdb, kstr_str(vname));
		fprintf(vn_out, "%s\t%d\t%lu\n", kstr_str(vname), read_start_id, PDB_NUM_SEQS(&pdb));
		free_kstring(vname);
		++vid;
	}
	
	FCLOSE(vn_out);
	destroy_packed_db(&pdb);
	free_kstring(vnn);
	
	dump_reads_info(wrk_dir, vid, num_reads);
	
	return 0;
}
