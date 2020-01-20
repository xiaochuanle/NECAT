#include "makedb_aux.h"

#include "../common/ontcns_aux.h"

static void
copy_wrk_dir_name(const char* wrk_dir, kstring_t* name)
{
	ksprintf(name, "%s", wrk_dir);
	if (wrk_dir[ strlen(wrk_dir) - 1 ] != '/') ksprintf(name, "%c", '/');
}

void
make_volume_name(const char* wrk_dir, const int vid, kstring_t* name)
{
	copy_wrk_dir_name(wrk_dir, name);
	ksprintf(name, "vol%d", vid);
}

void
make_reads_info_name(const char* wrk_dir, kstring_t* name)
{
	copy_wrk_dir_name(wrk_dir, name);
	ksprintf(name, "reads_info.txt");
}

void
make_volume_name_name(const char* wrk_dir, kstring_t* name)
{
	copy_wrk_dir_name(wrk_dir, name);
	ksprintf(name, "volume_names.txt");
}

void
dump_reads_info(const char* wrk_dir, const int num_volumes, const int num_reads)
{
	new_kstring(name);
	make_reads_info_name(wrk_dir, &name);
	DFOPEN(out, kstr_str(name), "w");
	fprintf(out, "%d\t%d\n", num_volumes, num_reads);
	FCLOSE(out);
	free_kstring(name);
}

int
load_num_volumes(const char* wrk_dir)
{
	new_kstring(name);
	make_reads_info_name(wrk_dir, &name);
	DFOPEN(in, kstr_str(name), "r");
	int num_volumes;
	SAFE_SCANF(fscanf, in, 1, "%d", &num_volumes);
	FCLOSE(in);
	free_kstring(name);
	return num_volumes;
}

int
load_num_reads(const char* wrk_dir)
{
	new_kstring(name);
	make_reads_info_name(wrk_dir, &name);
	DFOPEN(in, kstr_str(name), "r");
	int num_volumes, num_reads;
	SAFE_SCANF(fscanf, in, 2, "%d%d", &num_volumes, &num_reads);
	FCLOSE(in);
	free_kstring(name);
	return num_reads;
}

const char*
vi_volume_name(VolumesInfo* vi, int vid)
{
	size_t vn_start = kv_A(vi->volume_name_starts, vid);
	return kstr_str(vi->volume_names) + vn_start;
}

VolumesInfo*
load_volumes_info(const char* wrk_dir)
{
	VolumesInfo* vi = (VolumesInfo*)malloc(sizeof(VolumesInfo));
	kv_init(vi->read_start_id);
	kv_init(vi->read_count);
	kv_init(vi->volume_name_starts);
	kstr_init(vi->volume_names);
	vi->num_reads = load_num_reads(wrk_dir);
	vi->num_volumes = load_num_volumes(wrk_dir);
	
	new_kstring(vnn);
	char line[2048];
	size_t vn_start = 0;
	make_volume_name_name(wrk_dir, &vnn);
	DFOPEN(in, kstr_str(vnn), "r");
	int read_start_id, read_count;
	for (int i = 0; i < vi->num_volumes; ++i) {
		ontcns_getline(line, 2048, in);
		vn_start = kstr_size(vi->volume_names);
		kv_push(size_t, vi->volume_name_starts, vn_start);
		
		size_t k = 0, n = strlen(line);
		while (k < n && (!isspace(line[k]))) ++k;
		kputsn(line, k, &vi->volume_names);
		kputc('\0', &vi->volume_names);
		++k;
		read_start_id = atoi(line + k);
		while (k < n && (!isspace(line[k]))) ++k;
		++k;
		read_count = atoi(line + k);
		kv_push(int, vi->read_start_id, read_start_id);
		kv_push(int, vi->read_count, read_count);
	}
	FCLOSE(in);
	free_kstring(vnn);
	return vi;
}

VolumesInfo*
destroy_volumes_info(VolumesInfo* vi)
{
	kv_destroy(vi->read_start_id);
	kv_destroy(vi->read_count);
	kv_destroy(vi->volume_name_starts);
	free_kstring(vi->volume_names);
	free(vi);
	return 0;
}

void
print_volume_info(VolumesInfo* vi)
{
	fprintf(stdout, "number of reads: %d\n", vi->num_reads);
	fprintf(stdout, "number of volumes: %d\n", vi->num_volumes);
	for (int i = 0; i < vi->num_volumes; ++i) {
		fprintf(stdout, "%s\t%d\t%d\n", vi_volume_name(vi, i), kv_A(vi->read_start_id, i), kv_A(vi->read_count, i));
	}
}

PackedDB*
merge_volumes(const char* wrk_dir)
{
	VolumesInfo* vis = load_volumes_info(wrk_dir);
	PackedDB* reads = new_PackedDB();
	pdb_load(reads, vi_volume_name(vis, 0), TECH_NANOPORE);
	
	for (int i = 1; i < vis->num_volumes; ++i) {
		PackedDB* volume = new_PackedDB();
		pdb_load(volume, vi_volume_name(vis, i), TECH_NANOPORE);
		pdb_merge(volume, reads);
		free_PackedDB(volume);
	}
	
	destroy_volumes_info(vis);
	return reads;
}
