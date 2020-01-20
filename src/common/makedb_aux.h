#ifndef MAKEDB_AUX_H
#define MAKEDB_AUX_H

#include "../klib/kstring.h"
#include "../klib/kvec.h"
#include "../common/packed_db.h"

void
make_volume_name(const char* wk_dir, const int vid, kstring_t* name);

void
dump_reads_info(const char* wrk_dir, const int num_volumes, const int num_reads);

int
load_num_volumes(const char* wrk_dir);

int
load_num_reads(const char* wrk_dir);

void
make_volume_name_name(const char* wrk_dir, kstring_t* name);

typedef struct
{
	int num_reads;
	int num_volumes;
	vec_int	read_start_id;
	vec_int read_count;
	vec_size_type volume_name_starts;
	kstring_t volume_names;
} VolumesInfo;

const char*
vi_volume_name(VolumesInfo* vi, int vid);

VolumesInfo*
load_volumes_info(const char* wrk_dir);

VolumesInfo*
destroy_volumes_info(VolumesInfo* vi);

void
print_volume_info(VolumesInfo* vi);

PackedDB*
merge_volumes(const char* wrk_dir);

#endif // MAKEDB_AUX_H
