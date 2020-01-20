#ifndef PM4_AUX_H
#define PM4_AUX_H

#include "../klib/kstring.h"
#include "../common/m4_record.h"
#include "../common/record_writer.h"

void
make_partition_name(const char* prefix, const int pid, kstring_t* name);

void
make_partition_index_name(const char* prefix, kstring_t* name);

int 
load_num_partitions(const char* m4_path);

void
dump_num_partitions(const char* m4_path, const int np);

typedef struct {
	int nw;
	int mnw;
	FILE** file_list;
} AsmM4Writers;

void
pm4_main(const char* wrk_dir, 
		 const char* m4_path, 
		 const int num_dumpped_files, 
		 const int partition_size, 
		 const double min_ident_perc,
		 const int num_threads);

#endif // PCAN_AUX_H
