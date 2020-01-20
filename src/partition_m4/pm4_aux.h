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
	int max_output_files;
	int num_output_files;
	ResultsWriter** rlist;
} M4RecordPartitionWriter;

M4RecordPartitionWriter*
new_M4RecordPartitionWriter(const int max_output_files);

M4RecordPartitionWriter*
free_M4RecordPartitionWriter(M4RecordPartitionWriter* writer);

void
close_M4RecordPartitionWriter(M4RecordPartitionWriter* writer);

void
open_M4RecordPartitionWriter(const int num_output_files, const char* prefix, const int sfid, M4RecordPartitionWriter* writer);

void
dump_M4Record(M4RecordPartitionWriter* writer, const int fid, M4Record* m4);

#endif // PCAN_AUX_H
