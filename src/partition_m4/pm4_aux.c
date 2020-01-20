#include "pm4_aux.h"

#include <assert.h>

#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "../common/record_writer.h"

RESULT_WRITER_IMPL(M4Record, M4_RECORD, static)

void
make_partition_name(const char* prefix, const int pid, kstring_t* name)
{
	kstr_clear(*name);
	kputs(prefix, name);
	ksprintf(name, ".p%d", pid);
	kputc('\0', name);
}

void
make_partition_index_name(const char* prefix, kstring_t* name)
{
	kstr_clear(*name);
	kputs(prefix, name);
	kputs(".partitions", name);
	kputc('\0', name);
}

int 
load_num_partitions(const char* m4_path)
{
	new_kstring(path);
	make_partition_index_name(m4_path, &path);
	DFOPEN(in, kstr_str(path), "r");
	int n;
	fscanf(in, "%d", &n);
	FCLOSE(in);
	free_kstring(path);
	return n;
}

void
dump_num_partitions(const char* m4_path, const int np)
{
	new_kstring(path);
	make_partition_index_name(m4_path, &path);
	DFOPEN(out, kstr_str(path), "w");
	fprintf(out, "%d\n", np);
	FCLOSE(out);
	free_kstring(path);
}

M4RecordPartitionWriter*
new_M4RecordPartitionWriter(const int max_output_files)
{
	M4RecordPartitionWriter* writer = (M4RecordPartitionWriter*)malloc(sizeof(M4RecordPartitionWriter));
	writer->max_output_files = max_output_files;
	writer->num_output_files = 0;
	
	writer->rlist = (ResultsWriter**)malloc(sizeof(ResultsWriter*) * max_output_files);
	for (int i = 0; i < max_output_files; ++i) {
		writer->rlist[i] = new_ResultsWriter(NULL, NULL, TRUE, sizeof(M4Record), OcWriterBufferDefaultSize);
	}
	
	return writer;
}

M4RecordPartitionWriter*
free_M4RecordPartitionWriter(M4RecordPartitionWriter* writer)
{
	close_M4RecordPartitionWriter(writer);
	for (int i = 0; i < writer->max_output_files; ++i) {
		free_ResultsWriter(writer->rlist[i]);
	}
	free(writer);
	return 0;
}

void
close_M4RecordPartitionWriter(M4RecordPartitionWriter* writer)
{
	if (writer->num_output_files) {
		for (int i = 0; i < writer->num_output_files; ++i) {
			dump_ResultsWriter(writer->rlist[i]);
			FCLOSE(writer->rlist[i]->bin_writer->out);
			writer->rlist[i]->bin_writer->out = NULL;
			writer->rlist[i]->bin_writer->out_lock = NULL;
		}
		writer->num_output_files = 0;
	}
}

void
open_M4RecordPartitionWriter(const int num_output_files, const char* prefix, const int sfid, M4RecordPartitionWriter* writer)
{
	close_M4RecordPartitionWriter(writer);
	new_kstring(name);
	for (int i = 0; i < num_output_files; ++i) {
		make_partition_name(prefix, sfid + i, &name);
		const char* name_str = kstr_str(name);
		FOPEN(writer->rlist[i]->bin_writer->out, name_str, "w");
	}
	writer->num_output_files = num_output_files;
	free_kstring(name);
}

void
dump_M4Record(M4RecordPartitionWriter* writer, const int fid, M4Record* m4)
{
	oc_assert(fid < writer->num_output_files, "fid = %d, num_output_files = %d", fid, writer->num_output_files);
	RW_DUMP_ONE_DATA(M4_RECORD, writer->rlist[fid], m4) ;
}
