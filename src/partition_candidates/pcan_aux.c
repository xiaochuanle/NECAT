#include "pcan_aux.h"

#include <assert.h>

#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "../common/record_writer.h"

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
load_num_partitions(const char* candidates_path)
{
	new_kstring(path);
	make_partition_index_name(candidates_path, &path);
	DFOPEN(in, kstr_str(path), "r");
	int n;
	SAFE_SCANF(fscanf, in, 1, "%d", &n);
	FCLOSE(in);
	free_kstring(path);
	return n;
}

void
dump_num_partitions(const char* candidates_path, const int np)
{
	new_kstring(path);
	make_partition_index_name(candidates_path, &path);
	DFOPEN(out, kstr_str(path), "w");
	fprintf(out, "%d\n", np);
	FCLOSE(out);
	free_kstring(path);
}

CandidatePartitionWriter*
new_CandidatePartitionWriter(const int max_output_files)
{
	CandidatePartitionWriter* writer = (CandidatePartitionWriter*)malloc(sizeof(CandidatePartitionWriter));
	writer->max_output_files = max_output_files;
	writer->num_output_files = 0;
	
	writer->file_list = (FILE**)malloc(sizeof(FILE*) * max_output_files);
	for (int i = 0; i < max_output_files; ++i) {
		writer->file_list[i] = NULL;
	}
	
	return writer;
}

CandidatePartitionWriter*
free_CandidatePartitionWriter(CandidatePartitionWriter* writer)
{
	close_CandidatePartitionWriter(writer);
	free(writer->file_list);
	free(writer);
	return 0;
}

void
close_CandidatePartitionWriter(CandidatePartitionWriter* writer)
{
	if (writer->num_output_files) {
		for (int i = 0; i < writer->num_output_files; ++i) {
			FCLOSE(writer->file_list[i]);
		}
		writer->num_output_files = 0;
	}
}

void
open_CandidatePartitionWriter(const int num_output_files, const char* prefix, const int sfid, CandidatePartitionWriter* writer)
{
	close_CandidatePartitionWriter(writer);
	new_kstring(name);
	for (int i = 0; i < num_output_files; ++i) {
		make_partition_name(prefix, sfid + i, &name);
		const char* name_str = kstr_str(name);
		FOPEN(writer->file_list[i], name_str, "w");
	}
	writer->num_output_files = num_output_files;
	free_kstring(name);
}
