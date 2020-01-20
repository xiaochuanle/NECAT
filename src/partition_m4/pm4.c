#include "pm4.h"

#include <assert.h>

#include "pm4_aux.h"
#include "pm4_options.h"
#include "../common/makedb_aux.h"
#include "../common/ontcns_aux.h"
#include "../common/record_reader.h"

#define fix_m4_offsets(c, nc, query_is_target) \
do { \
	nc.vscore = c.vscore; \
    nc.ident_perc = c.ident_perc; \
	\
	if (query_is_target) { \
		nc.qid = c.sid; \
		nc.qsize = c.ssize; \
		nc.sid = c.qid; \
		nc.ssize = c.qsize; \
		\
		if (c.qdir == REV) { \
			nc.qdir = ReverseStrand(c.sdir); \
			nc.qend = c.ssize - c.soff; \
			nc.qoff = c.ssize - c.send; \
			nc.qext = c.ssize - 1 - c.sext; \
			nc.sdir = FWD; \
			nc.send = c.qsize - c.qoff; \
			nc.soff = c.qsize - c.qend; \
			nc.sext = c.qsize - 1 - c.qext; \
		} else { \
			assert(c.qdir == FWD); \
			nc.qdir = c.sdir; \
			nc.qend = c.send; \
			nc.qoff = c.soff; \
			nc.qext = c.sext; \
			nc.sdir = c.qdir; \
			nc.send = c.qend; \
			nc.soff = c.qoff; \
			nc.sext = c.qext; \
		} \
	} else { \
		nc.qid = c.qid; \
		nc.qsize = c.qsize; \
		nc.sid = c.sid; \
		nc.ssize = c.ssize;	\
		\
		if (c.sdir == REV) {	\
			nc.qdir = ReverseStrand(c.qdir);	\
			nc.qend = c.qsize - c.qoff; \
			nc.qoff = c.qsize - c.qend; \
			nc.qext = c.qsize - 1 - c.qext; \
			nc.sdir = FWD; \
			nc.send = c.ssize - c.soff; \
			nc.soff = c.ssize - c.send; \
			nc.sext = c.ssize - 1 - c.sext; \
		} else { \
			assert(c.sdir == FWD); \
			nc = c; \
		} \
	} \
} while(0)

void
pm4_main(PM4Options* options, const char* wrk_dir, const char* m4_path)
{
	int num_reads = load_num_reads(wrk_dir);
	int num_batches = (num_reads + options->batch_size - 1) / options->batch_size;
	dump_num_partitions(m4_path, num_batches);
	M4RecordPartitionWriter* writer = new_M4RecordPartitionWriter(options->num_output_files);
	M4Record m4, nm4;
	char job[1024];
	
	for (int fid = 0; fid < num_batches; fid += options->num_output_files) {
		int sfid = fid;
		int efid = OC_MIN(sfid + options->num_output_files, num_batches);
		int nfid = efid - sfid;
		int ssid = sfid * options->batch_size;
		int esid = efid * options->batch_size;
		sprintf(job, "dumping records for partitions [%d, %d)", sfid, efid);
		TIMING_START(job);
		open_M4RecordPartitionWriter(nfid, m4_path, sfid, writer);
		BinRecordReader* m4in = new_BinRecordReader(m4_path, sizeof(M4Record));
		while (get_BinRecordReader(m4in, &m4)) {
			if (m4.qid >= ssid && m4.qid < esid) {
				fix_m4_offsets(m4, nm4, TRUE);
				int i = (m4.qid - ssid) / options->batch_size;
				dump_M4Record(writer, i, &nm4);
			}
			if (m4.sid >= ssid && m4.sid < esid) {
				fix_m4_offsets(m4, nm4, FALSE);
				int i = (m4.sid - ssid) / options->batch_size;
				dump_M4Record(writer, i, &nm4);
			}
		}
		close_M4RecordPartitionWriter(writer);
		free_BinRecordReader(m4in);
		TIMING_END(job);
	}
	
	free_M4RecordPartitionWriter(writer);
}
