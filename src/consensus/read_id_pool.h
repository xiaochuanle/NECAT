#ifndef READ_ID_POOL_H
#define READ_ID_POOL_H

#include "../common/ontcns_defs.h"
#include "../common/packed_db.h"
#include "../common/gapped_candidate.h"
#include "../common/m4_record.h"

typedef struct {
	void* read_id_pool;
	pthread_mutex_t* id_lock;
} ReadIdPool;

ReadIdPool*
new_ReadIdPool(OcMutex* id_lock);

ReadIdPool*
free_ReadIdPool(ReadIdPool* pool);

void
clear_ReadIdPool(ReadIdPool* pool);

BOOL
read_id_exists(ReadIdPool* pool, int id);

void
add_read_id(ReadIdPool* pool, int id);

typedef struct {
	void* cns_read_ids;
	PackedDB* cns_reads;
} CnsReads;

CnsReads*
cns_reads_new(PackedGappedCandidate* cans, 
	const int ncan,
	const int max_can_per_template,
	const char* pac_reads_dir);

CnsReads*
cns_reads_free(CnsReads* cns_reads);

int
cns_reads_local_id(CnsReads* cns_reads, int global_id);

void
cns_reads_extract_sequence(CnsReads* cns_reads, int global_id, int strand, kstring_t* seq);

int
cns_reads_seq_size(CnsReads* cns_reads, int global_id);

const char*
cns_reads_seq_hdr(CnsReads* cns_reads, int global_id);

void 
truncate_m4_list(M4Record* m4v, int* nm4, const int max_m4_perc_read, const double ident_perc);

CnsReads*
cns_reads_new_m4(M4Record* m4s,
	const int nm4,
	const int max_m4_per_template,
	const double min_ident_perc,
	const char* pac_reads_dir);

#endif // READ_ID_POOL_H
