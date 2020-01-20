#include "map_aux.h"

BOOL
check_candidate_contain(vec_m4* m4list, GappedCandidate* can)
{
	for (size_t i = 0; i != kv_size(*m4list); ++i) {
		M4Record* m = kv_data(*m4list) + i;
		BOOL r = (can->qdir == m->qdir)
				 &&
				 (can->sid == m->sid)
				 &&
				 (can->qoff >= m->qoff && can->qoff <= m->qend)
				 &&
				 (can->soff >= m->soff && can->soff <= m->send);
		if (r) return TRUE;
	}
	return FALSE;
}

MappingThreadData*
new_MappingThreadData(int thread_id,
						MapOptions* options,
						PackedDB* reads,
						int read_start_id,
						PackedDB* reference,
						int reference_start_id,
						LookupTable* lktbl,
						FILE* out_file,
						pthread_mutex_t* out_lock,
						const int chunk_size,
						int* chunk_id,
						pthread_mutex_t* chunk_lock)
{
	MappingThreadData* mtd = (MappingThreadData*)malloc(sizeof(MappingThreadData));
	mtd->thread_id			= thread_id;
	mtd->options 			= options;
	mtd->reads 				= reads;
	mtd->read_start_id 		= read_start_id;
	mtd->reference 			= reference;
	mtd->reference_start_id = reference_start_id;
	mtd->lktbl 				= lktbl;
	mtd->chunk_id 			= chunk_id;
	mtd->chunk_size 		= chunk_size;
	mtd->chunk_lock 		= chunk_lock;
	mtd->output				= new_RecordWriter(out_file, out_lock, OcWriterBufferDefaultSize);
	
	return mtd;
}

MappingThreadData*
free_MappingThreadData(MappingThreadData* mtd)
{
	mtd->output = free_RecordWriter(mtd->output);
	free(mtd);
	return 0;
}

BOOL
get_next_read_chunk(const int chunk_size,
					int* chunk_id,
					pthread_mutex_t* chunk_lock,
					const int num_reads,
					int* sid,
					int* eid)
{
	int next_chunk_id;
	pthread_mutex_lock(chunk_lock);
	next_chunk_id = *chunk_id;
	++(*chunk_id);
	pthread_mutex_unlock(chunk_lock);
	int L = next_chunk_id * chunk_size;
	int R = L + chunk_size;
	if (L >= num_reads) return FALSE;
	*sid = L;
	*eid = (R > num_reads) ? num_reads : R;
	return TRUE;
}
