#ifndef MAP_AUX_H
#define MAP_AUX_H

#include "gapped_candidate.h"
#include "m4_record.h"
#include "map_options.h"
#include "packed_db.h"
#include "record_writer.h"
#include "../lookup_table/lookup_table.h"

typedef struct {
	int					thread_id;
	MapOptions*			options;
	PackedDB*			reads;
	int					read_start_id;
	PackedDB*			reference;
	int					reference_start_id;
	LookupTable*		lktbl;
	RecordWriter*		output;
	int					chunk_size;
	int*				chunk_id;
	pthread_mutex_t*	chunk_lock;
} MappingThreadData;

BOOL
get_next_read_chunk(const int chunk_size,
					int* chunk_id,
					pthread_mutex_t* chunk_lock,
					const int num_reads,
					int* sid,
					int* eid);

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
						pthread_mutex_t* chunk_lock);

MappingThreadData*
free_MappingThreadData(MappingThreadData* mtd);

BOOL
check_candidate_contain(vec_m4* m4list, GappedCandidate* can);

#endif // MAP_AUX_H
