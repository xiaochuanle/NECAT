#ifndef CONSENSUS_ONE_READ_H
#define CONSENSUS_ONE_READ_H

#include "consensus_aux.h"

void
consensus_one_read(CnsData* cns_data, size_t can_sid, size_t can_eid);

void
consensus_one_read_m4(PackedDB* reads, 
					  CnsReads* cns_reads,
                      kstring_t* target,
                      kstring_t* read,
                      M4Record* m4v,
                      int nm4,
					  CbCnsData* cbcns_data,
					  vec_int* cov_stats,
					  FullEdlibAlignData* align_data,
					  int trg_from,
					  int trg_to,
					  int* read_id,
					  FILE* out,
					  pthread_mutex_t* out_lock);

#endif // CONSENSUS_ONE_READ_H
