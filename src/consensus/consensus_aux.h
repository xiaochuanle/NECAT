#ifndef CONSENSUS_AUX_H
#define CONSENSUS_AUX_H

#include "../common/packed_db.h"
#include "../common/gapped_candidate.h"
#include "../common/record_writer.h"
#include "../gapped_align/oc_daligner.h"
#include "../gapped_align/oc_aligner.h"
#include "../edlib/edlib_wrapper.h"
#include "../tasc/cbcns.h"
#include "cns_options.h"
#include "overlaps_pool.h"
#include "read_id_pool.h"

#define MAX_EXAMINED_CAN 	300

typedef struct {
	kstring_t			query;
	kstring_t			target;
	kstring_t			qabuf;
	kstring_t			tabuf;
	CnsSeq				cns_seq;
	PackedDB* 			reads;
	CnsReads*			cns_reads;
	vec_pcan*			candidates;
	size_t*				next_can_id;
	pthread_mutex_t*	can_id_lock;
	CnsOptions*			options;
	RecordWriter*		cns_out;
	RecordWriter*		raw_out;
	CbCnsData* 			cns_data;
	OcAlignData* 		align_data;
	OcDalignData*		dalign_data;
	FullEdlibAlignData*	falign_data;
	vec_int				cov_stats;
	OverlapsPool* 		op;
	ReadIdPool*			extended_read_ids;
	ReadIdPool*			corrected_read_ids;
	vec_intpair			cov_ranges;
	int*				thread_id;
	pthread_mutex_t*	thread_id_lock;
	size_t 				can_sid;
	size_t 				can_eid;
	int 				template_id;
	int 				template_size;
} CnsData;

CnsData*
new_CnsData(CnsReads* cns_reads,
			PackedDB* reads,
			vec_pcan* candidates,
			size_t* next_can_id,
			pthread_mutex_t* can_id_lock,
			CnsOptions* options,
			FILE* cns_out,
			pthread_mutex_t* cns_out_lock,
			FILE* raw_out,
			pthread_mutex_t* raw_out_lock,
		    ReadIdPool* corrected_read_ids,
			int* thread_id,
			pthread_mutex_t* thread_id_lock);

CnsData*
free_CnsData(CnsData* data);

BOOL
extract_candidate_range(CnsData* cns_data, size_t* can_sid, size_t* can_eid);

BOOL
is_full_cov_ovlp(int ql, int qr, int qs, int tl, int tr, int ts, int ovlp_size, int tail_size);

BOOL
check_mapping_range(int ql, int qr, int qs, int tl, int tr, int ts, int min_ovlp_size, double mratio);

BOOL
cns_extension(GappedCandidate* can, 
			  OcAlignData* align_data,
			  OcDalignData* dalign_data,
			  FullEdlibAlignData* falign_data,
			  const char* query,
			  const char* target,
			  int min_align_size,
			  BOOL rescue_long_indels,
			  int* qoff,
			  int* qend,
			  int* toff,
			  int* tend,
			  double* ident_perc,
			  kstring_t** qaln,
			  kstring_t** taln);

#endif // CONSENSUS_AUX_H
