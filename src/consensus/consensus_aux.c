#include "consensus_aux.h"

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
			pthread_mutex_t* thread_id_lock)
{
	CnsData* data = (CnsData*)malloc(sizeof(CnsData));
	kstr_init(data->query);
	kstr_init(data->target);
	kstr_init(data->qabuf);
	kstr_init(data->tabuf);
	init_CnsSeq(data->cns_seq);
	data->cns_reads = cns_reads;
	data->reads = reads;
	data->candidates = candidates;
	data->next_can_id = next_can_id;
	data->can_id_lock = can_id_lock;
	data->options = options;
	data->cns_out = new_RecordWriter(cns_out, cns_out_lock, 64 * (1<<20));
	data->raw_out = new_RecordWriter(raw_out, raw_out_lock, 64 * (1<<20));
	data->cns_data = new_CbCnsData();
	data->align_data = new_OcAlignData(options->error);
	data->dalign_data = new_OcDalignData(options->error);
	data->falign_data = new_FullEdlibAlignData(options->error);
	kv_init(data->cov_stats);
	data->op = new_OverlapsPool();
	data->extended_read_ids = new_ReadIdPool(NULL);
	data->corrected_read_ids = corrected_read_ids;
	kv_init(data->cov_ranges);
	data->thread_id = thread_id;
	data->thread_id_lock = thread_id_lock;
	
	return data; 
}

CnsData*
free_CnsData(CnsData* data)
{
	free_kstring(data->query);
	free_kstring(data->target);
	free_kstring(data->qabuf);
	free_kstring(data->tabuf);
	free_CnsSeq(data->cns_seq);
	free_RecordWriter(data->cns_out);
	free_RecordWriter(data->raw_out);
	free_CbCnsData(data->cns_data);
	free_OcAlignData(data->align_data);
	free_OcDalignData(data->dalign_data);
	free_FullEdlibAlignData(data->falign_data);
	kv_destroy(data->cov_stats);
	free_OverlapsPool(data->op);
	free_ReadIdPool(data->extended_read_ids);
	kv_destroy(data->cov_ranges);
	free(data);
	return 0;
}

BOOL
extract_candidate_range(CnsData* cns_data, size_t* can_sid, size_t* can_eid)
{
	const size_t ncan = kv_size(*cns_data->candidates);
	size_t sid = -1, eid = -1;
	int template_id;
	PackedGappedCandidate* cans = kv_data(*cns_data->candidates);
	pthread_mutex_lock(cns_data->can_id_lock);
	sid = *cns_data->next_can_id;
	if (sid < ncan) {
		template_id = pcan_sid(cans[sid]);
		for (eid = sid + 1; eid < ncan; ++eid) {
			if (pcan_sid(cans[eid]) != template_id) break;
		}
		*cns_data->next_can_id = eid;
	}
	pthread_mutex_unlock(cns_data->can_id_lock);
	
	*can_sid = sid;
	*can_eid = eid;
	return sid < ncan;
}

BOOL
is_full_cov_ovlp(int ql, int qr, int qs, int tl, int tr, int ts, int ovlp_size, int tail_size)
{
  const int L = ovlp_size;
  const int M = tail_size;

  if (ql <= M && qs - qr <= M) return TRUE;
  if (tl <= M && ts - tr <= M) return TRUE;

  if (qs - qr <= M) {
    if (tl > M) return FALSE;
    if (qr - ql >= L) return TRUE;
  }

  if (ts - tr <= M) {
    if (ql > M) return FALSE;
    if (qr - ql >= L) return TRUE;
  }

  return FALSE;
}

BOOL
check_mapping_range(int ql, int qr, int qs, int tl, int tr, int ts, int min_ovlp_size, double mratio)
{
	BOOL r = (qr - ql >= min_ovlp_size) || (tr - tl >= min_ovlp_size);
	if (r) return r;
	r = (qr - ql >= qs * mratio) || (tr - tl >= ts * mratio);
	return r;
}

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
			  kstring_t** taln)
{
	const BOOL small_edlib = onc_align(query,
					   can->qoff,
					   can->qsize,
					   target,
					   can->soff,
					   can->ssize,
					   align_data,
					   kOcaBlockSize,
					   min_align_size,
                       ONC_TAIL_MATCH_LEN_LONG);
	BOOL r = small_edlib;
	if (r) {
		int qbeg = oca_query_start(*align_data);
		int qend = oca_query_end(*align_data);
		int lhang = (qbeg > can->qbeg) ? (qbeg - can->qbeg) : 0;
		int rhang = (qend < can->qend) ? (can->qend - qend) : 0;
		if (lhang + rhang > 200) r = FALSE;
	}
	if (r) {
		*qoff = oca_query_start(*align_data);
		*qend = oca_query_end(*align_data);
		*toff = oca_target_start(*align_data);
		*tend = oca_target_end(*align_data);
		*qaln = oca_query_mapped_string(*align_data);
		*taln = oca_target_mapped_string(*align_data);
		*ident_perc = oca_ident_perc(*align_data);
		return TRUE;
	}
	
	if (rescue_long_indels) {
		const BOOL dalign = ocda_go(query, 
									can->qoff,
									can->qsize,
									target,
									can->soff,
									can->ssize,
									dalign_data,
									min_align_size);
		if (dalign) {
			const BOOL large_edlib = edlib_go(query,
											  ocda_query_start(*dalign_data),
											  ocda_query_end(*dalign_data),
											  target,
											  ocda_target_start(*dalign_data),
											  ocda_target_end(*dalign_data),
											  falign_data,
											  ocda_distance(*dalign_data),
											  min_align_size,
											  TRUE,
											  4);
			if (large_edlib) {
				*qoff = falign_data->qoff;
				*qend = falign_data->qend;
				*toff = falign_data->toff;
				*tend = falign_data->tend;
				*qaln = &falign_data->query_align;
				*taln = &falign_data->target_align;
				*ident_perc = falign_data->ident_perc;
				return TRUE;
			}
		}
	}
	
	if (small_edlib) {
		*qoff = oca_query_start(*align_data);
		*qend = oca_query_end(*align_data);
		*toff = oca_target_start(*align_data);
		*tend = oca_target_end(*align_data);
		*qaln = oca_query_mapped_string(*align_data);
		*taln = oca_target_mapped_string(*align_data);
		*ident_perc = oca_ident_perc(*align_data);
		return TRUE;
	}
	
	return FALSE;
}
