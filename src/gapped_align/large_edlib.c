#include "large_edlib.h"
#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "edlib_ex.h"

static const int kMatchSize = 8;

LargeEdlibData*
new_LargeEdlibData(const double _error, int _maxNumBlock, int _maxSeqSize)
{
	LargeEdlibData* data = (LargeEdlibData*)malloc(sizeof(LargeEdlibData));
	data->edlib = new_EdlibAlignData(_maxNumBlock, _maxSeqSize);
	kstr_init(data->query_align);
	kstr_init(data->target_align);
	kstr_init(data->fqaln);
	kstr_init(data->rqaln);
	kstr_init(data->ftaln);
	kstr_init(data->rtaln);
	kstr_init(data->qfrag);
	kstr_init(data->tfrag);
	data->error = _error;
	
	return data;
}

void 
free_LargeEdlibData(LargeEdlibData* data)
{
	if (data) {
		if (data->edlib) data->edlib = free_EdlibAlignData(data->edlib);
		free_kstring(data->query_align);
		free_kstring(data->target_align);
		free_kstring(data->fqaln);
		free_kstring(data->rqaln);
		free_kstring(data->ftaln);
		free_kstring(data->rtaln);
		free_kstring(data->qfrag);
		free_kstring(data->tfrag);
		free(data);
	}
}

static void
validate_aligned_string(int qid,
						const char* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const char* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
			const char qc1 = DecodeDNA(right_extend ? query[x] : query[-x]);
			oc_assert(qc == qc1, "qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu",
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
			
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
			const char tc1 = DecodeDNA(right_extend ? target[y] : target[-y]);
			oc_assert(tc == tc1, "qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d",
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  tc,
					  tc1,
					  qoff,
					  qend,
					  toff,
					  tend);

			++y;
		}
	}
}

static void
calc_seq_size(int qsize_in,
			  int tsize_in,
			  int* qsize_out,
			  int* tsize_out)
{
	int A = tsize_in + 200;
	int B = tsize_in * 1.1;
	int C = OC_MIN(A, B);
	*qsize_out = OC_MIN(qsize_in, C);
	
	A = qsize_in + 200;
	B = qsize_in * 1.1;
	C = OC_MIN(A, B);
	*tsize_out = OC_MIN(tsize_in, C);
}

static void
trim_aligned_bases(const char* qaln,
				   const char* taln,
				   int* from,
				   int* to)
{
	*from = 0; 
	*to = 0;
	
	int alns = strlen(qaln);
	int k = alns - 1, m = 0, acnt = 0;
	while (k >= 0) {
		char qc = qaln[k];
		char tc = taln[k];
		if (qc == tc) {
			++m;
		} else {
			m = 0;
		}
		++acnt;
		if (m == kMatchSize) break;
		--k;
	}
	
	if (m != kMatchSize || k < 1) {
		acnt = 0;
		for (k = 0; k < alns; ++k) {
			if (qaln[k] != taln[k]) break;
			++acnt;
		}
		*to = acnt;
		return;
	} else {
		*to = alns - acnt + kMatchSize;
	}
	
	acnt = 0;
	m = 0;
	for (k = 0; k < alns; ++k) {
		if (qaln[k] == taln[k]) {
			++m;
		} else {
			m = 0;
		}
		++acnt;
		if (m == kMatchSize) break;
	}
	*from = acnt - kMatchSize;
}

static double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
	return 100.0 * n / align_size;
}

BOOL
ledlib_go(const char* query,
		  const int query_start,
		  const int query_size,
		  const char* target,
		  const int target_start,
		  const int target_size,
		  LargeEdlibData* oca_data,
		  const int min_align_size)
{
	alloc_space_EdlibAlignData(oca_data->edlib, query_size, target_size);
	
	int QS = query_start;
	int TS = target_start;
	kstr_clear(oca_data->query_align);
	kstr_clear(oca_data->target_align);

	int rqcnt = 0, rtcnt = 0;
	{
		int qblk, tblk;
		calc_seq_size(QS, TS, &qblk, &tblk);
		kstr_clear(oca_data->qfrag);
		for (int i = 0; i < qblk; ++i) {
			kputc(query[QS - 1 - i], &oca_data->qfrag);
		}
		kstr_clear(oca_data->tfrag);
		for (int i = 0; i < tblk; ++i) {
			kputc(target[TS - 1 - i], &oca_data->tfrag);
		}
		
		ks_reserve(&oca_data->rqaln, qblk * 2);
		ks_reserve(&oca_data->rtaln, tblk * 2);
		Edlib_align(kstr_str(oca_data->qfrag),
					qblk,
					kstr_str(oca_data->tfrag),
					tblk,
					oca_data->edlib,
					oca_data->error,
					kstr_str(oca_data->rqaln),
					kstr_str(oca_data->rtaln),
					&qblk,
					&tblk);
		
		int from, to;
		trim_aligned_bases(kstr_str(oca_data->rqaln), kstr_str(oca_data->rtaln), &from, &to);
		if (to - from > kMatchSize) {
			from += kMatchSize;
		} else {
			from = to;
		}
		for (int i = 0; i < from; ++i) {
			char qc = kstr_A(oca_data->rqaln, i);
			if (qc != GAP_CHAR) --QS;
			char tc = kstr_A(oca_data->rtaln, i);
			if (tc != GAP_CHAR) --TS;
		}
		
		for (int i = to; i > from; --i) {
			char qc = kstr_A(oca_data->rqaln, i - 1);
			kputc(qc, &oca_data->query_align);
			if (qc != GAP_CHAR) ++rqcnt;
			char tc = kstr_A(oca_data->rtaln, i - 1);
			kputc(tc, &oca_data->target_align);
			if (tc != GAP_CHAR) ++rtcnt;
		}
	}
	
	int fqcnt = 0, ftcnt = 0;
	{
		int qblk, tblk;
		calc_seq_size(query_size - QS, target_size - TS, &qblk, &tblk);
		ks_reserve(&oca_data->fqaln, qblk * 2);
		ks_reserve(&oca_data->ftaln, tblk * 2);
		Edlib_align(query + QS,
					qblk,
					target + TS,
					tblk,
					oca_data->edlib,
					oca_data->error,
					kstr_str(oca_data->fqaln),
					kstr_str(oca_data->ftaln),
					&qblk,
					&tblk);
		int from, to;
		trim_aligned_bases(kstr_str(oca_data->fqaln), kstr_str(oca_data->ftaln), &from, &to);
		for (int i = 0; i < to; ++i) {
			char qc = kstr_A(oca_data->fqaln, i);
			kputc(qc, &oca_data->query_align);
			if (qc != GAP_CHAR) ++fqcnt;
			char tc = kstr_A(oca_data->ftaln, i);
			kputc(tc, &oca_data->target_align);
			if (tc != GAP_CHAR) ++ftcnt;
		}
	}
	
	oca_data->qoff = QS - rqcnt;
	oca_data->qend = QS + fqcnt;
	oca_data->toff = TS - rtcnt;
	oca_data->tend = TS + ftcnt;
	
	validate_aligned_string(0,
							query,
							oca_data->qoff,
							oca_data->qend,
							kstr_str(oca_data->query_align),
							0,
							target,
							oca_data->toff,
							oca_data->tend,
							kstr_str(oca_data->target_align),
							kstr_size(oca_data->target_align),
						    TRUE);
	
	if (fqcnt + rqcnt < min_align_size || ftcnt + rtcnt < min_align_size) return FALSE;
	
	oca_data->ident_perc = calc_ident_perc(kstr_str(oca_data->query_align), kstr_str(oca_data->target_align), kstr_size(oca_data->target_align));
	return TRUE;
}
