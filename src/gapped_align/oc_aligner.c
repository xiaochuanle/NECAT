#include "oc_aligner.h"

#include <assert.h>

#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "edlib_ex.h"

static const int kOcaMatCnt = 8;

OcAlignData*
new_OcAlignData(double _error)
{
	OcAlignData* data = (OcAlignData*)malloc(sizeof(OcAlignData));
	data->edlib = new_EdlibAlignData(MaxNumBlocks, MaxSeqSize);
	kstr_init(data->query_align);
	kstr_init(data->target_align);
	kstr_init(data->fqaln);
	kstr_init(data->rqaln);
	kstr_init(data->ftaln);
	kstr_init(data->rtaln);
	kstr_init(data->qfrag);
	kstr_init(data->tfrag);
	data->qabuf = (char*)malloc(100000);
	data->tabuf = (char*)malloc(100000);
	data->error = _error;
	
	return data;
}

OcAlignData*
free_OcAlignData(OcAlignData* data)
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
		free(data->qabuf);
		free(data->tabuf);
		free(data);
	}
	return 0;
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

static BOOL
get_next_sequence_block(const char* query,
						int qidx,
						const int qsize,
						const char* target,
						int tidx,
						const int tsize,
						const int desired_block_size,
						const BOOL right_extend,
						kstring_t* qfrag,
						kstring_t* tfrag)
{
	BOOL last_block = FALSE;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	int qblk;
	int tblk;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = tleft * 1.3;
		qblk = OC_MIN(qblk, qleft);
		tblk = qleft * 1.3;
		tblk = OC_MIN(tblk, tleft);
		last_block = TRUE;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = FALSE;
	}
	
	kstr_clear(*qfrag);
	kstr_clear(*tfrag);
	if (right_extend) {
		const char* Q = query + qidx;
		for (int i = 0; i < qblk; ++i) kputc(Q[i], qfrag);
		const char* R = target + tidx;
		for (int i = 0; i < tblk; ++i) kputc(R[i], tfrag);
	} else {
		const char* Q = query - qidx;
		for (int i = 0; i < qblk; ++i) kputc(Q[-i], qfrag);
		const char* R = target - tidx;
		for (int i = 0; i < tblk; ++i) kputc(R[-i], tfrag);
	}
	
	return last_block;
}

static void
oca_extend(const char* query,
		   const int query_size,
		   const char* target,
		   const int target_size,
		   EdlibAlignData* align_data,
		   const int block_size,
		   const double error,
		   const BOOL right_extend,
		   kstring_t* qfrag,
		   kstring_t* tfrag,
		   char* qabuf,
		   char* tabuf,
		   kstring_t* qaln,
		   kstring_t* taln,
           const int tail_match_len)
{
	kstr_clear(*qaln);
	kstr_clear(*taln);
	int qidx = 0, tidx = 0;
	
	//OC_LOG("query_size = %d, target_size = %d", query_size, target_size);
	
	while (1) {
		int qfae, tfae, qfrag_size, tfrag_size;
		BOOL last_block = get_next_sequence_block(query,
												  qidx,
												  query_size,
												  target,
												  tidx,
												  target_size,
												  block_size,
												  right_extend,
												  qfrag,
												  tfrag);
		qfrag_size = kstr_size(*qfrag);
		tfrag_size = kstr_size(*tfrag);
		if (qfrag_size == 0 || tfrag_size == 0) break;
		
		Edlib_align(kstr_str(*qfrag),
					qfrag_size,
					kstr_str(*tfrag),
					tfrag_size,
					align_data,
					error,
					qabuf,
					tabuf,
					&qfae,
					&tfae);
		
		//OC_LOG("qsize = %d, tsize = %d, qfae = %d, tfae = %d", qfrag_size, tfrag_size, qfae, tfae);
		
		//OC_LOG("validating align string");
		validate_aligned_string(0,
								kstr_str(*qfrag),
								0,
								qfae,
								qabuf,
								0,
								kstr_str(*tfrag),
								0,
								tfae,
								tabuf,
								strlen(qabuf),
								TRUE);
		
		BOOL done = last_block;
		int acnt = 0, qcnt = 0, tcnt = 0;
		if (qfrag_size - qfae > 30 && tfrag_size - tfae > 30) done = 1;
		const int M = done ? tail_match_len : kOcaMatCnt;
		int align_size = strlen(qabuf);
		int k = align_size - 1, m = 0;
		while (k >= 0) {
			const char qc = qabuf[k];
			const char tc = tabuf[k];
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == M) break;
			--k;
		}
		
		if (m != M || k < 1) {
			align_size = 0;
			for (int i = 0; i < qfrag_size && i < tfrag_size; ++i) {
				const char qc = kstr_A(*qfrag, i);
				const char tc = kstr_A(*tfrag, i);
				if (qc != tc) break;
				qabuf[align_size] = DecodeDNA(qc);
				tabuf[align_size] = DecodeDNA(tc);
				++align_size;
			}
			done = TRUE;
		} else {
			align_size -= acnt;
			qidx += (qfae - qcnt);
			tidx += (tfae - tcnt);
			if (done) align_size += M;
		}
	
		kputsn(qabuf, align_size, qaln);
		kputsn(tabuf, align_size, taln);
		if (done) break;
	}
	
	if (1) {
		int qend = 0, tend = 0;
		for (size_t i = 0; i != kstr_size(*qaln); ++i) {
			if (kstr_A(*qaln, i) != GAP_CHAR) ++qend;
			if (kstr_A(*taln, i) != GAP_CHAR) ++tend;
		}
		//OC_LOG("validating align string");
		validate_aligned_string(0,
							query,
							0,
							qend,
							kstr_str(*qaln),
							0,
							target,
							0,
							tend,
							kstr_str(*taln),
							kstr_size(*qaln),
							right_extend);
	}
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
onc_align(const char* query,
		  const int query_start,
		  const int query_size,
		  const char* target,
		  const int target_start,
		  const int target_size,
		  OcAlignData* oca_data,
		  const int block_size,
		  const int min_align_size,
          const int tail_match_len)
{
	kstring_t* query_align = &oca_data->query_align;
	kstring_t* target_align = &oca_data->target_align;
	kstr_clear(*query_align);
	kstr_clear(*target_align);
	int QS = query_start;
	int TS = target_start;
	BOOL right_extend;
	
	right_extend = FALSE;
	oca_extend(query + QS - 1,
			   QS,
			   target + TS -1,
			   TS,
			   oca_data->edlib,
			   block_size,
			   oca_data->error,
			   right_extend,
			   &oca_data->qfrag,
			   &oca_data->tfrag,
			   oca_data->qabuf,
			   oca_data->tabuf,
			   &oca_data->rqaln,
			   &oca_data->rtaln,
               tail_match_len);
	int rqcnt = 0, rtcnt = 0;
	{
		int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
		int raln_size = kstr_size(oca_data->rqaln);
		for (int i = 0; i < raln_size; ++i) {
			const char qc = kstr_A(oca_data->rqaln, i);
			const char tc = kstr_A(oca_data->rtaln, i);
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == kOcaMatCnt) break;
		}
		if (m == kOcaMatCnt) {
			QS -= qcnt;
			TS -= tcnt;
			for (int i = raln_size; i > acnt; --i) {
				const char qc = kstr_A(oca_data->rqaln, i - 1);
				kputc(qc, query_align);
				if (qc != GAP_CHAR) ++rqcnt;
				const char tc = kstr_A(oca_data->rtaln, i - 1);
				kputc(tc, target_align);
				if (tc != GAP_CHAR) ++rtcnt;
			}
		}
	}
	
	right_extend = TRUE;
	oca_extend(query + QS,
			   query_size - QS,
			   target + TS,
			   target_size - TS,
			   oca_data->edlib,
			   block_size,
			   oca_data->error,
			   right_extend,
			   &oca_data->qfrag,
			   &oca_data->tfrag,
			   oca_data->qabuf,
			   oca_data->tabuf,
			   &oca_data->fqaln,
			   &oca_data->ftaln,
               tail_match_len);
	int fqcnt = 0, ftcnt = 0;
	if (kstr_size(*query_align) == 0) {
		int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
		int faln_size = kstr_size(oca_data->fqaln);
		for (int i = 0; i < faln_size; ++i) {
			const char qc = kstr_A(oca_data->fqaln, i);
			const char tc = kstr_A(oca_data->ftaln, i);
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == kOcaMatCnt) break;
		}
		if (m == kOcaMatCnt) {
			acnt -= kOcaMatCnt;
			qcnt -= kOcaMatCnt;
			tcnt -= kOcaMatCnt;
			QS += qcnt;
			TS += tcnt;
			for (int i = acnt; i < faln_size; ++i) {
				const char qc = kstr_A(oca_data->fqaln, i);
				if (qc != GAP_CHAR) ++fqcnt;
				kputc(qc, query_align);
				const char tc = kstr_A(oca_data->ftaln, i);
				if (tc != GAP_CHAR) ++ftcnt;
				kputc(tc, target_align);
			}
		}
	} else {
		int faln_size = kstr_size(oca_data->fqaln);
		for (int i = 0; i < faln_size; ++i) {
			const char qc = kstr_A(oca_data->fqaln, i);
			if (qc != GAP_CHAR) ++fqcnt;
			kputc(qc, query_align);
			const char tc = kstr_A(oca_data->ftaln, i);
			if (tc != GAP_CHAR) ++ftcnt;
			kputc(tc, target_align);
		}
	}
	
	oca_data->qoff = QS - rqcnt;
	oca_data->qend = QS + fqcnt;
	oca_data->toff = TS - rtcnt;
	oca_data->tend = TS + ftcnt;
	
	//OC_LOG("validating align string");
	validate_aligned_string(oca_data->qid,
							query,
							oca_data->qoff,
							oca_data->qend,
							kstr_str(*query_align),
							oca_data->tid,
							target,
							oca_data->toff,
							oca_data->tend,
							kstr_str(*target_align),
							kstr_size(*target_align),
						    TRUE);
	int align_size = kstr_size(*target_align);
	oca_data->ident_perc = calc_ident_perc(kstr_str(*query_align), kstr_str(*target_align), align_size);
	
	return align_size >= min_align_size;
}
