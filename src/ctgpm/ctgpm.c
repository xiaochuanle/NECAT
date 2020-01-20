#include "ctgpm.h"

#include "../common/oc_assert.h"

#include <assert.h>

static void
extract_subsequence(PackedDB* pctg, int ctg_id, int ctg_dir, idx from, idx size, const int right_extend, kstring_t* seq)
{
	kstr_clear(*seq);
	idx ctg_offset = PDB_SEQ_OFFSET(pctg, ctg_id);
	idx ctg_size = PDB_SEQ_SIZE(pctg, ctg_id);
	for (size_t i = 0; i < size; ++i) {
		idx k = right_extend ? (from + i) : (from - i);
		if (ctg_dir == FWD) {
			k = ctg_offset + k;
			char c = _get_pac(pctg->m_pac, k);
			kputc(c, seq);
		} else {
			assert(ctg_dir == REV);
			k = ctg_offset + (ctg_size - 1 - k);
			char c = _get_pac(pctg->m_pac, k);
			c = 3 - c;
			kputc(c, seq);
		}
	}
}

static BOOL
ctgpm_extension(OcDalignData* dalign_data,
			  FullEdlibAlignData* falign_data,
			  const char* query,
			  const int query_size,
			  const char* target,
			  const int target_size,
			  int min_align_size,
			  int* qend,
			  int* tend,
			  int* dist,
			  int* reach_end)
{
	*qend = 0;
	*tend = 0;
	*dist = 0;
	
	const BOOL dalign = ocda_go(query, 
								0,
								query_size,
								target,
								0,
								target_size,
								dalign_data,
								min_align_size);
	
	int dqend = ocda_query_end(*dalign_data);
	int dtend = ocda_target_end(*dalign_data);
	//printf("qe = %d, qs = %d, te = %d, ts = %d, dist = %d\n", dqend, query_size, dtend, target_size, ocda_distance(*dalign_data));
	*reach_end = (dqend == query_size) || (dtend == target_size);
	
	if (dalign) {
		const BOOL large_edlib = edlib_go(query,
										0,
										dqend,
										target,
										0,
										dtend,
										falign_data,
										ocda_distance(*dalign_data),
										min_align_size,
										TRUE,
									    8);
		if (large_edlib) {
			*qend = falign_data->qend;
			*tend = falign_data->tend;
			*dist = falign_data->dist;
			return TRUE;
		}
	}
	
	return FALSE;
}

static BOOL
calc_block_size(idx qleft, idx tleft, idx* qsize, idx* tsize)
{
	const idx block_size = 10000;
	const idx max_block_size = 15000;
	BOOL last_block = 0;
	if (qleft <= max_block_size || tleft <= max_block_size) {
		last_block = 1;
		idx x = tleft + 500;
		*qsize = OC_MIN(qleft, x);
		x = qleft + 500;
		*tsize = OC_MIN(x, tleft);
	} else {
		*qsize = block_size;
		*tsize = block_size;
	}
	return last_block;
}

static void
extend(PackedDB* pctg,
	   OcDalignData* dalign_data,
	   FullEdlibAlignData* falign_data,
	   kstring_t* query,
	   kstring_t* target,
	   int qid,
	   int qdir,
	   idx qstart,
	   idx qsize,
	   int tid,
	   int tdir,
	   idx tstart,
	   idx tsize,
	   idx* qend,
	   idx* tend,
	   idx* dist,
	   int right_extend)
{
	idx qidx = 0, tidx = 0;
	idx qleft = qsize, tleft = tsize;
	*dist = 0;
	while (1) {
		idx qblk, tblk;
		BOOL is_done = calc_block_size(qleft, tleft, &qblk, &tblk);
		if (qblk == 0 || tblk == 0) break;
		extract_subsequence(pctg, qid, qdir, qstart, qblk, right_extend, query);
		extract_subsequence(pctg, tid, tdir, tstart, tblk, right_extend, target);
		int qe = -1, te = -1, d = 0, full_align = 0;
		BOOL r = ctgpm_extension(dalign_data, falign_data, kstr_str(*query), qblk, kstr_str(*target), tblk, 1, &qe, &te, &d, &full_align);
		if (!r) break;
		if (right_extend) {
			qstart += qe;
			tstart += te;
		} else {
			qstart -= qe;
			tstart -= te;
		}
		qidx += qe;
		tidx += te;
		qleft -= qe;
		tleft -= te;
		*dist += d;
		is_done = is_done || (!full_align);
		if (is_done) break;
	}
	*qend = qidx;
	*tend = tidx;
}

BOOL
ctgmp_align(PackedDB* pctg,
			OcDalignData* dalign_data,
			FullEdlibAlignData* falign_data,
			GappedCandidate* _can,
			kstring_t* query,
			kstring_t* target,
			idx* qoff,
			idx* qend,
			idx* toff,
			idx* tend,
			double* ident_perc,
		    const int min_align_size)
{
	GappedCandidate can = *_can;
	
	idx ldist = 0;
	idx lqext, ltext;
	extend(pctg,
		   dalign_data,
		   falign_data,
		   query,
		   target,
		   can.qid,
		   can.qdir,
		   can.qoff - 1,
		   can.qoff,
		   can.sid,
		   can.sdir,
		   can.soff - 1,
		   can.soff,
		   &lqext,
		   &ltext,
		   &ldist,
		   0);
	
	idx rdist = 0;
	idx rqext, rtext;
	extend(pctg,
		   dalign_data,
		   falign_data,
		   query,
		   target,
		   can.qid,
		   can.qdir,
		   can.qoff,
		   can.qsize - can.qoff,
		   can.sid,
		   can.sdir,
		   can.soff,
		   can.ssize - can.soff,
		   &rqext,
		   &rtext,
		   &rdist,
		   1);
	
	*qoff = can.qoff - lqext;
	*qend = can.qoff + rqext;
	*toff = can.soff - ltext;
	*tend = can.soff + rtext;
	
	oc_assert(*qoff <= *qend);
	oc_assert(*qend <= can.qsize);
	oc_assert(*toff <= *tend);
	oc_assert(*tend <= can.ssize);
	
	idx qas = *qend - *qoff;
	idx tas = *tend - *toff;
	if (qas < min_align_size || tas < min_align_size) return 0;
	
	*ident_perc = 100.0 - 200.0 *(ldist + rdist) / (qas + tas);
	return 1;
}
