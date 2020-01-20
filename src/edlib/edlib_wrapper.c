#include "edlib_wrapper.h"

#include "../klib/kalloc.h"
#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"

FullEdlibAlignData*
new_FullEdlibAlignData(double error)
{
	FullEdlibAlignData* data = (FullEdlibAlignData*)malloc( sizeof(FullEdlibAlignData) );
	kstr_init(data->query);
	kstr_init(data->target);
	kstr_init(data->query_align);
	kstr_init(data->target_align);
	kstr_init(data->qAln);
	kstr_init(data->tAln);
	data->error = error;
	data->km = km_init();
	return data;
}

void
free_FullEdlibAlignData(FullEdlibAlignData* data)
{
	free_kstring(data->query);
	free_kstring(data->target);
	free_kstring(data->query_align);
	free_kstring(data->target_align);
	free_kstring(data->qAln);
	free_kstring(data->tAln);
	km_destroy(data->km);
	free(data);
}


static double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size,
			    int* dist)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
	*dist = align_size - n;
	return 100.0 * n / align_size;
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

int
edlib_go(const char* query,
		const int query_from,
		const int query_to,
		const char* target,
		const int target_from,
		const int target_to,
		FullEdlibAlignData* align_data,
		const int tolerance,
		const int min_align_size,
	    const BOOL find_path,
	    const int kMatchSize)
{
	kstr_clear(align_data->query);
	for ( int i = query_from; i < query_to; ++i) {
		char c = query[i];
		kputc("ACGT"[c], &align_data->query );
	}
	kstr_clear(align_data->target);
	for (int i = target_from; i < target_to; ++i) {
		char c = target[i];
		kputc("ACGT"[c], &align_data->target);
	}
	
	int query_size = query_to - query_from;
	int target_size = target_to - target_from;
	EdlibAlignTask task = find_path ? EDLIB_TASK_PATH : EDLIB_TASK_DISTANCE;
	EdlibAlignResult align = edlibAlign(kstr_str(align_data->query),            
										query_to - query_from,
                                        kstr_str(align_data->target),
										target_to - target_from,
                                        edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0),
									    align_data->km);
	if (align.numLocations == 0) return 0;
	if (!find_path) {
		double ident = 200.0 * align.editDistance / (query_size + target_size);
		ident = 100.0 - ident;
		align_data->qoff = query_from;
		align_data->qend = query_to;
		align_data->toff = target_from;
		align_data->tend = target_to;
		align_data->ident_perc = ident;
		edlibFreeAlignResult(align, align_data->km);
		return 1;
	}
	
	int  alignLen  = align.endLocations[0] - align.startLocations[0];
    double alignDiff = align.editDistance / (double)alignLen;
	if (alignLen < min_align_size) return 0;
	if (alignDiff > align_data->error) return 0;
	
    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;

	ks_resize(&align_data->qAln, align.alignmentLength + 1);
	ks_resize(&align_data->tAln, align.alignmentLength + 1);
    edlibAlignmentToStrings(align.alignment,
                            align.alignmentLength,
                            tBgn,
                            tEnd,
                            qBgn,
                            qEnd,
                            kstr_str(align_data->target),
                            kstr_str(align_data->query),
                            kstr_str(align_data->tAln),
                            kstr_str(align_data->qAln));
	
	const char* qAln = kstr_str(align_data->qAln);
	const char* tAln = kstr_str(align_data->tAln);
	int aln_from = 0, pqcnt = 0, ptcnt = 0;
	int m = 0;
	while (aln_from < align.alignmentLength) {
		if (qAln[aln_from] == tAln[aln_from]) {
			++m;
		} else {
			m = 0;
		}
		if (qAln[aln_from] != '-') ++pqcnt;
		if (tAln[aln_from] != '-') ++ptcnt;
		++aln_from;
		if (m == kMatchSize) break;
	}
	if (m == kMatchSize) {
		aln_from -= kMatchSize;
		pqcnt -= kMatchSize;
		ptcnt -= kMatchSize;
	} else {
		return 0;
	}
	
	int aln_to = align.alignmentLength, tqcnt = 0, ttcnt = 0;
	m = 0;
	while (aln_to) {
		if (qAln[aln_to - 1] == tAln[aln_to - 1]) {
			++m;
		} else {
			m = 0;
		}
		if (qAln[aln_to - 1] != '-') ++tqcnt;
		if (tAln[aln_to - 1] != '-') ++ttcnt;
		--aln_to;
		if (m == kMatchSize) break;
	}
	oc_assert(m == kMatchSize);
	aln_to += kMatchSize;
	tqcnt -= kMatchSize;
	ttcnt -= kMatchSize;
	qBgn += pqcnt;
	qEnd -= tqcnt;
	tBgn += ptcnt;
	tEnd -= ttcnt;

	kstr_clear(align_data->query_align);
	kstr_clear(align_data->target_align);
	kputsn(qAln + aln_from, aln_to - aln_from, &align_data->query_align);
	kputsn(tAln + aln_from, aln_to - aln_from, &align_data->target_align);
	
	align_data->qoff = query_from + qBgn;
	align_data->qend = query_from + qEnd;
	align_data->toff = target_from + tBgn;
	align_data->tend = target_from + tEnd;
	
	align_data->ident_perc = calc_ident_perc(kstr_str(align_data->query_align),
											 kstr_str(align_data->target_align),
											 kstr_size(align_data->query_align),
											 &align_data->dist);
	
	edlibFreeAlignResult(align, align_data->km);
	
	validate_aligned_string(0,
							query,
							align_data->qoff,
							align_data->qend,
							kstr_str(align_data->query_align),
							0,
							target,
							align_data->toff,
							align_data->tend,
							kstr_str(align_data->target_align),
							kstr_size(align_data->query_align),
							TRUE);
	
	return 1;
}
