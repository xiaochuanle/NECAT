#include "ksw2_aligner.h"

void
init_Ksw2Params(Ksw2Params* p)
{
	p->a = 2;
	p->b = 4;
	p->q = 4;
	p->e = 2;
	p->q2 = 24;
	p->e2 = 1;
	p->bw = 800;
	p->zdrop = 400;
}

Ksw2AlignData*
new_Ksw2AlignData()
{
	Ksw2AlignData* data = (Ksw2AlignData*)malloc(sizeof(Ksw2AlignData));
	init_Ksw2Params(&data->align_params);
	data->km = km_init();
	kstr_init(data->query_align);
	kstr_init(data->target_align);
	kstr_init(data->fqaln);
	kstr_init(data->rqaln);
	kstr_init(data->ftaln);
	kstr_init(data->rtaln);
	kv_init(data->qfrag);
	kv_init(data->tfrag);
	
	return data;
}

void
free_Ksw2AlignData(Ksw2AlignData* data)
{
	km_destroy(data->km);
	free_kstring(data->query_align);
	free_kstring(data->target_align);
	free_kstring(data->fqaln);
	free_kstring(data->rqaln);
	free_kstring(data->ftaln);
	free_kstring(data->rtaln);
	kv_destroy(data->qfrag);
	kv_destroy(data->tfrag);
	free(data);
}

static void
gen_simple_mat(int m, int8_t* mat, int8_t a, int8_t b)
{
	int i ,j;
	a = a < 0 ? -a : a;
	b = b > 0 ? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j ? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
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

void
gen_aligned_string(vec_u8* Q, vec_u8* T, ksw_extz_t* ez, kstring_t* qaln, kstring_t* taln)
{
	kstr_clear(*qaln);
	kstr_clear(*taln);
	const char gap = '-';
printf("n_cigar = %d\n", ez->n_cigar);
	
	int qidx = 0, tidx = 0;
	for (int i = 0; i < ez->n_cigar; ++i) {
		int num = ez->cigar[i] >> 4;
		char type = "MID"[ez->cigar[i] & 0xf];
		if (type == 'M') {
			for (int k = 0; k < num; ++k) {
				int c = kv_A(*Q, qidx);
				c = "ACGT"[c];
				kputc(c, qaln);
				++qidx;
				
				c = kv_A(*T, tidx);
				c = "ACGT"[c];
				kputc(c, taln);
				++tidx;
			}
		} else if (type == 'D') {
			for (int k = 0; k < num; ++k) {
				int c = kv_A(*Q, qidx);
				c = "ACGT"[c];
				kputc(c, qaln);
				++qidx;
				
				kputc(gap, taln);
			}
		} else {
			for (int k = 0; k < num; ++k) {
				kputc(gap, qaln);
				
				int c = kv_A(*T, tidx);
				c = "ACGT"[c];
				kputc(c, taln);
				++tidx;
			}
		}
	}
}

BOOL
ksw2_go(const char* query,
		  const int query_start,
		  const int query_size,
		  const char* target,
		  const int target_start,
		  const int target_size,
		  Ksw2AlignData* ksw_data,
		  const int min_align_size)
{
	kstr_clear(ksw_data->query_align);
	kstr_clear(ksw_data->target_align);
	int8_t mat[25];
	gen_simple_mat(5, mat, ksw_data->align_params.a, ksw_data->align_params.b);
	const int align_flag = KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR;
	ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
	int QS = query_start;
	int TS = target_start;
	
	// extend left 
	{
		kv_clear(ksw_data->qfrag);
		for (int i = query_start; i; --i) {
			uint8_t c = query[i - 1];
			kv_push(uint8_t, ksw_data->qfrag, c);
		}
		kv_clear(ksw_data->tfrag);
		for (int i = target_start; i; --i) {
			uint8_t c = target[i - 1];
			kv_push(uint8_t, ksw_data->tfrag, c);
		}
		ksw_extd2_sse(NULL, //ksw_data->km, 
					  QS, 
					  kv_data(ksw_data->qfrag), 
					  TS, 
					  kv_data(ksw_data->tfrag), 
					  5, 
					  mat, 
					  ksw_data->align_params.q, 
					  ksw_data->align_params.e, 
					  ksw_data->align_params.q2, 
					  ksw_data->align_params.e2, 
					  ksw_data->align_params.bw, 
					  ksw_data->align_params.zdrop, 
					  align_flag, 
					  &ez);
		gen_aligned_string(&ksw_data->qfrag, 
						   &ksw_data->tfrag, 
						   &ez, 
						   &ksw_data->rqaln,
						   &ksw_data->rtaln);
		kfree(NULL, ez.cigar);
		ez.cigar = 0;
	}
	
	// extend right
	{
		kv_clear(ksw_data->qfrag);
		for (int i = query_start; i < query_size; ++i) {
			uint8_t c = query[i];
			kv_push(uint8_t, ksw_data->qfrag, c);
		}
		kv_clear(ksw_data->tfrag);
		for (int i = target_start; i < target_size; ++i) {
			uint8_t c = target[i - 1];
			kv_push(uint8_t, ksw_data->tfrag, c);
		}
		ksw_extd2_sse(NULL, //ksw_data->km, 
					  QS, 
					  kv_data(ksw_data->qfrag), 
					  TS, 
					  kv_data(ksw_data->tfrag), 
					  5, 
					  mat, 
					  ksw_data->align_params.q, 
					  ksw_data->align_params.e, 
					  ksw_data->align_params.q2, 
					  ksw_data->align_params.e2, 
					  ksw_data->align_params.bw, 
					  ksw_data->align_params.zdrop, 
					  align_flag, 
					  &ez);
		gen_aligned_string(&ksw_data->qfrag, 
						   &ksw_data->tfrag, 
						   &ez, 
						   &ksw_data->fqaln,
						   &ksw_data->ftaln);
//		kfree(ksw_data->km, ez.cigar);
		kfree(NULL, ez.cigar);
	}
	
	int rqcnt = 0, rtcnt = 0;
	for (size_t i = kstr_size(ksw_data->rqaln); i; --i) {
		char c = kstr_A(ksw_data->rqaln, i - 1);
		if (c != GAP_CHAR) ++rqcnt;
		kputc(c, &ksw_data->query_align);
		
		c = kstr_A(ksw_data->rtaln, i - 1);
		if (c != GAP_CHAR) ++rtcnt;
		kputc(c, &ksw_data->target_align);
	}
	
	int fqcnt = 0, ftcnt = 0;
	for (size_t i = 0; i < kstr_size(ksw_data->fqaln); ++i) {
		char c = kstr_A(ksw_data->fqaln, i);
		if (c != GAP_CHAR) ++fqcnt;
		kputc(c, &ksw_data->query_align);
		
		c = kstr_A(ksw_data->ftaln, i);
		if (c != GAP_CHAR) ++ftcnt;
		kputc(c, &ksw_data->target_align);
	}
	
	ksw_data->qoff = QS - rqcnt;
	ksw_data->qend = QS + fqcnt;
	ksw_data->toff = TS - rtcnt;
	ksw_data->tend = TS + ftcnt;
	
	if (fqcnt + rqcnt < min_align_size || ftcnt + rtcnt < min_align_size) return FALSE;
	
	ksw_data->ident_perc = calc_ident_perc(kstr_str(ksw_data->query_align), 
										   kstr_str(ksw_data->target_align),
										   kstr_size(ksw_data->query_align));
	return TRUE;
}
