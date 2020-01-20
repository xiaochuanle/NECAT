#include "largest_cover_range.h"
#include "../common/m4_record.h"
#include "range_list.h"
#include "pm4_aux.h"
#include "../tasc/cbcns.h"
#include "../edlib/edlib_wrapper.h"

#include <assert.h>

static BOOL
is_qualified_m4(M4Record* m4)
{
	const int L = 2000;
	const int M = 20;
	idx qoff = m4->qoff, qend = m4->qend, qsize = m4->qsize;
	if (m4->qdir == REV) { qoff = qsize - m4->qend; qend = qsize - m4->qoff; }
	idx soff = m4->soff, send = m4->send, ssize = m4->ssize;
	if (m4->sdir == REV) { soff = ssize - m4->send; send = ssize - m4->soff; }

	if (qoff <= M && qsize - qend <= M) return TRUE;
	if (soff <= M && ssize - send <= M) return TRUE;

	if (qsize - qend <= M) {
		if (soff > M) return FALSE;
		if (qend - qoff >= L) return TRUE;
	}

	if (ssize - send <= M) {
		if (qoff > M) return FALSE;
		if (qend - qoff >= L) return TRUE;
	}

	return FALSE;
}

static void 
truncate_m4_list(M4Record* m4v, int* nm4, const double ident_perc)
{
	int n = *nm4;
	const int MaxNm4 = 300;
	if (n > MaxNm4) {
		int j = 0;
		for (int i = 0; i < n; ++i) {
			if (m4v[i].ident_perc >= ident_perc) m4v[j++] = m4v[i];
		}
		n = j;
	}
	if (n > MaxNm4) {
		ks_introsort_M4Record_IdentGT(n, m4v);
		n = MaxNm4;
	}
	*nm4 = n;
}

BOOL
largest_cover_range(M4Record* m4v,
                    int nm4,
                    int* fbgn,
                    int* fend,
                    double min_ident_perc,
                    int min_ovlp_size,
                    int min_cov)
{
  CovRangeList IL, ID;
  init_CovRangeList(&IL);
  init_CovRangeList(&ID);
  int nskip = 0, nused = 0;
  
/*
  const int MaxNm4 = 300;
  if (nm4 > MaxNm4)
  {
	  int j = 0;
	  for (int i = 0; i < nm4; ++i) 
		  if (m4v[i].ident_perc >= min_ident_perc) m4v[j++] = m4v[i];
		  else ++nskip;
	  nm4 = j;
  }
  if (nm4 > MaxNm4) {
	  ks_introsort_M4Record_IdentGT(nm4, m4v);
	  nm4 = MaxNm4;
  }
*/

  for (int i = 0; i < nm4; ++i) {
    int tbgn = m4v[i].soff;
    int tend = m4v[i].send;
    add_CovRangeList(&IL, tbgn, tend - tbgn, 0);
    ++nused;
  }

  if (min_cov > 0) {
    CovRangeList DE; init_CovRangeList(&DE);
    depth_from_CovRangeList(&DE, &IL);
    size_t it = 0;
    int ib = 0, ie = 0;

    while (it < kv_size(DE.list)) {
      if (CovRangeListDepth(DE, it) < min_cov) {
        if (ie > ib) add_CovRangeList(&ID, ib, ie - ib, 0);
        ib = 0;
        ie = 0;
      } else if (ib == 0 && ie == 0) {
        ib = CovRangeListLo(DE, it);
        ie = CovRangeListHi(DE, it);
      } else if (ie == CovRangeListLo(DE, it)) {
        ie = CovRangeListHi(DE, it);
      } else {
        if (ie > ib) add_CovRangeList(&ID, ib, ie - ib, 0);
        ib = CovRangeListLo(DE, it);
        ie = CovRangeListHi(DE, it);
      }
      ++it;
    }

    if (ie > ib) add_CovRangeList(&ID, ib, ie - ib, 0);

    destroy_CovRangeList(&DE);
  }

  merge_CovRangeList(&IL, min_ovlp_size);

  if (min_cov > 0) {
    CovRangeList FI;
    init_CovRangeList(&FI);
    size_t li = 0, di = 0;
    while (li < kv_size(IL.list) && di < kv_size(ID.list)) {
      int ll = CovRangeListLo(IL, li);
      int lh = CovRangeListHi(IL, li);
      int dl = CovRangeListLo(ID, di);
      int dh = CovRangeListHi(ID, di);
      int nl = 0;
      int nh = 0;

      if (ll <= dl && dl < lh) {
        nl = dl;
        nh = OC_MIN(lh, dh);
      }

      if (dl <= ll && ll < dh) {
        nl = ll;
        nh = OC_MIN(lh, dh);
      }

      if (nl < nh) add_CovRangeList(&FI, nl, nh - nl, 0);

      if (lh <= dh) ++li;

      if (dh <= lh) ++di;
    }

    copy_CovRangeList(&IL, &FI);
    destroy_CovRangeList(&FI);
  }

  if(kv_size(IL.list) == 0) return FALSE;

  int max_l = 0, max_r = 0;
  for (size_t i = 0; i < kv_size(IL.list); ++i) {
    int l = kv_A(IL.list, i).lo;
    int h = kv_A(IL.list, i).hi;
    if (h - l > max_r - max_l) {
      max_l = l;
      max_r = h;
    }
  }

  *fbgn = max_l;
  *fend = max_r;
  
  destroy_CovRangeList(&IL);
  destroy_CovRangeList(&ID);
  return TRUE;
}

static void
load_partition_m4(const char* m4_path, const int pid, vec_m4* m4v, vec_int* idx_range)
{
	new_kstring(pname);
	make_partition_name(m4_path, pid, &pname);
	size_t file_bytes = FILE_SIZE(kstr_data(pname));
	size_t num_m4 = file_bytes / sizeof(M4Record);
	if (num_m4 == 0) return;
	kv_resize(M4Record, *m4v, num_m4);
	DFOPEN(m4_in, kstr_str(pname), "rb");
	FREAD(kv_data(*m4v), sizeof(M4Record), num_m4, m4_in);
	FCLOSE(m4_in);
	free_kstring(pname);
	ks_introsort_M4Record_SidLT(kv_size(*m4v), kv_data(*m4v));
	
	M4Record* m4s = kv_data(*m4v);
	int nm4 = kv_size(*m4v);
	int i = 0;
	kv_clear(*idx_range);
	while (i < nm4) {
		int j = i + 1;
		while (j < nm4 && m4s[j].sid == m4s[i].sid) ++j;
		kv_push(int, *idx_range, i);
		i = j;
	}
	kv_push(int, *idx_range, nm4);
}

typedef struct {
	PackedDB* reads;
	M4Record* m4v;
	int nm4;
	int* idx_range;
	int nrange;
	int next_range_id;
	pthread_mutex_t range_get_lock;
	int* read_id;
	FILE* out;
	pthread_mutex_t* out_lock;
	double min_ident_perc;
	int min_ovlp_size;
	int min_cov;
} LcrData;

LcrData*
new_LcrData(M4Record* m4v,
			int nm4,
			PackedDB* reads,
			int* idx_range,
			int nrange,
			int* read_id,
			FILE* out,
			pthread_mutex_t* out_lock,
		   	double min_ident_perc,
			int min_ovlp_size,
			int min_cov)
{
	LcrData* data = (LcrData*)malloc(sizeof(LcrData));
	data->m4v = m4v;
	data->reads = reads;
	data->nm4 = nm4;
	data->idx_range = idx_range;
	data->nrange = nrange;
	data->next_range_id = 0;
	pthread_mutex_init(&data->range_get_lock, NULL);
	data->read_id = read_id;
	data->out = out;
	data->out_lock = out_lock;
	data->min_ident_perc = min_ident_perc;
	data->min_ovlp_size = min_ovlp_size;
	data->min_cov = min_cov;
	return data;
}

LcrData*
free_LcrData(LcrData* data)
{
	free(data);
	return NULL;
}

static M4Record*
get_next_range(LcrData* data, int* nm4)
{
	M4Record* m4v = NULL;
	int i = 0;
	pthread_mutex_lock(&data->range_get_lock);
	i = data->next_range_id;
	++data->next_range_id;
	pthread_mutex_unlock(&data->range_get_lock);
	
	if (i < data->nrange) {
		m4v = data->m4v + data->idx_range[i];
		*nm4 = data->idx_range[i+1] - data->idx_range[i];
	}
	return m4v;
}

void*
lcr_worker(void* arg)
{
	LcrData* data = (LcrData*)(arg);
	M4Record* m4v = NULL;
	int nm4;
	int read_id, left, right, size;
	new_kstring(target);
	new_kstring(read);
	new_kvec(vec_int, cov_stats);
	CbCnsData* cbcns_data = new_CbCnsData();
	FullEdlibAlignData* align_data = new_FullEdlibAlignData(0.5);
	while ((m4v = get_next_range(data, &nm4))) {
		truncate_m4_list(m4v, &nm4, data->min_ident_perc);
		BOOL r = largest_cover_range(m4v,
									 nm4,
									 &left,
									 &right,
									 data->min_ident_perc,
									 data-> min_ovlp_size,
									 data->min_cov);
									 
		if (r) {
			consensus_one_read_m4(data->reads, 
					  NULL,
                      &target,
                      &read,
                      m4v,
                      nm4,
					  cbcns_data,
					  &cov_stats,
					  align_data,
					  left,
					  right,
					  data->read_id,
					  data->out,
					  data->out_lock);
		}
	}

	free_FullEdlibAlignData(align_data);
	free_CbCnsData(cbcns_data);
	free_kvec(cov_stats);
	free_kstring(read);
	free_kstring(target);
	return NULL;
}

void
get_largest_cover_range_for_one_partition(const char* m4_path, 
		const int pid, 
		PackedDB* reads,
		int* read_id,
		FILE* out, 
		pthread_mutex_t* out_lock,
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		const int num_threads)
{
	new_kvec(vec_m4, m4v);
	new_kvec(vec_int, idx_range);
	load_partition_m4(m4_path, pid, &m4v, &idx_range);
	if (kv_size(m4v) == 0) {
		free_kvec(m4v);
		free_kvec(idx_range);
		return;
	}

	LcrData* lcr_data = new_LcrData(kv_data(m4v),
									kv_size(m4v),
									reads,
									kv_data(idx_range),
									kv_size(idx_range) - 1,
									read_id,
									out,
									out_lock,
									min_ident_perc,
									min_ovlp_size,
									min_cov);
	pthread_t jobs[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		pthread_create(jobs + i, NULL, lcr_worker, (void*)(lcr_data));
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(jobs[i], NULL);
	}
	free_LcrData(lcr_data);
	free_kvec(m4v);
	free_kvec(idx_range);
}
