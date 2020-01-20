#include "cbcns.h"

#include "../common/oc_assert.h"

#include "../klib/ksort.h"

#define TrimReadsIntP_SoffLT(a, b) \
(((a).first < (b).first) || ((a).first == (b).first && (a).second > (b).second))

KSORT_INIT(TrimReadsIntP_SoffLT, IntPair, TrimReadsIntP_SoffLT)

CbCnsData*
new_CbCnsData()
{
	CbCnsData* data = (CbCnsData*)malloc(sizeof(CbCnsData));
	kv_init(data->tags);
	kv_init(data->backbone);
	kv_init(data->coverage);
	data->li_alloc = new_OcObjectAllocator(sizeof(LinkInfo));
	data->dci_alloc = new_OcObjectAllocator(sizeof(DeltaCovInfo));
	return data;
}

CbCnsData*
free_CbCnsData(CbCnsData* data)
{
	kv_destroy(data->tags);
	kv_destroy(data->backbone);
	kv_destroy(data->coverage);
	data->li_alloc = free_OcObjectAllocator(data->li_alloc);
	data->dci_alloc = free_OcObjectAllocator(data->dci_alloc);
	free(data);
	return 0;
}

void
clear_CbCnsData(CbCnsData* data)
{
	kv_clear(data->tags);
	kv_clear(data->backbone);
	kv_clear(data->coverage);
	clear_OcObjectAllocator(data->li_alloc);
	clear_OcObjectAllocator(data->dci_alloc);
}

void
add_one_align(CbCnsData* cns_data,
			  const char* qaln,
			  const char* taln,
			  const size_t aln_size,
			  kstring_t* target,
			  int toff,
			  int tend,
			  double weight)
{
	get_cns_tags(qaln,
				 taln,
				 aln_size,
				target,
				 toff,
				 tend,
				 weight,
				 &cns_data->tags);
}

static BOOL
calc_max_consensus_range(IntPair* ranges, const int nranges, const int min_length, int* _left, int* _right)
{
	if (nranges == 0) return FALSE;
	ks_introsort_TrimReadsIntP_SoffLT(nranges, ranges);
	int i = 0, j;
	int left = ranges[0].first, right;
	int max_left = 0, max_right = 0;

	while (i < nranges) {
		j = i + 1;
		while (j < nranges && ranges[j].second <= ranges[i].second) ++j;
		if (j == nranges) {
			right = ranges[i].second;
			if (right - left > max_right - max_left) {
				max_left = left;
				max_right = right;
			}
			break;
		}

		if (ranges[j].first > ranges[i].second || ranges[i].second - ranges[j].first < 100) {
			right = ranges[i].second;
			if (right - left > max_right - max_left) {
				max_left = left;
				max_right = right;
			}
			left = ranges[j].first;
		}

		i = j;
	}

	if (max_right - max_left >= min_length) {
		*_left = max_left;
		*_right = max_right;
		return TRUE;
	}
	return FALSE;
}

void
consensus_broken(CbCnsData* cns_data,
				  const int min_cov,
				  const int min_size,
				  const int template_id,
				  const int template_size,
				  vec_intpair* cns_intvs,
				  CnsSeq* cns_seq,
				  RecordWriter* out,
				  const int check_chimeric_read,
				  vec_intpair* cov_ranges)
{
	build_backbone(kv_data(cns_data->tags), 
				   kv_size(cns_data->tags), 
				   template_size, 
				   cns_data->dci_alloc, 
				   cns_data->li_alloc,
				   &cns_data->backbone,
				   &cns_data->coverage);
	
	int cns_from = 0, cns_to = template_size;
	
	int i = cns_from;
	int* coverage = kv_data(cns_data->coverage);
	while (i < cns_to) {
		while (i < cns_to && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < cns_to && coverage[j] >= min_cov) ++j;
		if (j - i >= min_size * 0.85) {
//printf("i = %d, j = %d, min_cov = %d\n", i, j, min_cov);
			clear_CnsSeq(*cns_seq);
			consensus_backbone_segment(kv_data(cns_data->backbone),
									   i,
									   j,
									   coverage,
									   &cns_seq->cns_seq,
									   NULL,
									   NULL);
			if (kstr_size(cns_seq->cns_seq) >= min_size) {
				cns_seq->id = template_id;
				cns_seq->left = i;
				cns_seq->right = j;
				cns_seq->org_seq_size = template_size;
				for (size_t k = 0; k != kstr_size(cns_seq->cns_seq); ++k) {
					int c = kstr_A(cns_seq->cns_seq, k);
					kstr_A(cns_seq->cns_seq, k) = DecodeDNA(c);
				}
				RW_DUMP_ONE_DATA(CnsSeq, DUMP_CNS_SEQ, FALSE, out, cns_seq);
				IntPair ip;
				ip.first = i;
				ip.second = j;
				kv_push(IntPair, *cns_intvs, ip);
			}
		}
		i = j;
	}
}

typedef struct {
	int raw_from, raw_to;
	int cns_from, cns_to;
} CnsIntv;
typedef kvec_t(CnsIntv) vec_cns_intv;

int
consensus_unbroken(CbCnsData* cns_data,
			   const int min_cov,
			   const int min_size,
			   const char* raw_read,
			   const int template_size,
			   kstring_t* cns_seq)
{
	build_backbone(kv_data(cns_data->tags), 
				   kv_size(cns_data->tags), 
				   template_size, 
				   cns_data->dci_alloc, 
				   cns_data->li_alloc,
				   &cns_data->backbone,
				   &cns_data->coverage);
	int i = 0;
	int raw_from, raw_to, cns_from, cns_to;
	new_kstring(cns_frag);
	new_kstring(cns_all);
	new_kvec(vec_cns_intv, cns_intvs);
	CnsIntv intv;
	int* coverage = kv_data(cns_data->coverage);
	kstr_clear(*cns_seq);
	int num_cns_intvs = 0;
	while (i < template_size) {
		while (i < template_size && coverage[i] < min_cov) ++i;
		int j = i + 1;
		while (j < template_size && coverage[j] >= min_cov) ++j;
		
		if (j - i >= min_size * 0.85) {
			consensus_backbone_segment(kv_data(cns_data->backbone),
									   i,
									   j,
									   coverage,
									   &cns_frag,
									   &raw_from,
									   &raw_to);
			if (kstr_size(cns_frag) >= min_size) {
				cns_from = kstr_size(cns_all);
				cns_to = cns_from + kstr_size(cns_frag);
				kputsn(kstr_str(cns_frag), kstr_size(cns_frag), &cns_all);
				intv.raw_from = raw_from;
				intv.raw_to = raw_to;
				intv.cns_from = cns_from;
				intv.cns_to = cns_to;
				kv_push(CnsIntv, cns_intvs, intv);
				++num_cns_intvs;
			}
		}
		i = j;
	}
	
	if (kv_size(cns_intvs) == 0) {
		free_kstring(cns_frag);
		free_kstring(cns_all);
		free_kvec(cns_intvs);
		return num_cns_intvs;
	}
	
	int last_raw_to = 0;
	int last_cns_to = 0;
	for (size_t k = 0; k != kv_size(cns_intvs); ++k) {
		raw_from = kv_A(cns_intvs, k).raw_from;
		raw_to = kv_A(cns_intvs, k).raw_to;
		cns_from = kv_A(cns_intvs, k).cns_from;
		cns_to = kv_A(cns_intvs, k).cns_to;
		oc_assert(cns_from == last_cns_to);
		if (raw_from > last_raw_to) {
			oc_assert(raw_from < template_size);
			kputsn(raw_read + last_raw_to, raw_from - last_raw_to, cns_seq);
		}
		kputsn(kstr_str(cns_all) + cns_from, cns_to - cns_from, cns_seq);
		last_raw_to = raw_to;
		last_cns_to = cns_to;
		oc_assert(last_raw_to <= template_size);
	}
	
	if (last_raw_to != template_size) {
		oc_assert(last_raw_to < template_size);
		kputsn(raw_read + last_raw_to, template_size - last_raw_to, cns_seq);
	}
	
	for (size_t k = 0; k != kstr_size(*cns_seq); ++k) {
		int c = kstr_A(*cns_seq, k);
		oc_assert(c >= 0 && c < 4);
		kstr_A(*cns_seq, k) = DecodeDNA(c);
	}

	free_kstring(cns_frag);
	free_kstring(cns_all);
	free_kvec(cns_intvs);
	
	return num_cns_intvs;
}
