#include "cns_ctg_subseq.h"

#include "chain_dp.h"
#include "mem_finder.h"
#include "fc_correct_one_read.h"
#include "../common/makedb_aux.h"
#include "../common/packed_db.h"
#include "../common/m4_record.h"
#include "../gapped_align/oc_daligner.h"
#include "../gapped_align/oc_aligner.h"
#include "../edlib/edlib_wrapper.h"
#include "../klib/ksort.h"

#define IS_LOWER(c) ((c) >= 'a' && (c) <= 'z')
#define IS_UPPER(c) ((c) >= 'A' && (c) <= 'Z')

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

static int
cov_stats_is_full(u8* cov_stats, int soff, int send)
{
	int n = 0;
	for (int i = soff; i < send; ++i) if (cov_stats[i] < MAX_CNS_COV) ++n;
    return n < 200;
}

static void
extend_ctg_m4s(PackedDB* reads,
    const u8* ctg_subseq,
    const idx ctg_subseq_offset,
    const int ctg_subseq_size,
    M4Record* m4_array,
    const int m4_count,
	CnsData* fc_cns_data)
{
    MaximalExactMatchWorkData* mem_data = MaximalExactMatchWorkDataNew(CNS_MEM_KMER_SIZE, CNS_MEM_WINDOW_SIZE, CNS_MEM_MEM_SIZE);
    MaximalExactMatchWorkData_Init(mem_data, ctg_subseq, ctg_subseq_size);
    ChainDpWorkData* chain_data = ChainDpWorkDataNew(1, 200);
    OcAlignData* align_data = new_OcAlignData(0.5);
    OcDalignData* dalign_data = new_OcDalignData(0.5);
    FullEdlibAlignData* falign_data = new_FullEdlibAlignData(0.5);
    new_kstring(read);
    int qb, qe, sb, se;
    double ident_perc;
    kstring_t* qaln = NULL;
    kstring_t* saln = NULL;
    u8* cov_stats = (u8*)calloc(ctg_subseq_size, sizeof(u8));
    int num_extended_m4 = 0;

    for (int i = 0; i < m4_count; ++i) {
        M4Record* m4 = m4_array + i;
        oc_assert(m4->sdir == FWD);
        sb = (m4->soff >= ctg_subseq_offset) ? (m4->soff - ctg_subseq_offset) : 0;
        oc_assert(m4->send > ctg_subseq_offset);
        se = m4->send - ctg_subseq_offset;
        se = OC_MIN(se, ctg_subseq_size);
        if (cov_stats_is_full(cov_stats, sb, se)) {
            //OC_LOG("%d --- %d is full", sb, se);
            continue;
        }
        pdb_extract_sequence(reads, m4->qid, m4->qdir, &read);
        const u8* rs = (u8*)(kstr_data(read));
        const int rl = kstr_size(read);
        oc_assert(rl == m4->qsize, "rl = %d", rl);
        if (!MaximalExactMatchWorkData_FindCandidates(mem_data, chain_data, rs, rl)) continue;
        oc_assert(kv_size(mem_data->can_list) > 0);
        GappedCandidate can = kv_A(mem_data->can_list, 0);
        can.qsize = rl;
        can.ssize = ctg_subseq_size;
        //DUMP_GAPPED_CANDIDATE(fprintf, stderr, can);
        int r = cns_extension(&can,
                    align_data,
                    dalign_data,
                    falign_data,
                    (const char*)rs,
                    (const char*)ctg_subseq,
                    1000,
                    1,
                    &qb,
                    &qe,
                    &sb,
                    &se,
                    &ident_perc,
                    &qaln,
                    &saln);
        ++num_extended_m4;
        if (!r) continue;
        r = (ident_perc >= 90.0) && (qe - qb >= rl * 0.6 || sb - se >= ctg_subseq_size * 0.6);
        if (!r) continue;
        //OC_LOG("[%d, %d, %d] x [%d, %d, %d], %g", qb, qe, rl, sb, se, ctg_subseq_size, ident_perc);
        for (int p = sb; p < se; ++p) ++cov_stats[p];
		add_align_tags(kstr_data(*qaln), kstr_data(*saln), kstr_size(*qaln),
			qb, qe, sb, se, 1.0, &fc_cns_data->tags);
    }

    MaximalExactMatchWorkDataFree(mem_data);
    ChainDpWorkDataFree(chain_data);
    free_OcAlignData(align_data);
    free_OcDalignData(dalign_data);
    free_FullEdlibAlignData(falign_data);
    free_kstring(read);
    free(cov_stats);

    OC_LOG("extended m4: %d", num_extended_m4);
}

void
cns_ctg_subseq(PackedDB* reads,
    u8* ctg_subseq,
    const idx ctg_subseq_offset,
    const int ctg_subseq_size,
    M4Record* m4_array,
    const int m4_count,
	kstring_t* cns_ctg)
{
	CnsData* fc_cns_data = cns_data_new();
    extend_ctg_m4s(reads, ctg_subseq, ctg_subseq_offset, ctg_subseq_size, m4_array, m4_count, fc_cns_data);
	get_cns_from_align_tags(fc_cns_data, ctg_subseq_size, MIN_CNS_COV);

	for (int i = 0; i < ctg_subseq_size; ++i) ctg_subseq[i] = DecodeDNA(ctg_subseq[i]);
    int last_raw_idx = 0;
	const char* cns_seq = kstr_data(fc_cns_data->cns_seq);
	const int n = kstr_size(fc_cns_data->cns_seq);
	const int* t_pos = kv_data(fc_cns_data->t_pos);
	int i = 0;
    while (i < n) {
        while (i < n && IS_LOWER(cns_seq[i])) ++i;
        if (i >= n) break;
        size_t j = i + 1;
        while (j < n && IS_UPPER(cns_seq[j])) ++j;
        if (j - i >= 1000) {
            int curr_raw_idx = t_pos[i];
            const char* tmp = (const char*)(ctg_subseq + last_raw_idx);
            if (curr_raw_idx > last_raw_idx) kputsn(tmp, curr_raw_idx - last_raw_idx, cns_ctg);
            kputsn(cns_seq + i, j - i, cns_ctg);
            last_raw_idx = t_pos[j - 1] + 1;
        }
        i = j;
    }
    const char* tmp = (const char*)(ctg_subseq + last_raw_idx);
    if (ctg_subseq_size > last_raw_idx) kputsn(tmp, ctg_subseq_size - last_raw_idx, cns_ctg);
	cns_data_free(fc_cns_data);
}
