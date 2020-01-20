#include "consensus_one_read.h"

#include "../common/oc_assert.h"
#include "../common/record_writer.h"
#include "../klib/kvec.h"
#include "error_estimate.h"
#include "overlaps_pool.h"

#include <assert.h>

#define calc_cns_weight(ident_perc, weight) \
        do { \
                double __error = (100.0 - (ident_perc)) / 100.0 / 2.0; \
                weight = (1.0 - __error) * (1.0 - __error) + __error * __error / 3.0; \
                if (100.0 - ident_perc <= 1.0e-6) weight = 1.0; \
        } while(0)


static void
get_raw_intvs(int read_size, vec_intpair* cns_intvs, vec_intpair* raw_intvs)
{
	BOOL first = TRUE;
	BOOL last_interval_open = TRUE;
	int start_offset, end_offset, filter_start, filter_end;
	int left = 0, right;
	IntPair intv;
	
	start_offset = 0;
	end_offset = read_size - 1;
	for (size_t i = 0; i != kv_size(*cns_intvs); ++i) {
		filter_start = start_offset + kv_A(*cns_intvs, i).first;
		filter_end = start_offset + kv_A(*cns_intvs, i).second;
		if (first) {
			last_interval_open = 0;
			first = 0;
			if (filter_start > start_offset) {
				left = start_offset;
			} else {
				left = filter_end + 1;
				continue;
			}
		}
		
		right = filter_start - 1;
		intv.first = left;
		intv.second = right;
		if (right - left + 1 >= 1000) kv_push(IntPair, *raw_intvs, intv);
		
		if (filter_end >= end_offset) {
			last_interval_open = 0;
			break;
		} else {
			left = filter_end + 1;
		}
	}
	
	if (last_interval_open) {
		right = end_offset;
		if (right - left + 1 >= 1000) {
			intv.first = left;
			intv.second = right;
			kv_push(IntPair, *raw_intvs, intv);
		}
	}
}

static void
consensus_one_read_broken(CbCnsData* cbcns_data,
						 kstring_t* target,
						 const int min_cov,
						 const int min_size,
						 const int template_id,
						 const int template_size,
						 CnsSeq* cns_seq,
						 RecordWriter* cns_out,
						 RecordWriter* raw_out,
						 const int check_chimeric_read,
						 vec_intpair* cov_ranges)
{
	new_kvec(vec_intpair, cns_intvs);
	new_kvec(vec_intpair, raw_intvs);
	consensus_broken(cbcns_data,
					  min_cov,
					  min_size,
					  template_id,
					  template_size,
					  &cns_intvs,
					  cns_seq,
					  cns_out,
					  check_chimeric_read,
					  cov_ranges);
	get_raw_intvs(template_size, &cns_intvs, &raw_intvs);
	for (size_t i = 0; i != kv_size(raw_intvs); ++i) {
		int from = kv_A(raw_intvs, i).first;
		int to = kv_A(raw_intvs, i).second + 1;
		cns_seq->id = template_id;
		cns_seq->left = from;
		cns_seq->right = to;
		cns_seq->org_seq_size = template_size;
		kstr_clear(cns_seq->cns_seq);
		for (int k = from; k < to; ++k) {
			char c = DecodeDNA( kstr_A(*target, k) );
			kputc(c, &cns_seq->cns_seq);
		}
		RW_DUMP_ONE_DATA(CnsSeq, DUMP_CNS_SEQ, FALSE, raw_out, cns_seq);
	}
	free_kvec(cns_intvs);
	free_kvec(raw_intvs);
}

static int 
consensus_one_read_unbroken(CbCnsData* cbcns_data,
							  kstring_t* target,
							  const int min_cov,
							  const int min_size,
							  const int template_id,
							  const int template_size,
							  CnsSeq* cns_seq,
							  RecordWriter* cns_out,
							  RecordWriter* raw_out)
{
	int num_cns_intvs = consensus_unbroken(cbcns_data, min_cov, min_size, kstr_str(*target), template_size, &cns_seq->cns_seq);
	if (num_cns_intvs) {
		cns_seq->id = template_id;
		cns_seq->left = 0;
		cns_seq->right = template_size;
		cns_seq->org_seq_size = template_size;
		RW_DUMP_ONE_DATA(CnsSeq, DUMP_CNS_SEQ, FALSE, cns_out, cns_seq);
	}
	return num_cns_intvs;
}

void
update_cov_stats(int* cov_stats, int toff, int tend)
{
	for (int i = toff; i < tend; ++i) ++cov_stats[i];
}

BOOL
region_coverage_is_full(const int* cov_stats, int toff, int tend, int max_cov)
{
	int n = 0;
	for (int i = toff; i < tend; ++i) if (cov_stats[i] >= max_cov) ++n;
	return n == tend - toff;
}

int
add_extended_overlaps(OverlapsPool* op,
					  vec_int* cov_stats,
					  ReadIdPool* extended_read_ids,
					  const double weight,
					  const double ident_cutoff,
					  kstring_t* qabuf,
					  kstring_t* tabuf,
					  CbCnsData* cbcns_data,
					  kstring_t* target,
					  const int target_size,
					  vec_intpair* cov_ranges,
					  const int min_align_size,
					  const double mapping_ratio)
{
	int num_ovlps = 0;
	for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
		OverlapIndex* oip = oc_op_oip(*op, i);
		if (oip->ident_perc < ident_cutoff) continue;
		if (!check_mapping_range(oip->qoff, oip->qend, oip->qsize, 
								 oip->toff, oip->tend, target_size,
								 min_align_size,
								 mapping_ratio))
			continue;
		++num_ovlps;
        double w;
        calc_cns_weight(oip->ident_perc, w);
		add_one_align(cbcns_data, oc_op_qstr(*op, i), oc_op_tstr(*op, i), oip->align_size, target, oip->toff, oip->tend, w);
		update_cov_stats(kv_data(*cov_stats), oip->toff, oip->tend);
		add_read_id(extended_read_ids, oip->qid);
		if (is_full_cov_ovlp(oip->qoff, oip->qend, oip->qsize, oip->toff, oip->tend, target_size, 1000, 200)) {
			IntPair ip; ip.first = oip->toff; ip.second = oip->tend;
			kv_push(IntPair, *cov_ranges, ip);
		}
	}
	return num_ovlps;
}

static int 
PackedGappedCandidate_CnsQidLT(PackedGappedCandidate a, PackedGappedCandidate b)
{
	int a_item, b_item;
	
	a_item = pcan_qid(a);
	b_item = pcan_qid(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_score(a);
	b_item = pcan_score(b);
	if (a_item != b_item) return a_item > b_item;
	
	a_item = pcan_qdir(a);
	b_item = pcan_qdir(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_qbeg(a);
	b_item = pcan_qbeg(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_sbeg(a);
	b_item = pcan_sbeg(b);
	if (a_item != b_item) return a_item < b_item;
	
	return 0;
}

KSORT_INIT(PackedGappedCandidate_CnsQidLT, PackedGappedCandidate, PackedGappedCandidate_CnsQidLT)

void
consensus_one_read(CnsData* cns_data, size_t can_sid, size_t can_eid)
{
	if (cns_data->options->min_cov > can_eid - can_sid) return;
	
	CnsReads* cns_reads = cns_data->cns_reads;
	PackedDB* reads = cns_data->reads;
	PackedGappedCandidate* candidates = kv_data(*cns_data->candidates);
	CnsOptions* options = cns_data->options;
	RecordWriter* cns_out = cns_data->cns_out;
	RecordWriter* raw_out = cns_data->raw_out;
	ReadIdPool* corrected_read_ids = cns_data->corrected_read_ids;
	ReadIdPool* extended_read_ids = cns_data->extended_read_ids;
	CbCnsData* cbcns_data = cns_data->cns_data;
	vec_int* cov_stats = &cns_data->cov_stats;
	OverlapsPool* op = cns_data->op;
	const int min_align_size = options->min_align_size;
	const int min_size = options->min_size;
	const int min_cov = options->min_cov;
	const int max_cov = options->max_cov;
	const int template_id = pcan_sid(candidates[can_sid]);
	const int template_size = reads ? PDB_SEQ_SIZE(reads, template_id) : cns_reads_seq_size(cns_reads, template_id);
	kstring_t* query = &cns_data->query;
	kstring_t* target = &cns_data->target;
	kstring_t* qabuf = &cns_data->qabuf;
	kstring_t* tabuf = &cns_data->tabuf;
	CnsSeq* cns_seq = &cns_data->cns_seq;
	cns_seq->hdr = reads ? PDB_SEQ_NAME(reads, template_id) : cns_reads_seq_hdr(cns_reads, template_id);
	if (reads) pdb_extract_sequence(reads, template_id, FWD, target);
	else cns_reads_extract_sequence(cns_reads, template_id, FWD, target);
	oc_assert(template_size == kstr_size(*target));
	clear_ReadIdPool(extended_read_ids);
	const size_t kMaxExaminedCan = MAX_EXAMINED_CAN;

	if (cns_reads) {
		if (can_eid - can_sid > kMaxExaminedCan) {
			can_eid = can_sid + kMaxExaminedCan;
		} else {
			ks_introsort_PackedGappedCandidate_CnsScoreGT(can_eid - can_sid, candidates + can_sid);
		}
	} else {
		ks_introsort_PackedGappedCandidate_CnsScoreGT(can_eid - can_sid, candidates + can_sid);
		if (can_eid - can_sid > kMaxExaminedCan) {
			can_eid = can_sid + kMaxExaminedCan;
		}
	}

	clear_CbCnsData(cbcns_data);
	cbcns_data->template_size = template_size;
	kv_resize(int, *cov_stats, template_size);
	kv_zero(int, *cov_stats);
	kv_clear(cns_data->cov_ranges);
	
	size_t last_extended_can_id;
	double ident_cutoff;
	const double cns_weight = 1.0;
	size_t next_can_id = can_sid;
	
	if (options->use_fixed_ident_cutoff) {
		next_can_id = can_sid;
		ident_cutoff = 100.0 * (1.0 - options->error);
		cns_seq->num_ovlps = 0;
		cns_seq->num_can = 0;
		cns_seq->ident_cutoff = ident_cutoff;
	} else {
		if (!get_good_overlaps(candidates,
						   can_sid,
						   can_eid,
						   &last_extended_can_id,
						   cns_data->align_data,
							cns_data->dalign_data,
							cns_data->falign_data,
							options->rescue_long_indels,
						   op,
						   query,
						   target,
						   reads,
						   cns_reads,
						   extended_read_ids,
						   min_align_size,
						   &ident_cutoff)) {
			return;
		}
		cns_seq->ident_cutoff = ident_cutoff;
		cns_seq->num_ovlps = add_extended_overlaps(op, 
							 cov_stats, 
							 extended_read_ids, 
							 cns_weight, 
							 ident_cutoff, 
							 qabuf, 
							 tabuf, 
							 cbcns_data, 
							 target, 
							 template_size, 
							 &cns_data->cov_ranges,
							 min_align_size,
							 options->mapping_ratio);
		next_can_id = last_extended_can_id + 1;
		cns_seq->num_can = next_can_id - can_sid;
	}
	
	GappedCandidate _can, *can = &_can;
	while (1) {
		if (next_can_id >= can_eid) { break; }
		if (region_coverage_is_full(kv_data(*cov_stats), 0, template_size, max_cov)) { break; }
		size_t from = next_can_id;
		size_t to = from + 50; to = OC_MIN(to, can_eid);
		next_can_id = to;
		for (size_t i = from; i < to; ++i) {
			unpack_candidate(can, candidates + i);
			can->qsize = reads ? PDB_SEQ_SIZE(reads, can->qid) : cns_reads_seq_size(cns_reads, can->qid);
			can->ssize = template_size;
			oc_assert(can->sdir == FWD);
			oc_assert(can->sid == template_id);
			if (read_id_exists(extended_read_ids, can->qid)) { continue; }
			if (region_coverage_is_full(kv_data(*cov_stats), can->sbeg, can->send, max_cov)) { continue; }
			if (reads) pdb_extract_sequence(reads, can->qid, can->qdir, query);
			else cns_reads_extract_sequence(cns_reads, can->qid, can->qdir, query);
			oc_assert(can->qsize == kstr_size(*query));
			
			int qoff, qend, toff, tend;
			kstring_t* qalign;
			kstring_t* talign;
			double ident_perc;
			BOOL r = cns_extension(can, 
								   cns_data->align_data,
								   cns_data->dalign_data,
								   cns_data->falign_data,
								   kstr_str(*query),
								   kstr_str(*target),
								   min_align_size,
								   options->rescue_long_indels,
								   &qoff,
								   &qend,
								   &toff,
								   &tend,
								   &ident_perc,
								   &qalign,
								   &talign);
			++cns_seq->num_can;
			if (!r) { continue; }
			if (ident_perc < ident_cutoff && (!is_full_cov_ovlp(qoff, qend, can->qsize, toff, tend, can->ssize, 5000, 100))) continue;
			if (!check_mapping_range(qoff, qend, can->qsize, toff, tend, can->ssize, options->min_align_size, options->mapping_ratio)) continue;
            double w;
            calc_cns_weight(ident_perc, w);
			++cns_seq->num_ovlps;
			add_one_align(cbcns_data, 
						  kstr_data(*qalign), 
						  kstr_data(*talign), 
						  kstr_size(*qalign), 
						  target, 
						  toff, 
						  tend, 
						  w);
			update_cov_stats(kv_data(*cov_stats), toff, tend);
			add_read_id(extended_read_ids, can->qid);
		}
	}
	
	int r = 1;
	if (options->full_consensus) {
		r = consensus_one_read_unbroken(cbcns_data, target, min_cov, min_size, template_id, template_size, cns_seq, cns_out, raw_out);
	} else {
		consensus_one_read_broken(cbcns_data, 
								  target, 
								  min_cov, 
								  min_size, 
								  template_id, 
								  template_size, 
								  cns_seq, 
								  cns_out, 
								  raw_out,
								  FALSE,
								  &cns_data->cov_ranges);
	}
	
	if (r) add_read_id(corrected_read_ids, template_id);
}

BOOL
region_coverage_is_full_m4(const int* cov_stats, 
						   int trg_from,
						   int trg_to,
						   int toff, 
						   int tend, 
						   int max_cov)
{
	int n = 0;
	int from = OC_MAX(trg_from, toff);
	int to = OC_MIN(trg_to, tend);
	for (int i = from; i < to; ++i) if (cov_stats[i] >= max_cov) ++n;
	return n == to - from;
}

void 
consensus_one_read_m4(PackedDB* reads, 
					  CnsReads* cns_reads,
                      kstring_t* target,
                      kstring_t* read,
                      M4Record* m4v,
                      int nm4,
					  CbCnsData* cbcns_data,
					  vec_int* cov_stats,
					  FullEdlibAlignData* align_data,
					  int trg_from,
					  int trg_to,
					  int* read_id,
					  FILE* out,
					  pthread_mutex_t* out_lock)
{
	for (int i = 0; i < nm4; ++i) {
		assert(m4v[i].qid > 0);
		assert(m4v[i].sid > 0);
		--m4v[i].qid;
		--m4v[i].sid;
	}
	const int target_id = m4v[0].sid;
	const int target_size = m4v[0].ssize;
	kv_resize(int, *cov_stats, target_size);
	kv_fill(*cov_stats, 0);
	clear_CbCnsData(cbcns_data);
	if (reads) pdb_extract_sequence(reads, target_id, FWD, target);
	else cns_reads_extract_sequence(cns_reads, target_id, FWD, target);
	assert(kstr_size(*target) == target_size);
	const int min_cov = 1;
	const int max_cov = 12;
	const int min_size = 500;
	int next_m4_id = 0;
	int qoff, qend, toff, tend;
	double ident_perc;
	kstring_t* qalign = NULL;
	kstring_t* talign = NULL;
	int num_extended_m4 = 0;

	while (1) {
		if (next_m4_id >= nm4) break;
		if (region_coverage_is_full_m4(kv_data(*cov_stats), trg_from, trg_to, trg_from, trg_to, max_cov)) { break; }
		int from = next_m4_id;
		int to = OC_MIN(from + 50, nm4);
		next_m4_id = to;
		for (int i = from; i < to; ++i) {
			if (m4v[i].send <= trg_from || m4v[i].soff >= trg_to) continue;
			if (region_coverage_is_full_m4(kv_data(*cov_stats), trg_from, trg_to, m4v[i].soff, m4v[i].send, max_cov)) continue;
			if (reads) pdb_extract_sequence(reads, m4v[i].qid, m4v[i].qdir, read);
			else (cns_reads_extract_sequence(cns_reads, m4v[i].qid, m4v[i].qdir, read));
			assert(kstr_size(*read) == m4v[i].qsize);
			int tolerence = (100.0 - m4v[i].ident_perc) / 100.0 * (m4v[i].send - m4v[i].soff) * 1.5;
			if (m4v[i].qdir == FWD) {
				qoff = m4v[i].qoff;
				qend = m4v[i].qend;
			} else {
				assert(m4v[i].qdir == REV);
				qoff = m4v[i].qsize - m4v[i].qend;
				qend = m4v[i].qsize - m4v[i].qoff;
			}
			int r = edlib_go(kstr_str(*read),
							qoff,
							qend,
							kstr_str(*target),
							m4v[i].soff,
							m4v[i].send,
							align_data,
							tolerence,
							400,
							TRUE,
							4);
			++num_extended_m4;
			if (!r) continue;
			qoff = align_data->qoff;
			qend = align_data->qend;
			toff = align_data->toff;
			tend = align_data->tend;
			qalign = &align_data->query_align;
			talign = &align_data->target_align;
			ident_perc = align_data->ident_perc;
			if (ident_perc < 90.0) continue;
			//r = ((qend - qoff) >= m4v[i].qsize*0.6) || ((tend - toff) >= m4v[i].ssize * 0.6);
			//if (!r) continue;
			add_one_align(cbcns_data, 
						  kstr_data(*qalign), 
						  kstr_data(*talign), 
						  kstr_size(*qalign), 
						  target, 
						  toff, 
						  tend, 
						  1.0);
			update_cov_stats(kv_data(*cov_stats), toff, tend);
		}
	}

	kstring_t cns_seq;
	kstr_init(cns_seq);
	int max_from = 0, max_to = 0, max_size = 0;
	int i = trg_from;
	while (i < trg_to) {
		while (i < trg_to && kv_A(*cov_stats, i) < min_cov) ++i;
		int j = i + 1;
		while (j < trg_to && kv_A(*cov_stats, j) >= min_cov) ++j;
		if (j - i > max_size) {
			max_from = i;
			max_to = j;
			max_size = j - i;
		}
		i = j;
	}
	int cns_from, cns_to;
	if (max_size >= min_size) {
		build_backbone(kv_data(cbcns_data->tags), 
				   	   kv_size(cbcns_data->tags), 
				       target_size, 
				   	   cbcns_data->dci_alloc, 
				   	   cbcns_data->li_alloc,
				       &cbcns_data->backbone,
				       &cbcns_data->coverage);
		consensus_backbone_segment(kv_data(cbcns_data->backbone), max_from, max_to, kv_data(*cov_stats), &cns_seq, &cns_from, &cns_to);
		if (kstr_size(cns_seq) >= min_size) {
			pthread_mutex_lock(out_lock);
			fprintf(out, ">%d\n", *read_id);
			++(*read_id);
			for (size_t i = 0; i < kstr_size(cns_seq); ++i) {
				int c = kstr_A(cns_seq, i);
				c = "ACGT"[c];
				fprintf(out, "%c", c);
			}
			fprintf(out, "\n");
			pthread_mutex_unlock(out_lock);
		}
	}
	free_kstring(cns_seq);
}
