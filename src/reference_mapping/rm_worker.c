#include "rm_worker.h"

#include "../common/m4_record.h"
#include "../common/map_aux.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../lookup_table/lookup_table.h"
#include "../word_finder/word_finder.h"
#include "../gapped_align/oc_aligner.h"
#include "../gapped_align/oc_daligner.h"
#include "../edlib/edlib_wrapper.h"

#include <assert.h>

static int
GappedCandidate_RmScoreGT(GappedCandidate a, GappedCandidate b)
{
	if (a.score != b.score) return a.score > b.score;
	if (a.qdir != b.qdir) return a.qdir < b.qdir;
	if (a.sid != b.sid) return a.sid < b.sid;
	if (a.qoff != b.qoff) return a.qoff < b.qoff;
	if (a.soff != b.soff) return a.soff < b.soff;
	return 0;
}

KSORT_INIT(GappedCandidate_RmScoreGT, GappedCandidate, GappedCandidate_RmScoreGT)

PackedDB*
load_mapping_reads(kseq_t* read)
{
	if (kseq_read(read) < 0) return 0;
	const idx VolumeSize = 2000000000;
	
	PackedDB* reads = new_PackedDB();
	pdb_enlarge_size(reads, VolumeSize + 2000000);
	while (1) {
		pdb_add_one_seq(reads, read, TECH_NANOPORE);
		if (PDB_SIZE(reads) >= VolumeSize) break;
		if (kseq_read(read) < 0) break;
	}
	return reads;
}

static void
calc_reference_range(GappedCandidate* can, idx* from, idx* to, int* soff)
{
	idx n = (can->qoff < can->soff) ? (can->qoff * 1.3) : (can->soff);
	n = OC_MIN(n, can->soff);
	*from = can->soff - n;
	*soff = n;
	
	idx sr = can->ssize - can->soff;
	idx qr = can->qsize - can->qoff;
	n = (qr < sr) ? (qr * 1.3) : (sr);
	n = OC_MIN(n, sr);
	*to = can->soff + n;
}

static BOOL
rm_extend_candidate(GappedCandidate* can,
					OcAlignData* align_data,
					OcDalignData* dalign_data,
					FullEdlibAlignData* falign_data,
					const char* read,
					PackedDB* reference,
					kstring_t* subject,
					const int min_align_size,
					M4Record* m4)
{
	idx sfrom, sto;
	int ssize, soff;
	BOOL check_mapped_range = (can->qoff >= can->qbeg) && (can->qoff < can->qend) && (can->soff >= can->sbeg) && (can->soff < can->send);
	calc_reference_range(can, &sfrom, &sto, &soff);
	ssize = sto - sfrom;
	pdb_extract_subsequence(reference, can->sid, sfrom, sto, FWD, subject);
	int qbeg = INT_MIN, qend = INT_MIN, sbeg = INT_MIN, send = INT_MIN;
	double ident_perc = 0.0;
	
	const BOOL small_edlib = onc_align(read,
									   can->qoff,
									   can->qsize,
									   kstr_str(*subject),
									   soff,
									   ssize,
									   align_data,
									   kOcaBlockSize,
									   min_align_size,
                                       ONC_TAIL_MATCH_LEN_SHORT);
	if (!small_edlib) return FALSE;
	
	BOOL need_rescue = FALSE;
	if (check_mapped_range) {
		qbeg = oca_query_start(*align_data);
		qend = oca_query_end(*align_data);
		int lhang = (qbeg > can->qbeg) ? (qbeg - can->qbeg) : 0;
		int rhang = (qend < can->qend) ? (can->qend - qend) : 0;
		if (lhang + rhang > 500) need_rescue = TRUE;
	}
	if (!need_rescue) {
		qbeg = oca_query_start(*align_data);
		qend = oca_query_end(*align_data);
		sbeg = oca_target_start(*align_data);
		send = oca_target_end(*align_data);
		ident_perc = oca_ident_perc(*align_data);
	} else {
		const BOOL dalign = ocda_go(read, 
									can->qoff,
									can->qsize,
									kstr_str(*subject),
									soff,
									ssize,
									dalign_data,
									min_align_size);
		if (!dalign) return FALSE;
		const BOOL large_edlib = edlib_go(read,
										  ocda_query_start(*dalign_data),
										  ocda_query_end(*dalign_data),
										  kstr_str(*subject),
										  ocda_target_start(*dalign_data),
										  ocda_target_end(*dalign_data),
										  falign_data,
										  ocda_distance(*dalign_data),
										  min_align_size,
										  TRUE,
										  4);
		if (large_edlib) {
			qbeg = falign_data->qoff;
			qend = falign_data->qend;
			sbeg = falign_data->toff;
			send = falign_data->tend;
			ident_perc = falign_data->ident_perc;
		} else {
			return FALSE;
		}
	}
	
	m4->qid = can->qid;
	m4->sid = can->sid;
	m4->ident_perc = ident_perc;
	m4->vscore = can->score;
	m4->qdir = can->qdir;
	m4->qoff = qbeg;
	m4->qend = qend;
	m4->qext = can->qoff;
	m4->qsize = can->qsize;
	m4->sdir = FWD;
	m4->soff = sbeg + sfrom;
	m4->send = send + sfrom;
	m4->sext = can->soff;
	m4->ssize = can->ssize;
	
	if (can->qdir == REV) {
		qbeg = m4->qsize - m4->qend;
		qend = m4->qsize - m4->qoff;
		m4->qoff = qbeg;
		m4->qend = qend;
	}
	
	return TRUE;
}

static void
rm_extend_candidates(GappedCandidate* cans,
				  const int ncan,
				  OcAlignData* align_data,
				  OcDalignData* dalign_data,
				  FullEdlibAlignData* falign_data,
				  const char* fwd_read,
				  const char* rev_read,
				  kstring_t* subject,
				  PackedDB* reads,
				  PackedDB* reference,
				  const int min_align_size,
				  vec_m4* m4list)
{
	ks_introsort_GappedCandidate_RmScoreGT(ncan, cans);
	kv_clear(*m4list);
	M4Record m4;
	for (int i = 0; i < ncan; ++i) {
		if (check_candidate_contain(m4list, cans + i)) continue;
		const char* read = (cans[i].qdir == FWD) ? fwd_read : rev_read;
		BOOL r = rm_extend_candidate(cans + i,
									 align_data,
									 dalign_data,
									 falign_data,
									 read,
									 reference,
									 subject,
									 min_align_size,
									 &m4);
		if (r) kv_push(M4Record, *m4list, m4);
	}
}

void*
rm_search_one_volume(void* arg)
{
	MappingThreadData* mtd 	= (MappingThreadData*)(arg);
	RecordWriter* out 		= mtd->output;
	PackedDB* reads 		= mtd->reads;
	const int num_reads		= PDB_NUM_SEQS(reads);
	PackedDB* reference 	= mtd->reference;
	int read_start_id		= mtd->read_start_id;
	int reference_start_id	= mtd->reference_start_id;
	LookupTable* lktbl		= mtd->lktbl;
	MapOptions* options		= mtd->options;
	
	new_kstring(fwd_read);
	new_kstring(rev_read);
	new_kstring(subject);
	WordFindData* wfdata = new_WordFindData(PDB_SIZE(reference), options->block_size, options->kmer_size, options->block_score_cutoff);
	new_kvec(vec_can, candidates);
	new_kvec(vec_m4, m4list);
	OcAlignData* align_data = new_OcAlignData(options->error);
	OcDalignData* dalign_data = new_OcDalignData(options->error);
	FullEdlibAlignData* falign_data = new_FullEdlibAlignData(options->error);
	
	int sid, eid;
	while (get_next_read_chunk(mtd->chunk_size, mtd->chunk_id, mtd->chunk_lock, num_reads, &sid, &eid)) {
		OC_LOG("mapping read %d --- %d", sid, eid);
		for (int i = sid; i < eid; ++i) {
			//OC_LOG("mapping read %d", i);
			kv_clear(candidates);
			pdb_extract_sequence(reads, i, FWD, &fwd_read);
			find_candidates(kstr_str(fwd_read),
							kstr_size(fwd_read),
							i,
							FWD,
							read_start_id,
							reference_start_id,
							FALSE,
							reference,
							lktbl,
							options,
							wfdata,
							&candidates);
			pdb_extract_sequence(reads, i, REV, &rev_read);
			find_candidates(kstr_str(rev_read),
							kstr_size(rev_read),
							i,
							REV,
							read_start_id,
							reference_start_id,
							FALSE,
							reference,
							lktbl,
							options,
							wfdata,
							&candidates);
			
			ks_introsort(GappedCandidate_RmScoreGT, kv_size(candidates), kv_data(candidates));
			if (kv_size(candidates) > options->num_candidates) kv_resize(GappedCandidate, candidates, options->num_candidates);
				//OC_LOG("number of candidates: %lu", kv_size(candidates));
			rm_extend_candidates(kv_data(candidates),
								  kv_size(candidates),
								  align_data,
								  dalign_data,
								  falign_data,
								  kstr_str(fwd_read),
								  kstr_str(rev_read),
								  &subject,
								  reads,
								  reference,
								  options->align_size_cutoff,
								  &m4list);
			for (size_t k = 0; k != kv_size(m4list); ++k) {
				kv_A(m4list, k).qid += read_start_id;
				kv_A(m4list, k).sid += reference_start_id;
				M4Record* m = &kv_A(m4list, k);
				if (options->use_hdr_as_id) {
					const char* qhdr = PDB_SEQ_NAME(reads, m->qid - read_start_id);
					const char* shdr = PDB_SEQ_NAME(reference, m->sid - reference_start_id);
					RW_DUMP_ONE_DATA(M4Record, DUMP_ASM_M4_HDR_ID, options->binary_output, out, m);
				} else {
					RW_DUMP_ONE_DATA(M4Record, DUMP_ASM_M4, options->binary_output, out, m);
				}
			}
		}
	}
	
	free_kstring(fwd_read);
	free_kstring(rev_read);
	free_kstring(subject);
	free_WordFindData(wfdata);
	kv_destroy(candidates);
	kv_destroy(m4list);
	free_OcAlignData(align_data);
	free_OcDalignData(dalign_data);
	free_FullEdlibAlignData(falign_data);
	return NULL;
}

void
rm_main(MapOptions* options, const char* reads_path, const char* reference_path, const char* output)
{
	PackedDB* reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_NANOPORE);
	LookupTable* lktbl = build_lookup_table(reference, options->kmer_size, options->kmer_cnt_cutoff, options->num_threads);
	const int reference_start_id = 0;
	
	OcMutex out_lock;
	pthread_mutex_init(&out_lock, NULL);
	OcMutex chunk_lock;
	pthread_mutex_init(&chunk_lock, NULL);
	const int num_threads = options->num_threads;
	pthread_t tids[num_threads];
	int read_start_id = 0;
	int chunk_id;
	DFOPEN(out, output, "w");
	MappingThreadData* mdata[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		mdata[i] = new_MappingThreadData(i,
										   options,
										   NULL,
										   read_start_id,
										   reference,
										   reference_start_id,
										   lktbl,
										   out,
										   &out_lock,
										   500,
										   &chunk_id,
										   &chunk_lock);
	}
	
	DGZ_OPEN(reads_in, reads_path, "r");
	kseq_t* read = kseq_init(reads_in);
	char job[1024];
	while (1) {
		PackedDB* reads = load_mapping_reads(read);
		if (reads == NULL) break;
		sprintf(job, "mapping %lu reads", PDB_NUM_SEQS(reads));
		TIMING_START(job);
		chunk_id = 0;
		for (int i = 0; i < num_threads; ++i) {
			mdata[i]->reads = reads;
			mdata[i]->read_start_id = read_start_id;
			pthread_create(tids + i, NULL, rm_search_one_volume, mdata[i]);
		}
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(tids[i], NULL);
		}
		read_start_id += PDB_NUM_SEQS(reads);
		free_PackedDB(reads);
		TIMING_END(job);
	}
	GZ_CLOSE(reads_in);
	kseq_destroy(read);
	
	for (int i = 0; i < num_threads; ++i) {
		free_MappingThreadData(mdata[i]);
	}
	
	free_PackedDB(reference);
	destroy_lookup_table(lktbl);
	FCLOSE(out);
}
