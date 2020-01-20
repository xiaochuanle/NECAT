#include "ctgpm.h"

#include "../common/m4_record.h"
#include "../common/map_aux.h"
#include "../common/map_options.h"
#include "../common/oc_assert.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

static MapOptions		options;
static vec_can	 		candidates;
static pthread_mutex_t 	can_lock;
static size_t 			next_can_id;
static PackedDB* 		pctg = NULL;
static FILE* 			out = NULL;
static pthread_mutex_t	out_lock;
static const char* 		can_path = NULL;
static const char* 		ctg_path = NULL;
static const char* 		output = NULL;

static int
GappedCandidate_PmScoreGT(GappedCandidate a, GappedCandidate b)
{
	if (a.score != b.score) return a.score > b.score;
	if (a.qdir != b.qdir) return a.qdir < b.qdir;
	if (a.sid != b.sid) return a.sid < b.sid;
	if (a.qoff != b.qoff) return a.qoff < b.qoff;
	if (a.soff != b.soff) return a.soff < b.soff;
	return 0;
}

KSORT_INIT(GappedCandidate_PmScoreGT, GappedCandidate, GappedCandidate_PmScoreGT)

#define GappedCandidate_QidLt(a, b) ((a).qid < (b).qid)

KSORT_INIT(GappedCandidate_QidLt, GappedCandidate, GappedCandidate_QidLt)

void
load_candidates();

static void
init_global_data()
{
	oc_assert(can_path);
	oc_assert(ctg_path);
	oc_assert(output);
	
	kv_init(candidates);
	load_candidates();
	pthread_mutex_init(&can_lock, NULL);
	next_can_id = 0;
	pctg = new_PackedDB();
	pdb_load(pctg, ctg_path, TECH_NANOPORE);
	FOPEN(out, output, "w");
	pthread_mutex_init(&out_lock, NULL);
}

static void
destroy_global_data()
{
	kv_destroy(candidates);
	free_PackedDB(pctg);
	FCLOSE(out);
}

static BOOL
extract_next_can_ids(size_t* can_sid, size_t* can_eid)
{
	BOOL r = 1;
	size_t ncan = kv_size(candidates);
	size_t sid, eid;
	
	pthread_mutex_lock(&can_lock);
	if (next_can_id < ncan) {
		sid = next_can_id;
		eid = sid + 1;
		while (eid < ncan && kv_A(candidates, sid).qid == kv_A(candidates, eid).qid) ++eid;
		next_can_id = eid;
		*can_sid = sid;
		*can_eid = eid;
	} else {
		r = 0;
	}
	pthread_mutex_unlock(&can_lock);
	
	return r;
}

static void
dump_candidates(M4Record* m4v, const size_t n)
{
	pthread_mutex_lock(&out_lock);
	for (size_t i = 0; i < n; ++i) {
		int qid = m4v[i].qid;
		int sid = m4v[i].sid;
		const char* qhdr = PDB_SEQ_NAME(pctg, qid);
		const char* shdr = PDB_SEQ_NAME(pctg, sid);
		DUMP_ASM_M4_HDR_ID(fprintf, out, m4v[i]);
	}
	pthread_mutex_unlock(&out_lock);
}

void
load_candidates()
{
	DGZ_OPEN(in, can_path, "r");
	char line[1024];
	GappedCandidate can;
	while (gzgets(in, line, 1024)) {
		LOAD_GAPPED_CANDIDATE(sscanf, line, can);
		kv_push(GappedCandidate, candidates, can);
	}
	GZ_CLOSE(in);
	
	ks_introsort_GappedCandidate_QidLt(kv_size(candidates), kv_data(candidates));
}

static void
validate_candidate(PackedDB* pctg, 
				   GappedCandidate* can,
				   int qid,
				   int qdir,
				   idx qoff,
				   int sid,
				   int sdir,
				   idx soff)
{
	int m = 10;//options.kmer_size;
	idx qb = PDB_SEQ_OFFSET(pctg, qid);
	idx qs = PDB_SEQ_SIZE(pctg, qid);
	idx tb = PDB_SEQ_OFFSET(pctg, sid);
	idx ts = PDB_SEQ_SIZE(pctg, sid);
	
	char qc, tc;
	for (int i = 0; i < m; ++i) {
		idx k = qoff + i;
		if (qdir == REV) {
			k = qs - 1 - k + qb;
			qc = _get_pac(pctg->m_pac, k);
			qc = 3 - qc;
		} else {
			k = k + qb;
			qc = _get_pac(pctg->m_pac, k);
		}
		
		k = soff + i;
		if (sdir == REV) {
			k = ts - 1 - k + tb;
			tc = _get_pac(pctg->m_pac, k);
			tc = 3 - tc;
		} else {
			k = k + tb;
			tc = _get_pac(pctg->m_pac, k);
		}
		if (qc != tc) {
			DUMP_GAPPED_CANDIDATE(fprintf, stdout, *can);
			printf("m = %d, i = %d, qc = %d, tc = %d\n", m, i, (int)qc, (int)tc);
			exit (0);
		}
	}
}

static void
extend_candidates(GappedCandidate* cans,
				  const int ncan,
				  OcDalignData* dalign_data,
				  FullEdlibAlignData* falign_data,
				  kstring_t* query,
				  kstring_t* subject,
				  const int min_align_size,
				  vec_m4* m4list)
{
	ks_introsort_GappedCandidate_PmScoreGT(ncan, cans);
	kv_clear(*m4list);
	M4Record m4;
	for (int i = 0; i < ncan; ++i) {
		if (1) {
			GappedCandidate* can = cans + i;
			validate_candidate(pctg, can, can->qid, can->qdir, can->qoff, can->sid, can->sdir, can->soff);
			validate_candidate(pctg, can, can->qid, can->qdir, can->qend - 13, can->sid, can->sdir, can->send - 13);
		}
		if (check_candidate_contain(m4list, cans + i)) continue;
		idx qoff, qend, toff, tend;
		double ident_perc;
		BOOL r = ctgmp_align(pctg,
							 dalign_data,
							 falign_data,
							 cans + i,
							 query,
							 subject,
							 &qoff,
							 &qend,
							 &toff,
							 &tend,
							 &ident_perc,
							 min_align_size);
		if (r) {
			m4.qid = cans[i].qid;
			m4.sid = cans[i].sid;
			m4.ident_perc = ident_perc;
			m4.vscore = cans[i].score;
			m4.qdir = cans[i].qdir;
			m4.qoff = qoff;
			m4.qend = qend;
			m4.qext = cans[i].qoff;
			m4.qsize = cans[i].qsize;
			m4.sdir = cans[i].sdir;
			m4.soff = toff;
			m4.send = tend;
			m4.sext = cans[i].soff;
			m4.ssize = cans[i].ssize;
			if (m4.qdir == REV) {
				qoff = m4.qsize - m4.qend;
				qend = m4.qsize - m4.qoff;
				m4.qoff = qoff;
				m4.qend = qend;
				m4.qext = m4.qsize - 1 - m4.qext;
			}
			if (m4.sdir == REV) {
				toff = m4.ssize - m4.send;
				tend = m4.ssize - m4.soff;
				m4.soff = toff;
				m4.send = tend;
				m4.sext = m4.ssize - 1 - m4.sext;
			}
			kv_push(M4Record, *m4list, m4);
		}
	}
}

static void*
ctgpm_worker()
{
	new_kstring(query);
	new_kstring(target);
	OcDalignData* dalign_data = new_OcDalignData(0.3);
	FullEdlibAlignData* falign_data = new_FullEdlibAlignData(0.3);
	new_kvec(vec_m4, m4list);
	size_t can_sid, can_eid;
	while (extract_next_can_ids(&can_sid, &can_eid)) {
		GappedCandidate* cans = kv_data(candidates) + can_sid;
		size_t ncan = can_eid - can_sid;
		ncan = OC_MIN(ncan, 500);
		extend_candidates(cans, ncan, dalign_data, falign_data, &query, &target, 1000, &m4list);
		dump_candidates(kv_data(m4list), kv_size(m4list));
	}
	
	free_kvec(m4list);
	free_OcDalignData(dalign_data);
	free_FullEdlibAlignData(falign_data);
	free_kstring(query);
	free_kstring(target);
	
	return NULL;
}

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE\n");
	fprintf(out, "%s [options] ctg-path can-path output\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc < 4) {
		print_usage(argv[0]);
		return 1;
	}

	options = sDefaultPairwiseMapingOptions;
	if (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS) {
		print_usage(argv[0]);
		return 1;
	}
	
	print_MapOptions(&options);
	
	ctg_path 	= argv[argc - 3];
	can_path 	= argv[argc - 2];
	output 		= argv[argc - 1];
	init_global_data();
	
	pthread_t tids[options.num_threads];
	for (int i = 0; i < options.num_threads; ++i) {
		pthread_create(tids + i, NULL, ctgpm_worker, NULL);
	}
	for (int i = 0; i < options.num_threads; ++i) {
		pthread_join(tids[i], NULL);
	}
	
	destroy_global_data();
}
