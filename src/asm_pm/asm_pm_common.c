#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <limits.h>

#include "../common/m4_record.h"
#include "../common/makedb_aux.h"
#include "../common/map_aux.h"
#include "../common/oc_assert.h"
#include "../common/record_writer.h"
#include "../common/packed_db.h"
#include "../lookup_table/lookup_table.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"
#include "find_mem.h"
#include "hbn_align.h"

#define ZV 1000
static int BC = 10;
#define SM 60
#define SI 61
#define MAXC 100
#define PLL 50
#define ErrorRate 0.5

#define KMER_SIZE	13
#define MAX_KMER_OCC	500
#define KMER_CNT_CUTOFF	5
#define BLOCK_SCORE_CUTOFF	4

pthread_mutex_t mutilock; 
int runnumber=0, readcount, terminalnum;
int seed_len = KMER_SIZE;
PackedDB* reference = 0;
int reference_start_id;
PackedDB* reads = 0;
int read_start_id;
LookupTable* lktbl = 0;

static FILE* out = 0;

static pthread_mutex_t m4_sink_lock;
static M4Record* m4_sink;
static int max_num_m4;
static int num_m4;
static int output_ident_perc = 1;
static int binary_output = 0;
static int num_extended_can = MAXC;

static void
alloc_m4_sink()
{
	const size_t mem = 64 * (1<<20);
	max_num_m4 = mem / sizeof(M4Record);
	m4_sink = (M4Record*)malloc(sizeof(M4Record) * max_num_m4);
	num_m4 = 0;
	pthread_mutex_init(&m4_sink_lock, NULL);
}

static void
dump_m4_sink()
{
	if (binary_output) {
		for (int i = 0; i < num_m4; ++i) {
			if (!output_ident_perc) m4_sink[i].ident_perc = (100.0 - m4_sink[i].ident_perc) / 100.0;
			++m4_sink[i].qid;
			++m4_sink[i].sid;
		}
		FWRITE(m4_sink, sizeof(M4Record), num_m4, out);
	} else {
		for (int i = 0; i < num_m4; ++i) {
			if (!output_ident_perc) m4_sink[i].ident_perc = (100.0 - m4_sink[i].ident_perc) / 100.0;
			M4Record* m = m4_sink + i;
			oc_assert(m->qid >= read_start_id);
			const char* qhdr = PDB_SEQ_NAME(reads, m->qid - read_start_id);
			const char* shdr = PDB_SEQ_NAME(reference, m->sid - reference_start_id);
			DUMP_ASM_M4_HDR_ID(fprintf, out, *m); 
		}
	}
	num_m4 = 0;
}

static void
dump_m4_list(vec_m4* m4list)
{
	pthread_mutex_lock(&m4_sink_lock);
	
	int n = kv_size(*m4list);
	if (n + num_m4 > max_num_m4) dump_m4_sink();
	memcpy(m4_sink + num_m4, kv_data(*m4list), sizeof(M4Record) * n);
	num_m4 += n;
	
	pthread_mutex_unlock(&m4_sink_lock);
}

static void
free_m4_sink()
{
	dump_m4_sink();
	free(m4_sink);
	m4_sink = 0;
	max_num_m4 = 0;
	num_m4 = 0;
}

typedef struct{
	int readno;
	int score;
	char chain;
	int target_start, target_size, target_id;
	int query_start;
}AsmGappedCandidate;

#define DUMP_ASM_CAN(can) fprintf(stderr, "qid = %d, sid = %d, qstart = %d,sstart = %d, score = %d\n", \
    (can).readno, (can).target_id, (can).query_start, (can).target_start, (can).score)

typedef kvec_t(AsmGappedCandidate) vec_asm_can;

static void
fill_gapped_candidate(AsmGappedCandidate* asm_can, GappedCandidate* can)
{
	can->qid = asm_can->readno;
	can->sid = asm_can->target_id;
	can->qdir = asm_can->chain;
	can->sdir = FWD;
	can->score = asm_can->score;
	can->qoff = asm_can->query_start;
	can->soff = asm_can->target_start;
	can->ssize = asm_can->target_size;
}

static void
normalise_candidate_qoff(GappedCandidate* can)
{
    if (can->qdir == REV) {
        can->qdir = FWD;
        can->sdir = 1 - can->sdir;
    }
}

static int AsmGappedCandidate_ScoreGT(AsmGappedCandidate a, AsmGappedCandidate b)
{
	if (a.score != b.score) return a.score > b.score;
	if (a.chain != b.chain) return a.chain < b.chain;
	if (a.target_id != b.target_id) return a.target_id < b.target_id;
	if (a.query_start != b.query_start) return a.query_start < b.query_start;
	if (a.target_start != b.target_start) return a.target_start < b.target_start;
	return 0;
}

KSORT_INIT(AsmGappedCandidate_ScoreGT, AsmGappedCandidate, AsmGappedCandidate_ScoreGT)

static BOOL
check_asm_candidate_contain(vec_m4* m4list, GappedCandidate* can)
{
	for (size_t i = 0; i != kv_size(*m4list); ++i) {
		M4Record* m = kv_data(*m4list) + i;
		BOOL r = (can->qdir == m->qdir)
				 &&
				 (can->sid == m->sid)
				 &&
				 (can->sdir == m->sdir);
		if (r) return TRUE;
	}
	return FALSE;
}

#if 0
static void
extend_candidates(ChainWorkData* chain_data,
				  OcDalignData* dalign_data,
				  OcDalignData* sensitive_dalign_data,
				  AsmGappedCandidate* cans,
				  const int ncan,
				  const char* fwd_read,
				  const int read_id,
				  const int read_size,
				  kstring_t* target,
				  vec_m4* m4list)
{
	kv_clear(*m4list);
	if (ncan == 0) return;
	ks_introsort(AsmGappedCandidate_ScoreGT, ncan, cans);
	GappedCandidate can, tmp_can;
	M4Record m4;
	int num_m4 = 0;
	new_kvec(vec_u8, read);
	new_kvec(vec_u8, subject);
	kv_resize(u8, read, read_size);
	memcpy(kv_data(read), fwd_read, read_size);
	const int kmer_size = 10;
	const int window_size = 6;
	const int mem_size = 15;
	new_kvec(vec_kmif, kmif_list);
	build_kmif_list(&read, kmer_size, 1, &kmif_list);
	sort_kmif_list(&kmif_list);
	int r;
	new_kvec(vec_mem, chain_mems);
	int min_align_size = 400;

	for (int i = 0; i < ncan && i < num_extended_can; ++i) {
		fill_gapped_candidate(cans + i, &can);
		can.qsize = read_size;
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
		normalise_candidate_qoff(&can);
		if (check_asm_candidate_contain(m4list, &can)) continue;
		pdb_extract_sequence(reference, can.sid, can.sdir, target);
		kv_resize(u8, subject, kstr_size(*target));
		memcpy(kv_data(subject), kstr_str(*target), kstr_size(*target));
		r = compute_align_range_1(chain_data, 
				&subject, 
				&read, 
				&kmif_list, 
				kmer_size, 
				window_size, 
				mem_size, 
				&tmp_can, 
				&chain_mems);
		if (!r) continue;
		can.qbeg = tmp_can.sbeg;
		can.qend = tmp_can.send;
		can.qoff = tmp_can.soff;
		can.sbeg = tmp_can.qbeg;
		can.send = tmp_can.qend;
		can.soff = tmp_can.qoff;
		can.score = tmp_can.score;
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
//OC_LOG("compute align");
		r = ocda_go(fwd_read,
				can.qoff,
				can.qsize,
				kstr_str(*target),
				can.soff,
				can.ssize,
				dalign_data,
				min_align_size);
		if (!r) continue;
		int qb = ocda_query_start(*dalign_data);
        int qe = ocda_query_end(*dalign_data);
        int tb = ocda_target_start(*dalign_data);
        int te = ocda_target_end(*dalign_data);
		double ident_perc = ocda_ident_perc(*dalign_data);
        int q_delta_l = (qb > tmp_can.sbeg) ? (qb - tmp_can.sbeg) : 0;
        int q_delta_r = (qe < tmp_can.send) ? (tmp_can.send - qe) : 0;
        int t_delta_l = (tb > tmp_can.qbeg) ? (tb - tmp_can.qbeg) : 0;
        int t_delta_r = (te < tmp_can.qend) ? (tmp_can.qend - te) : 0;
        if ((q_delta_l + q_delta_r > 300) && (t_delta_l + t_delta_r > 300)) r = 0;
		//const int E = 200;
		//r = !((q_delta_l > E) || (q_delta_r > E) || (t_delta_l > E) || (t_delta_r > E));

		if (!r) {
			r = ocda_go(fwd_read,
					can.qoff,
					can.qsize,
					kstr_str(*target),
					can.soff,
					can.ssize,
					sensitive_dalign_data,
					min_align_size);		
			if (r) {
				int print = 0;
				//if (ident_perc > 90.0) print = 1;
				if (print) printf("[%d, %d, %d] x [%d, %d, %d], %g\n", qb, qe, can.qsize, tb, te, can.ssize, ident_perc);
				qb = ocda_query_start(*sensitive_dalign_data);
        		qe = ocda_query_end(*sensitive_dalign_data);
        		tb = ocda_target_start(*sensitive_dalign_data);
        		te = ocda_target_end(*sensitive_dalign_data);	
				ident_perc = ocda_ident_perc(*sensitive_dalign_data);	
				if (print) DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
				if (print) printf("[%d, %d, %d] x [%d, %d, %d], %g\n", qb, qe, can.qsize, tb, te, can.ssize, ident_perc);		
			}	
		}
		if (!r) continue;
		
		++num_m4;
		m4.qid = read_id;
		m4.sid = can.sid;
		m4.ident_perc = ident_perc;
		m4.vscore = can.score;
		m4.qdir = can.qdir;
		m4.qoff = qb;
		m4.qend = qe;
		m4.qext = can.qoff;
		m4.qsize = read_size;
		m4.sdir = can.sdir;
		m4.soff = tb;
		m4.send = te;
		m4.sext = can.soff;
		m4.ssize = can.ssize;
		kv_push(M4Record, *m4list, m4);
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
		//DUMP_M4_RECORD(fprintf, stdout, m4);
	}
	
	for (int i = 0; i < num_m4; ++i) {
		M4Record* pm4 = &kv_A(*m4list, i);
		oc_assert(pm4->qid >= 0);
		pm4->qid += read_start_id;
		pm4->sid += reference_start_id;
		if (pm4->qdir == REV) {
			int qoff = read_size - pm4->qend;
			int qend = read_size - pm4->qoff;
			pm4->qoff = qoff;
			pm4->qend = qend;
		}
        if (pm4->sdir == REV) {
            int soff = pm4->ssize - pm4->send;
            int send = pm4->ssize - pm4->soff;
            pm4->soff = soff;
            pm4->send = send;
        }
	}

	kv_destroy(read);
	kv_destroy(subject);
	kv_destroy(kmif_list);
	kv_destroy(chain_mems);
}
#endif

static void
extend_candidates(ChainWorkData* chain_data,
				  HbnMapData* map_data,
				  AsmGappedCandidate* cans,
				  const int ncan,
				  const char* fwd_read,
				  const int read_id,
				  const int read_size,
				  kstring_t* target,
				  vec_m4* m4list)
{
	kv_clear(*m4list);
	if (ncan == 0) return;
	ks_introsort(AsmGappedCandidate_ScoreGT, ncan, cans);
	GappedCandidate can, tmp_can;
	M4Record m4;
	int num_m4 = 0;
	new_kvec(vec_u8, read);
	new_kvec(vec_u8, subject);
	kv_resize(u8, read, read_size);
	memcpy(kv_data(read), fwd_read, read_size);
	const int kmer_size = 10;
	const int window_size = 6;
	const int mem_size = 15;
	new_kvec(vec_kmif, kmif_list);
	build_kmif_list(&read, kmer_size, 1, &kmif_list);
	sort_kmif_list(&kmif_list);
	int r;
	new_kvec(vec_mem, chain_mems);
	int min_align_size = 400;
	double min_ident_perc = 65.0;
	new_kstring(qaln);
	new_kstring(taln);

	for (int i = 0; i < ncan && i < num_extended_can; ++i) {
		fill_gapped_candidate(cans + i, &can);
		can.qsize = read_size;
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
		normalise_candidate_qoff(&can);
		if (check_asm_candidate_contain(m4list, &can)) continue;
		pdb_extract_sequence(reference, can.sid, can.sdir, target);
		kv_resize(u8, subject, kstr_size(*target));
		memcpy(kv_data(subject), kstr_str(*target), kstr_size(*target));
		r = compute_align_range_1(chain_data, 
				&subject, 
				&read, 
				&kmif_list, 
				kmer_size, 
				window_size, 
				mem_size, 
				&tmp_can, 
				&chain_mems);
		if (!r) continue;
		can.qbeg = tmp_can.sbeg;
		can.qend = tmp_can.send;
		can.qoff = tmp_can.soff;
		can.qsize = read_size;
		can.sbeg = tmp_can.qbeg;
		can.send = tmp_can.qend;
		can.soff = tmp_can.qoff;
		can.ssize = kv_size(subject);
		can.score = tmp_can.score;
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
		int qb, qe, tb, te;
		double ident_perc;
		r = hbn_map_extend(map_data,
				kv_data(read),
				(const u8*)kstr_str(*target),
				&can,
				min_align_size,
				min_ident_perc,
				&qb,
				&qe,
				&tb,
				&te,
				&ident_perc,
				&qaln,
				&taln);
		if (!r) continue;
		
		++num_m4;
		m4.qid = read_id;
		m4.sid = can.sid;
		m4.ident_perc = ident_perc;
		m4.vscore = can.score;
		m4.qdir = can.qdir;
		m4.qoff = qb;
		m4.qend = qe;
		m4.qext = can.qoff;
		m4.qsize = read_size;
		m4.sdir = can.sdir;
		m4.soff = tb;
		m4.send = te;
		m4.sext = can.soff;
		m4.ssize = can.ssize;
		kv_push(M4Record, *m4list, m4);
		//DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
		//fprintf(stderr, "%d\t", num_m4); DUMP_M4_RECORD(fprintf, stdout, m4);
	}
	
	for (int i = 0; i < num_m4; ++i) {
		M4Record* pm4 = &kv_A(*m4list, i);
		oc_assert(pm4->qid >= 0);
		pm4->qid += read_start_id;
		pm4->sid += reference_start_id;
		if (pm4->qdir == REV) {
			int qoff = read_size - pm4->qend;
			int qend = read_size - pm4->qoff;
			pm4->qoff = qoff;
			pm4->qend = qend;
		}
        if (pm4->sdir == REV) {
            int soff = pm4->ssize - pm4->send;
            int send = pm4->ssize - pm4->soff;
            pm4->soff = soff;
            pm4->send = send;
        }
	}

	kv_destroy(read);
	kv_destroy(subject);
	kv_destroy(kmif_list);
	kv_destroy(chain_mems);
	free_kstring(qaln);
	free_kstring(taln);
}

struct Back_List{
	int score,loczhi[SM],seedno[SM],seednum;
	int index;
};

int 
extract_hash_values(char *seqm,vec_int* hash_list, int *endn,int len_str,int readnum)
{
	int eit=0,temp;
	int i,j,start,num;
	num=(len_str-readnum)/BC+1;
	*endn=(len_str-readnum)%BC;
	kv_resize(int, *hash_list, num);
	int* value = kv_data(*hash_list);
	for(i=0;i<num;i++){
		eit=0;start=i*BC;
		for(j=0;j<readnum;j++){
            temp=seqm[start+j];
			if(temp==4){eit=-1;break;}
			eit=eit<<2;
			eit=eit+temp;
		}
		value[i]=eit;
	}
	return(num);
}

int find_location(int *t_loc,int *t_seedn,int *t_score,int *loc,int k,int *rep_loc,float len,int read_len1){
*rep_loc = INT_MAX;
	int i,j,maxval=0,maxi = 0,rep=0,lasti = 0,tempi = 0;
	for (i = 0; i < k; ++i) oc_assert(t_seedn[i] > 0);
	for(i=0;i<k;i++)t_score[i]=0;
	for(i=0;i<k-1;i++)for(j=i+1,tempi=t_seedn[i];j<k;j++)if(tempi!=t_seedn[j]&&t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<0.10){t_score[i]++;t_score[j]++;tempi=t_seedn[j];}
	
	for(i=0;i<k;i++){
		if(maxval<t_score[i]){maxval=t_score[i];maxi=i;rep=0;}
		else if(maxval==t_score[i]){rep++;lasti=i;}
	}
	for(i=0;i<4;i++)loc[i]=0;
	if(maxval>=5&&rep==maxval){loc[0]=t_loc[maxi],loc[1]=t_seedn[maxi];*rep_loc=maxi;loc[2]=t_loc[lasti],loc[3]=t_seedn[lasti];return(1);}
	else if(maxval>=5&&rep!=maxval){
		for(j=0;j<maxi;j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_len1&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*len)-1)<0.10){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
			}
		j=maxi;
		if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
		else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
		for(j=maxi+1;j<k;j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_len1&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*len)-1)<0.10){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
			}
		return(1);
	}
	else return(0);
}

void* pairwise_mapping(){
	int cleave_num,read_len,s_k;
	new_kvec(vec_int, hash_list);
	kv_reserve(int, hash_list, 1000000);
	int *mvalue,*leadarray,flag_end,u_k;
	int count1=0,i,j,k;
	struct Back_List *database,*temp_spr,*temp_spr1;
	int location_loc[4],repeat_loc,*index_list,*index_spr;
	int *index_score,*index_ss;
	int temp_list[150],temp_seedn[150],temp_score[150],start_loc;
	int endnum,ii;
	int loc,templong,localnum,read_i,read_end;
	int cc1,readno,loc_seed,loc_list;
	int length1,num1,num2;
	int seedcount;
	AsmGappedCandidate canidate_temp;
	new_kvec(vec_asm_can, candidates);
	new_kvec(vec_m4, m4list);
	int readstart1,readend1,left_length1,right_length1,left_length2,right_length2;
	
	ChainWorkData* chain_data = chain_data_new();
	HbnMapData* map_data = hbn_map_data_new();
	new_kstring(fwd_read);
	new_kstring(rev_read);
	new_kstring(target);
	kstring_t* read;
	
	j = PDB_SIZE(reference) / ZV+5;
	index_list=(int *)malloc(j*sizeof(int));
	index_score=(int *)malloc(j*sizeof(int));
	database=(struct Back_List *)malloc(j*sizeof(struct Back_List));
	for(i=0,temp_spr=database;i<j;temp_spr++,i++){temp_spr->score=0;temp_spr->index=-1;}
	
	while(1){
        pthread_mutex_lock(&mutilock);
        localnum=runnumber;
        runnumber++;
        pthread_mutex_unlock(&mutilock);
        if(localnum>=terminalnum) break;
		if(localnum==terminalnum-1)read_end=readcount;
		else read_end = (localnum + 1)*PLL;
		OC_LOG("mapping %d --- %d (%d)", localnum*PLL, read_end, readcount);
		for(read_i=localnum*PLL;read_i<read_end;read_i++) {
			//if (read_i + read_start_id != 1002) continue;
			//OC_LOG("mapping read %d", read_i);
			int max_tid = reference_start_id + PDB_NUM_SEQS(reference);
			int soff_max = INT_MAX;
			if (read_i + read_start_id >= reference_start_id && read_i + read_start_id < max_tid) {
				soff_max = PDB_SEQ_OFFSET(reference, read_i);	
			}

			pdb_extract_sequence(reads, read_i, FWD, &fwd_read);
			pdb_extract_sequence(reads, read_i, REV, &rev_read);
			read_len = PDB_SEQ_SIZE(reads, read_i);
			//OC_LOG("len = %d", read_len);
			kv_clear(candidates);
			for(ii=0;ii<2;ii++) {
				read = ii ? (&rev_read) : (&fwd_read);
				cleave_num = extract_hash_values(kstr_str(*read),&hash_list,&endnum,read_len,seed_len);
				mvalue = kv_data(hash_list);
				j=0;
				index_spr=index_list;
				index_ss=index_score;
				endnum=0;
				for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
					{
						u64* match_list;
						u64 match_list_size;
						match_list = extract_kmer_list(lktbl, mvalue[k], &match_list_size);
						assert(match_list_size <= MAX_KMER_OCC);
						if (0) for (u64 np = 0; np != match_list_size; ++np) {
								int qidx = k * BC;
								int tidx = match_list[np];
								for (int kp = 0; kp < seed_len; ++kp) {
									int qc = kstr_A(*read, qidx); 
									int tc = _get_pac(reference->m_pac, tidx);
									oc_assert(qc == tc, "read %d, size = %d, strand %d, k = %d, np = %d, kp = %d, qc = %d, tc = %d",
											  read_i, read_len, ii, k, np, kp, qc, tc);
									++qidx;
									++tidx;
								}
							}
						count1 = match_list_size;
						for(i=0; i<count1; i++,leadarray++)
						{
							if (match_list[i] >= soff_max) continue;
							templong=(match_list[i])/ZV;
							u_k=(match_list[i])%ZV;
							if(templong>=0)
							{
								temp_spr=database+templong;
								if(temp_spr->score==0||temp_spr->seednum<k+1)
								{
									loc=++(temp_spr->score);
									if (temp_spr->score > SM) temp_spr->score = SM;
									if(loc<=SM)
									{
										temp_spr->loczhi[loc-1]=u_k;
										temp_spr->seedno[loc-1]=k+1;
									}
									if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
									else s_k=temp_spr->score;
									if(endnum<s_k)endnum=s_k;
									if(temp_spr->index==-1)
									{
										*(index_spr++)=templong;
										*(index_ss++)=s_k;
										temp_spr->index=j;
										j++;
									}
									else index_score[temp_spr->index]=s_k;
								}
								temp_spr->seednum=k+1;
							}
						}
					}
				
				cc1=j;
				for(i=0,index_spr=index_list,index_ss=index_score;i<cc1;i++,index_spr++,index_ss++)if(*index_ss>KMER_CNT_CUTOFF){
						temp_spr=database+*index_spr;
						if(temp_spr->score==0)continue;
						s_k=temp_spr->score;
						start_loc=(*index_spr)*ZV;
						if(*index_spr>0){loc=(temp_spr-1)->score;if(loc>0)start_loc=(*index_spr-1)*ZV;}
						else loc=0;
						if(loc==0)for(j=0,u_k=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j];temp_seedn[u_k]=temp_spr->seedno[j];
oc_assert(temp_spr->seedno[j]>0, "i = %d, j = %d, score = %d, seedno = %d, bid = %d", i, j, temp_spr->score, temp_spr->seedno[j], *index_spr); u_k++;}
						else {
							k=loc;
							u_k=0;
							temp_spr1=temp_spr-1;
							for(j=0;j<k&&j<SM;j++){temp_list[u_k]=temp_spr1->loczhi[j];temp_seedn[u_k]=temp_spr1->seedno[j];oc_assert(temp_spr1->seedno[j] > 0, "j = %d, k = %d, seedno = %d, bid = %d\n", j, k, temp_spr1->seedno[j], *index_spr); u_k++;}
							for(j=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j]+ZV;temp_seedn[u_k]=temp_spr->seedno[j];oc_assert(temp_spr->seedno[j] > 0, "j = %d, k = %d, seedno = %d, bid = %d", j, k, temp_spr->seedno[j], *index_spr); u_k++;}
						}
						flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len);
						if(flag_end==0)continue;
						if(temp_score[repeat_loc]<BLOCK_SCORE_CUTOFF)continue;
						canidate_temp.score=temp_score[repeat_loc];
						loc_seed=temp_seedn[repeat_loc];
						location_loc[0]=start_loc+location_loc[0];
						loc_list=location_loc[0];
						readno = pdb_offset_to_id(reference, location_loc[0]);
						readstart1 = PDB_SEQ_OFFSET(reference, readno);
						length1 = PDB_SEQ_SIZE(reference, readno);
						readend1 = readstart1 + length1;
						
						//if(readno + reference_start_id > read_i + read_start_id)continue;
						if(readno + reference_start_id == read_i + read_start_id){
							u_k=readstart1/ZV;temp_spr=database+u_k;s_k=readstart1%ZV;
							for(j=0,k=0;j<temp_spr->score&&j<SM;j++)if(temp_spr->loczhi[j]<s_k){temp_spr->loczhi[k]=temp_spr->loczhi[j];k++;}
							temp_spr->score=k;
							for(temp_spr++,u_k++,k=readend1/ZV;u_k<k;u_k++,temp_spr++)temp_spr->score=0;
							for(j=0,k=0,s_k=readend1%ZV;j<temp_spr->score&&j<SM;j++)if(temp_spr->loczhi[j]>s_k){temp_spr->loczhi[k]=temp_spr->loczhi[j];k++;}
							temp_spr->score=k;
							continue;
						}

						canidate_temp.readno=read_i;
						//canidate_temp.readstart=readstart1;
						location_loc[1]=(location_loc[1]-1)*BC;
						canidate_temp.target_id = readno;
						canidate_temp.target_start = location_loc[0] - readstart1;
						canidate_temp.target_size = length1;
						canidate_temp.query_start = location_loc[1];
						
						left_length1=location_loc[0]-readstart1+seed_len-1;right_length1=readend1-location_loc[0];
						left_length2=location_loc[1]+seed_len-1;right_length2=read_len-location_loc[1];
						if(left_length1>=left_length2)num1=left_length2;
						else num1=left_length1;
						if(right_length1>=right_length2)num2=right_length2;
						else num2=right_length1;
						if(num1+num2<400)continue;
						seedcount=0;
						//find all left seed 
						for(u_k=*index_spr-2,k=num1/ZV,temp_spr1=temp_spr-2;u_k>=0&&k>=0;temp_spr1--,k--,u_k--)if(temp_spr1->score>0){
								start_loc=u_k*ZV;
								for(j=0,s_k=0;j<temp_spr1->score && j < SM;j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<0.10){
										seedcount++;s_k++;
									}
								if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
							}
						//find all right seed
						for(u_k=*index_spr+1,k=num2/ZV,temp_spr1=temp_spr+1;k>0;temp_spr1++,k--,u_k++)if(temp_spr1->score>0){
								start_loc=u_k*ZV;
								for(j=0,s_k=0;j<temp_spr1->score && j < SM;j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<0.10){
										seedcount++;s_k++;
									}
								if(s_k*1.0/temp_spr1->score>0.4)temp_spr1->score=0;
							}
						canidate_temp.score=canidate_temp.score+seedcount;
						canidate_temp.chain = ii;
						kv_push(AsmGappedCandidate, candidates, canidate_temp);
					}
				for(i=0,index_spr=index_list;i<cc1;i++,index_spr++){database[*index_spr].score=0;database[*index_spr].index=-1;}
			} // for (ii = 0; ii < 2; ++ii)
			
			extend_candidates(chain_data,
							  map_data,
							  kv_data(candidates),
							  kv_size(candidates),
							  kstr_str(fwd_read),
							  read_i,
							  read_len,
							  &target,
							  &m4list);
			
			dump_m4_list(&m4list);			
		} // for(read_i=localnum*PLL;read_i<read_end;read_i++)
	} // while(1)
	
	free_kstring(fwd_read);
	free_kstring(rev_read);
	free_kstring(target);
	free_kvec(candidates);
	free_kvec(m4list);
	hbn_map_data_free(map_data);
	chain_data_free(chain_data);
	free(index_list);
	free(index_score);
	free(database);
	kv_destroy(hash_list);
	
	return NULL;
}
