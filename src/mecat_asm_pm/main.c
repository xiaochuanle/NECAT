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
#include "../gapped_align/oc_aligner.h"
#include "../lookup_table/lookup_table.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

#define ZV 1000
#define BC 10
#define SM 60
#define SI 61
#define MAXC 100
#define PLL 500
#define ErrorRate 0.2

#define KMER_SIZE	15
#define MAX_KMER_OCC	1000
#define KMER_CNT_CUTOFF	6
#define BLOCK_SCORE_CUTOFF	5

pthread_t *thread; 
int threadnum;
pthread_mutex_t mutilock; 
int runnumber=0, readcount, terminalnum;
char *savework,workpath[300];
int seed_len;
PackedDB* reference;
int reference_start_id;
PackedDB* reads;
int read_start_id;
LookupTable* lktbl;

static FILE* out;

static pthread_mutex_t m4_sink_lock;
static M4Record* m4_sink;
static int max_num_m4;
static int num_m4;

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
	for (int i = 0; i < num_m4; ++i) DUMP_ASM_M4(fprintf, out, m4_sink[i]); //{ dump_asm_m4(out, m4_sink + i); }
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

typedef kvec_t(AsmGappedCandidate) vec_asm_can;
#define AsmGappedCandidate_ScoreGT(a, b) ((a).score > (b).score)
KSORT_INIT(AsmGappedCandidate_ScoreGT, AsmGappedCandidate, AsmGappedCandidate_ScoreGT)

static BOOL
check_asm_candidate_contain(vec_m4* m4list, AsmGappedCandidate* can)
{
	for (size_t i = 0; i != kv_size(*m4list); ++i) {
		M4Record* m = kv_data(*m4list) + i;
		BOOL r = (can->chain == m->qdir)
				 &&
				 (can->target_id == m->sid)
				 &&
				 (can->query_start >= m->qoff && can->query_start <= m->qend)
				 &&
				 (can->target_start >= m->soff && can->target_start <= m->send);
		if (r) return TRUE;
	}
	return FALSE;
}

static void
extend_candidates(AsmGappedCandidate* cans,
				  const int ncan,
				  OcAlignData* align_data,
				  const char* fwd_read,
				  const char* rev_read,
				  const int read_id,
				  const int read_size,
				  kstring_t* target,
				  vec_m4* m4list)
{
	ks_introsort(AsmGappedCandidate_ScoreGT, ncan, cans);
	kv_clear(*m4list);
	M4Record m4;
	int num_m4 = 0;
	for (int i = 0; i < ncan && i < MAXC; ++i) {
		if (check_asm_candidate_contain(m4list, cans + i)) continue;
		const char* read = (cans[i].chain == FWD) ? fwd_read : rev_read;
		pdb_extract_sequence(reference, cans[i].target_id, FWD, target);
		BOOL r = onc_align(read, 
						   cans[i].query_start, 
						   read_size, 
						   kstr_str(*target), 
						   cans[i].target_start, 
						   cans[i].target_size, 
						   align_data, 
						   400);
		if (!r) continue;
		++num_m4;
		m4.qid = read_id;
		m4.sid = cans[i].target_id;
		m4.ident_perc = align_data->ident_perc;
		m4.vscore = cans[i].score;
		m4.qdir = cans[i].chain;
		m4.qoff = align_data->qoff;
		m4.qend = align_data->qend;
		m4.qext = cans[i].query_start;
		m4.qsize = read_size;
		m4.sdir = FWD;
		m4.soff = align_data->toff;
		m4.send = align_data->tend;
		m4.sext = cans[i].target_start;
		m4.ssize = cans[i].target_size;
		kv_push(M4Record, *m4list, m4);
	}
	
	for (int i = 0; i < num_m4; ++i) {
		M4Record* pm4 = &kv_A(*m4list, i);
		pm4->qid += read_start_id + 1;
		pm4->sid += reference_start_id + 1;
		if (pm4->qdir == REV) {
			int qoff = read_size - pm4->qend;
			int qend = read_size - pm4->qoff;
			pm4->qoff = qoff;
			pm4->qend = qend;
		}
		
	}
}

struct Back_List{
	short int score,loczhi[SM],seedno[SM],seednum;
	int index;
};

int 
extract_hash_values(char *seqm,int *value,int *endn,int len_str,int readnum)
{
	int eit=0,temp;
	int i,j,start,num;
	num=(len_str-readnum)/BC+1;
	*endn=(len_str-readnum)%BC;
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
	int i,j,maxval=0,maxi = 0,rep=0,lasti = 0,tempi = 0;
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
	int mvalue[50000],*leadarray,flag_end,u_k;
	int count1=0,i,j,k;
	struct Back_List *database,*temp_spr,*temp_spr1;
	int location_loc[4],repeat_loc,*index_list,*index_spr;
	short int *index_score,*index_ss;
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
	
	OcAlignData* align_data = new_OcAlignData(ErrorRate);
	new_kstring(fwd_read);
	new_kstring(rev_read);
	new_kstring(target);
	kstring_t* read;
	
	j = PDB_SIZE(reference) / ZV+5;
	index_list=(int *)malloc(j*sizeof(int));
	index_score=(short int *)malloc(j*sizeof(short int));
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
fprintf(stdout, "mapping read %d - %d\n", localnum * PLL, read_end);
		for(read_i=localnum*PLL;read_i<read_end;read_i++) {
int max_tid = reference_start_id + PDB_NUM_SEQS(reference);
int soff_max = INT_MAX;
if (read_i + read_start_id >= reference_start_id && read_i + read_start_id < max_tid) {
	soff_max = PDB_SEQ_OFFSET(reference, read_i);	
}

			pdb_extract_sequence(reads, read_i, FWD, &fwd_read);
			pdb_extract_sequence(reads, read_i, REV, &rev_read);
			read_len = PDB_SEQ_SIZE(reads, read_i);
			kv_clear(candidates);
			for(ii=0;ii<2;ii++) {
				read = ii ? (&rev_read) : (&fwd_read);
				cleave_num = extract_hash_values(kstr_str(*read),mvalue,&endnum,read_len,seed_len);
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
						if(loc==0)for(j=0,u_k=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j];temp_seedn[u_k]=temp_spr->seedno[j];u_k++;}
						else {
							k=loc;
							u_k=0;
							temp_spr1=temp_spr-1;
							for(j=0;j<k&&j<SM;j++){temp_list[u_k]=temp_spr1->loczhi[j];temp_seedn[u_k]=temp_spr1->seedno[j];u_k++;}
							for(j=0;j<s_k&&j<SM;j++){temp_list[u_k]=temp_spr->loczhi[j]+ZV;temp_seedn[u_k]=temp_spr->seedno[j];u_k++;}
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

						canidate_temp.readno=readno;
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
			
			extend_candidates(kv_data(candidates),
							  kv_size(candidates),
							  align_data,
							  kstr_str(fwd_read),
							  kstr_str(rev_read),
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
	free_OcAlignData(align_data);
	free(index_list);
	free(index_score);
	free(database);
	
	return NULL;
}

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s [map options] wrk_dir volume_id output\n", prog);
}

int main(int argc, char* argv[])
{
	MapOptions options;
	if (argc < 4 || (parse_MapOptions(argc - 3, argv, &options) != ARG_PARSE_SUCCESS)) {
		print_usage(argv[0]);
		return 1;
	}
print_MapOptions(&options);
	
	const char* wrk_dir = argv[argc - 3];
	const int vid = atoi(argv[argc - 2]);
	const char* output = argv[argc - 1];
	out = fopen(output, "w");
	alloc_m4_sink();
	pthread_t jobs[options.num_threads];
	char job_name[1024];
	
	VolumesInfo* volumes = load_volumes_info(wrk_dir);
	const char* reference_path = vi_volume_name(volumes, vid);
	reference = new_PackedDB();
	pdb_load(reference, reference_path, TECH_NANOPORE);
	lktbl = build_lookup_table(reference, options.kmer_size, options.kmer_cnt_cutoff, options.num_threads);
	seed_len = options.kmer_size;
	reference_start_id = kv_A(volumes->read_start_id, vid);
	
	for (int i = vid; i < volumes->num_volumes; ++i) {
		sprintf(job_name, "pairwise mapping v%d vs v%d", vid, i);
		TIMING_START(job_name);
		const char* reads_path = vi_volume_name(volumes, i);
		reads = new_PackedDB();
		pdb_load(reads, reads_path, TECH_NANOPORE);
		read_start_id = kv_A(volumes->read_start_id, i);
		runnumber=0;
		readcount = PDB_NUM_SEQS(reads);
		terminalnum = (readcount + PLL - 1) / PLL;
		
		for (int k = 0; k < options.num_threads; ++k) {
			pthread_create(jobs + k, NULL, pairwise_mapping, NULL);
		}
		for (int k = 0; k < options.num_threads; ++k) {
			pthread_join(jobs[k], NULL);
		}
		reads = free_PackedDB(reads);
		TIMING_END(job_name);
	}
	
	free_m4_sink();
	fclose(out);
	lktbl = destroy_lookup_table(lktbl);
	reference = free_PackedDB(reference);
	volumes = destroy_volumes_info(volumes);
}
