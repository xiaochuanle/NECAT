#include "word_finder.h"

#include <math.h>
#include <assert.h>

#include "../klib/ksort.h"
#include "chain_dp.h"
#include "../gapped_align/oc_aligner.h"

#define DDFS_CUTOFF 0.25

KSORT_INIT(ChainSeed, ChainSeed, ChainSeedLT)

WordFindData*
new_WordFindData(idx reference_size, 
				 int block_size,
				 int kmer_size,
				 int block_score_cutoff)
{
	WordFindData* wfd = (WordFindData*)malloc(sizeof(WordFindData));
	int nblk = reference_size / block_size + 5;
	wfd->_blk_list = (ScoringBlock*)malloc(sizeof(ScoringBlock) * (nblk + 1));
	wfd->blk_list = wfd->_blk_list + 1;
	wfd->blk_idx_list = (ScoringBlockIndex*)malloc(sizeof(ScoringBlockIndex) * nblk);
	for (int i = 0; i < nblk; ++i) {
		int blk_id = i;
		wfd->blk_list[blk_id].index = -1;
		wfd->blk_list[blk_id].score = 0;
		wfd->blk_list[blk_id].last_kmer_id = -1;
	}
	wfd->_blk_list[0].score = 0;
	wfd->nblk = 0;
	kv_init(wfd->hash_list);
	kv_init(wfd->chain_seed_list);
	wfd->chain_data = new_ChainDpData(kmer_size, block_score_cutoff);
	return wfd;
}

void
clear_WordFindData(WordFindData* wfd)
{
	for (int i = 0; i < wfd->nblk; ++i) {
		int blk_id = wfd->blk_idx_list[i].block_idx;
		wfd->blk_list[blk_id].index = -1;
		wfd->blk_list[blk_id].score = 0;
		wfd->blk_list[blk_id].last_kmer_id = -1;
	}
	kv_clear(wfd->hash_list);
	kv_clear(wfd->chain_seed_list);
	wfd->nblk = 0;
}

WordFindData*
free_WordFindData(WordFindData* wfd)
{
	free(wfd->_blk_list);
	free(wfd->blk_idx_list);
	kv_destroy(wfd->hash_list);
	kv_destroy(wfd->chain_seed_list);
	free_ChainDpData(wfd->chain_data);
	free(wfd);
	return 0;
}

static int
extract_hash_values(const char* read,
					const int read_size,
					const int kmer_size,
					const int scan_window,
					vec_u64* hash_list)
{
    kv_clear(*hash_list);
    for (int i = 0; i <= read_size - kmer_size; i += scan_window) {
        u64 hash = 0;
        for (int j = 0; j < kmer_size; ++j) {
            const u8 c = (u8)read[i + j];
            hash = (hash << 2) | c;
        }
        kv_push(u64, *hash_list, hash);
    }
    return (int)kv_size(*hash_list);
}

static void
fill_one_seed(const int kmer_id,
			  const int blk_id,
			  const short blk_offset,
			  WordFindData* wfd)
{
	ScoringBlock* sblk = wfd->blk_list + blk_id;
	if (sblk->last_kmer_id >= kmer_id + 1) return;
	if (sblk->score >= BLK_SEEDS) return;
	int seed_id = sblk->score;
	++sblk->score;
	sblk->blk_offset[seed_id] = blk_offset;
	sblk->kmer_id[seed_id] = kmer_id + 1;
	sblk->last_kmer_id = kmer_id + 1;
	if (sblk->index == -1) {
		sblk->index = wfd->nblk;
		++wfd->nblk;
		wfd->blk_idx_list[sblk->index].block_idx = blk_id;
	}
	wfd->blk_idx_list[sblk->index].score = sblk->score + (sblk - 1)->score;
}

static void
collect_seeds(vec_u64* hash_list,
			  const char* read,
			  const int read_size,
			  const int read_id,
			  const int read_start_id,
			  const int reference_start_id,
			  PackedDB* reference,
			  LookupTable* lktbl,
			  const int block_size,
			  const int kmer_size,
			  const int scan_window,
			  const BOOL pairwise,
			  WordFindData* wfd)
{
	idx soff_max = IDX_MAX;
    if (pairwise) {
        int max_rid = reference_start_id + (int)PDB_NUM_SEQS(reference);
        if (read_id + read_start_id >= reference_start_id && read_id + read_start_id < max_rid) {
            soff_max = PDB_SEQ_OFFSET(reference, read_id);
        }
    }
	int n_kmer = extract_hash_values(read, read_size, kmer_size, scan_window, hash_list);
	
	for (int i = 0; i < n_kmer; ++i) {
		u64 n_match;
		u64* match_list = extract_kmer_list(lktbl, kv_A(*hash_list, i), &n_match);
		for (u64 k = 0; k < n_match; ++k) {
			if (match_list[k] >= soff_max) continue;
			fill_one_seed(i, match_list[k] / block_size, match_list[k] % block_size, wfd);
		}
	}
}

static int 
scoring_seeds(int *t_loc, int *t_seedn, int *t_score, int *loc, int k, int *rep_loc, float scan_window, int read_size)
{
	int i,j,maxval=0,maxi = 0,rep=0,lasti,tempi;
	for(i=0;i<k;i++)t_score[i]=0;
	for(i=0;i<k-1;i++)for(j=i+1,tempi=t_seedn[i];j<k;j++)if(tempi!=t_seedn[j]&&t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_size&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*scan_window)-1.0)<DDFS_CUTOFF){t_score[i]++;t_score[j]++;tempi=t_seedn[j];}
	
	for(i=0;i<k;i++){
		if(maxval<t_score[i]){maxval=t_score[i];maxi=i;rep=0;}
		else if(maxval==t_score[i]){rep++;lasti=i;}
	}
	for(i=0;i<4;i++)loc[i]=0;
	if(maxval>=5&&rep==maxval){loc[0]=t_loc[maxi],loc[1]=t_seedn[maxi];*rep_loc=maxi;loc[2]=t_loc[lasti],loc[3]=t_seedn[lasti];return(1);}
	else if(maxval>=5&&rep!=maxval){
		for(j=0;j<maxi;j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_size&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*scan_window)-1.0)<DDFS_CUTOFF){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
			}
		j=maxi;
		if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
		else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
		for(j=maxi+1;j<k;j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_size&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*scan_window)-1.0)<DDFS_CUTOFF){
				if(loc[0]==0){loc[0]=t_loc[j];loc[1]=t_seedn[j];*rep_loc=j;}
				else { loc[2]=t_loc[j];loc[3]=t_seedn[j];}
			}
		return(1);
	}
	else return(0);
}

static void
clear_block_scores(GappedCandidate* can, const size_t subject_start, const size_t block_size, WordFindData* wfd)
{
	idx sbeg = can->sbeg + subject_start;
	idx send = can->send + subject_start;
	idx sblk = sbeg / block_size;
	idx eblk = send / block_size;
	
	for (size_t i = sblk; i <= eblk; ++i) {
		wfd->blk_list[i].score = 0;
	}
}

static int
find_candidate_for_one_block(WordFindData* wfd,
							 const int block_id,
							 PackedDB* reference,
							 const double ddfs_cutoff,
							 const int block_score_cutoff,
							 const size_t block_size,
							 const int align_size_cutoff,
							 vec_chain_seed* chain_seed_list,
							 ChainDpData* chain_data,
							 const int scan_window,
							 const int qid,
							 const int qdir,
							 const idx qsize,
							 vec_can* candidates)
{
	const int list_size = BLK_SEEDS * 2;
	int kmer_id_list[list_size];
	int blk_offset_list[list_size];
	int score_list[list_size];
	int n_seeds = 0;
	int A = 0;
	idx blk_start = (idx)block_size * block_id;;

	if (wfd->blk_list[block_id - 1].score) {
		ScoringBlock* sblk = wfd->blk_list + block_id - 1;
		for (int i = 0; i < sblk->score; ++i) {
			kmer_id_list[n_seeds] = sblk->kmer_id[i];
			blk_offset_list[n_seeds] = sblk->blk_offset[i];
			++n_seeds;
		}
		A = block_size;
		blk_start = (idx)block_size * (block_id - 1);
	}
	ScoringBlock* sblk = wfd->blk_list + block_id;
	for (int i = 0; i < sblk->score; ++i) {
		kmer_id_list[n_seeds] = sblk->kmer_id[i];
		blk_offset_list[n_seeds] = sblk->blk_offset[i] + A;
		++n_seeds;
	}

	int max_score_id = -1;
	int seed_candidates[4];
	int r = scoring_seeds(blk_offset_list, kmer_id_list, score_list, seed_candidates, n_seeds, &max_score_id, scan_window, qsize);
	if (!r) return 0;
	if (score_list[max_score_id] < block_score_cutoff) return 0;
	
	idx seed_toff = seed_candidates[0] + blk_start;
	idx seed_qoff = (seed_candidates[1] - 1) * scan_window;
	int seed_bid = seed_toff / block_size;
	idx seed_tid = pdb_offset_to_id(reference, seed_toff);
	idx seed_tsize = PDB_SEQ_SIZE(reference, seed_tid);
	idx seed_tstart = PDB_SEQ_OFFSET(reference, seed_tid);
	idx seed_tend = seed_tstart + seed_tsize;
	seed_toff -= seed_tstart;
	idx L = OC_MIN(seed_toff, seed_qoff);
	int bid_start = seed_bid - L / block_size - 1;
	if (bid_start < 0) bid_start = 0;
	idx tr = seed_tsize - seed_toff;
	idx qr = qsize - seed_qoff;
	L = OC_MIN(tr, qr);
	int bid_end = seed_bid + (L + block_size - 1) / block_size;
	kv_clear(*chain_seed_list);
	ChainSeed seed;
	
	int seed_score = 0;
	for (int i = bid_start; i <= seed_bid; ++i) {
		ScoringBlock* sblk = wfd->blk_list + i;
		if (!sblk->score) continue;
		blk_start = (idx)i * block_size;
		int relevant = 0;
		for (int k = 0; k < sblk->score; ++k) {
			idx toff = blk_start + sblk->blk_offset[k];
			idx qoff = (sblk->kmer_id[k] - 1) * scan_window;
			if (toff < seed_tstart) {
				continue;
			}
			toff -= seed_tstart;
			if (toff < seed_toff && qoff < seed_qoff) {
				double s = 1.0 * (seed_toff - toff) / (seed_qoff - qoff) - 1.0;
				r = fabs(s) < DDFS_CUTOFF;
				if (!r) {
					continue;
				}
				++relevant;
				seed.qoff = qoff;
				seed.soff = toff;
				kv_push(ChainSeed, *chain_seed_list, seed);
			}
		}
		if (i != seed_bid && 1.0 * relevant / sblk->score >= 0.4) sblk->score = 0;
		seed_score += relevant;
	}

	seed.qoff = seed_qoff;
	seed.soff = seed_toff;
	kv_push(ChainSeed, *chain_seed_list, seed);

	for (int i = seed_bid; i <= bid_end; ++i) {
		ScoringBlock* sblk = wfd->blk_list + i;
		if (!sblk->score) continue;
		blk_start = (idx)i * block_size;
		int relevant = 0;
		for (int k = 0; k < sblk->score; ++k) {
			idx toff = blk_start + sblk->blk_offset[k];
			idx qoff = (sblk->kmer_id[k] - 1) * scan_window;
			if (toff >= seed_tend) {
				continue;
			}
			toff -= seed_tstart;
			if (toff > seed_toff && qoff > seed_qoff) {
				double s = 1.0 * (toff - seed_toff) / (qoff - seed_qoff) - 1.0;
				r = fabs(s) < DDFS_CUTOFF;
				if (!r) { 
					continue;
				}
				++relevant;
				seed.qoff = qoff;
				seed.soff = toff;
				kv_push(ChainSeed, *chain_seed_list, seed);
			}
		}
		if (i != seed_bid && 1.0 * relevant / sblk->score >= 0.4) sblk->score = 0;
		seed_score += relevant;
	}

	size_t n_chain_seeds = kv_size(*chain_seed_list);
//printf("============\n");
//for (size_t i = 0; i < n_chain_seeds; ++i) printf("%u\t%d\n", kv_A(*chain_seed_list,i).qoff, kv_A(*chain_seed_list, i).soff);
	ks_introsort(ChainSeed, n_chain_seeds, kv_data(*chain_seed_list));
	chain_data->chain_seeds = chain_seed_list;
	chain_dp(chain_data, qid, qdir, qsize, seed_tid, seed_tsize);
	size_t ncan = kv_size(chain_data->lcanv);
	if (!ncan) return 0;

	GappedCandidate can = kv_A(chain_data->lcanv, 0);
	if (seed_qoff >= can.qbeg && seed_qoff < can.qend && seed_toff >= can.sbeg && seed_toff < can.send) {
		can.score = seed_score;
		can.qoff = seed_qoff;
		can.soff = seed_toff;
		clear_block_scores(&can, seed_tstart, block_size, wfd);
		BOOL r = (can.send - can.sbeg >= align_size_cutoff) || (can.qend - can.qbeg >= align_size_cutoff);
		if (r) kv_push(GappedCandidate, *candidates, can);
		return r;
	}

	size_t max_i = ncan;
	int max_cov = 0;
	for (size_t i = 0; i < ncan; ++i) {
		can = kv_A(chain_data->lcanv, i);
		if (seed_qoff >= can.qbeg && seed_qoff < can.qend && seed_toff >= can.sbeg && seed_toff < can.send) {
			int cov = can.qend - can.qbeg;
			if (cov > max_cov) { max_cov = cov; max_i = i; }
		}
	}
	if (max_i < ncan) {
		can.score = seed_score;
		can.qoff = seed_qoff;
		can.soff = seed_toff;
		clear_block_scores(&can, seed_tstart, block_size, wfd);
		BOOL r = (can.send - can.sbeg >= align_size_cutoff) || (can.qend - can.qbeg >= align_size_cutoff);
		if (r) kv_push(GappedCandidate, *candidates, can);
		return r;
	}

	can = kv_A(chain_data->lcanv, 0);
	if (can.qend - can.qbeg >= 5000) {
		can.score = seed_score;
		can.qoff = seed_qoff;
		can.soff = seed_toff;
		clear_block_scores(&can, seed_tstart, block_size, wfd);
		BOOL r = (can.send - can.sbeg >= align_size_cutoff) || (can.qend - can.qbeg >= align_size_cutoff);
		if (r) kv_push(GappedCandidate, *candidates, can);
		return r;
	}
	return 0;
}

KSORT_INIT(sbi_gt, ScoringBlockIndex, sbi_gt)

void
find_candidates(const char* read,
				const int read_size,
				const int qid,
				const int qdir,
				const int read_start_id,
				const int reference_start_id,
				BOOL pairwise,
				PackedDB* reference,
				LookupTable* lktbl,
				MapOptions* options,
				WordFindData* wfd,
				vec_can* candidates)
{
	clear_WordFindData(wfd);
	collect_seeds(&wfd->hash_list,
				  read,
				  read_size,
				  qid,
				  read_start_id,
				  reference_start_id,
				  reference,
				  lktbl,
				  options->block_size,
				  options->kmer_size,
				  options->scan_window,
				  pairwise,
				  wfd);
	
	for (int i = 0; i < wfd->nblk; ++i) {
		int blk_id = wfd->blk_idx_list[i].block_idx;
		if (wfd->blk_list[blk_id].score >= options->block_score_cutoff)
		if (wfd->blk_idx_list[i].score >= 2 * options->block_score_cutoff) {
			find_candidate_for_one_block(wfd,
										 wfd->blk_idx_list[i].block_idx,
										 reference,
										 options->ddfs_cutoff,
										 2 * options->block_score_cutoff,
										 options->block_size,
										 options->align_size_cutoff,
										 &wfd->chain_seed_list,
										 wfd->chain_data,
										 options->scan_window,
										 qid,
										 qdir,
										 read_size,
										 candidates);
		}
	}
}
