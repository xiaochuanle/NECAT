#ifndef CHAIN_H
#define CHAIN_H

#include "../common/gapped_candidate.h"
#include "../klib/kvec.h"
#include "word_finder_aux.h"

typedef struct
{
	vec_int			f;
	vec_int 		p;
	vec_int 		t;
	vec_int 		v;
	vec_intpair 	u;
	vec_can			lcanv;
	vec_chain_seed*	chain_seeds;
	int				kmer_size;
	int 			max_dist;
	int 			bw;
	int 			max_skip;
	int 			min_cnt;
	int				min_sc;
} ChainDpData;

ChainDpData*
new_ChainDpData(int kmer_size, int block_score_cutoff);

ChainDpData*
free_ChainDpData(ChainDpData* chain_data);

int
chain_dp_forward(ChainDpData* data, ChainSeed* rm_seed);

int
chain_dp_backward(ChainDpData* data, ChainSeed* lm_seed);

#endif // CHAIN_H
