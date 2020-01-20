#ifndef WORD_FINDER_AUX_H
#define WORD_FINDER_AUX_H

#include <limits.h>

#include "../common/gapped_candidate.h"
#include "../common/ontcns_defs.h"

#define BLK_SEEDS	40

typedef struct {
	idx qoff;
	idx soff;
} ChainSeed;

typedef kvec_t(ChainSeed) vec_chain_seed;
#define ChainSeedLT(a, b) ((a).soff < (b).soff || ((a).soff == (b).soff && (a).qoff < (b).qoff))

typedef struct {
	short 	score;
	short 	blk_offset[BLK_SEEDS];
	int 	kmer_id[BLK_SEEDS];
	int 	last_kmer_id;
	int 	index;
} ScoringBlock;

typedef struct {
	int score;
	int	block_idx;
} ScoringBlockIndex;

#define sbi_gt(a, b) (((a).score > (b).score) || ((a).score == (b).score && (a).block_idx < (b).block_idx))

#endif // WORD_FINDER_AUX_H
