#ifndef WORD_FINDER_H
#define WORD_FINDER_H

#include "../common/map_options.h"
#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"
#include "../lookup_table/lookup_table.h"
#include "chain_dp.h"
#include "word_finder_aux.h"

typedef struct {
	ScoringBlock*		_blk_list;
	ScoringBlock*		blk_list;
	ScoringBlockIndex*	blk_idx_list;
	int					nblk;
	ChainDpData*		chain_data;
	vec_u64				hash_list;
	vec_chain_seed		chain_seed_list;
} WordFindData;

WordFindData*
new_WordFindData(idx reference_size, 
				 int block_size,
				 int kmer_size,
				 int block_score_cutoff);

void
clear_WordFindData(WordFindData* wfd);

WordFindData*
free_WordFindData(WordFindData* wfd);

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
				WordFindData* wfdata,
				vec_can* candidates);

#endif // WORD_FINDER_H
