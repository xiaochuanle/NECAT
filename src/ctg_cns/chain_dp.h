#ifndef CHAIN_DP_H
#define CHAIN_DP_H

#include "../common/ontcns_aux.h"
#include "../common/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int match_size;
    idx qoff;
    idx soff;
} ChainDpSeed;

typedef kvec_t(ChainDpSeed) vec_chain_seed;

void ks_introsort_chain_seed_soff_lt(size_t n, ChainDpSeed* a);

typedef vec_intpair vec_int_pair;

typedef struct {
    int max_dist_ref;
    int max_dist_qry;
    int max_band_width;
    int max_skip;
    int min_cnt;
    int min_score;
    int dump_info;

    vec_int     f;
    vec_int     p;
    vec_int     t;
    vec_int     v;
    vec_int_pair u;
    vec_chain_seed seeds;
} ChainDpWorkData;

ChainDpWorkData*
ChainDpWorkDataNew(int min_seed_cnt, int min_can_score);

ChainDpWorkData*
ChainDpWorkDataFree(ChainDpWorkData* data);

void
ChainDpWorkDataSetup(ChainDpWorkData* data, int n);

int chaining_seeds(ChainDpWorkData* data,
        int* best_seed_index,
        int* best_seed_score,
        vec_chain_seed* chain_seeds);

int chaining_find_candidates(ChainDpWorkData* data,
        vec_can* can_list,
        vec_chain_seed* chain_seed_list,
        vec_int_pair* chain_seed_info_list);

#ifdef __cplusplus
}
#endif
#endif // CHAIN_DP_H