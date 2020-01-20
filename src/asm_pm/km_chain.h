#ifndef KM_CHAIN_H
#define KM_CHAIN_H

#include "../common/gapped_candidate.h"
#include "align_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int match_size;
    int query_offset;
    int reference_offset;
} MaximalExactMatch;

typedef kvec_t(MaximalExactMatch) vec_mem;

void ks_introsort_mem_soff_lt(size_t n, MaximalExactMatch* a);

#define DUMP_MEM(output_func, out, mem) output_func(out, "%d\t%d\t%d\n", (mem).query_offset, (mem).reference_offset, (mem).match_size)
#define DUMP_MEM_STD(mem) DUMP_MEM(fprintf, stdout, mem)

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
} ChainWorkData;

ChainWorkData*
chain_data_new();

ChainWorkData*
chain_data_free(ChainWorkData* data);

void
chain_data_reset(ChainWorkData* data, int n);

void mem_chain_dp(ChainWorkData* data,
        MaximalExactMatch* mems,
        const int n,
        const int read_id,
        const int read_dir,
        const int read_size,
        const int ref_id,
        const idx ref_size,
        vec_can* cans);

int mem_find_best_can(ChainWorkData* data,
        MaximalExactMatch* mems,
        const int n,
        const int read_id,
        const int read_dir,
        const int read_size,
        const int ref_id,
        const idx ref_size,
        GappedCandidate* best_can,
        vec_mem* chain_mems);

int find_max_score_mem(ChainWorkData* data,
        MaximalExactMatch* v_mem,
        const int n,
        const int min_score);

#ifdef __cplusplus
}
#endif
#endif // KM_CHAIN_H