#ifndef FIND_MEM_H
#define FIND_MEM_H

#include "../common/ontcns_aux.h"
#include "km_chain.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    u64 hash;
    int offset;
    int occ;
} KmerInfo;

typedef kvec_t(KmerInfo) vec_kmif;

void
sort_kmif_list(vec_kmif* kmif_list);

int
build_kmif_list(vec_u8* sequence, int kmer_size, int window_size, vec_kmif* kmif_list);

int
compute_align_range_1(ChainWorkData* chain_data,
    vec_u8* query,
    vec_u8* target,
    vec_kmif* target_kmif_list,
    const int kmer_size,
    const int window_size,
    const int min_mem_size,
    GappedCandidate* can,
    vec_mem* chain_mems);

#ifdef __cplusplus
}
#endif
#endif // FIND_MEM_H