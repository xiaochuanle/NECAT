#ifndef MEM_FINDER_H
#define MEM_FINDER_H

#include "../common/gapped_candidate.h"
#include "../klib/kvec.h"
#include "chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int query_offset;
    int reference_offset;
    int match_size;
} MaximalExactMatch;

typedef kvec_t(MaximalExactMatch) vec_mem;

#define DUMP_MEM(output_func, out, mem) \
    output_func(out, "%d\t%d\t%d\n", (mem).query_offset, (mem).reference_offset, (mem).match_size)

#define DUMP_MEM_STD(mem) DUMP_MEM(fprintf, stderr, mem)

void ks_introsort_mem_soff_lt(size_t n, MaximalExactMatch* a);

typedef struct {
    u64 hash;
    int offset;
    int occ;
} KmerInfo;

typedef kvec_t(KmerInfo) vec_kmif;

typedef struct {
    vec_kmif qry_kmif_list;
    vec_kmif ref_kmif_list;
    vec_chain_seed chain_seed_list;
    vec_can can_list;
    vec_int_pair chain_seed_info_list;
    const u8* reference;
    int ref_size;
    int kmer_size;
    int window_size;
    int mem_size;
} MaximalExactMatchWorkData;

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataNew(int kmer_size, int window_size, int mem_size);

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataFree(MaximalExactMatchWorkData* data);

void 
MaximalExactMatchWorkData_Init(MaximalExactMatchWorkData* data, const u8* reference, const int ref_size);

int
MaximalExactMatchWorkData_FindCandidate(
    MaximalExactMatchWorkData* data, 
    ChainDpWorkData* chain_data,
    const u8* read, const int read_size,
    GappedCandidate* can);

int
MaximalExactMatchWorkData_FindCandidates(
    MaximalExactMatchWorkData* data, 
    ChainDpWorkData* chain_data,
    const u8* read, 
    const int read_size);

#ifdef __cplusplus
}
#endif
#endif // MEM_FINDER_H