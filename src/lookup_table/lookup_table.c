#include "lookup_table.h"

#include <assert.h>

#include "../common/ontcns_aux.h"
#include "hash_list_bucket_sort.h"

#define get_next_hash { \
	u64 __k = offset + j; \
	const u8 c = PDB_GET_CHAR(pdb, __k); \
    hash = ((hash << 2) | c) & hash_mask; \
}

u64*
get_kmer_counts(PackedDB* pdb,
				const int kmer_size,
				const int kmer_cnt_cutoff,
				u64* num_kmers)
{
	TIMING_START(__func__);
	
	u64* cnt_table = (u64*)malloc(sizeof(u64) * (kmer_cnt_cutoff + 2));
	for (u64 i = 0; i <= kmer_cnt_cutoff; ++i) cnt_table[i] = i + 1;
	cnt_table[kmer_cnt_cutoff + 1] = kmer_cnt_cutoff + 1;
	u64 max_hash = U64_ONE << (kmer_size * 2);
	u64 hash_mask = max_hash - 1;
	u64* kmer_cnts = (u64*)calloc(max_hash, sizeof(u64));
	const u64 ns = PDB_NUM_SEQS(pdb);
	
	for (u64 i = 0; i != ns; ++i) {
		const u64 offset = PDB_SEQ_OFFSET(pdb, i);
		const u64 size = PDB_SEQ_SIZE(pdb, i);
		u64 hash = 0;
		for (u64 j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (u64 j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			kmer_cnts[hash] = cnt_table[ kmer_cnts[hash] ];
		}
	}
	
	u64 n = 0;
	for (u64 i = 0; i != max_hash; ++i) {
		if (kmer_cnts[i] > kmer_cnt_cutoff) kmer_cnts[i] = 0;
		n += kmer_cnts[i];
	}
	free(cnt_table);
	for (u64 i = 0; i != max_hash; ++i) {
		kmer_cnts[i] = kmer_cnts[i] << OffsetBits;
	}
	
	TIMING_END(__func__);
	
	*num_kmers = n;
	return kmer_cnts;
}

u64*
get_offset_list(PackedDB* pdb,
				u64* kmer_stats,
				const int kmer_size,
				const u64 num_kmers)
{
	TIMING_START(__func__);
	
	u64* offset_list = (u64*)malloc( sizeof(u64) * num_kmers );
	u64 cnt = 0;
	u64 max_hash = U64_ONE << (kmer_size * 2);
	u64 hash_mask = max_hash - 1;
	u64 ns = PDB_NUM_SEQS(pdb);
	
	for (u64 i = 0; i != ns; ++i) {
		const u64 offset = PDB_SEQ_OFFSET(pdb, i);
		const u64 size = PDB_SEQ_SIZE(pdb, i);
		u64 hash = 0;
		for (u64 j = 0; j < kmer_size - 1; ++j) {
			get_next_hash;
		}
		for (u64 j = kmer_size - 1; j < size; ++j) {
			get_next_hash;
			u64 k = offset + j + 1 - kmer_size;
			u64 n = KmerStats_Cnt(kmer_stats[hash]);
			if (n) offset_list[cnt++] = PACK_OFFSET(hash, k);
		}
	}
	
	assert(cnt == num_kmers);
	TIMING_END(__func__);
	return offset_list;
}

void
build_kmer_starts(u64* offset_list, const u64 num_kmers, const int kmer_size, u64* kmer_stats)
{
	TIMING_START(__func__);
	
	u64 i = 0;
	while (i < num_kmers) {
		u64 hash = OFFSET_HASH(offset_list[i]);
		u64 j = i + 1;
		while (j < num_kmers) {
			u64 h = OFFSET_HASH(offset_list[j]);
			if (h != hash) break;
			++j;
		}
		kmer_stats[hash] |= i;
		i = j;
	}
	
	TIMING_END(__func__);
}

void
clear_hash_in_offset_list(u64* offset_list, const u64 num_kmers, const u64 max_offset)
{
	TIMING_START(__func__);
	
	for (u64 i = 0; i != num_kmers; ++i) {
		offset_list[i] = offset_list[i] & OffsetMask;
		assert(offset_list[i] < max_offset);
	}
	
	TIMING_END(__func__);
}

u64 hash_extractor(void* list, const u64 i)
{
	u64* ul = list;
	u64 u = ul[i];
	return OFFSET_HASH(u);
}

u64 offset_extractor(void* list, const u64 i)
{
	u64* ul = list;
	u64 u = ul[i];
	return OFFSET_OFFSET(u);
}

void set_offset_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
	u64* su = src;
	u64* du = dst;
	du[dst_idx] = su[src_idx];
}

LookupTable*
build_lookup_table(PackedDB* pdb, 
				   const int kmer_size, 
				   const int kmer_cnt_cutoff, 
				   const int num_threads)
{
	u64 num_kmers;
	u64* kmer_stats = get_kmer_counts(pdb, kmer_size, kmer_cnt_cutoff, &num_kmers);
	u64* offset_list = get_offset_list(pdb, kmer_stats, kmer_size, num_kmers);
	radix_sort(offset_list, sizeof(u64), num_kmers, num_threads, offset_extractor, hash_extractor, set_offset_value);
	build_kmer_starts(offset_list, num_kmers, kmer_size, kmer_stats);
	clear_hash_in_offset_list(offset_list, num_kmers, PDB_SIZE(pdb));
	LookupTable* lktbl = (LookupTable*)malloc(sizeof(LookupTable));
	lktbl->offset_list = offset_list;
	lktbl->kmer_stats = kmer_stats;
	return lktbl;
}

LookupTable*
destroy_lookup_table(LookupTable* lktbl)
{
	if (lktbl->offset_list) free(lktbl->offset_list);
	if (lktbl->kmer_stats) free(lktbl->kmer_stats);
	free(lktbl);
	return 0;
}

u64*
extract_kmer_list(LookupTable* lktbl, const u64 hash, u64* n)
{
	u64* list = 0;
	*n = 0;
	
	const u64 u = lktbl->kmer_stats[hash];
	u64 cnt = KmerStats_Cnt(u);
	u64 start = KmerStats_Offset(u);
	if (cnt) {
		list = lktbl->offset_list + start;
		*n = cnt;
	}
	
	return list;
}
