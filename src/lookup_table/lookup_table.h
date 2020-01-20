#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H

#include "../common/packed_db.h"

typedef struct 
{
	u64*	offset_list;
	u64*	kmer_stats;
} LookupTable;

#define OffsetBits	((u64)34)
#define HashBits	((u64)30)
#define CntBits		((u64)30)
#define OffsetMask	((U64_ONE << OffsetBits)-1)

#define KmerStats_Offset(u) ((u)&OffsetMask)
#define KmerStats_Cnt(u)	((u)>>OffsetBits)
#define PACK_OFFSET(h, o)	(((h) << OffsetBits) | (o))
#define OFFSET_HASH(u)		((u)>>OffsetBits)
#define OFFSET_OFFSET(u)	((u) & OffsetMask)

LookupTable*
build_lookup_table(PackedDB* pdb, 
				   const int kmer_size, 
				   const int kmer_cnt_cutoff, 
				   const int num_threads);

LookupTable*
destroy_lookup_table(LookupTable* lktbl);

u64*
extract_kmer_list(LookupTable* lktbl, const u64 hash, u64* n);

#endif // LOOKUP_TABLE_H
