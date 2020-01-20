#ifndef HASH_LIST_BUCKET_SORT_H
#define HASH_LIST_BUCKET_SORT_H

#include "../common/ontcns_defs.h"

typedef u64 (*ValueExtractor)(void* list, const u64 i);
typedef void (*SetListValue)(void* src, const u64 src_idx, void* dst, const u64 dst_idx);

void
radix_sort(void* src, 
		   const u64 item_size,
		   const u64 item_count, 
		   const int num_threads,
		   ValueExtractor offset_extractor, 
		   ValueExtractor hash_extractor,
		   SetListValue set_list_value);

#endif // HASH_LIST_BUCKET_SORT_H
