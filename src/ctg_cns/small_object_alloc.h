#ifndef SMALL_OBJECT_ALLOC_H
#define SMALL_OBJECT_ALLOC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../common/ontcns_aux.h"
#include "../klib/kvec.h"

typedef struct {
    char* data;
    u32 avail_data;
    u32 block_size;
} CnsMemoryChunk;

typedef kvec_t(CnsMemoryChunk*) vec_cns_mem_chunk;

typedef struct {
    vec_cns_mem_chunk   chunk_list;
    u32                 first_avail_chunk;
    u32                 block_size;
} CnsMemoryAlloc;

CnsMemoryAlloc*
cns_mem_alloc_new(u32 block_size);

CnsMemoryAlloc*
cns_mem_alloc_free(CnsMemoryAlloc* alloc);

CnsMemoryAlloc*
cns_mem_alloc_clear(CnsMemoryAlloc* alloc);

void*
cns_mem_alloc_alloc(CnsMemoryAlloc* alloc, const u32 num_blocks);

#ifdef __cplusplus
}
#endif
#endif // SMALL_OBJECT_ALLOC_H