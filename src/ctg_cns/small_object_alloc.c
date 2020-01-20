#include "small_object_alloc.h"

#include "../common/oc_assert.h"

static const u32 kCnsMemChunkSize = 8 * 1024 * 1024;

CnsMemoryChunk*
cns_mem_chunk_new(const u32 block_size)
{
    CnsMemoryChunk* chunk = (CnsMemoryChunk*)malloc(sizeof(CnsMemoryChunk));
    chunk->data = (char*)malloc(kCnsMemChunkSize);
    chunk->avail_data = 0;
    chunk->block_size = block_size;
    return chunk;
}

CnsMemoryChunk* 
cns_mem_chunk_free(CnsMemoryChunk* chunk)
{
    free(chunk->data);
    free(chunk);
    return NULL;
}

void
cns_mem_chunk_clear(CnsMemoryChunk* chunk)
{
    chunk->avail_data = 0;
}

void*
cns_mem_chunk_alloc(CnsMemoryChunk* chunk, u32 num_blocks)
{
    u32 alloc_bytes = num_blocks * chunk->block_size;
    if (chunk->avail_data + alloc_bytes > kCnsMemChunkSize) return NULL;
    char* p = chunk->data + chunk->avail_data;
    chunk->avail_data += alloc_bytes;
    return (void*)p;
}

static u32
round_up_16(const u32 s)
{
    const u32 mask = 15;
    u32 r = 16 - (s & mask);
    return s + r;
}

CnsMemoryAlloc*
cns_mem_alloc_new(u32 block_size)
{
    CnsMemoryAlloc* alloc = (CnsMemoryAlloc*)malloc(sizeof(CnsMemoryAlloc));
    kv_init(alloc->chunk_list);
    alloc->first_avail_chunk = 0;
    alloc->block_size = round_up_16(block_size);
    CnsMemoryChunk* chunk = cns_mem_chunk_new(alloc->block_size);
    kv_push(CnsMemoryChunk*, alloc->chunk_list, chunk);

    return alloc;
}

CnsMemoryAlloc*
cns_mem_alloc_free(CnsMemoryAlloc* alloc)
{
    size_t n_chunk = kv_size(alloc->chunk_list);
    for (size_t i = 0; i < n_chunk; ++i) {
        cns_mem_chunk_free(kv_A(alloc->chunk_list, i));
    }
    kv_destroy(alloc->chunk_list);
    free(alloc);

    return NULL;
}

CnsMemoryAlloc*
cns_mem_alloc_clear(CnsMemoryAlloc* alloc)
{
    size_t n_chunk = kv_size(alloc->chunk_list);
    for (size_t i = 0; i < n_chunk; ++i) {
        cns_mem_chunk_clear(kv_A(alloc->chunk_list, i));
    }
    alloc->first_avail_chunk = 0;

    return NULL;
}

void*
cns_mem_alloc_alloc(CnsMemoryAlloc* alloc, const u32 num_blocks)
{
    void* p = cns_mem_chunk_alloc(kv_A(alloc->chunk_list, alloc->first_avail_chunk), num_blocks);
    if (!p) {
        if (alloc->first_avail_chunk == kv_size(alloc->chunk_list) - 1) {
            CnsMemoryChunk* chunk = cns_mem_chunk_new(alloc->block_size);
            kv_push(CnsMemoryChunk*, alloc->chunk_list, chunk);
        }
        alloc->first_avail_chunk++;
        p = cns_mem_chunk_alloc(kv_A(alloc->chunk_list, alloc->first_avail_chunk), num_blocks);
        oc_assert(p != NULL);
    }
    return p;
}