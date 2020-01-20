#include "soa.h"

#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"

static const u32 kOcChunkSize = 8 * 1024 * 1024;

OcChunk*
new_OcChunk(const u32 blockSize)
{
	OcChunk* chunk = (OcChunk*)malloc(sizeof(OcChunk));
	chunk->data = (char*)malloc(kOcChunkSize);
	chunk->nextAvailableData = 0;
	chunk->blockSize = blockSize;
	return chunk;
}

OcChunk*
free_OcChunk(OcChunk* chunk)
{
	free(chunk->data);
	free(chunk);
	return 0;
}

void
clear_OcChunk(OcChunk* chunk)
{
	chunk->nextAvailableData = 0;
}

void*
alloc_OcChunk(OcChunk* chunk, u32 numBlocks)
{
	u32 allocBytes = numBlocks * chunk->blockSize;
	if (chunk->nextAvailableData + allocBytes > kOcChunkSize) return 0;
	char* p = chunk->data + chunk->nextAvailableData;
	chunk->nextAvailableData += allocBytes;
	return (void*)p;
}

static u32
RoundupTo16(const u32 s)
{
	const u32 mask = 15;
	u32 r = 16 - (s & mask);
	return s + r;
}

OcObjectAllocator*
new_OcObjectAllocator(u32 blockSize)
{
	OcObjectAllocator* alloc = (OcObjectAllocator*)malloc(sizeof(OcObjectAllocator));
	kv_init(alloc->chunkList);
	alloc->firstAvailableChunkIdx = 0;
	alloc->blockSize = RoundupTo16(blockSize);
	OcChunk* chunk = new_OcChunk(alloc->blockSize);
	kv_push(OcChunk*, alloc->chunkList, chunk);

    return alloc;
}

OcObjectAllocator*
free_OcObjectAllocator(OcObjectAllocator* alloc)
{
	size_t n_chunk = kv_size(alloc->chunkList);
	for (size_t i = 0; i != n_chunk; ++i) {
		free_OcChunk(kv_A(alloc->chunkList, i));
	}
	kv_destroy(alloc->chunkList);
	free(alloc);
	return 0;
}

void
clear_OcObjectAllocator(OcObjectAllocator* alloc)
{
	size_t n_chunk = kv_size(alloc->chunkList);
	for (size_t i = 0; i != n_chunk; ++i) {
		clear_OcChunk(kv_A(alloc->chunkList, i));
	}
	alloc->firstAvailableChunkIdx = 0;
}

void*
alloc_OcObjectAllocator(OcObjectAllocator* alloc, const u32 numBlocks)
{
	void* p = alloc_OcChunk(kv_A(alloc->chunkList, alloc->firstAvailableChunkIdx), numBlocks);
	if (!p) {
		if (alloc->firstAvailableChunkIdx == kv_size(alloc->chunkList) - 1) {
			OcChunk* chunk = new_OcChunk(alloc->blockSize);
			kv_push(OcChunk*, alloc->chunkList, chunk);
		}
		++alloc->firstAvailableChunkIdx;
		p = alloc_OcChunk(kv_A(alloc->chunkList, alloc->firstAvailableChunkIdx), numBlocks);
		oc_assert(p, "alloc failed");
	}
	return p;
}
