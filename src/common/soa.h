#ifndef SOA_H
#define SOA_H

#include "../klib/kvec.h"
#include "../common/ontcns_defs.h"

typedef struct {
	char* data;
	u32	nextAvailableData;
	u32 blockSize;
} OcChunk;
typedef kvec_t(OcChunk*) vec_occhunk;

typedef struct {
	vec_occhunk 	chunkList;
	u32				firstAvailableChunkIdx;
	u32				blockSize;
} OcObjectAllocator;

OcObjectAllocator*
new_OcObjectAllocator(u32 blockSize);

OcObjectAllocator*
free_OcObjectAllocator(OcObjectAllocator* alloc);

void
clear_OcObjectAllocator(OcObjectAllocator* alloc);

void*
alloc_OcObjectAllocator(OcObjectAllocator* alloc, const u32 numBlocks);

#endif // SOA_H
