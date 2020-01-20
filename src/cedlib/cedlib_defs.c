#include "cedlib_defs.h"

EdlibData*
new_EdlibData(const int maxQueryLength, const int maxTargetLength)
{
    EdlibData* data = (EdlibData*)malloc( sizeof(EdlibData) );
    kv_init(data->Ps);
    kv_init(data->Ms);
    kv_init(data->scores);
    kv_init(data->firstBlocks);
    kv_init(data->lastBlocks);
    kv_init(data->blocks);
    kv_init(data->peq);
    kv_init(data->startLocations);
	kv_init(data->endLocations);
	kv_init(data->alignment);

    int maxNumBlocks = CALC_NUM_BLOCKS(maxQueryLength);
    int space = maxNumBlocks * maxTargetLength;
    kv_reserve(Word, data->Ps, space);
    kv_reserve(Word, data->Ms, space);
    kv_reserve(int, data->scores, space);
    kv_reserve(int, data->firstBlocks, maxTargetLength);
    kv_reserve(int, data->lastBlocks, maxTargetLength);
    kv_reserve(EdlibBlock, data->blocks, maxNumBlocks);
    space = (AlphabetSize + 1) * maxNumBlocks;
    kv_reserve(Word, data->peq, space);
	
	data->status = EDLIB_STATUS_OK;
	data->editDistance = -1;

    return data;
}

EdlibData*
free_EdlibData(EdlibData* data)
{
    free_kvec(data->Ps);
    free_kvec(data->Ms);
    free_kvec(data->scores);
    free_kvec(data->firstBlocks);
    free_kvec(data->lastBlocks);
    free_kvec(data->blocks);
    free_kvec(data->peq);
	free_kvec(data->startLocations);
	free_kvec(data->endLocations);
	free_kvec(data->alignment);
    free(data);
    return 0;
}

void
reset_EdlibData(EdlibData* data, const int queryLength, const int targetLength)
{
    int maxNumBlocks = CALC_NUM_BLOCKS(queryLength);
    int space = maxNumBlocks * targetLength;
    kv_resize(Word, data->Ps, space);
    kv_resize(Word, data->Ms, space);
    kv_resize(int, data->scores, space);
    kv_resize(int, data->firstBlocks, targetLength);
    kv_resize(int, data->lastBlocks, targetLength);
    kv_resize(EdlibBlock, data->blocks, maxNumBlocks);
    space = (AlphabetSize + 1) * maxNumBlocks;
    kv_resize(Word, data->peq, space);
	kv_clear(data->startLocations);
	kv_clear(data->endLocations);
	kv_clear(data->alignment);
	data->status = EDLIB_STATUS_OK;
	data->editDistance = -1;
}
