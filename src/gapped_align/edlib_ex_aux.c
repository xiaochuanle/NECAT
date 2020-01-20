#include "edlib_ex_aux.h"

EdlibAlignMatrix*
new_EdlibAlignMatrix(int maxNumBlocks, int targetLength)
{
    EdlibAlignMatrix* data = (EdlibAlignMatrix*)malloc(sizeof(EdlibAlignMatrix));
    data->Ps = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->Ms = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->scores = (int*)malloc(sizeof(int) * maxNumBlocks * targetLength);
    data->first_blocks = (int*)malloc(sizeof(int) * targetLength);
    data->last_blocks = (int*)malloc(sizeof(int) * targetLength);

    return data;
}

EdlibAlignMatrix*
free_EdlibAlignMatrix(EdlibAlignMatrix* data)
{
    free(data->Ps);
    free(data->Ms);
    free(data->scores);
    free(data->first_blocks);
    free(data->last_blocks);
    free(data);
    return 0;
}

SmallEdlibAlignResult*
new_SmallEdlibAlignResult()
{
    SmallEdlibAlignResult* result = (SmallEdlibAlignResult*)malloc(sizeof(SmallEdlibAlignResult));
    kv_init(result->end_locations);
    kv_init(result->start_locations);
    kv_init(result->alignment);
    return result;
}

SmallEdlibAlignResult*
free_SmallEdlibAlignResult(SmallEdlibAlignResult* result)
{
    kv_destroy(result->end_locations);
    kv_destroy(result->start_locations);
    kv_destroy(result->alignment);
    free(result);
    return 0;
}

void
clear_SmallEdlibAlignResult(SmallEdlibAlignResult* result)
{
    result->edit_distance = -1;
    kv_clear(result->end_locations);
    kv_clear(result->start_locations);
    kv_clear(result->alignment);
}

EdlibAlignData*
new_EdlibAlignData(int maxNumBlocks, int targetLength)
{
    EdlibAlignData* data = (EdlibAlignData*)malloc(sizeof(EdlibAlignData));
    data->align_matrix = new_EdlibAlignMatrix(maxNumBlocks, targetLength);
    data->result = new_SmallEdlibAlignResult();
    data->blocks = (EdlibBlock*)malloc(sizeof(EdlibBlock) * maxNumBlocks);
	data->peq = (Word*)malloc(sizeof(Word) * (AlphabetSize + 1) * maxNumBlocks);
    kv_init(data->cigar);
	data->maxNumBlocks = maxNumBlocks;
	data->maxTargetLength = targetLength;
    return data;
}

EdlibAlignData*
free_EdlibAlignData(EdlibAlignData* data)
{
    data->align_matrix = free_EdlibAlignMatrix(data->align_matrix);
    data->result = free_SmallEdlibAlignResult(data->result);
    free(data->blocks);
	free(data->peq);
    kv_destroy(data->cigar);
    free(data);
    return 0;
}

void
alloc_space_EdlibAlignData(EdlibAlignData* data, int query_length, int target_length)
{
	int newMaxNumBlocks = (query_length + WORD_SIZE - 1) / WORD_SIZE;
	int newMaxTargetLength = target_length;
	int newMaxSpace = sizeof(Word) * newMaxNumBlocks * newMaxTargetLength;
	
	if (newMaxNumBlocks > data->maxNumBlocks || newMaxTargetLength > data->maxTargetLength) {
		data->align_matrix->Ps = (Word*)realloc(data->align_matrix->Ps, newMaxSpace);
		data->align_matrix->Ms = (Word*)realloc(data->align_matrix->Ms, newMaxSpace);
		data->align_matrix->scores = (int*)realloc(data->align_matrix->scores, sizeof(int) * newMaxNumBlocks * newMaxTargetLength);
		data->align_matrix->first_blocks = (int*)realloc(data->align_matrix->first_blocks, sizeof(int) * newMaxTargetLength);
		data->align_matrix->last_blocks = (int*)realloc(data->align_matrix->last_blocks, sizeof(int) * newMaxTargetLength);
		
		data->blocks = (EdlibBlock*)realloc(data->blocks, sizeof(EdlibBlock) * newMaxNumBlocks);
		data->peq = (Word*)realloc(data->peq, sizeof(Word) * (AlphabetSize + 1) * newMaxNumBlocks);
		
		data->maxNumBlocks = newMaxNumBlocks;
		data->maxTargetLength = newMaxTargetLength;
	}
}
