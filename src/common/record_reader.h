#ifndef RECORD_READER_H
#define RECORD_READER_H

#include "../common/ontcns_aux.h"

#define OcReaderBufferDefaultSize (64*(1<<20))

typedef struct {
	FILE* in;
	char* buffer;
	char* curr;
	size_t record_size;
	size_t l, m;
	size_t max_num_records;
} BinRecordReader;

BinRecordReader*
new_BinRecordReader(const char* path, size_t record_size);

BinRecordReader*
free_BinRecordReader(BinRecordReader* reader);

BOOL
get_BinRecordReader(BinRecordReader* reader, void* rptr);

#endif // RECORD_READER_H
