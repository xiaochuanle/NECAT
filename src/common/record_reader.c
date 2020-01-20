#include "record_reader.h"

BinRecordReader*
new_BinRecordReader(const char* path, size_t record_size)
{
	BinRecordReader* brr = (BinRecordReader*)malloc( sizeof(BinRecordReader) );
	FOPEN(brr->in, path, "rb");
	brr->record_size = record_size;
	size_t n = OcReaderBufferDefaultSize / record_size;
	size_t bytes = n * record_size;
	brr->buffer = (char*)malloc(bytes);
	brr->curr = brr->buffer;
	brr->l = 0;
	brr->m = 0;
	brr->max_num_records = n;
	
	return brr;
}

BinRecordReader*
free_BinRecordReader(BinRecordReader* reader)
{
	FCLOSE(reader->in);
	free(reader->buffer);
	free(reader);
	return 0;
}

BOOL
get_BinRecordReader(BinRecordReader* reader, void* rptr)
{
	if (reader->l >= reader->m) {
		size_t n = fread(reader->buffer, reader->record_size, reader->max_num_records, reader->in);
		if (n == 0) return FALSE;
		reader->l = 0;
		reader->m = n;
		reader->curr = reader->buffer;
	}
	memcpy(rptr, reader->curr, reader->record_size);
	reader->curr += reader->record_size;
	++reader->l;
	return TRUE;
}
