#ifndef RECORD_WRITER_H
#define RECORD_WRITER_H

#include <pthread.h>

#include "ontcns_aux.h"
#include "../klib/kstring.h"

#define OcWriterBufferDefaultSize (8*(1<<20))

typedef struct {
    FILE*       out;
    OcMutex*    out_lock;
    kstring_t   buffer;
    size_t      buffer_size;
} RecordWriter;

RecordWriter*
new_RecordWriter(FILE* out, OcMutex* out_lock, size_t buffer_size);

RecordWriter*
free_RecordWriter(RecordWriter* writer);

void
dump_RecordWriter(RecordWriter* writer);

void
check_RecordWriter(RecordWriter* writer);

#define RW_DUMP_ONE_DATA(structName, outputFunc, bin, out, rptr) do { \
	if (bin) { \
		kputsn_((const void*)rptr, sizeof(structName), &(out)->buffer); \
	} else { \
		outputFunc(ksprintf, &(out)->buffer, *rptr); \
	} \
	check_RecordWriter(out); \
} while(0)

#endif // RECORD_WRITER_H
