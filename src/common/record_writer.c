#include "record_writer.h"

RecordWriter*
new_RecordWriter(FILE* out, OcMutex* out_lock, size_t buffer_size)
{
	RecordWriter* w = (RecordWriter*)malloc( sizeof(RecordWriter) );
	kstr_init(w->buffer);
	ks_reserve(&w->buffer, buffer_size * 1.5);
	w->buffer_size = 0;
	w->out = out;
	w->out_lock = out_lock;
	return w;
}

RecordWriter*
free_RecordWriter(RecordWriter* writer)
{
	dump_RecordWriter(writer);
	free_kstring(writer->buffer);
	free(writer);
	return 0;
}

void
dump_RecordWriter(RecordWriter* writer)
{
	size_t n = kstr_size(writer->buffer);
	if (!n) return;
	if (writer->out_lock) {
		pthread_mutex_lock(writer->out_lock);
	}
	FWRITE(kstr_str(writer->buffer), 1, n, writer->out);
	if (writer->out_lock) {
		pthread_mutex_unlock(writer->out_lock);
	}
	kstr_clear(writer->buffer);
}

void
check_RecordWriter(RecordWriter* writer)
{
	size_t n = kstr_size(writer->buffer);
	if (n >= writer->buffer_size) {
		dump_RecordWriter(writer);
	}
}
