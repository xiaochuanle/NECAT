#include "pm4_aux.h"

#include <assert.h>

#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "../common/record_writer.h"

static FILE* m4_in = NULL;
static pthread_mutex_t m4_read_lock;
static AsmM4Writers* m4_out = NULL;
static pthread_mutex_t m4_write_lock;
static int min_read_id = 0;
static int max_read_id = 0;
static int batch_size = 0;
static double ident_perc_cutoff = 0.0;

void
make_partition_name(const char* prefix, const int pid, kstring_t* name)
{
	kstr_clear(*name);
	kputs(prefix, name);
	ksprintf(name, ".p%d", pid);
	kputc('\0', name);
}

void
make_partition_index_name(const char* prefix, kstring_t* name)
{
	kstr_clear(*name);
	kputs(prefix, name);
	kputs(".partitions", name);
	kputc('\0', name);
}

int 
load_num_partitions(const char* m4_path)
{
	new_kstring(path);
	make_partition_index_name(m4_path, &path);
	DFOPEN(in, kstr_str(path), "r");
	int n;
	SAFE_SCANF(fscanf, in, 1, "%d", &n);
	FCLOSE(in);
	free_kstring(path);
	return n;
}

void
dump_num_partitions(const char* m4_path, const int np)
{
	new_kstring(path);
	make_partition_index_name(m4_path, &path);
	DFOPEN(out, kstr_str(path), "w");
	fprintf(out, "%d\n", np);
	FCLOSE(out);
	free_kstring(path);
}

AsmM4Writers*
new_AsmM4Writers(const int num_dumpped_files)
{
	AsmM4Writers* w = (AsmM4Writers*)malloc( sizeof(AsmM4Writers) );
	w->nw = 0;
	w->mnw = num_dumpped_files;
	w->file_list = (FILE**)malloc( sizeof(FILE*) * num_dumpped_files );
	for (int i = 0; i < w->mnw; ++i) w->file_list[i] = NULL;
	return w;
}

AsmM4Writers*
free_AsmM4Writers(AsmM4Writers* w)
{
	free(w->file_list);
	free(w);
	return 0;
}

void
open_AsmM4Writers(AsmM4Writers* w, const char* m4_name, int sfid, int efid)
{
	int nfid = efid - sfid;
	w->nw = nfid;
	new_kstring(name);
	for (int i = 0; i < nfid; ++i) {
		make_partition_name(m4_name, i + sfid, &name);
		FOPEN(w->file_list[i], kstr_str(name), "wb");
	}
	free_kstring(name);
}

void
close_AsmM4Writers(AsmM4Writers* w)
{
	for (int i = 0; i < w->nw; ++i) FCLOSE(w->file_list[i]);
}

static size_t
load_m4(vec_m4* m4v)
{
	const size_t buff_size = U64_ONE * 256 * 1024 * 1024;
	const size_t N = buff_size / sizeof(M4Record);
	const size_t M = 2 * N;
	kv_clear(*m4v);
	kv_reserve(M4Record, *m4v, M);
	size_t n = 0;
	pthread_mutex_lock(&m4_read_lock);
	n = fread(kv_data(*m4v), sizeof(M4Record), N, m4_in);
	pthread_mutex_unlock(&m4_read_lock);
	kv_resize(M4Record, *m4v, n);
	return n;
}

#define id_in_range(id) ((id) >= min_read_id && (id) < max_read_id)

#define fix_asm_m4_offsets(c, nc, query_is_target) \
	do { \
	nc = c; \
		if (query_is_target) { \
		nc.qid = c.sid; \
		nc.qdir = c.sdir; \
		nc.qoff = c.soff; \
		nc.qend = c.send; \
		nc.qext = c.sext; \
		nc.qsize = c.ssize; \
		nc.sid = c.qid; \
		nc.sdir = c.qdir; \
		nc.soff = c.qoff; \
		nc.send = c.qend; \
		nc.sext = c.qext; \
		nc.ssize = c.qsize; \
	} \
		if (nc.sdir == REV) { \
		nc.sdir = ReverseStrand(nc.sdir); \
		nc.qdir = ReverseStrand(nc.qdir); \
	} \
} while(0)

static void*
pm4_thread_func(void* arg)
{
	new_kvec(vec_m4, m4v);
	size_t n;
	size_t m;
	M4Record m4;
	new_kvec(vec_size_type, idx_range);
	while ((n = load_m4(&m4v))) {
		m = 0;
		M4Record* m4s = kv_data(m4v);
		for (size_t i = 0; i < n; ++i) {
			if (m4s[i].ident_perc < ident_perc_cutoff) continue;
			int qid = m4s[i].qid;
			int sid = m4s[i].sid;
			BOOL r = id_in_range(qid) || id_in_range(sid);
			if (r) m4s[m++] = m4s[i];
		}
		if (m == 0) continue;
		n = m;
		
		for (size_t i = 0; i < m; ++i) {
			BOOL sid_is_in = FALSE;
			int sid = m4s[i].sid;
			if (id_in_range(sid)) {
				sid_is_in = TRUE;
			}
			int qid = m4s[i].qid;
			if (id_in_range(qid)) {
				fix_asm_m4_offsets(m4s[i], m4, 1);
				if (sid_is_in) m4s[n++] = m4;
				else m4s[i] = m4;
			}
		}
		
		ks_introsort_M4Record_SidLT(n, m4s);
		kv_clear(idx_range);
		kv_push(size_t, idx_range, 0);
		size_t i = 0;
		while (i < n) {
			const int sid = m4s[i].sid;
			const int bid = sid / batch_size;
			const int sid_from = bid * batch_size;
			const int sid_to = sid_from + batch_size;
			size_t j = i + 1;
			while (j < n && m4s[j].sid < sid_to) ++j;
			kv_push(size_t, idx_range, i);
			i = j;
		}
		kv_push(size_t, idx_range, n);
		
		pthread_mutex_lock(&m4_write_lock);
		oc_assert(kv_size(idx_range) > 1);
		for (i = 0; i < kv_size(idx_range) - 1; ++i) {
			size_t from = kv_A(idx_range, i);
			size_t to = kv_A(idx_range, i + 1);
			m = to - from;
			int bid = (m4s[from].sid - min_read_id) / batch_size;
			oc_assert(bid >= 0 && bid < m4_out->nw);
			FWRITE(m4s + from, sizeof(M4Record), m, m4_out->file_list[bid]);
		}
		pthread_mutex_unlock(&m4_write_lock);
	}
	
	return NULL;
}

void
pm4_main(const char* wrk_dir, 
		 const char* m4_path, 
		 const int num_dumpped_files, 
		 const int partition_size, 
		 const double min_ident_perc,
		 const int num_threads)
{
	int num_reads = load_num_reads(wrk_dir);
	int num_batches = (num_reads + partition_size - 1) / partition_size;
	dump_num_partitions(m4_path, num_batches);
	m4_out = new_AsmM4Writers(num_dumpped_files);
	batch_size = partition_size;
	ident_perc_cutoff = min_ident_perc;
	char job[1024];
	pthread_t job_ids[num_threads];
	
	for (int fid = 0; fid < num_batches; fid += num_dumpped_files) {
		int sfid = fid;
		int efid = OC_MIN(sfid + num_dumpped_files, num_batches);
		int ssid = sfid * partition_size;
		int esid = efid * partition_size;
		min_read_id = ssid;
		max_read_id = esid;
		sprintf(job, "dumping records for partitions [%d, %d)", sfid, efid);
		TIMING_START(job);
		open_AsmM4Writers(m4_out, m4_path, sfid, efid);
		FOPEN(m4_in, m4_path, "rb");
		for (int i = 0; i < num_threads; ++i) {
			pthread_create(job_ids + i, NULL, pm4_thread_func, NULL);
		}
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(job_ids[i], NULL);
		}
		FCLOSE(m4_in);
		close_AsmM4Writers(m4_out);
		TIMING_END(job);
	}
	free_AsmM4Writers(m4_out);
}
