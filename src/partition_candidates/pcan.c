#include "pcan.h"

#include <assert.h>

#include "pcan_aux.h"
#include "pcan_options.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../common/ontcns_aux.h"
#include "../common/record_reader.h"

static FILE* pcan_in = NULL;
static pthread_mutex_t pcan_read_lock;
static CandidatePartitionWriter* pcan_writer = NULL;
static pthread_mutex_t pcan_write_lock;
static int min_read_id = 0;
static int max_read_id = 0;
static int batch_size = 0;

static size_t
load_pcan(vec_pcan* pcanv)
{
	const size_t buff_size = U64_ONE * 256 * 1024 * 1024;
	const size_t N = buff_size / sizeof(PackedGappedCandidate);
	const size_t M = 2 * N;
	kv_clear(*pcanv);
	kv_reserve(PackedGappedCandidate, *pcanv, M);
	size_t n = 0;
	pthread_mutex_lock(&pcan_read_lock);
	n = fread(kv_data(*pcanv), sizeof(PackedGappedCandidate), N, pcan_in);
	pthread_mutex_unlock(&pcan_read_lock);
	kv_resize(PackedGappedCandidate, *pcanv, n);
	return n;
}

#define id_in_range(id) ((id) >= min_read_id && (id) < max_read_id)

static void*
pcan_func(void* arg)
{
	new_kvec(vec_pcan, pcanv);
	size_t n;
	size_t m;
	PackedGappedCandidate pcan;
	new_kvec(vec_size_type, idx_range);
	while ((n = load_pcan(&pcanv))) {
		m = 0;
		PackedGappedCandidate* cans = kv_data(pcanv);
		for (size_t i = 0; i < n; ++i) {
			int qid = pcan_qid(cans[i]);
			int sid = pcan_sid(cans[i]);
			BOOL r = id_in_range(qid) || id_in_range(sid);
			if (r) cans[m++] = cans[i];
		}
		if (m == 0) continue;
		n = m;
		
		for (size_t i = 0; i < m; ++i) {
			BOOL sid_is_in = FALSE;
			int sid = pcan_sid(cans[i]);
			if (id_in_range(sid)) {
				sid_is_in = TRUE;
			}
			int qid = pcan_qid(cans[i]);
			if (id_in_range(qid)) {
				change_pcan_roles(&cans[i], &pcan);
				if (sid_is_in) cans[n++] = pcan;
				else cans[i] = pcan;
			}
		}
		
		ks_introsort_PackedGappedCandidate_SidLT(n, cans);
		kv_clear(idx_range);
		kv_push(size_t, idx_range, 0);
		size_t i = 0;
		while (i < n) {
			const int sid = pcan_sid(cans[i]);
			const int bid = sid / batch_size;
			const int sid_from = bid * batch_size;
			const int sid_to = sid_from + batch_size;
			size_t j = i + 1;
			while (j < n && pcan_sid(cans[j]) < sid_to) ++j;
			kv_push(size_t, idx_range, j);
			i = j;
		}
		oc_assert(kv_back(idx_range) == n);
		
		pthread_mutex_lock(&pcan_write_lock);
		oc_assert(kv_size(idx_range) > 1, "n = %d, size = %d", n, kv_size(idx_range));
		for (i = 0; i < kv_size(idx_range) - 1; ++i) {
			size_t from = kv_A(idx_range, i);
			size_t to = kv_A(idx_range, i + 1);
			m = to - from;
			int sid = pcan_sid(cans[from]);
			int bid = (sid - min_read_id) / batch_size;
			oc_assert(bid < pcan_writer->num_output_files);
			FWRITE(cans + from, sizeof(PackedGappedCandidate), m, pcan_writer->file_list[bid]);
		}
		pthread_mutex_unlock(&pcan_write_lock);
	}
	
	return NULL;
}

void
pcan_main(PcanOptions* options,
		  const char* wrk_dir,
		  const char* can_path)
{
	int num_reads = load_num_reads(wrk_dir);
	int num_batches = (num_reads + options->batch_size - 1) / options->batch_size;
	dump_num_partitions(can_path, num_batches);
	pcan_writer = new_CandidatePartitionWriter(options->num_output_files);
	pthread_mutex_init(&pcan_read_lock, NULL);
	pthread_mutex_init(&pcan_write_lock, NULL);
	batch_size = options->batch_size;
	char job[1024];
	const int num_threads = OC_MIN(8, options->num_threads);
	pthread_t job_ids[num_threads];
	
	for (int fid = 0; fid < num_batches; fid += options->num_output_files) {
		int sfid = fid;
		int efid = OC_MIN(sfid + options->num_output_files, num_batches);
		int nfid = efid - sfid;
		min_read_id = sfid * batch_size;
		max_read_id = efid * batch_size;
		sprintf(job, "dumping candidates for partitions [%d, %d)", sfid, efid);
		TIMING_START(job);
		open_CandidatePartitionWriter(nfid, can_path, sfid, pcan_writer);
		FOPEN(pcan_in, can_path, "rb");
		
		for (int i = 0; i < num_threads; ++i) {
			pthread_create(job_ids + i, NULL, pcan_func, NULL);
		}
		
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(job_ids[i], NULL);
		}

		FCLOSE(pcan_in);
		close_CandidatePartitionWriter(pcan_writer);
		TIMING_END(job);
	}
	
	free_CandidatePartitionWriter(pcan_writer);
}
