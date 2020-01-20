#include "consensus_one_partition.h"

#include "consensus_aux.h"
#include "consensus_one_read.h"
#include "../common/ontcns_aux.h"
#include "../common/record_reader.h"
#include "../partition_candidates/pcan_aux.h"

static void
load_partition_candidates(PackedDB* reads,
						  const char* can_path, 
						  const int pid, 
						  int* min_read_id, 
						  int* max_read_id, 
						  vec_pcan* candidates)
{
	kv_clear(*candidates);
	size_t file_bytes = 0;
	size_t num_can = 0;
	
	{
		char job[1024];
		sprintf(job, "loading candidates");
		TIMING_START(job);

		new_kstring(pname);
		make_partition_name(can_path, pid, &pname);
		file_bytes = FILE_SIZE(kstr_data(pname));
		num_can = file_bytes / sizeof(PackedGappedCandidate);
		if (num_can == 0) return;
		kv_resize(PackedGappedCandidate, *candidates, num_can);
		DFOPEN(can_in, kstr_str(pname), "rb");
		FREAD(kv_data(*candidates), sizeof(PackedGappedCandidate), num_can, can_in);
		FCLOSE(can_in);
		TIMING_END(job);
	}
	
	{
		char job[1024];
		sprintf(job, "sorting candidates");
		TIMING_START(job);
		ks_introsort_PackedGappedCandidate_SidLT(kv_size(*candidates), kv_data(*candidates));
		TIMING_END(job);
	}
	
	{
		PackedGappedCandidate* pcan = kv_data(*candidates);
		*min_read_id = pcan_sid(*pcan);
		pcan += (num_can - 1);
		*max_read_id = pcan_sid(*pcan);
	}
}

static void
normalise_partition_candidates(PackedDB* reads,
                          CnsReads* cns_reads,
						  vec_pcan* candidates)
{
    size_t num_can = kv_size(*candidates);
    if (reads) {
	    for (size_t i = 0; i < num_can; ++i) {
			PackedGappedCandidate* pcan = kv_data(*candidates) + i;
            if (pcan_sdir(*pcan) == FWD) continue;
			const int qid = pcan_qid(*pcan);
			const int sid = pcan_sid(*pcan);
			const u32 qsize = PDB_SEQ_SIZE(reads, qid);
			const u32 ssize = PDB_SEQ_SIZE(reads, sid);
			normalise_pcan_sdir(pcan, qsize, ssize);
	    }
    } else {
        size_t i = 0;
        while (i < num_can) {
            PackedGappedCandidate* pcan = kv_data(*candidates) + i;
            int sid = pcan_sid(*pcan);
            size_t j = i + 1;
            while (j < num_can) {
                pcan = kv_data(*candidates) + j;
                int ssid = pcan_sid(*pcan);
                if (ssid != sid) break;
                ++j;
            }
            int n = j - i;
            if (n > MAX_EXAMINED_CAN) n = MAX_EXAMINED_CAN;
            const u32 ssize = cns_reads_seq_size(cns_reads, sid);
            for (int k = 0; k < n; ++k) {
                pcan = kv_data(*candidates) + i + k;
                if (pcan_sdir(*pcan) == FWD) continue;
                const int qid = pcan_qid(*pcan);
                const u32 qsize = cns_reads_seq_size(cns_reads, qid);
                normalise_pcan_sdir(pcan, qsize, ssize);
            }
            i = j;
        }
    }
}

static void*
consensus_thread(void* arg)
{
	CnsData* cns_data = (CnsData*)(arg);
	size_t can_sid;
	size_t can_eid;
	while (extract_candidate_range(cns_data, &can_sid, &can_eid)) {
		consensus_one_read(cns_data, can_sid, can_eid);
	}
	return NULL;
}

void
consensus_one_partition(const char* pac_reads_dir,
						PackedDB* reads,
						const char* can_path,
						CnsOptions* options,
						FILE* cns_out,
						FILE* raw_out,
						const int pid)
{
	char job_name[1024];
	sprintf(job_name, "consensus partition %d", pid);
	TIMING_START(job_name);
	
	size_t next_can_id = 0;
	OcMutex can_id_lock;
	pthread_mutex_init(&can_id_lock, NULL);
	OcMutex cns_out_lock;
	pthread_mutex_init(&cns_out_lock, NULL);
	OcMutex raw_out_lock;
	pthread_mutex_init(&raw_out_lock, NULL);
	OcMutex read_id_lock;
	pthread_mutex_init(&read_id_lock, NULL);
	ReadIdPool* corrected_read_ids = new_ReadIdPool(&read_id_lock);
	new_kvec(vec_pcan, candidates);
	int max_read_id = 0, min_read_id = 0;
	load_partition_candidates(reads, can_path, pid, &min_read_id, &max_read_id, &candidates);
	CnsReads* cns_reads = NULL;
	if (reads == NULL) cns_reads = cns_reads_new(kv_data(candidates), kv_size(candidates), MAX_EXAMINED_CAN, pac_reads_dir);
    normalise_partition_candidates(reads, cns_reads, &candidates);
	
	CnsData** cns_data_pool = (CnsData**)malloc(options->num_threads * sizeof(CnsData*));
	int thread_id = 0;
	pthread_mutex_t thread_id_lock;
	pthread_mutex_init(&thread_id_lock, NULL);
	for (int i = 0; i < options->num_threads; ++i) {
		cns_data_pool[i] = new_CnsData(cns_reads,
									reads,
									&candidates,
									&next_can_id,
									&can_id_lock,
									options,
									cns_out,
									&cns_out_lock,
									raw_out,
									&raw_out_lock,
									corrected_read_ids, 
									&thread_id, 
									&thread_id_lock);
	}
	
	pthread_t* tids = (pthread_t*)malloc(sizeof(pthread_t) * options->num_threads);
	for (int i = 0; i < options->num_threads; ++i) {
		pthread_create(&tids[i], NULL, consensus_thread, cns_data_pool[i]);
	}
	for (int i = 0; i < options->num_threads; ++i) {
		pthread_join(tids[i], NULL);
	}
	free(tids);
	
	for (int i = 0; i < options->num_threads; ++i) {
		cns_data_pool[i] = free_CnsData(cns_data_pool[i]);
	}
	free(cns_data_pool);
	
if (reads) {
	new_CnsSeq(cns_seq);
	cns_seq.num_ovlps = 0;
	cns_seq.ident_cutoff = 0.0;
	cns_seq.num_can = 0;
	for (int i = min_read_id; i < max_read_id; ++i) {
		if (!read_id_exists(corrected_read_ids, i)) {
			cns_seq.id = i;
			cns_seq.hdr = reads ? PDB_SEQ_NAME(reads, i) : cns_reads_seq_hdr(cns_reads, i);
			cns_seq.left = 0;
			cns_seq.right = reads ? PDB_SEQ_SIZE(reads, i) : cns_reads_seq_size(cns_reads, i);
			cns_seq.org_seq_size = cns_seq.right;
			if (reads) pdb_extract_sequence(reads, i, FWD, &cns_seq.cns_seq);
			else cns_reads_extract_sequence(cns_reads, i, FWD, &cns_seq.cns_seq);
			for (size_t k = 0; k != kstr_size(cns_seq.cns_seq); ++k) {
				kstr_A(cns_seq.cns_seq, k) = DecodeDNA(kstr_A(cns_seq.cns_seq, k));
			}
			DUMP_CNS_SEQ(fprintf, raw_out, cns_seq);
		}
	}
	free_CnsSeq(cns_seq);
}
	
	if (cns_reads) cns_reads_free(cns_reads);
	free_kvec(candidates);
	free_ReadIdPool(corrected_read_ids);
	
	TIMING_END(job_name);
}
