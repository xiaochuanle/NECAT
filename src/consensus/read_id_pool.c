#include "read_id_pool.h"

#include "../klib/khash.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"

KHASH_SET_INIT_INT(32)
typedef khash_t(32) khash_int;

ReadIdPool*
new_ReadIdPool(OcMutex* id_lock)
{
	ReadIdPool* pool = (ReadIdPool*)malloc(sizeof(ReadIdPool));
	pool->read_id_pool = (void*)kh_init(32);
	pool->id_lock = id_lock;
	return pool;
}

ReadIdPool*
free_ReadIdPool(ReadIdPool* pool)
{
	kh_destroy(32, pool->read_id_pool);
	free(pool);
	return 0;
}

void
clear_ReadIdPool(ReadIdPool* pool)
{
	kh_clear_32(pool->read_id_pool);
}

BOOL
read_id_exists(ReadIdPool* pool, int id)
{
	khash_int* hid = (khash_int*)pool->read_id_pool;
	return kh_get(32, hid, id) != kh_end(hid);
}

void
add_read_id(ReadIdPool* pool, int id)
{
	khash_int* hid = (khash_int*)pool->read_id_pool;
	int ret;
	if (pool->id_lock) pthread_mutex_lock(pool->id_lock);
	kh_put(32, hid, id, &ret);
	if (pool->id_lock) pthread_mutex_unlock(pool->id_lock);
}

//////////////////////

KHASH_MAP_INIT_INT(map_int_int, int)
typedef khash_t(map_int_int) kmap_int_int;

static void
cns_reads_add_id(CnsReads* cns_reads, int global_id, int local_id)
{
	kmap_int_int* cns_read_ids = (kmap_int_int*)(cns_reads->cns_read_ids);
	khiter_t iter = kh_get(map_int_int, cns_read_ids, global_id);
	if (iter == kh_end(cns_read_ids)) {
		int ret;
		iter = kh_put(map_int_int, cns_read_ids, global_id, &ret);
		kh_value(cns_read_ids, iter) = local_id;
	}	
}

static int
cns_reads_change_id(CnsReads* cns_reads, int global_id, int local_id)
{
	kmap_int_int* cns_read_ids = (kmap_int_int*)(cns_reads->cns_read_ids);
	khiter_t iter = kh_get(map_int_int, cns_read_ids, global_id);
	if (iter != kh_end(cns_read_ids)) {
		kh_value(cns_read_ids, iter) = local_id;
		return 1;
	} else {
		return 0;
	}
}

static void
decode_sequence(kstring_t* seq)
{
    char*s = kstr_str(*seq);
    size_t size = kstr_size(*seq);
    for (size_t i = 0; i < size; ++i) {
        int c = s[i];
        oc_assert(c >= 0 && c < 4);
        c = "ACGT"[c];
        s[i] = c;
    }
}

CnsReads*
cns_reads_new(PackedGappedCandidate* cans, 
	const int ncan,
	const int max_can_per_template,
	const char* pac_reads_dir)
{
	CnsReads* cns_reads = (CnsReads*)calloc(1, sizeof(CnsReads));
	kmap_int_int* cns_read_ids = kh_init(map_int_int);
	cns_reads->cns_read_ids = (void*)(cns_read_ids);
	int i = 0;
	while (i < ncan) {
		int sid = pcan_sid(cans[i]);
		int j = i + 1;
		while (j < ncan && pcan_sid(cans[j]) == sid) ++j;
		int n = j - i;
		if (n > max_can_per_template) {
			ks_introsort_PackedGappedCandidate_CnsScoreGT(n, cans + i);
			n = max_can_per_template;
		}
		cns_reads_add_id(cns_reads, sid, 1);
		for (int k = 0; k < n; ++k) {
			int qid = pcan_qid(cans[i+k]);
			cns_reads_add_id(cns_reads, qid, 1);
		}
		i = j;
	}

	cns_reads->cns_reads = new_PackedDB();
	VolumesInfo* vis = load_volumes_info(pac_reads_dir);
	int global_id = 0;
	int local_id = 0;
	kseq_t* read = kseq_init(NULL);
	for (int i = 0; i < vis->num_volumes; ++i) {
		PackedDB* volume = new_PackedDB();
		pdb_load(volume, vi_volume_name(vis, i), TECH_NANOPORE);
		for (idx lid = 0; lid < PDB_NUM_SEQS(volume); ++lid) {
			int r = cns_reads_change_id(cns_reads, global_id, local_id);
			if (r) {
				pdb_extract_sequence(volume, lid, FWD, &read->seq);
                decode_sequence(&read->seq);
				const char* hdr = PDB_SEQ_NAME(volume, lid);
				kstr_clear(read->name);
				kputsn(hdr, strlen(hdr), &read->name);
				pdb_add_one_seq(cns_reads->cns_reads, read, TECH_NANOPORE);
				++local_id;
			}
			++global_id;
		}
		free_PackedDB(volume);
	}
	destroy_volumes_info(vis);

	OC_LOG("Load %d reads, total %zu bps", local_id, PDB_SIZE(cns_reads->cns_reads));
	return cns_reads;
}

CnsReads*
cns_reads_free(CnsReads* cns_reads)
{
	kh_destroy(map_int_int, cns_reads->cns_read_ids);
	free_PackedDB(cns_reads->cns_reads);
	free(cns_reads);
	return NULL;
}

int
cns_reads_local_id(CnsReads* cns_reads, int global_id)
{
	kmap_int_int* cns_read_ids = (kmap_int_int*)(cns_reads->cns_read_ids);
	khiter_t iter = kh_get(map_int_int, cns_read_ids, global_id);
	oc_assert (iter != kh_end(cns_read_ids), "global_id = %d", global_id);
	return kh_value(cns_read_ids, iter);
}

void
cns_reads_extract_sequence(CnsReads* cns_reads, int global_id, int strand, kstring_t* seq)
{
	int local_id = cns_reads_local_id(cns_reads, global_id);
	pdb_extract_sequence(cns_reads->cns_reads, local_id, strand, seq);
}

int
cns_reads_seq_size(CnsReads* cns_reads, int global_id)
{
	int local_id = cns_reads_local_id(cns_reads, global_id);
	return PDB_SEQ_SIZE(cns_reads->cns_reads, local_id);
}

const char*
cns_reads_seq_hdr(CnsReads* cns_reads, int global_id)
{
	int local_id = cns_reads_local_id(cns_reads, global_id);
	return PDB_SEQ_NAME(cns_reads->cns_reads, local_id);
}

void 
truncate_m4_list(M4Record* m4v, int* nm4, const int max_m4_perc_read, const double ident_perc)
{
	int n = *nm4;
	const int MaxNm4 = max_m4_perc_read;
	if (n > MaxNm4) {
		int j = 0;
		for (int i = 0; i < n; ++i) {
			if (m4v[i].ident_perc >= ident_perc) m4v[j++] = m4v[i];
		}
		n = j;
	}
	if (n > MaxNm4) {
		ks_introsort_M4Record_IdentGT(n, m4v);
		n = MaxNm4;
	}
	*nm4 = n;
}

CnsReads*
cns_reads_new_m4(M4Record* m4s,
	const int nm4,
	const int max_m4_per_template,
	const double min_ident_perc,
	const char* pac_reads_dir)
{
	CnsReads* cns_reads = (CnsReads*)calloc(1, sizeof(CnsReads));
	kmap_int_int* cns_read_ids = kh_init(map_int_int);
	cns_reads->cns_read_ids = (void*)(cns_read_ids);
	int i = 0;
	while (i < nm4) {
		int sid = m4s[i].sid;
		int j = i + 1;
		while (j < nm4 && m4s[j].sid == sid) ++j;
		int n = j - i;
		if (n > max_m4_per_template) {
			truncate_m4_list(m4s + i, &n, max_m4_per_template, min_ident_perc);
		}
		oc_assert(sid > 0);
		cns_reads_add_id(cns_reads, sid - 1, 1);
		for (int k = 0; k < n; ++k) {
			int qid = m4s[i+k].qid;
			oc_assert(qid > 0);
			cns_reads_add_id(cns_reads, qid-1, 1);
		}
		i = j;
	}

	cns_reads->cns_reads = new_PackedDB();
	VolumesInfo* vis = load_volumes_info(pac_reads_dir);
	int global_id = 0;
	int local_id = 0;
	kseq_t* read = kseq_init(NULL);
	for (int i = 0; i < vis->num_volumes; ++i) {
		PackedDB* volume = new_PackedDB();
		pdb_load(volume, vi_volume_name(vis, i), TECH_NANOPORE);
		for (idx lid = 0; lid < PDB_NUM_SEQS(volume); ++lid) {
			int r = cns_reads_change_id(cns_reads, global_id, local_id);
			if (r) {
				pdb_extract_sequence(volume, lid, FWD, &read->seq);
				const char* hdr = PDB_SEQ_NAME(volume, lid);
				kstr_clear(read->name);
				kputsn(hdr, strlen(hdr), &read->name);
				pdb_add_one_seq(cns_reads->cns_reads, read, TECH_NANOPORE);
				++local_id;
			}
			++global_id;
		}
		free_PackedDB(volume);
	}
	destroy_volumes_info(vis);

	OC_LOG("Load %d reads, total %ll bps", local_id, PDB_SIZE(cns_reads->cns_reads));
	return cns_reads;
}
