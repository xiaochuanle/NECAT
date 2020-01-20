#include "mem_finder.h"

#include "../common/oc_assert.h"
#include "../klib/ksort.h"

#define mem_soff_lt(a, b) ( \
    ((a).reference_offset < (b).reference_offset) \
    || \
    ((a).reference_offset == (b).reference_offset && (a).query_offset < (b).query_offset) \
)
KSORT_INIT(mem_soff_lt, MaximalExactMatch, mem_soff_lt);

#define kmif_hash_lt(a, b) ( \
    ((a).hash < (b).hash) \
    || \
    ((a).hash == (b).hash && (a).offset < (b).offset) \
)

KSORT_INIT(kmif_hash_lt, KmerInfo, kmif_hash_lt);

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataNew(int kmer_size, int window_size, int mem_size)
{
    MaximalExactMatchWorkData* data = 
        (MaximalExactMatchWorkData*)calloc(1, sizeof(MaximalExactMatchWorkData));
    kv_init(data->qry_kmif_list);
    kv_reserve(KmerInfo, data->qry_kmif_list, 100000);
    kv_init(data->ref_kmif_list);
    kv_reserve(KmerInfo, data->ref_kmif_list, 100000);
    kv_init(data->chain_seed_list);
    kv_reserve(ChainDpSeed, data->chain_seed_list, 100000);
    kv_init(data->can_list);
    kv_init(data->chain_seed_info_list);
    data->kmer_size = kmer_size;
    data->window_size = window_size;
    data->mem_size = mem_size;
    data->reference = NULL;
    data->ref_size = 0;
    return data;
}

MaximalExactMatchWorkData*
MaximalExactMatchWorkDataFree(MaximalExactMatchWorkData* data)
{
    kv_destroy(data->qry_kmif_list);
    kv_destroy(data->ref_kmif_list);
    kv_destroy(data->chain_seed_list);
    kv_destroy(data->can_list);
    kv_destroy(data->chain_seed_info_list);
    free(data);
    return NULL;
}

int
build_kmif_list(const u8* sequence, const size_t seq_len, int kmer_size, int window_size, vec_kmif* kmif_list)
{
    kv_clear(*kmif_list);
    KmerInfo kmif = { 0, 0, 0 };
    if (seq_len < kmer_size) return 0;
    const u8* s = sequence;
    const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;

	if (!intersect) {
		for (size_t j = 0; j <= seq_len - kmer_size; j += window_size) {
			u64 hash = 0;
			for (int k = 0; k < kmer_size; ++k) {
				size_t pos = j + k;
				u8 c = s[pos];
                oc_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
            kmif.hash = hash;
            kmif.offset = j;
            kv_push(KmerInfo, *kmif_list, kmif);
		}
	} else {
		u64 hash = 0;
		for (int j = 0; j < kmer_size; ++j) {
			size_t pos = j;
			u8 c = s[pos];
            oc_assert(c >= 0 && c < 4);
			hash = (hash << 2) | c;
		}
		kmif.offset = hash;
        kmif.offset = 0;
        kv_push(KmerInfo, *kmif_list, kmif);
		for (u64 j = window_size; j <= seq_len - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
				size_t pos = j + k;
                oc_assert(pos < seq_len, "p = %d, read_size = %d, j = %d, k = %d, stride = %d", 
                    pos, seq_len, j, k, stride);
				u8 c = s[pos];
                oc_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
			kmif.hash = hash;
            kmif.offset = j;
            kv_push(KmerInfo, *kmif_list, kmif);
		}
	}

    return kv_size(*kmif_list);
}

static void
sort_kmif_list(vec_kmif* kmif_list)
{
    size_t n = kv_size(*kmif_list);
    KmerInfo* a = kv_data(*kmif_list);
    ks_introsort_kmif_hash_lt(n, a);
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && a[j].hash == a[i].hash) ++j;
        a[i].occ = j - i;
        i = j;
    }
}

void 
MaximalExactMatchWorkData_Init(MaximalExactMatchWorkData* data, const u8* reference, const int ref_size)
{
    build_kmif_list(reference, ref_size, data->kmer_size, 1, &data->ref_kmif_list);
    sort_kmif_list(&data->ref_kmif_list);
    kv_clear(data->chain_seed_list);
    kv_clear(data->chain_seed_info_list);
    kv_clear(data->can_list);
    data->reference = reference;
    data->ref_size = ref_size;
}

static void
find_kmer_match(vec_kmif* query_kmif_list,
    vec_kmif* target_kmif_list,
    const int kmer_size,
    const int max_occ,
    vec_chain_seed* mem_list)
{
    kv_clear(*mem_list);
    ChainDpSeed m;
    size_t qn = kv_size(*query_kmif_list);
    KmerInfo* qki = kv_data(*query_kmif_list);
    size_t tn = kv_size(*target_kmif_list);
    KmerInfo* tki = kv_data(*target_kmif_list);
    size_t qi = 0, ti = 0;
    while (1) {
        while (qi < qn && qki[qi].hash < tki[ti].hash) qi += qki[qi].occ;
        if (qi >= qn) break;
        while (ti < tn && tki[ti].hash < qki[qi].hash) ti += tki[ti].occ;
        if (ti >= tn) break;
        if (qki[qi].hash == tki[ti].hash) {
            int n_occ = qki[qi].occ * tki[ti].occ;
            if (n_occ <= max_occ) {
                for (int i = 0; i < qki[qi].occ; ++i) {
                    m.match_size = kmer_size;
                    m.qoff = qki[qi + i].offset;
                    for (int j = 0; j < tki[ti].occ; ++j) {
                        m.soff = tki[ti + j].offset;
                        kv_push(ChainDpSeed, *mem_list, m);
                    }
                }
            }
            qi += qki[qi].occ;
            ti += tki[ti].occ;
        }
        if (qi >= qn || ti >= tn) break;
    }
}

static int
left_extend(const u8* query,
    const u8* target,
    int qoff,
    int toff)
{
    int n = 0;
    while (qoff && toff) {
        --qoff;
        --toff;
        if (query[qoff] == target[toff]) {
            ++n;
        } else {
            break;
        }
    }
    return n;
}

static int
right_extend(const u8* query,
    const u8* target,
    int qoff,
    int toff,
    int qsize,
    int tsize)
{
    int n = 0;
    while (qoff < qsize && toff < tsize) {
        if (query[qoff] == target[toff]) {
            ++n;
        } else {
            break;
        }
        ++qoff;
        ++toff;
    }
    return n;
}

static void
extend_kmer_match(vec_chain_seed* mem_list,
    const u8* query,
    const int qsize,
    const u8* target,
    const int tsize,
    const int min_mem_size)
{
    size_t n = kv_size(*mem_list);
    ChainDpSeed* a = kv_data(*mem_list);
    if (n == 0) return;
    ks_introsort_chain_seed_soff_lt(n, a);
    const u8* q = query;
    const u8* t = target;
    int ql, qr, tl, tr;
    for (size_t i = 0; i < n; ++i) {
        if (a[i].match_size == 0) continue;
        ql = a[i].qoff;
        tl = a[i].soff;
        qr = ql + a[i].match_size;
        tr = tl + a[i].match_size;
        int e = left_extend(q, t, ql, tl);
        oc_assert(ql >= e);
        ql -= e;
        oc_assert(tl >= e);
        tl -= e;
        e = right_extend(q, t, qr, tr, qsize, tsize);
        oc_assert(qr + e <= qsize);
        qr += e;
        oc_assert(tr + e <= tsize);
        tr += e;
        a[i].qoff = ql;
        a[i].soff = tl;
        a[i].match_size = qr - ql;

        size_t j = i + 1;
        while (j < n && a[j].soff < tr) {
            if (a[j].qoff < qr && a[j].soff < tr
                &&
                qr - a[j].qoff == tr - a[j].soff) {
                a[j].match_size = 0;
            }
            ++j;
        }

        if (a[i].match_size < min_mem_size) a[i].match_size = 0;
    }

    size_t i = 0;
    for (size_t k = 0; k < n; ++k) {
        if (a[k].match_size) a[i++] = a[k];
    }
    kv_resize(ChainDpSeed, *mem_list, i);
}

int
MaximalExactMatchWorkData_FindBestCandidate(
    MaximalExactMatchWorkData* data, 
    ChainDpWorkData* chain_data,
    const u8* read, const int read_size,
    GappedCandidate* can)
{
    const int verbose = 0;
    const int kMaxKmOcc = 20;
    build_kmif_list(read, read_size, data->kmer_size, data->window_size, &data->qry_kmif_list);
    sort_kmif_list(&data->qry_kmif_list);
    if (verbose) OC_LOG("number of km: %d", kv_size(data->qry_kmif_list));
    kv_clear(chain_data->seeds);
    find_kmer_match(&data->qry_kmif_list,
        &data->ref_kmif_list,
        data->kmer_size,
        kMaxKmOcc,
        &chain_data->seeds);
    if (verbose) OC_LOG("number of match: %d", kv_size(chain_data->seeds));
    extend_kmer_match(&chain_data->seeds, read, read_size, data->reference, data->ref_size, data->mem_size);
    if (verbose) OC_LOG("number of match: %d", kv_size(chain_data->seeds));
    if (verbose) {
        for (size_t i = 0; i < kv_size(chain_data->seeds); ++i) {
            OC_LOG("%d:\t%d\t%d\t%d\n", i, 
                kv_A(chain_data->seeds, i).qoff,
                kv_A(chain_data->seeds, i).soff,
                kv_A(chain_data->seeds, i).match_size);
        }
    }
    ks_introsort_chain_seed_soff_lt(kv_size(chain_data->seeds), kv_data(chain_data->seeds));
    int best_mem_index = -1;
    int best_mem_score = 0;
    vec_chain_seed* chain_seeds = &data->chain_seed_list;
    //HBN_LOG("min_cnt = %d, min_score = %d", chain_data->min_cnt, chain_data->min_score);
    int r = chaining_seeds(chain_data, &best_mem_index, &best_mem_score, chain_seeds);
    if (!r) return r;
    can->qbeg = kv_A(*chain_seeds, 0).qoff;
    can->sbeg = kv_A(*chain_seeds, 0).soff;
    can->qend = kv_back(*chain_seeds).qoff + kv_back(*chain_seeds).match_size;
    can->send = kv_back(*chain_seeds).soff + kv_back(*chain_seeds).match_size;
    oc_assert(best_mem_index >= 0);
    oc_assert(best_mem_index < kv_size(chain_data->seeds));
    can->qoff = kv_A(chain_data->seeds, best_mem_index).qoff;
    can->qoff = kv_A(chain_data->seeds, best_mem_index).soff;
    can->qsize = read_size;
    can->ssize = data->ref_size;
    can->score = best_mem_score;
    //DUMP_GAPPED_CANDIDATE(fprintf, stderr, *can);
    return r;
}

int
MaximalExactMatchWorkData_FindCandidates(
    MaximalExactMatchWorkData* data, 
    ChainDpWorkData* chain_data,
    const u8* read, 
    const int read_size)
{
    const int verbose = 0;
    const int kMaxKmOcc = 20;
    build_kmif_list(read, read_size, data->kmer_size, data->window_size, &data->qry_kmif_list);
    sort_kmif_list(&data->qry_kmif_list);
    if (verbose) OC_LOG("number of km: %d", kv_size(data->qry_kmif_list));
    kv_clear(chain_data->seeds);
    find_kmer_match(&data->qry_kmif_list,
        &data->ref_kmif_list,
        data->kmer_size,
        kMaxKmOcc,
        &chain_data->seeds);
    if (verbose) OC_LOG("number of match: %d", kv_size(chain_data->seeds));
    extend_kmer_match(&chain_data->seeds, read, read_size, data->reference, data->ref_size, data->mem_size);
    if (verbose) OC_LOG("number of match: %d", kv_size(chain_data->seeds));
    if (verbose) {
        for (size_t i = 0; i < kv_size(chain_data->seeds); ++i) {
            OC_LOG("%d:\t%d\t%d\t%d\n", i, 
                kv_A(chain_data->seeds, i).qoff,
                kv_A(chain_data->seeds, i).soff,
                kv_A(chain_data->seeds, i).match_size);
        }
    }
    ks_introsort_chain_seed_soff_lt(kv_size(chain_data->seeds), kv_data(chain_data->seeds));
    chaining_find_candidates(chain_data, &data->can_list, &data->chain_seed_list, &data->chain_seed_info_list);
    return kv_size(data->can_list) > 0;
}
