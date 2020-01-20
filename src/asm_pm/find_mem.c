#include "find_mem.h"

#include "../common/oc_assert.h"

#define kmif_hash_lt(a, b) ( \
    ((a).hash < (b).hash) \
    || \
    ((a).hash == (b).hash && (a).offset < (b).offset) \
)

KSORT_INIT(kmif_hash_lt, KmerInfo, kmif_hash_lt);

#define mem_soff_lt(a, b) ( \
    ((a).reference_offset < (b).reference_offset) \
    || \
    ((a).reference_offset == (b).reference_offset && (a).query_offset < (b).query_offset) \
)

KSORT_INIT(mem_soff_lt, MaximalExactMatch, mem_soff_lt);

int
build_kmif_list(vec_u8* sequence, int kmer_size, int window_size, vec_kmif* kmif_list)
{
    kv_clear(*kmif_list);
    KmerInfo kmif = { 0, 0, 0 };
    const size_t seq_len = kv_size(*sequence);
    if (seq_len < kmer_size) return 0;
    const u8* s = kv_data(*sequence);
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

void
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
find_kmer_match(vec_kmif* query_kmif_list,
    vec_kmif* target_kmif_list,
    const int kmer_size,
    const int max_occ,
    vec_mem* mem_list)
{
    kv_clear(*mem_list);
    MaximalExactMatch m;
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
                    m.query_offset = qki[qi + i].offset;
                    for (int j = 0; j < tki[ti].occ; ++j) {
                        m.reference_offset = tki[ti + j].offset;
                        kv_push(MaximalExactMatch, *mem_list, m);
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

void
extend_kmer_match(vec_mem* mem_list,
    vec_u8* query,
    vec_u8* target,
    const int min_mem_size)
{
    size_t n = kv_size(*mem_list);
    MaximalExactMatch* a = kv_data(*mem_list);
    if (n == 0) return;
    ks_introsort_mem_soff_lt(n, a);
    const u8* q = kv_data(*query);
    const int qsize = kv_size(*query);
    const u8* t = kv_data(*target);
    const int tsize = kv_size(*target);
    int ql, qr, tl, tr;
    for (size_t i = 0; i < n; ++i) {
        if (a[i].match_size == 0) continue;
        ql = a[i].query_offset;
        tl = a[i].reference_offset;
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
        a[i].query_offset = ql;
        a[i].reference_offset = tl;
        a[i].match_size = qr - ql;

        size_t j = i + 1;
        while (j < n && a[j].reference_offset < tr) {
            if (a[j].query_offset < qr && a[j].reference_offset < tr
                &&
                qr - a[j].query_offset == tr - a[j].reference_offset) {
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
    kv_resize(MaximalExactMatch, *mem_list, i);
}

int
compute_align_range_1(ChainWorkData* chain_data,
    vec_u8* query,
    vec_u8* target,
    vec_kmif* target_kmif_list,
    const int kmer_size,
    const int window_size,
    const int min_mem_size,
    GappedCandidate* can,
    vec_mem* chain_mems)
{
    const int verbose = 0;
    const int kMaxKmOcc = 20;
    new_kvec(vec_kmif, query_kmif_list);
    new_kvec(vec_mem, mem_list);
    build_kmif_list(query, kmer_size, window_size, &query_kmif_list);
    sort_kmif_list(&query_kmif_list);
    if (verbose) OC_LOG("number of km: %d", kv_size(query_kmif_list));
    find_kmer_match(&query_kmif_list,
        target_kmif_list,
        kmer_size,
        kMaxKmOcc,
        &mem_list);
    if (verbose) OC_LOG("number of match: %d", kv_size(mem_list));
    extend_kmer_match(&mem_list, query, target, min_mem_size);
    if (verbose) OC_LOG("number of match: %d", kv_size(mem_list));
    ks_introsort_mem_soff_lt(kv_size(mem_list), kv_data(mem_list));
    int r = mem_find_best_can(chain_data, 
                kv_data(mem_list),
                kv_size(mem_list),
                0,
                0,
                0,
                0,
                0,
                can,
                chain_mems);
    kv_destroy(query_kmif_list);
    kv_destroy(mem_list);
    return r;
}