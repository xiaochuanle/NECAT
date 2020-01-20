#include "km_chain.h"

#include "../klib/ksort.h"

ChainWorkData*
chain_data_new()
{
    ChainWorkData* data = (ChainWorkData*)malloc(sizeof(ChainWorkData));
    kv_init(data->f);
    kv_init(data->p);
    kv_init(data->t);
    kv_init(data->v);
    kv_init(data->u);
    data->max_dist_qry = 3000;
    data->max_dist_ref = 3000;
    data->max_band_width = 500;
    data->max_skip = 25;
    data->min_cnt = 1;
    data->min_score = 100;
    data->dump_info = 0;
    return data;
}

ChainWorkData*
chain_data_free(ChainWorkData* data)
{
    kv_destroy(data->f);
    kv_destroy(data->p);
    kv_destroy(data->t);
    kv_destroy(data->v);
    kv_destroy(data->u);
    free(data);
    return data;
}

void
chain_data_reset(ChainWorkData* data, int n)
{
    kv_resize(int, data->f, n);
    kv_resize(int, data->p, n);
    kv_resize(int, data->t, n);
    kv_resize(int, data->v, n);
    kv_zero(int, data->t);
    kv_zero(int, data->f);
    kv_fill(data->p, -1);
    kv_resize(IntPair, data->u, n);
}

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

static void
remove_block_repetitive_kms(MaximalExactMatch* v_mem,
    int from,
    int to,
    int max_kms,
    int* last_i_,
    int min_ql,
    int max_qr,
    int min_sl,
    int max_sr)
{
    int last_i = *last_i_;
    int next_i = last_i;
    for (int k = from; k < to; ++k) {
        int qoff = v_mem[k].query_offset;
        int soff = v_mem[k].reference_offset;
        int r = (qoff >= min_ql) && (qoff < max_qr) && (soff >= min_sl) && (soff < max_sr);
        if (r) v_mem[next_i++] = v_mem[k];
    }
    int m = next_i - last_i;
    if (m > max_kms) {
        int stride = m / max_kms;
        int p = last_i;
        for (int k = last_i; k < next_i; k += stride) {
            v_mem[p++] = v_mem[k];
        }
        next_i = p;
    }
    *last_i_ = next_i;
}

static void
remove_repetitive_kms(MaximalExactMatch* v_mem, int* n_mem_)
{
   const int n_mem = *n_mem_;
   ks_introsort_mem_soff_lt(n_mem, v_mem);
   const idx kBlockSize = 2000;
   const int kMaxKm = 400;
   int i = 0;
   int last_i = 0;
   int min_ql = 738, max_qr = 43023;
   int min_sl = 3888, max_sr = 41511;
   
   while (i < n_mem) {
       int block_id = v_mem[i].reference_offset / kBlockSize;
       idx block_end = (block_id + 1) * kBlockSize;
       int j = i + 1;
       while (j < n_mem && v_mem[j].reference_offset < block_end) ++j;
       int m = j - i;
       if (m > kMaxKm) {
#if 0
           const int stride = m / kMaxKm;
           for (int k = i; k < j; k += stride) v_mem[last_i++] = v_mem[k];
           HBN_LOG("i = %d, j = %d, m = %d", i, j, m);
#endif
        remove_block_repetitive_kms(v_mem, i, j, kMaxKm, &last_i, min_ql, max_qr, min_sl, max_sr);
       } else {
           last_i += m;
       }
       i = j;
   }
   OC_LOG("before: %d, after: %d", n_mem, last_i);
   *n_mem_ = last_i;
}

static void
scoring_mems(ChainWorkData* data,
    MaximalExactMatch* v_mem,
    int n)
{
    const int max_dist_ref = data->max_dist_ref;
    const int max_dist_qry = data->max_dist_qry;
    const int band_width = data->max_band_width;
    const int max_skip = data->max_skip;
    int sum_cov = 0;
    for (int i = 0; i < n; ++i) sum_cov += v_mem[i].match_size;
    const int avg_cov = sum_cov / n;
    int st = 0;
    chain_data_reset(data, n);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    //remove_repetitive_kms(v_mem, &n);

    // fill the score and backtrack arrays
    for (int i = 0; i < n; ++i) {
        //if (!kmer_match_is_valid(v_mem[i])) continue;
        idx ri = v_mem[i].reference_offset;
        int max_j = -1;
        int qi = v_mem[i].query_offset;
        int cov = v_mem[i].match_size;
        int max_f = cov, n_skip = 0, min_d;
        while (st < i && ri > v_mem[st].reference_offset + max_dist_ref) ++st;
        for (int j = i - 1; j >= st; --j) {
            //if (!kmer_match_is_valid(v_mem[j])) continue;
            if (v_mem[j].query_offset + v_mem[j].match_size >= qi ||
                v_mem[j].reference_offset + v_mem[j].match_size >= ri) continue;
            idx dr = ri - v_mem[j].reference_offset;
            int dq = qi - v_mem[j].query_offset;
            int dd, sc, log_dd;
            if (dr == 0 || dq <= 0) continue;
            if (dq > max_dist_qry || dr > max_dist_ref) continue;
            dd = (dr > dq) ? (dr - dq) : (dq - dr);
            if (dd > band_width) continue;
            min_d = OC_MIN(dq, dr);
            sc = (min_d > cov) ? cov : OC_MIN(dq, dr);
            log_dd = dd ? ilog2_32(dd) : 0;
            sc -= (int)(dd * .01 * avg_cov) + (log_dd>>1);
            sc += f[j];
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                if (n_skip) --n_skip;
            } else if (t[j] == i) {
                if (++n_skip > max_skip) { break; }
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        f[i] = max_f;
        p[i] = max_j;
        // v[i] keeps the peak score up to i;
        // f[i] is the score ending at i, not always the peak score
        v[i] = (max_j >= 0 && v[max_j] > max_f) ? v[max_j] : max_f;
    }
}

void mem_chain_dp(ChainWorkData* data,
        MaximalExactMatch* mems,
        const int n,
        const int read_id,
        const int read_dir,
        const int read_size,
        const int ref_id,
        const idx ref_size,
        vec_can* cans)
{
    if (n == 0) return;
    GappedCandidate can;
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    scoring_mems(data, mems, n);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    IntPair* u = kv_data(data->u);

    // find the ending position of chains
    memset(t, 0, sizeof(int) * n);
    for (int i = 0; i < n; ++i) {
        if (p[i] >= 0) t[p[i]] = 1;
    }
    int n_u;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) ++n_u;
    }
    if (n_u == 0) return;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) {
            int j = i;
            while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maxizes f
            if (j < 0) {
                j = i;
            }
            u[n_u].first = f[j];
            u[n_u].second = j;
            ++n_u;
        }
    }
    ks_introsort_IntPair(n_u, u);
    // reverse u, such that highest scoring chains appear first
    for (int i = 0; i < n_u>>1; ++i) {
        IntPair tmp = u[i];
        u[i] = u[n_u - i - 1];
        u[n_u - i - 1] = tmp;
    }

    // backtrack
    memset(t, 0, sizeof(int) * n);
    int n_v, k;
    for (int i = n_v = k = 0; i < n_u; ++i) { // start from the highest score
        int n_v0 = n_v, k0 = k;
        int j = u[i].second;
        //HBN_LOG("score = %d", u[i].first);
        do {
            v[n_v++] = j;
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first;
                u[k].second = n_v - n_v0;
                ++k;
//for (int fg = n_v0; fg < n_v; ++fg) DUMP_MEM_STD(mems[v[fg]]);
////HBN_LOG("number of mems: %d", n_v - n_v0);
                can.qend = mems[v[n_v0]].query_offset + mems[v[n_v0]].match_size;
                can.send = mems[v[n_v0]].reference_offset + mems[v[n_v0]].match_size;
                can.qbeg = mems[v[n_v-1]].query_offset;
                can.sbeg = mems[v[n_v-1]].reference_offset;
                can.score = u[k-1].first;
                can.qdir = mems[v[n_v0]].match_size;
                can.qoff = mems[v[n_v0]].query_offset;
                can.soff = mems[v[n_v0]].reference_offset;
                can.qid = read_id;
                can.qdir = read_dir;
                can.qsize = read_size;
                can.sid = ref_id;
                can.sdir = FWD;
                can.ssize = ref_size;
                //HBN_LOG("dump_info = %d", data->dump_info);
                if (data->dump_info == 1) {
                    OC_LOG("find candidate:\t");
                    DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
                    for (int x = n_v0; x < n_v; ++x) DUMP_MEM(fprintf, stdout, mems[v[x]]);
                    data->dump_info = -1;
                }
                kv_push(GappedCandidate, *cans, can);
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
//for (int fg = n_v0; fg < n_v; ++fg) DUMP_MEM_STD(mems[v[fg]]);
//HBN_LOG("number of mems: %d", n_v - n_v0);
                can.qend = mems[v[n_v0]].query_offset + mems[v[n_v0]].match_size;
                can.send = mems[v[n_v0]].reference_offset + mems[v[n_v0]].match_size;
                can.qbeg = mems[v[n_v-1]].query_offset;
                can.sbeg = mems[v[n_v-1]].reference_offset;
                can.score = u[k-1].first;
                can.qdir = mems[v[n_v0]].match_size;
                can.qoff = mems[v[n_v0]].query_offset;
                can.soff = mems[v[n_v0]].reference_offset;
                can.qid = read_id;
                can.qdir = read_dir;
                can.qsize = read_size;
                can.sid = ref_id;
                can.sdir = FWD;
                can.ssize = ref_size;
                if (data->dump_info == 1) {
                    OC_LOG("find candidate: score1 = %d, score2 = %d, min_score = %d\t", u[i].first, f[j], min_score);
                    DUMP_GAPPED_CANDIDATE(fprintf, stdout, can);
                    for (int x = n_v0; x < n_v; ++x) DUMP_MEM(fprintf, stdout, mems[v[x]]);
                    data->dump_info = -1;
                }
                kv_push(GappedCandidate, *cans, can);
            }
        }
        if (k0 == k) n_v = n_v0;
    }
}

int mem_find_best_can(ChainWorkData* data,
        MaximalExactMatch* mems,
        const int n,
        const int read_id,
        const int read_dir,
        const int read_size,
        const int ref_id,
        const idx ref_size,
        GappedCandidate* best_can,
        vec_mem* chain_mems)
{
    if (n == 0) return 0;
    GappedCandidate can;
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    if (0) {
        OC_LOG("********************, number of mems: %d", n);
        for (int i = 0; i < n; ++i) {
            printf("%d\t", i);
            DUMP_MEM_STD(mems[i]);
        }
    }
    scoring_mems(data, mems, n);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    IntPair* u = kv_data(data->u);

    // find the ending position of chains
    memset(t, 0, sizeof(int) * n);
    for (int i = 0; i < n; ++i) {
        if (p[i] >= 0) t[p[i]] = 1;
    }
    int n_u;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) ++n_u;
    }
    if (n_u == 0) return 0;
    for (int i = n_u = 0; i < n; ++i) {
        if (t[i] == 0 && v[i] >= min_score) {
            int j = i;
            while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maxizes f
            if (j < 0) {
                j = i;
            }
            u[n_u].first = f[j];
            u[n_u].second = j;
            ++n_u;
        }
    }
    ks_introsort_IntPair(n_u, u);
    // reverse u, such that highest scoring chains appear first
    for (int i = 0; i < n_u>>1; ++i) {
        IntPair tmp = u[i];
        u[i] = u[n_u - i - 1];
        u[n_u - i - 1] = tmp;
    }

    // backtrack
    memset(t, 0, sizeof(int) * n);
    int n_v, k;
    int find_can = 0;
    for (int i = n_v = k = 0; i < n_u; ++i) { // start from the highest score
        int n_v0 = n_v, k0 = k;
        int j = u[i].second;
        //HBN_LOG("score = %d", u[i].first);
        do {
            v[n_v++] = j;
            t[j] = 1;
            j = p[j];
        } while (j >= 0 && t[j] == 0);
        if (j < 0) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first;
                u[k].second = n_v - n_v0;
                ++k;
                find_can = 1;
                //HBN_LOG("find candidate");
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
                find_can = 1;
               // HBN_LOG("find candidate");
            }
        }
        if (find_can) {
            can.qend = mems[v[n_v0]].query_offset + mems[v[n_v0]].match_size;
            can.send = mems[v[n_v0]].reference_offset + mems[v[n_v0]].match_size;
            can.qbeg = mems[v[n_v-1]].query_offset;
            can.sbeg = mems[v[n_v-1]].reference_offset;
            can.score = u[k-1].first;
            can.qdir = mems[v[n_v0]].match_size;
            can.qoff = mems[v[n_v0]].query_offset + mems[v[n_v0]].match_size / 2;
            can.soff = mems[v[n_v0]].reference_offset + mems[v[n_v0]].match_size / 2;
            can.qid = read_id;
            can.qdir = read_dir;
            can.qsize = read_size;
            can.sid = ref_id;
            can.sdir = FWD;
            can.ssize = ref_size;
            //DUMP_GAPPED_CANDIDATE_STD(can);
            *best_can = can;
            kv_clear(*chain_mems);
            for (int x = n_v; x > n_v0; --x) {
                int y = v[x-1];
                if (0) {
                    if (mems[y].match_size >= 20 || x == n_v || x == n_v0+1) {
                        kv_push(MaximalExactMatch, *chain_mems, mems[y]);
                    }
                } else {
                    kv_push(MaximalExactMatch, *chain_mems, mems[y]);
                }
            }
            size_t n = kv_size(*chain_mems) / 2;
            //OC_LOG("number of mems: %d", n);
            MaximalExactMatch m = kv_A(*chain_mems, n);
            best_can->qoff = m.query_offset + m.match_size / 2;
            best_can->soff = m.reference_offset + m.match_size / 2;
            return 1;
        }

        if (k0 == k) n_v = n_v0;
    }    
    return 0;
}
