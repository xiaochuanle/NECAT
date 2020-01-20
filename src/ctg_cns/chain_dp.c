#include "chain_dp.h"

#include "../klib/ksort.h"

#define chain_seed_soff_lt(a, b) ( \
    ((a).soff < (b).soff) \
    || \
    ((a).soff == (b).soff && (a).qoff < (b).qoff) \
)
KSORT_INIT(chain_seed_soff_lt, ChainDpSeed, chain_seed_soff_lt);

ChainDpWorkData*
ChainDpWorkDataNew(int min_seed_cnt, int min_can_score)
{
    ChainDpWorkData* data = (ChainDpWorkData*)calloc(1, sizeof(ChainDpWorkData));
    kv_init(data->f);
    kv_init(data->p);
    kv_init(data->t);
    kv_init(data->v);
    kv_init(data->u);
    kv_init(data->seeds);
    data->max_dist_qry = 3000;
    data->max_dist_ref = 3000;
    data->max_band_width = 1500;
    data->max_skip = 25;
    data->min_cnt = min_seed_cnt;
    data->min_score = min_can_score;
    return data;
}

ChainDpWorkData*
ChainDpWorkDataFree(ChainDpWorkData* data)
{
    kv_destroy(data->f);
    kv_destroy(data->p);
    kv_destroy(data->t);
    kv_destroy(data->v);
    kv_destroy(data->u);
    kv_destroy(data->seeds);
    free(data); 
    return NULL;   
}

void
ChainDpWorkDataSetup(ChainDpWorkData* data, int n)
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
scoring_seeds(ChainDpWorkData* data)
{
    const int max_dist_ref = data->max_dist_ref;
    const int max_dist_qry = data->max_dist_qry;
    const int band_width = data->max_band_width;
    const int max_skip = data->max_skip;
    const ChainDpSeed* seeds = kv_data(data->seeds);
    const int n = kv_size(data->seeds);
    int sum_cov = 0;
    for (int i = 0; i < n; ++i) sum_cov += seeds[i].match_size;
    const int avg_cov = sum_cov / n;
    int st = 0;
    ChainDpWorkDataSetup(data, n);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);

    //for (int i = 0; i < n; ++i) fprintf(stderr, "%d:\t%d\t%d\t%d\n", i, seeds[i].qoff, seeds[i].soff, seeds[i].match_size);

    // fill the score and backtrack arrays
    for (int i = 0; i < n; ++i) {
        idx ri = seeds[i].soff;
        idx qi = seeds[i].qoff;
        int max_j = -1;
        int cov = seeds[i].match_size;
        int max_f = cov, n_skip = 0, min_d;
        while (st < i && ri > seeds[st].soff + max_dist_ref) ++st;
        for (int j = i - 1; j >= st; --j) {
            if (seeds[j].qoff >= qi || seeds[j].soff >= ri) continue;
            idx dr = ri - seeds[j].soff;
            idx dq = qi - seeds[j].qoff;
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
        //HBN_LOG("i = %d, st = %d, score = %d, p[i] = %d", i, st, f[i], p[i]);
    }
}

int chaining_seeds(ChainDpWorkData* data,
        int* best_seed_index,
        int* best_seed_score,
        vec_chain_seed* chain_seeds)
{
    const ChainDpSeed* seeds = kv_data(data->seeds);
    const int n = kv_size(data->seeds);
    if (n == 0) return 0;
    //HBN_LOG("min_score = %d, min_cnt = %d", data->min_score, data->min_cnt);
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    scoring_seeds(data);
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
    int find_best_seed = 0;
    for (int i = n_v = k = 0; i < n_u; ++i) { // start from the highest score
        int n_v0 = n_v, k0 = k;
        int j = u[i].second;
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
                find_best_seed = 1;
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
                find_best_seed = 1;
            }
        }
        //HBN_LOG("find_best_seed = %d, n_v = %d, n_v0 = %d, min_cnt = %d, j = %d, f[j] = %d, score = %d",
        //    find_best_seed, n_v, n_v0, min_cnt, j, f[j], u[i].first);

        if (find_best_seed) {
            if (best_seed_index) *best_seed_index = v[n_v0];
            if (best_seed_score) *best_seed_score = u[k-1].first;
            if (chain_seeds) {
                kv_clear(*chain_seeds);
                
                for (int x = n_v; x > n_v0; --x) {
                    int y = v[x-1];
                    
                    kv_push(ChainDpSeed, *chain_seeds, seeds[y]);
                }                
            }
            return 1;
        }

        if (k0 == k) n_v = n_v0;
    }    
    return 0;
}

int chaining_find_candidates(ChainDpWorkData* data,
        vec_can* can_list,
        vec_chain_seed* chain_seed_list,
        vec_int_pair* chain_seed_info_list)
{
    kv_clear(*can_list);
    kv_clear(*chain_seed_list);
    kv_clear(*chain_seed_info_list);

    const ChainDpSeed* seeds = kv_data(data->seeds);
    const int n = kv_size(data->seeds);
    if (n == 0) return 0;
    //HBN_LOG("min_score = %d, min_cnt = %d", data->min_score, data->min_cnt);
    const int min_cnt = data->min_cnt;
    const int min_score = data->min_score;
    scoring_seeds(data);
    int* f = kv_data(data->f);
    int* p = kv_data(data->p);
    int* t = kv_data(data->t);
    int* v = kv_data(data->v);
    IntPair* u = kv_data(data->u);
    GappedCandidate can;
    IntPair ip;

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
        find_can = 0;
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
            }
        } else if (u[i].first - f[j] >= min_score) {
            if (n_v - n_v0 >= min_cnt) {
                u[k].first = u[i].first - f[j];
                u[k].second = n_v - n_v0;
                ++k;
                find_can = 1;
            }
        }
        //HBN_LOG("find_best_seed = %d, n_v = %d, n_v0 = %d, min_cnt = %d, j = %d, f[j] = %d, score = %d",
        //    find_best_seed, n_v, n_v0, min_cnt, j, f[j], u[i].first);

        if (find_can) {
            ip.first = kv_size(*chain_seed_list);
            ip.second = n_v - n_v0;
            kv_push(IntPair, *chain_seed_info_list, ip);

            can.qend = seeds[v[n_v0]].qoff + seeds[v[n_v0]].match_size;
            can.send = seeds[v[n_v0]].soff + seeds[v[n_v0]].match_size;
            can.qbeg = seeds[v[n_v-1]].qoff;
            can.sbeg = seeds[v[n_v-1]].soff;
            can.score = u[k-1].first;
            can.qdir = FWD;
            can.qid = 0;
            can.qdir = FWD;
            can.qsize = 0;
            can.sid = 0;
            can.sdir = FWD;
            can.ssize = 0;
            int max_size = 0;
            //fprintf(stderr, "========================\n");
            for (int x = n_v; x > n_v0; --x) {
                int y = v[x-1];
                //fprintf(stderr, "%d\t[%d, %d, %d]\n", x, seeds[y].qoff, seeds[y].soff, seeds[y].match_size);
                kv_push(ChainDpSeed, *chain_seed_list, seeds[y]);
                if (seeds[y].match_size > max_size) {
                    max_size = seeds[y].match_size;
                    can.qoff = seeds[y].qoff + max_size / 2;
                    can.soff = seeds[y].soff + max_size / 2;
                }
            }
            //DUMP_GAPPED_CANDIDATE_STD(can);
            kv_push(GappedCandidate, *can_list, can);
        }

        if (k0 == k) n_v = n_v0;
    }    
    return 0;
}