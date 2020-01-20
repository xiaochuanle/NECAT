#include "chain_dp.h"

#include "../common/ontcns_defs.h"
#include "../klib/ksort.h"

#include <assert.h>

#define IntPair_ChainDpGT(a,b) (((a).first > (b).first) || ((a).first == (b).first && (a).second < (b).second))
KSORT_INIT(IntPair_ChainDpGT, IntPair, IntPair_ChainDpGT)

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	register uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

static int
GappedCandidate_cdpScoreGT(GappedCandidate a, GappedCandidate b)
{
	if (a.score != b.score) return a.score > b.score;
	if (a.qoff != b.qoff) return a.qoff < b.qoff;
	if (a.soff != b.soff) return a.soff < b.soff;
	return 0;
}

KSORT_INIT(GappedCandidate_cdpScoreGT, GappedCandidate, GappedCandidate_cdpScoreGT)

void
chain_dp(ChainDpData* data, const int qid, const int qdir, const idx qsize, const int sid, const idx ssize)
{
	vec_intpair* u = &data->u;
	vec_can* lcanv = &data->lcanv;
	ChainSeed* seeds = kv_data(*data->chain_seeds);
	const size_t n_seeds = kv_size(*data->chain_seeds);
	const int kmer_size = data->kmer_size;
	const int max_dist = data->max_dist;
	const int bw = data->bw;
	const int max_skip = data->max_skip;
	const int min_cnt = data->min_cnt;
	const int min_sc = data->min_sc;
	
	kv_clear(*lcanv);
	kv_resize(int, data->f, n_seeds); kv_zero(int, data->f);  	int* f = kv_data(data->f);
	kv_resize(int, data->p, n_seeds); kv_fill(data->p, -1);		int* p = kv_data(data->p);
	kv_resize(int, data->t, n_seeds); kv_zero(int, data->t);	int* t = kv_data(data->t);
	kv_resize(int, data->v, n_seeds); kv_zero(int, data->v);	int* v = kv_data(data->v);
	int i, j, k, st = 0;

	for (i = 0; i < n_seeds; ++i) {
		idx ri = seeds[i].soff;
		idx qi = seeds[i].qoff;
		int max_j = -1;
		int max_f = kmer_size;
		int n_skip = 0;
		while (st < i && ri - seeds[st].soff > max_dist) ++st;
		for (j = i - 1; j >= st; --j) {
			if (ri <= seeds[j].soff || qi <= seeds[j].qoff || qi - seeds[j].qoff > max_dist) continue;
			idx dr = ri - seeds[j].soff;
			idx dq = qi - seeds[j].qoff;
			idx dd = dr > dq ? dr - dq : dq - dr;
			if (dd > bw) continue;
			idx min_d = OC_MIN(dq, dr);
			int sc = OC_MIN(min_d, kmer_size);
			int log_dd = dd ? ilog2_32(dd) : 0;
			sc -= (int)(dd * 0.01 * kmer_size) + (log_dd >> 1);
			sc += f[j];
			if (sc > max_f) {
				max_f = sc;
				max_j = j;
				if (n_skip > 0) --n_skip;
			} else if (t[j] == i) {
				if (++n_skip > max_skip) break;
			}
			if (p[j] >= 0) t[p[j]] = i;
		}
		f[i] = max_f;
		p[i] = max_j;
		v[i] = max_j >= 0 && v[max_j] > max_f ? v[max_j] : max_f;
	}

	memset(t, 0, sizeof(int) * n_seeds);
	for (i = 0; i < n_seeds; ++i)
		if (p[i] >= 0) t[p[i]] = 1;
	int n_u;
	for (i = n_u = 0; i < n_seeds; ++i) {
		if (t[i] == 0 && v[i] >= min_sc) ++n_u;
	}
	if (n_u == 0) return;
	
	kv_clear(*u);
	IntPair intp;
	for (i = n_u = 0; i < n_seeds; ++i) {
		if (t[i] == 0 && v[i] >= min_sc) {
			j = i;
			while (j >= 0 && f[j] < v[j]) j = p[j];
			if (j < 0) j = i;
			intp.first = f[j];
			intp.second = j;
			kv_push(IntPair, *u,intp);
			++n_u;
		}
	}
	ks_introsort(IntPair_ChainDpGT, n_u, kv_data(*u));
	
	GappedCandidate can;
	can.qid = qid;
	can.qdir = qdir;
	can.qsize = qsize;
	can.sid = sid;
	can.sdir = FWD;
	can.ssize = ssize;
	
	memset(t, 0, sizeof(int) * n_seeds);
	int n_v = 0;
	for (i = n_v = k = 0; i < n_u; ++i) {
		int n_v0 = n_v, k0 = k;
		j = kv_A(*u, i).second;
		can.qend = seeds[j].qoff + kmer_size;
		can.send = seeds[j].soff + kmer_size;
		can.qoff = can.qend;
		can.soff = can.send;
		int last_j = j;
		do {
			last_j = j;
			n_v++;
			t[j] = 1;
			j = p[j];
		} while (j >= 0 && t[j] == 0);
		if (j < 0) {
			if (n_v - n_v0 >= min_cnt) {
				can.qbeg = seeds[last_j].qoff;
				can.sbeg = seeds[last_j].soff;
				can.score = kv_A(*u, i).first;
				kv_push(GappedCandidate, *lcanv, can);
				++k;
			}
		}  else if (kv_A(*u, i).first - f[j] >= min_sc) {
			if (n_v - n_v0 >= min_cnt) {
				can.qbeg = seeds[last_j].qoff;
				can.sbeg = seeds[last_j].soff;
				can.score = kv_A(*u, i).first - f[j];
				kv_push(GappedCandidate, *lcanv, can);
				++k;
			}
		}
		if (k0 == k) n_v = n_v0;
	}
	if (kv_size(*lcanv) == 0) return;
	
	ks_introsort(GappedCandidate_cdpScoreGT, kv_size(*lcanv), kv_data(*lcanv));
}

ChainDpData*
new_ChainDpData(int kmer_size, int block_score_cutoff)
{
	ChainDpData* chain_data = (ChainDpData*)malloc(sizeof(ChainDpData));
	kv_init(chain_data->f);
	kv_init(chain_data->p);
	kv_init(chain_data->t);
	kv_init(chain_data->v);
	kv_init(chain_data->u);
	kv_init(chain_data->lcanv);
	chain_data->chain_seeds = NULL;
	
	chain_data->kmer_size = kmer_size;
	chain_data->max_dist = 5000;
	chain_data->bw = 500;
	chain_data->max_skip = 25;
	chain_data->min_cnt = block_score_cutoff;
	chain_data->min_sc = 30;
	
	return chain_data;
}

ChainDpData*
free_ChainDpData(ChainDpData* chain_data)
{
	kv_destroy(chain_data->f);
	kv_destroy(chain_data->p);
	kv_destroy(chain_data->t);
	kv_destroy(chain_data->v);
	kv_destroy(chain_data->u);
	kv_destroy(chain_data->lcanv);

	free(chain_data);
	return 0;
}
