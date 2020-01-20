#include "chain.h"
#include "word_finder_aux.h"
#include "../common/ontcns_aux.h"

#include <limits.h>

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

int
chain_dp_backward(ChainDpData* data, ChainSeed* lm_seed)
{
  const int kmer_size = data->kmer_size;
  const int max_dist = data->max_dist;
  const int bw = data->bw;
  const int max_skip = data->max_skip;
  const int min_cnt = data->min_cnt;
  const int min_sc = data->min_sc;

  vec_intpair* u = &data->u;
  ChainSeed* seeds = kv_data(*data->chain_seeds);
  const int n_seeds = kv_size(*data->chain_seeds);

  kv_resize(int, data->f, n_seeds); kv_zero(int, data->f);  int* f = kv_data(data->f);
  kv_resize(int, data->p, n_seeds); kv_fill(data->p, -1);   int* p = kv_data(data->p);
  kv_resize(int, data->t, n_seeds); kv_zero(int, data->t);  int* t = kv_data(data->t);
  kv_resize(int, data->v, n_seeds); kv_zero(int, data->v);  int* v = kv_data(data->v);

  int i , j, k;
  int st = 0; 

  for (i = 0; i < n_seeds; ++i) {
    idx roff = seeds[i].soff;
    idx qoff = seeds[i].qoff;
    int max_f = kmer_size, max_j = -1, n_skip = 0, min_d, max_f_past = INT_MIN;
    while (st < i && roff - seeds[st].soff > max_dist) ++st;
    for (j = i - 1; j >= st; --j) {
      if (seeds[j].soff >= roff || seeds[j].qoff >= qoff || qoff - seeds[j].qoff > max_dist) continue;
      int dr = roff - seeds[j].soff;
      int dq = qoff - seeds[j].qoff;
      int dd = (dr > dq) ? (dr - dq) : (dq - dr);
      if (dd > bw) continue;
      max_f_past = OC_MAX(max_f_past, f[j]);
      int min_d = OC_MIN(dq, dr);
      int sc = OC_MIN(min_d, kmer_size);
      sc -= (int)(dd * 0.01 * kmer_size) + (ilog2_32(dd) >> 1);
      sc += f[j];
      if (sc > max_f) {
        max_f = sc;
        max_j = j;
        if (n_skip) --n_skip;
      } else if (t[j] == i) {
        if (++n_skip > max_skip) break;
      }
      if (p[j] >= 0) t[p[j]] = i;
    }

    f[i] = max_f; p[i] = max_j; v[i] = max_f_past;
  }


  int score = 0;
  int next_i = n_seeds - 1;
  while (1) {
    if (p[next_i] == -1) break;
    next_i = p[next_i];
    *lm_seed = seeds[next_i];
    ++score;
  }
  return score;
}

int
chain_dp_forward(ChainDpData* data, ChainSeed* rm_seed)
{
  const int kmer_size = data->kmer_size;
  const int max_dist = data->max_dist;
  const int bw = data->bw;
  const int max_skip = data->max_skip;
  const int min_cnt = data->min_cnt;
  const int min_sc = data->min_sc;

  vec_intpair* u = &data->u;
  ChainSeed* seeds = kv_data(*data->chain_seeds);
  const int n_seeds = kv_size(*data->chain_seeds);

  kv_resize(int, data->f, n_seeds); kv_zero(int, data->f);  int* f = kv_data(data->f);
  kv_resize(int, data->p, n_seeds); kv_fill(data->p, -1);   int* p = kv_data(data->p);
  kv_resize(int, data->t, n_seeds); kv_zero(int, data->t);  int* t = kv_data(data->t);
  kv_resize(int, data->v, n_seeds); kv_zero(int, data->v);  int* v = kv_data(data->v);

  int i, j, k;
  int st = n_seeds - 1;
  
  for (i = n_seeds - 1; i >= 0; --i) {
    idx roff = seeds[i].soff;
    idx qoff = seeds[i].qoff;
    int max_f = kmer_size, max_j = -1, n_skip = 0, min_d, max_f_past = INT_MIN;
    while (st > i && seeds[st].soff - roff > max_dist) --st;
    for (j = i + 1; j <= st; ++j) {
      if (seeds[j].soff <= roff || seeds[j].qoff <= qoff || seeds[j].qoff - qoff > max_dist) continue;
      int dr = seeds[j].soff - roff;
      int dq = seeds[j].qoff - qoff;
      int dd = dr > dq ? (dr - dq) : (dq - dr);
      if (dd > bw) continue;
      max_f_past = OC_MAX(max_f_past, f[j]);
      int min_d = OC_MIN(dq, dr);
      int sc = OC_MIN(kmer_size, min_d);
      sc -= (int)(dd * 0.01 * kmer_size) + (ilog2_32(dd) >> 1);
      sc += f[j];
      if (sc > max_f) {
        max_f = sc;
        max_j = j;
        if (n_skip) --n_skip;
      } else if (t[j] == i) {
        if (++n_skip > max_skip) break;
      }
      if (p[j] >= 0) t[p[j]] = i;
    }
    f[i] = max_f, p[i] = max_j, v[i] = max_f_past;
  }

  int score = 0;
  int next_i = 0;
  while(1) {
    if (p[next_i] == -1) break;
    next_i = p[next_i];
    *rm_seed = seeds[next_i];
    ++score;
  }

  return score;
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
	chain_data->max_dist = 200;
	chain_data->bw = 500;
	chain_data->max_skip = 25;
	chain_data->min_cnt = block_score_cutoff;
	chain_data->min_sc = 40;
	
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
