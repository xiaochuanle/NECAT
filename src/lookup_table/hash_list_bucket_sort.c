#include "hash_list_bucket_sort.h"
#include "../common/ontcns_aux.h"

#include <assert.h>
#include <pthread.h>

#define BS_B			((u64)8)
#define BS_PASS			((u64)4)
#define BS_MASK			((u64)255)
#define BS_BucketSize	((u64)256)

typedef struct
{
	void*			src;
	void* 			trg;
	u64				item_count;
	u64 			pass;
	int				thread_id;
	int				num_threads;
	u64**			buckets;
	u64***			next_buckets;
	ValueExtractor 	offset_extractor;
	ValueExtractor 	hash_extractor;
	SetListValue	set_list_value;
} BucketSortThreadData;

void
init_buckets(void* src,
			 ValueExtractor hash_extractor,
			 const u64 item_count,
			 const int num_threads,
			 u64** buckets)
{
	TIMING_START(__func__);
	
	const u64 part = (item_count + num_threads - 1) / num_threads;
	for (int i = 0; i < num_threads; ++i) {
		u64 from = i * part;
		u64 to = from + part; 
		if (to > item_count) to = item_count;
		for (u64 k = from; k < to; ++k) {
			u64 hash = hash_extractor(src, k);
			u64 bucket = hash & BS_MASK;
			++buckets[i][bucket];
		}
	}
	
	TIMING_END(__func__);
}

void
fill_bucket_from_next(u64** buckets, u64*** next_buckets, const int num_threads)
{
    for (int i = 0; i < num_threads; ++i) {
        for (u64 j = 0; j < BS_BucketSize; ++j) {
            u64 c = 0;
            for (int k = 0; k < num_threads; ++k) {
                c += next_buckets[k][i][j];
            }
            buckets[i][j] = c;
        }
    }
}

void
update_buckets(u64** buckets, u64** tmp_buckets, const u64 src_size, const int num_threads)
{
    for (int i = 0; i < num_threads; ++i) {
		memcpy(tmp_buckets[i], buckets[i], sizeof(u64) * BS_BucketSize);
    }

    for (int t = 0; t < num_threads; ++t) {
        for (u64 b = 0; b < BS_BucketSize; ++b) {
            u64 cnt = 0;
            for (int u = 0; u < num_threads; ++u) {
                for (u64 c = 0; c < b; ++c) {
                    cnt += tmp_buckets[u][c];
                }
            }
            for (int u = 0; u < t; ++u) {
                cnt += tmp_buckets[u][b];
            }
            buckets[t][b] = cnt;
        }
    }
}

void*
bucket_sort_thread(void* arg)
{
	BucketSortThreadData* bstd = (BucketSortThreadData*)(arg);
	const u64 shift = bstd->pass * BS_B;
	const u64 item_count = bstd->item_count;
	const int num_threads = bstd->num_threads;
	const int thread_id = bstd->thread_id;
	void* src = bstd->src;
	void* trg = bstd->trg;
	u64** buckets = bstd->buckets;
	u64*** next_buckets = bstd->next_buckets;
	ValueExtractor hash_extractor = bstd->hash_extractor;
	SetListValue set_list_value = bstd->set_list_value;
	const u64 part = (item_count + num_threads - 1) / num_threads;
	const u64 from = part * thread_id;
	u64 to = from + part; if (to > item_count) to = item_count;
	
	for (u64 i = from; i < to; ++i) {
		u64 c = hash_extractor(src, i);
		u64 b = c >> shift;
		u64 b_idx = b & BS_MASK;
		u64 x = buckets[thread_id][b_idx]++;
		assert(x < item_count);
		set_list_value(src, i, trg, x);
		++next_buckets[thread_id][x/part][(b>>BS_B)&BS_MASK];
	}
	
	return NULL;
}

void
validate_ordered_list(void* src,
					  const u64 item_count,
					  ValueExtractor hash_extractor)
{
	for (u64 i = 0; i != item_count - 1; ++i) {
		u64 hi = hash_extractor(src, i);
		u64 hi1 = hash_extractor(src, i + 1);
		if (hi > hi1) {
			OC_ERROR("i = %u, order invalid", i);
		}
	}
}

void
radix_sort(void* src, 
		   const u64 item_size,
		   const u64 item_count, 
		   const int num_threads,
		   ValueExtractor offset_extractor, 
		   ValueExtractor hash_extractor,
		   SetListValue set_list_value)
{
	/// buckets
	u64** buckets = (u64**)malloc(sizeof(u64*) * num_threads);
	u64** tmp_buckets = (u64**)malloc(sizeof(u64*) * num_threads);
	for (int i = 0; i < num_threads; ++i) {
		buckets[i] = (u64*)malloc(sizeof(u64) * BS_BucketSize);
		memset(buckets[i], 0, sizeof(u64) * BS_BucketSize);
		tmp_buckets[i] = (u64*)malloc(sizeof(u64) * BS_BucketSize);
	}
	/// next buckets
	u64*** next_buckets = (u64***)malloc(sizeof(u64**) * num_threads);
	for (int t = 0; t < num_threads; ++t) {
		next_buckets[t] = (u64**)malloc(sizeof(u64*) * num_threads);
		for (int i = 0; i < num_threads; ++i) {
			next_buckets[t][i] = (u64*)malloc(sizeof(u64) * BS_BucketSize);
		}
	}
	
	BucketSortThreadData* bstds = (BucketSortThreadData*)malloc(sizeof(BucketSortThreadData) * num_threads);
	pthread_t jobs[num_threads];
	void* trg = malloc(item_size * item_count);
	char job_name[1024];
	
	for (u64 p = 0; p != BS_PASS; ++p) {
		sprintf(job_name, "bucket sort pass %lu", p);
		TIMING_START(job_name);
		
		if (p == 0) {
			init_buckets(src, hash_extractor, item_count, num_threads, buckets);
		} else {
			fill_bucket_from_next(buckets, next_buckets, num_threads);
		}
		
		update_buckets(buckets, tmp_buckets, item_count, num_threads);
		for (int i = 0; i < num_threads; ++i)
			for (int j = 0; j < num_threads; ++j)
				memset(next_buckets[i][j], 0, sizeof(u64) * BS_BucketSize);
		
		for (int i = 0; i < num_threads; ++i) {
			bstds[i].src = 				src;
			bstds[i].trg = 				trg;
			bstds[i].item_count = 		item_count;
			bstds[i].pass = 			p;
			bstds[i].thread_id = 		i;
			bstds[i].num_threads = 		num_threads;
			bstds[i].buckets = 			buckets;
			bstds[i].next_buckets = 	next_buckets;
			bstds[i].offset_extractor =	offset_extractor;
			bstds[i].hash_extractor = 	hash_extractor;
			bstds[i].set_list_value = 	set_list_value;
			
			pthread_create(jobs + i, NULL, bucket_sort_thread, bstds + i);
		}
		
		for (int i = 0; i < num_threads; ++i) {
			pthread_join(jobs[i], NULL);
		}
		
		void* tmp = src;
		src = trg;
		trg = tmp;
		
		TIMING_END(job_name);
	}
	
	for (int i = 0; i < num_threads; ++i) {
		free(buckets[i]);
		free(tmp_buckets[i]);
	}
	free(buckets);
	free(tmp_buckets);
	
	for (int i = 0; i < num_threads; ++i) {
		for (int j = 0; j < num_threads; ++j) {
			free(next_buckets[i][j]);
		}
		free(next_buckets[i]);
	}
	free(next_buckets);
	
	free(bstds);
	free(trg);
	
	validate_ordered_list(src, item_count, hash_extractor);
}
