#ifndef LARGEST_COVER_RANGE_H
#define LARGEST_COVER_RANGE_H

#include "../common/m4_record.h"
#include "range_list.h"
#include "../common/packed_db.h"

void
get_largest_cover_range_for_one_partition(const char* m4_path, 
		const int pid, 
		PackedDB* reads,
		int* read_id,
		FILE* out, 
		pthread_mutex_t* out_lock,
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		const int num_threads);

#endif // LARGEST_COVER_RANGE_H
