#ifndef LARGEST_COVER_RANGE_H
#define LARGEST_COVER_RANGE_H

#include "../common/m4_record.h"
#include "range_list.h"
#include "../common/packed_db.h"

typedef struct {
	int left;
	int right;
	int size;
} LargestCoverRange;

typedef kvec_t(LargestCoverRange) vec_lcr;

#define invalidate_lcr(lcr) 	((lcr).left = -1)
#define lcr_is_valid(lcr) 		((lcr).left >= 0)

int lcr_is_complete(LargestCoverRange* lcr);
int lcr_is_trimmed(LargestCoverRange* lcr);


void load_lcrs(const char* path, vec_lcr* lcrv);

void
get_largest_cover_range_for_one_partition(const char* m4_path, 
		const int pid, 
		LargestCoverRange* lcrv,
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		const int num_threads);

#endif // LARGEST_COVER_RANGE_H
