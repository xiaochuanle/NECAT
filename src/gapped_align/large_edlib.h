#ifndef LARGE_EDLIB_H
#define LARGE_EDLIB_H

#include "edlib_ex_aux.h"

typedef struct {
	EdlibAlignData*		edlib;
	int					qoff;
	int					qend;
	int					toff;
	int					tend;
	double				ident_perc;
	kstring_t			query_align;
	kstring_t			target_align;
	kstring_t			fqaln;
	kstring_t			rqaln;
	kstring_t			ftaln;
	kstring_t			rtaln;
	kstring_t			qfrag;
	kstring_t			tfrag;
	double				error;
} LargeEdlibData;

LargeEdlibData*
new_LargeEdlibData(const double _error, int _maxNumBlock, int _maxSeqSize);

void 
free_LargeEdlibData(LargeEdlibData* data);

BOOL
ledlib_go(const char* query,
		  const int query_start,
		  const int query_size,
		  const char* target,
		  const int target_start,
		  const int target_size,
		  LargeEdlibData* oca_data,
		  const int min_align_size);

#endif // LARGE_EDLIB_H
