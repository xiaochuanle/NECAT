#ifndef EDLIB_WRAPPER_H
#define EDLIB_WRAPPER_H

#include "edlib.h"
#include "../klib/kstring.h"
#include "../common/ontcns_defs.h"

typedef struct {
	kstring_t query;
	kstring_t target;
	kstring_t query_align;
	kstring_t target_align;
	kstring_t qAln;
	kstring_t tAln;
	int qoff;
	int qend;
	int toff;
	int tend;
	int dist;
	double ident_perc;
	double error;
	void* km;
} FullEdlibAlignData;

FullEdlibAlignData*
new_FullEdlibAlignData(double error);

void
free_FullEdlibAlignData(FullEdlibAlignData* data);

int
edlib_go(const char* query,
		 const int query_from,
		 const int query_to,
		 const char* target,
		 const int target_from,
		 const int target_to,
		 FullEdlibAlignData* align_data,
		 const int tolerance,
		 const int min_align_size,
		 const BOOL find_path,
		 const int kMatchSize);

#endif // EDLIB_WRAPPER_H
