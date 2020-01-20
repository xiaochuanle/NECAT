#ifndef CEDLIB_ALIGNER_H
#define CEDLIB_ALIGNER_H

#include "cedlib_defs.h"

#include "../klib/kstring.h"

typedef struct {
	void* km;
	EdlibData*	edlibData1;
	EdlibData*	edlibData2;
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
	double ident_perc;
} CEdlibAlignData;

CEdlibAlignData*
new_CEdlibAlignData();

CEdlibAlignData*
free_CEdlibAlignData(CEdlibAlignData* data);

BOOL
cedlib_align(const char* query,
			 const int query_from,
			 const int query_to,
			 const char* target,
			 const int target_from,
			 const int target_to,
			 CEdlibAlignData* align_data,
			 const int tolerance,
			 const int min_align_size);

#endif // CEDLIB_ALIGNER_H
