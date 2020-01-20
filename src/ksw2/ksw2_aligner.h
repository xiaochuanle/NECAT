#ifndef KSW2_ALIGNER_H
#define KSW2_ALIGNER_H

#include "kalloc.h"
#include "ksw2.h"

#include "../common/ontcns_defs.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef kvec_t(uint8_t) vec_u8;

typedef struct {
	int a;
	int b;
	int q;
	int e;
	int q2;
	int e2;
	int bw;
	int zdrop;
} Ksw2Params;

void
init_Ksw2Params(Ksw2Params* p);

typedef struct {
	Ksw2Params 			align_params;
	void* 				km;
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
	vec_u8				qfrag;
	vec_u8				tfrag;
} Ksw2AlignData;

Ksw2AlignData*
new_Ksw2AlignData();

void
free_Ksw2AlignData(Ksw2AlignData* data);

BOOL
ksw2_go(const char* query,
		const int query_start,
		const int query_size,
		const char* target,
		const int target_start,
		const int target_size,
		Ksw2AlignData* ksw_data,
		const int min_align_size);

#endif // KSW2_ALIGNER_H
