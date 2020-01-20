#ifndef ALIGN_TAGS_H
#define ALIGN_TAGS_H

#include "../klib/kstring.h"
#include "../klib/kvec.h"
#include "../common/ontcns_defs.h"

typedef struct {
	double weight;
	int t_pos;
	int p_t_pos;
	u8 delta;
	u8 p_delta;
	char q_base;
	char p_q_base;
} AlignTag;

typedef kvec_t(AlignTag) vec_align_tag;

void
ks_introsort_AlignTag(size_t n, AlignTag* tags);

BOOL
get_cns_tags(const char* qaln,
			 const char* taln,
			 const size_t aln_size,
			 kstring_t* target,
			 int toff,
			 int tend,
			 double weight,
			 vec_align_tag* tags);

#endif // ALIGN_TAGS_H
