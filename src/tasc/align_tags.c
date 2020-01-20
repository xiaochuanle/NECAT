#include "align_tags.h"

#include "../common/oc_assert.h"
#include "../klib/ksort.h"

#define AlignTag_LT(a, b) ( \
	((a).t_pos < (b).t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta < (b).delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base < (b).q_base)\
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos < (b).p_t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta < (b).p_delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base < (b).p_q_base) \
	)

KSORT_INIT(AlignTag, AlignTag, AlignTag_LT)

BOOL
get_cns_tags(const char* qaln,
			 const char* taln,
			 const size_t aln_size,
			 kstring_t* target,
			 int toff,
			 int tend,
			 double weight,
			 vec_align_tag* tags)
{
	AlignTag tag;
    tag.weight = weight;
    int jj = 0;
    int j = toff - 1;
    int p_j = -1;
    int p_jj = 0;
    char p_q_base = GAP_CHAR;
	
	for (size_t i = 0; i != aln_size; ++i) {
		if (qaln[i] != GAP_CHAR) ++jj;
		if (taln[i] != GAP_CHAR) jj = 0;
		if (jj >= U8_MAX || p_jj >= U8_MAX) return FALSE;
	}

    jj = 0;
    p_jj = 0;
    for (size_t i = 0; i != aln_size; ++i) {
        if (qaln[i] != GAP_CHAR) ++jj;
        if (taln[i] != GAP_CHAR) {
            ++j;
            jj = 0;
        }

        tag.t_pos = j;
		oc_assert(j < tend);
        tag.p_t_pos = p_j;
        tag.delta = (u8)jj;
        tag.p_delta = (u8)p_jj;
        tag.q_base = qaln[i];
        tag.p_q_base = p_q_base;

        p_j = j;
        p_jj = jj;
        p_q_base = qaln[i];

        kv_push(AlignTag, *tags, tag);
    }

    return TRUE;
}
