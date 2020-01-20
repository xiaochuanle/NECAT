#ifndef ALIGN_DEFS_H
#define ALIGN_DEFS_H

#include "../common/ontcns_aux.h"
#include "../common/oc_assert.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

#define hbn_min OC_MIN
#define hbn_max OC_MAX
#define HBN_ERR OC_ERROR
#define ks_init kstr_init
#define ks_destroy free_kstring
#define hbn_assert oc_assert
#define ks_clear kstr_clear
#define ks_size kstr_size
#define ks_s kstr_str
#define ks_A kstr_A

typedef kvec_t(u8) vec_u8;
typedef vec_intpair vec_int_pair;

static inline int ks_set_size(kstring_t* s, size_t size)
{
	ks_reserve(s, size);
	s->l = size;
	return 0;
}

#endif // ALIGN_DEFS_H