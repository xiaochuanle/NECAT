#ifndef CTGPM_H
#define CTGPM_H

#include "../common/gapped_candidate.h"
#include "../common/packed_db.h"
#include "../klib/kstring.h"

#include "../gapped_align/oc_daligner.h"
#include "../edlib/edlib_wrapper.h"

BOOL
ctgmp_align(PackedDB* pctg,
	OcDalignData* dalign_data,
	FullEdlibAlignData* falign_data,
	GappedCandidate* can,
	kstring_t* query,
	kstring_t* target,
	idx* qoff,
	idx* qend,
	idx* toff,
	idx* tend,
	double* ident_perc,
    const int min_align_size);

#endif // CTGPM_H
