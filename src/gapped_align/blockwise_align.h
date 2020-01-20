#ifndef BLOCKWISE_ALIGN_H
#define BLOCKWISE_ALIGN_H

#include "../common/gapped_candidate.h"
#include "../common/packed_db.h"
#include "../klib/kstring.h"

#include "../gapped_align/oc_daligner.h"
#include "../edlib/edlib_wrapper.h"

BOOL
blockwise_align(PackedDB* reads,
				PackedDB* reference,
				OcDalignData* dalign_data,
				FullEdlibAlignData* falign_data,
				const idx block_size,
				GappedCandidate* _can,
				kstring_t* query,
				kstring_t* target,
				idx* qoff,
				idx* qend,
				idx* toff,
				idx* tend,
				double* ident_perc,
				const int min_align_size);

#endif // BLOCKWISE_ALIGN_H
