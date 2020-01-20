#ifndef ERROR_ESTIMATE_H
#define ERROR_ESTIMATE_H

#include <math.h>

#include "../common/gapped_candidate.h"
#include "../common/m4_record.h"
#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"
#include "overlaps_pool.h"
#include "read_id_pool.h"
#include "../gapped_align/oc_daligner.h"
#include "../edlib/edlib_wrapper.h"

BOOL
get_good_overlaps(PackedGappedCandidate* candidates,
				  size_t can_sid,
				  size_t can_eid,
				  size_t* last_extended_can_id,
				  OcAlignData* align_data,
				  OcDalignData* dalign_data,
				  FullEdlibAlignData*	falign_data,
				  BOOL rescue_long_indels,
				  OverlapsPool* op,
				  kstring_t* query,
				  kstring_t* target,
				  PackedDB* reads,
				  CnsReads* cns_reads,
				  ReadIdPool* extended_reads,
				  int min_align_size,
				  double* ident_cutoff);

#endif // ERROR_ESTIMATE_H
