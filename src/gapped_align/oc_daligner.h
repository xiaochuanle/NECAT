#ifndef OC_DALIGNER_H
#define OC_DALIGNER_H

#include "align.h"
#include "../common/ontcns_defs.h"
#include "../klib/kstring.h"

typedef struct {
  Align_Spec*   align_spec;
  Work_Data*    work;
  float         basis_freq[4];
  Path          path;
  Alignment     align;
  kstring_t     aseq;
  kstring_t     bseq;
  double        ident_perc;
  kstring_t		query_align;
  kstring_t		target_align;
} OcDalignData;

OcDalignData*
new_OcDalignData(double error);

OcDalignData*
free_OcDalignData(OcDalignData* data);

#define ocda_ident_perc(ocda)     ((ocda).ident_perc)
#define ocda_query_start(ocda)    ((ocda).path.abpos)
#define ocda_query_end(ocda)      ((ocda).path.aepos)
#define ocda_target_start(ocda)   ((ocda).path.bbpos)
#define ocda_target_end(ocda)     ((ocda).path.bepos)
#define ocda_distance(ocda)		  ((ocda).path.diffs)

BOOL
ocda_go(const char* query,
    const int query_start,
    const int query_size,
    const char* target,
    const int target_start,
    const int target_size,
    OcDalignData* ocda,
    const int min_align_size);

#endif // OC_DALIGNER_H
