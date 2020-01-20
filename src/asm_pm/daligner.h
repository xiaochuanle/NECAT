#ifndef OC_DALIGNER_H
#define OC_DALIGNER_H

#include "align.h"

#include "align_defs.h"
#include "../common/gapped_candidate.h"

typedef struct {
  Align_Spec*   align_spec;
  Work_Data*    work;
  float         basis_freq[4];
  Path          path;
  Alignment     align;
  kstring_t     aseq;
  kstring_t     bseq;
  double        sequencing_error;
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

int
ocda_go(OcDalignData* ocda_data,
	const u8* query,
	const u8* target,
	GappedCandidate* can,
	const int min_align_size,
	const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align,
	kstring_t* target_align);

#endif // OC_DALIGNER_H
