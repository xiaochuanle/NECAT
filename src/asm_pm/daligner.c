#include "daligner.h"

#include <assert.h>

OcDalignData*
new_OcDalignData(double error)
{
	OcDalignData* data = (OcDalignData*)malloc( sizeof(OcDalignData) );
	data->basis_freq[0] = .25;
	data->basis_freq[1] = .25;
	data->basis_freq[2] = .25;
	data->basis_freq[3] = .25;
	data->align_spec = New_Align_Spec(1.0 - error, 100, data->basis_freq);
  data->sequencing_error = error;
	data->work = New_Work_Data();
	ks_init(data->aseq);
	ks_init(data->bseq);
	ks_init(data->query_align);
	ks_init(data->target_align);
	return data;
}

OcDalignData*
free_OcDalignData(OcDalignData* data)
{
	if (data->align_spec) Free_Align_Spec(data->align_spec);
	Free_Work_Data(data->work);
	ks_destroy(data->aseq);
	ks_destroy(data->bseq);
	ks_destroy(data->query_align);
	ks_destroy(data->target_align);
	free(data);
	return 0;
}

int
ocda_go(OcDalignData* ocda,
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
	kstring_t* target_align)
{
  const int query_start = can->qoff;
  const int query_size = can->qsize;
  const int target_start = can->soff;
  const int target_size = can->ssize;
  ks_set_size(&ocda->aseq, query_size + 2);
  ks_set_size(&ocda->bseq, target_size + 2);
  char* aseq = ks_s(ocda->aseq);
  char* bseq = ks_s(ocda->bseq);
  *aseq = 4;
  memcpy(aseq + 1, query, query_size);
  aseq[query_size + 1] = 4;
  *bseq = 4;
  memcpy(bseq + 1, target, target_size);
  bseq[target_size + 1] = 4;

  for (int i = 0; i < query_size; ++i) {
    assert(query[i] >= 0 && query[i] < 4);
  }
  for (int i = 0; i < target_size; ++i) {
    assert(target[i] >= 0 && target[i] < 4);
  }

  Alignment* align = &ocda->align;
  align->flags = 0;
  align->path = &ocda->path;
  align->aseq = aseq + 1;
  align->bseq = bseq + 1;
  align->alen = query_size;
  align->blen = target_size;

  Local_Alignment(align, 
      ocda->work, 
      ocda->align_spec, 
      query_start - target_start, 
      query_start - target_start,
      query_start + target_start,
      -1,
      -1);

  int asize = ocda->path.aepos - ocda->path.abpos;
  int bsize = ocda->path.bepos - ocda->path.bbpos;
  BOOL r = (asize >= min_align_size) && (bsize >= min_align_size);
  if (!r) return FALSE;

  double error = 200.0 * ocda->path.diffs / (asize + bsize);
  ocda->ident_perc = 100.0 - error;
  r = ocda->ident_perc >= min_ident_perc;
  if (!r) return FALSE;
  *qbeg = ocda_query_start(*ocda);
  *qend = ocda_query_end(*ocda);
  *tbeg = ocda_target_start(*ocda);
  *tend = ocda_target_end(*ocda);
  *ident_perc = ocda_ident_perc(*ocda);
  return r;
}
