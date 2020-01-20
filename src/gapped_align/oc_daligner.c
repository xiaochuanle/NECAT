#include "oc_daligner.h"

OcDalignData*
new_OcDalignData(double error)
{
	OcDalignData* data = (OcDalignData*)malloc( sizeof(OcDalignData) );
	data->basis_freq[0] = .25;
	data->basis_freq[1] = .25;
	data->basis_freq[2] = .25;
	data->basis_freq[3] = .25;
	data->align_spec = New_Align_Spec(1.0 - error, 100, data->basis_freq);
	data->work = New_Work_Data();
	kstr_init(data->aseq);
	kstr_init(data->bseq);
	kstr_init(data->query_align);
	kstr_init(data->target_align);
	return data;
}

OcDalignData*
free_OcDalignData(OcDalignData* data)
{
	Free_Align_Spec(data->align_spec);
	Free_Work_Data(data->work);
	free_kstring(data->aseq);
	free_kstring(data->bseq);
	free_kstring(data->query_align);
	free_kstring(data->target_align);
	free(data);
	return 0;
}

BOOL
ocda_go(const char* query,
    const int query_start,
    const int query_size,
    const char* target,
    const int target_start,
    const int target_size,
    OcDalignData* ocda,
    const int min_align_size)
{
  ks_resize(&ocda->aseq, query_size + 2);
  ks_resize(&ocda->bseq, target_size + 2);
  char* aseq = kstr_data(ocda->aseq);
  char* bseq = kstr_data(ocda->bseq);
  *aseq = 4;
  memcpy(aseq + 1, query, query_size);
  aseq[query_size + 1] = 4;
  *bseq = 4;
  memcpy(bseq + 1, target, target_size);
  bseq[target_size + 1] = 4;

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
  return TRUE;
}
