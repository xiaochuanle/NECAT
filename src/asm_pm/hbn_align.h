#ifndef HBN_ALIGN_H
#define HBN_ALIGN_H

#include "blockwise_edlib.h"
#include "daligner.h"

typedef struct {
    BlockwiseEdlibData* edlib_data;
    OcDalignData* dalign_data;
    vec_u8 qstr;
    vec_u8 tstr;
} HbnMapData;

#define hbn_map_query_align(hbnm)   ks_s((hbnm).edlib_data->query_align)
#define hbn_map_query_start(hbnm)   ((hbnm).edlib_data->qoff)
#define hbn_map_query_end(hbnm)     ((hbnm).edlib_data->qend)
#define hbn_map_target_align(hbnm)  ks_s(hbnm).edlib_data->target_align)
#define hbn_map_target_start(hbnm)  ((hbnm).edlib_data->toff)
#define hbn_map_target_end(hbnm)    ((hbnm).edlib_data->tend)
#define hbn_map_ident_perc(hbnm)    ((hbnm).edlib_data->ident_perc)

#define hbn_map_info(hbnm, qoff, qend, toff, tend, ident_perc, qalign, talign) do { \
    (qoff) = hbn_map_query_start(hbnm); \
    (qend) = hbn_map_query_end(hbnm); \
    (toff) = hbn_map_target_start(hbnm); \
    (tend) = hbn_map_target_end(hbnm); \
    (ident_perc) = hbn_map_ident_perc(hbnm); \
    if (qalign) qalign = hbn_map_query_align(hbnm); \
    if (talign) talign = hbn_map_target_align(hbnm); \
} while(0)

HbnMapData*
hbn_map_data_new();

HbnMapData*
hbn_map_data_free(HbnMapData* data);

BOOL 
hbn_cns_extend(HbnMapData* hbn_map_data,
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

BOOL 
hbn_map_extend(HbnMapData* hbn_map_data,
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

#endif // HBN_ALIGN_H