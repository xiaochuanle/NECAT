#ifndef OC_ALIGNER_H
#define OC_ALIGNER_H

#include "edlib_ex_aux.h"

typedef struct {
	EdlibAlignData*		edlib;
	int					qoff;
	int					qend;
	int					toff;
	int					tend;
	double				ident_perc;
	kstring_t			query_align;
	kstring_t			target_align;
	kstring_t			fqaln;
	kstring_t			rqaln;
	kstring_t			ftaln;
	kstring_t			rtaln;
	kstring_t			qfrag;
	kstring_t			tfrag;
	char*				qabuf;
	char*				tabuf;
	double				error;
	int					qid;
	int					tid;
} OcAlignData;

#define oca_ident_perc(oca) 			((oca).ident_perc)
#define oca_query_start(oca) 			((oca).qoff)
#define oca_query_end(oca)				((oca).qend)
#define oca_target_start(oca)			((oca).toff)
#define oca_target_end(oca)				((oca).tend)
#define oca_query_mapped_string(oca)	(&((oca).query_align))
#define oca_target_mapped_string(oca)	(&((oca).target_align))

OcAlignData*
new_OcAlignData(double _error);

OcAlignData*
free_OcAlignData(OcAlignData* data);

#define ONC_TAIL_MATCH_LEN_LONG     4
#define ONC_TAIL_MATCH_LEN_SHORT    1

BOOL
onc_align(const char* query,
		  const int query_start,
		  const int query_size,
		  const char* target,
		  const int target_start,
		  const int target_size,
		  OcAlignData* oca_data,
		  const int block_size,
		  const int min_align_size,
          const int tail_match_len);

#endif // OC_ALIGNER_H
