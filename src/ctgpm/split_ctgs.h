#ifndef CTG_PM_H
#define CTG_PM_H

#include "../common/ontcns_aux.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef struct {
	int ctg_id;
	int seq_size;
	idx ctg_offset;
	idx ctg_size;
} CtgSeqHdrInfo;

typedef kvec_t(CtgSeqHdrInfo) vec_hdr_info;

void
make_ctg_seq_hdr(int ctg_id, int seq_size, idx ctg_offset, idx ctg_size, kstring_t* hdr);

void
extract_ctg_seq_hdr(const char* hdr, int* ctg_id, int* seq_size, idx* ctg_offset, idx* ctg_size);

void
split_contig(int ctg_id, kstring_t* ctg, FILE* out);

#endif // CTG_PM_H
