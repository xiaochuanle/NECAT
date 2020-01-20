#include "split_ctgs.h"

#define CtgSeqSize 50000
#define MaxCtgLeftSize 10000
#define MaxCtgSeqSize 60000

void
make_ctg_seq_hdr(int ctg_id, int seq_size, idx ctg_offset, idx ctg_size, kstring_t* hdr)
{
	kstr_clear(*hdr);
	ksprintf(hdr, ">%d_%d_%" PRIdx "_%" PRIdx, ctg_id, seq_size, ctg_offset, ctg_size);
}

void
extract_ctg_seq_hdr(const char* hdr, int* ctg_id, int* seq_size, idx* ctg_offset, idx* ctg_size)
{
	sscanf(hdr, "%d_%d_%" PRIdx "_%" PRIdx, ctg_id, seq_size, ctg_offset, ctg_size);
}

void
split_contig(int ctg_id, kstring_t* ctg, FILE* out)
{
	new_kstring(hdr);
	size_t ctg_size = kstr_size(*ctg);
	size_t left = ctg_size;
	size_t from = 0;
	while (left) {
		size_t to = from + CtgSeqSize;
		to = OC_MIN(to, ctg_size);
		left = ctg_size - to;
		if (left <= MaxCtgLeftSize) {
			to = ctg_size;
			left = 0;
		}
		make_ctg_seq_hdr(ctg_id, to - from, from, ctg_size, &hdr);
		fprintf(out, "%s\n", kstr_str(hdr));
		for (size_t i = from; i < to; ++i) fprintf(out, "%c", kstr_A(*ctg, i));
		fprintf(out, "\n");
		
		from = to;
	}
	free_kstring(hdr);
}
