#ifndef CNS_SEQ_H
#define CNS_SEQ_H

#include "../common/ontcns_defs.h"
#include "../klib/kstring.h"

typedef struct {
	const char* hdr;
	int id;
	int left, right;
	int org_seq_size;
	int num_can;
	int num_ovlps;
	double ident_cutoff;
	kstring_t cns_seq;
} CnsSeq;

extern const char* ontsa_hdr_prefix;

BOOL is_ontsa_hdr(const char* hdr);
int extract_ontsa_id(const char* hdr);
   
#define DUMP_CNS_SEQ(output_func, out, cns) \
	do { \
		if (is_ontsa_hdr((cns).hdr)) { \
			output_func(out, ">%s_", (cns).hdr); \
		} else { \
			output_func(out, ">%s_%d_", ontsa_hdr_prefix, (cns).id); \
		} \
		output_func(out, "(%d_%d_%d_%d_%d_%d_%lf)", \
					(cns).left, \
					(cns).right, \
					(int)kstr_size((cns).cns_seq), \
					(cns).org_seq_size, \
					(cns).num_can, \
					(cns).num_ovlps, \
					(cns).ident_cutoff); \
		output_func(out, "\n"); \
		DUMP_KSTRING(output_func, out, (cns).cns_seq); \
		output_func(out, "\n"); \
	} while(0)
   
#define clear_CnsSeq(cns) (kstr_clear((cns).cns_seq))
#define new_CnsSeq(cns) CnsSeq cns; kstr_init((cns).cns_seq)
#define free_CnsSeq(cns) free_kstring((cns).cns_seq)
#define init_CnsSeq(cns) kstr_init((cns).cns_seq)

#endif // CNS_SEQ_H
