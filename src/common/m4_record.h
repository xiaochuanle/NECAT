#ifndef M4_RECORD_H
#define M4_RECORD_H

#include "../common/ontcns_aux.h"
#include "../klib/kvec.h"
#include "../common/oc_assert.h"

#include <stdio.h>

typedef struct {
	int qid;
    int qdir;
    idx qoff;
    idx qend;
    idx qext;
    idx qsize;
    int sid;
    int sdir;
    idx soff;
    idx send;
    idx sext;
    idx ssize;
    double ident_perc;
    int vscore;
} M4Record;

typedef kvec_t(M4Record) vec_m4;

#define M4Record_SidLT(a, b) ((a).sid < (b).sid)
void ks_introsort_M4Record_SidLT(size_t n, M4Record* a);

#define M4Record_ScoreGT(a, b) ((a).vscore > (b).vscore)
void ks_introsort_M4Record_ScoreGT(size_t n, M4Record* a);

#define M4Record_IdentGT(a, b) ((a).ident_perc > (b).ident_perc)
void ks_introsort_M4Record_IdentGT(size_t n, M4Record* a);

#define DUMP_M4_RECORD(output_func, out, m) \
	output_func(out, "%d\t%d\t%.2f\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
				(m).qid, \
				(m).sid, \
				(m).ident_perc, \
				(m).vscore, \
				(m).qdir, \
				(m).qoff, \
				(m).qend, \
				(m).qext, \
				(m).qsize, \
				(m).sdir, \
				(m).soff, \
				(m).send, \
				(m).sext, \
				(m).ssize)
   
#define LOAD_M4_RECORD(input_func, in, m) \
	input_func(in, "%d%d%lf%d%d%lu%lu%lu%lu%d%lu%lu%lu%lu", \
			   &(m).qid, \
			   &(m).sid, \
			   &(m).ident_perc, \
			   &(m).vscore, \
			   &(m).qdir, \
			   &(m).qoff, \
			   &(m).qend, \
			   &(m).qext, \
			   &(m).qsize, \
			   &(m).sdir, \
			   &(m).soff, \
			   &(m).send, \
			   &(m).sext, \
			   &(m).ssize)

#define DUMP_ASM_M4(output_func, out, m4)  do { \
	output_func(out, 	"%d\t" \
						"%d\t" \
						"%.2f\t" \
						"%d\t" \
						"%d\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%d\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\n", \
						(m4).qid, \
						(m4).sid, \
						(m4).ident_perc, \
						(m4).vscore, \
						(m4).qdir, \
						(m4).qoff, \
						(m4).qend, \
						(m4).qsize, \
						(m4).sdir, \
						(m4).soff, \
						(m4).send, \
						(m4).ssize); \
} while(0)

#define DUMP_ASM_M4_HDR_ID(output_func, out, m4)  do { \
	output_func(out, 	"%s\t" \
						"%s\t" \
						"%.2f\t" \
						"%d\t" \
						"%d\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%d\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\t" \
						"%" PRIdx "\n", \
						qhdr, \
						shdr, \
						(m4).ident_perc, \
						(m4).vscore, \
						(m4).qdir, \
						(m4).qoff, \
						(m4).qend, \
						(m4).qsize, \
						(m4).sdir, \
						(m4).soff, \
						(m4).send, \
						(m4).ssize); \
} while(0)

#define LOAD_ASM_M4(input_func, stream, m4) \
	do { \
		SAFE_SCANF(input_func, stream, 12, \
				   "%d" \
				   "%d" \
				   "%lf" \
				   "%d" \
				   "%d" \
				   "%" PRIdx \
				   "%" PRIdx \
				   "%" PRIdx \
				   "%d" \
				   "%" PRIdx \
				   "%" PRIdx \
				   "%" PRIdx \
					, \
				   &(m4).qid, \
				   &(m4).sid, \
				   &(m4).ident_perc, \
				   &(m4).vscore, \
					&(m4).qdir, \
					&(m4).qoff, \
					&(m4).qend, \
					&(m4).qsize, \
				   &(m4).sdir, \
				   &(m4).soff, \
				   &(m4).send, \
				   &(m4).ssize); \
	} while(0)

typedef enum {
	eM4_Contained,
	eM4_Contains,
	eM4_Overlap,
	eM4_None
} EM4OverlapType;

EM4OverlapType
detect_m4_type(const M4Record* m, const idx unmapped_boundary);

#endif // M4_RECORD_H
