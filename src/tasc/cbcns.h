#ifndef CBCNS_H
#define CBCNS_H

#include "cns_aux.h"
#include "../common/cns_seq.h"
#include "../common/record_writer.h"

typedef struct {
	vec_align_tag			tags;
	vec_backbone_item		backbone;
	vec_int					coverage;
	LinkInfoAllocator*		li_alloc;
	DeltaCovInfoAllocator*	dci_alloc;
	int 					template_size;
} CbCnsData;

CbCnsData*
new_CbCnsData();

CbCnsData*
free_CbCnsData(CbCnsData* data);

void
clear_CbCnsData(CbCnsData* data);

void
add_one_align(CbCnsData* cns_data,
			  const char* qaln,
			  const char* taln,
			  const size_t aln_size,
			  kstring_t* target,
			  int toff,
			  int tend,
			  double weight);

void
consensus_broken(CbCnsData* cns_data,
				 const int min_cov,
				 const int min_size,
				 const int template_id,
				 const int template_size,
				 vec_intpair* cns_intvs,
				 CnsSeq* cns_seq,
				 RecordWriter* out,
				 const int check_chimeric_read,
				 vec_intpair* cov_ranges);

int
consensus_unbroken(CbCnsData* cns_data,
				   const int min_cov,
				   const int min_size,
				   const char* raw_read,
				   const int template_size,
				   kstring_t* cns_seq);

#endif // CBCNS_H
