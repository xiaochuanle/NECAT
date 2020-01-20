#ifndef OVERLAPS_POOL_H
#define OVERLAPS_POOL_H

#include "../common/m4_record.h"
#include "../gapped_align/oc_aligner.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef struct {
	int qid, qdir;
	int qoff, qend, qsize;
	int toff, tend;
	double ident_perc;
	size_t can_id;
	size_t align_size;
	size_t qalign_start;
	size_t talign_start;
} OverlapIndex;

typedef kvec_t(OverlapIndex) vec_overlap_index;

typedef struct {
	vec_overlap_index oiv;
	kstring_t	align_string;
} OverlapsPool;

#define oc_op_num_align(op)	(kv_size((op).oiv))
#define oc_op_init(op) 		do { kv_init((op).oiv); kstr_init((op).align_string); } while(0)
#define oc_op_clear(op) 	do { kv_clear((op).oiv); kstr_clear((op).align_string); } while(0)
#define oc_op_oip(op, i) 	(kv_data((op).oiv) + (i))
#define oc_op_qstr(op, i)	(kstr_str((op).align_string) + kv_A((op).oiv, i).qalign_start)
#define oc_op_tstr(op, i)  	(kstr_str((op).align_string) + kv_A((op).oiv, i).talign_start)

OverlapsPool*
new_OverlapsPool();

OverlapsPool*
free_OverlapsPool(OverlapsPool* op);

/*
void
op_add_align(OcAlignData* align_data,
			 const int qid,
			 const int qdir,
			 const int qsize,
			 size_t can_id,
			 OverlapsPool* op);
			 */

void
op_add_align(int qid,
			 int qdir,
			 int qsize,
			 size_t can_id,
			 int qoff,
			 int qend,
			 int toff,
			 int tend,
			 double ident_perc,
			 kstring_t* qaln,
			 kstring_t* taln,
			 OverlapsPool* op);

#endif // OVERLAPS_POOL_H
