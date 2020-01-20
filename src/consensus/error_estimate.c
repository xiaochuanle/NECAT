#include "error_estimate.h"

#include "../common/oc_assert.h"
#include "../klib/kstring.h"
#include "consensus_aux.h"

static BOOL
is_good_overlap(int qoff,
				int qend,
				int qsize,
				int soff,
				int send,
				int ssize)
{
	int qlh = qoff;
	int qrh = qsize - qend;
	int slh = soff;
	int srh = ssize - send;
	const int M = 200;
	
	BOOL r =  (qlh <= M && qrh <= M)
			  ||
			  (slh <= M && srh <= M)
			  ||
			  (qrh <= M && slh <= M)
			  ||
			  (srh <= M && qlh <= M);
	
	return r;
}

static BOOL
estimate_ident_lower_bound(double ident[],
						   const int n,
						   double* avg_error,
						   double* error_cutoff)
{
	//if (n < 5) return FALSE;
	if (n < 5) {
		if (avg_error) *avg_error = 100.0;
		*error_cutoff = 0.0;
		return TRUE;
	}
	
	double sum = 0.0;
	for (int i = 0; i < n; ++i) sum += ident[i];
	double avg = sum / n;
	double se = 0.0;
	for (int i = 0; i < n; ++i) se += (avg - ident[i]) * (avg - ident[i]);
	se /= n;
	se = sqrt(se);
	double s = avg - se * 5;
	
	if (avg_error) *avg_error = avg;
	if (error_cutoff) *error_cutoff = s;

if (0) {
	for (int i = 0; i < n; ++i) printf("%lf\t", ident[i]);	
	printf("\n");
	printf("avg = %lf, se = %lf, cutoff = %lf\n\n", avg, se, s);
}

	return TRUE;
}

static void
get_idents(OverlapsPool* op,
		   double* ident,
		   const int NIdent,
		   int* n_ident,
		   const int tsize)
{
	int n = 0;
	for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
		OverlapIndex* oip = oc_op_oip(*op, i);
		if (is_good_overlap(oip->qoff, oip->qend, oip->qsize, oip->toff, oip->tend, tsize)) {
			ident[n++] = oip->ident_perc;
			if (n == NIdent) break;
		}
	}
	if (n < NIdent) {
		n = 0;
		for (size_t i = 0; i != oc_op_num_align(*op); ++i) {
			OverlapIndex* oip = oc_op_oip(*op, i);
			BOOL r = (oip->qend - oip->qoff >= oip->qsize * 0.6) || (oip->tend - oip->toff >= tsize * 0.6);
			if (r) {
				ident[n++] = oip->ident_perc;
				if (n == NIdent) break;
			}
		}
	}
	*n_ident = n;
}

BOOL
get_good_overlaps(PackedGappedCandidate* candidates,
				  size_t can_sid,
				  size_t can_eid,
				  size_t* last_extended_can_id,
				  OcAlignData* align_data,
				  OcDalignData* dalign_data,
				  FullEdlibAlignData*	falign_data,
				  BOOL rescue_long_indels,
				  OverlapsPool* op,
				  kstring_t* query,
				  kstring_t* target,
				  PackedDB* reads,
				  CnsReads* cns_reads,
				  ReadIdPool* extended_reads,
				  int min_align_size,
				  double* ident_cutoff)
{
	oc_op_clear(*op);
	const int NIdent = 15;
	int n_ident = 0;
	double ident[NIdent];
	size_t i;
	GappedCandidate _can, *can = &_can;;
	for (i = can_sid; i < can_eid && i < can_sid + 50; ++i) {
		unpack_candidate(can, candidates + i);
		can->qsize = reads ? PDB_SEQ_SIZE(reads, can->qid) : cns_reads_seq_size(cns_reads, can->qid);
		can->ssize = reads ? PDB_SEQ_SIZE(reads, can->sid) : cns_reads_seq_size(cns_reads, can->sid);
		if (read_id_exists(extended_reads, can->qid)) continue;
		oc_assert(can->sdir == FWD);
		if (reads) pdb_extract_sequence(reads, can->qid, can->qdir, query);
		else cns_reads_extract_sequence(cns_reads, can->qid, can->qdir, query);
		oc_assert(can->qsize == kstr_size(*query));
		oc_assert(can->ssize == kstr_size(*target));
		align_data->qid = can->qid;
		align_data->tid = can->sid;

		int qoff, qend, toff, tend;
		double ident_perc;
		kstring_t* qaln;
		kstring_t* taln;
		BOOL r = cns_extension(can, 
							   align_data,
							   dalign_data,
							   falign_data,
							   kstr_str(*query),
							   kstr_str(*target),
							   min_align_size,
							   rescue_long_indels,
							   &qoff,
							   &qend,
							   &toff,
							   &tend,
							   &ident_perc,
							   &qaln,
							   &taln);

		if (!r) { continue; }
		op_add_align(can->qid,
					 can->qdir,
					 can->qsize,
					 i,
					 qoff,
					 qend,
					 toff,
					 tend,
					 ident_perc,
					 qaln,
					 taln,
					 op);

		add_read_id(extended_reads, can->qid);
		r = is_good_overlap(qoff, 
							qend, 
							can->qsize, 
							toff, 
							tend,
							can->ssize);
		if (r) {
			ident[n_ident++] = ident_perc;
			if (n_ident == NIdent) break;
		}
	}
	*last_extended_can_id = i - 1;
	
	if (n_ident < NIdent) get_idents(op, ident, NIdent, &n_ident, can->ssize);
	ks_introsort_double_gt(n_ident, ident);
	if (n_ident >= 8) n_ident = n_ident * 0.7;
	return estimate_ident_lower_bound(ident, n_ident, NULL, ident_cutoff);
}
