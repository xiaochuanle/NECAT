#ifndef GAPPED_CANDIDATE_H
#define GAPPED_CANDIDATE_H

#include "ontcns_aux.h"
#include "ontcns_defs.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"

typedef struct 
{
	int qid;
	int sid;
	int qdir;
	int sdir;
	int score;
	idx qbeg, qend, qsize;
	idx sbeg, send, ssize;
	idx qoff, soff;
} GappedCandidate;

#define GappedCandidate_SidLT(a, b) ((a).sid < (b).sid)
void ks_introsort_GappedCandidate_SidLT(size_t n, GappedCandidate* a);

typedef kvec_t(GappedCandidate) vec_can;

#define DUMP_GAPPED_CANDIDATE(output_func, out, can) \
	output_func(out, \
			"%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
			(can).qid, \
			(can).sid, \
			(can).score, \
			(can).qdir, \
			(can).qbeg, \
			(can).qend, \
			(can).qoff, \
			(can).qsize, \
			(can).sdir, \
			(can).sbeg, \
			(can).send, \
			(can).soff, \
			(can).ssize) \

#define LOAD_GAPPED_CANDIDATE(input_func, in, can) \
	input_func(in, "%d\t%d\t%d\t%d\t%lu\t%lu\t%lu\t%lu\t%d\t%lu\t%lu\t%lu\t%lu\n", \
					&(can).qid, \
					&(can).sid, \
					&(can).score, \
					&(can).qdir, \
					&(can).qbeg, \
					&(can).qend, \
					&(can).qoff, \
					&(can).qsize, \
					&(can).sdir, \
					&(can).sbeg, \
					&(can).send, \
					&(can).soff, \
					&(can).ssize)
   
typedef struct {
	u32 item[7];
} PackedGappedCandidate;

#define PcanSdirShift 	31
#define PcanSdirMask 	(U32_ONE << PcanSdirShift)
#define PcanQdirShift	30
#define PcanQdirMask	(U32_ONE << PcanQdirShift)
#define PcanOffShift	29
#define PcanOffMask		(U32_ONE << PcanOffShift)
#define PcanScoreMask	((U32_ONE << PcanOffShift) - 1)

#define pcan_sdir(pcan)		(((pcan).item[0]&PcanSdirMask)>>PcanSdirShift)
#define pcan_qdir(pcan)		(((pcan).item[0]&PcanQdirMask)>>PcanQdirShift)
#define pcan_score(pcan)	((pcan).item[0]&PcanScoreMask)
#define pcan_off_flag(pcan)	((pcan).item[0]&PcanOffMask)
#define pcan_sid(pcan)		((pcan).item[1])
#define pcan_sbeg(pcan)		((pcan).item[2])
#define pcan_send(pcan)		((pcan).item[3])
#define pcan_qid(pcan)		((pcan).item[4])
#define pcan_qbeg(pcan)		((pcan).item[5])
#define pcan_qend(pcan)		((pcan).item[6])

#define DUMP_PACKED_GAPPED_CANDIDATE(fprintf, out, pcan) \
	do { \
	GappedCandidate __can; \
	unpack_candidate(&__can, &(pcan)); \
	DUMP_GAPPED_CANDIDATE(fprintf, out, __can); \
} while(0)

typedef kvec_t(PackedGappedCandidate) vec_pcan;

#define PackedGappedCandidate_SidLT(a, b) (pcan_sid(a) < pcan_sid(b))
void ks_introsort_PackedGappedCandidate_SidLT(size_t n, PackedGappedCandidate* a);
void ks_introsort_PackedGappedCandidate_CnsScoreGT(size_t n, PackedGappedCandidate* a);

void
pack_candidate(GappedCandidate* can, PackedGappedCandidate* pcan);

void
unpack_candidate(GappedCandidate* can, PackedGappedCandidate* pcan);
   
void
change_pcan_roles(PackedGappedCandidate* src, PackedGappedCandidate* dst);

void
normalise_pcan_sdir(PackedGappedCandidate* pcan, const u32 qsize, const u32 ssize);

#endif // GAPPED_CANDIDATE_H
