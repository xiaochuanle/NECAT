#include "gapped_candidate.h"

#include "../common/oc_assert.h"
#include "../klib/ksort.h"

#include <assert.h>

KSORT_INIT(GappedCandidate_SidLT, GappedCandidate, GappedCandidate_SidLT)

KSORT_INIT(PackedGappedCandidate_SidLT, PackedGappedCandidate, PackedGappedCandidate_SidLT)

void
pack_candidate(GappedCandidate* can, PackedGappedCandidate* pcan)
{
	memset(pcan, 0, sizeof(PackedGappedCandidate));
	if (can->sdir == REV) pcan->item[0] |= PcanSdirMask;
	if (can->qdir == REV) pcan->item[0] |= PcanQdirMask;
	if (can->qoff == can->qbeg) pcan->item[0] |= PcanOffMask;
	u32 score = OC_MIN(1000000, can->score);
	pcan->item[0] |= score;
	
	pcan->item[1] = can->sid;
	pcan->item[2] = can->sbeg;
	pcan->item[3] = can->send;
	pcan->item[4] = can->qid;
	pcan->item[5] = can->qbeg;
	pcan->item[6] = can->qend;
}

void
unpack_candidate(GappedCandidate* can, PackedGappedCandidate* pcan)
{
	can->sdir = pcan_sdir(*pcan);
	can->qdir = pcan_qdir(*pcan);
	can->score = pcan_score(*pcan);
	
	can->sid = pcan->item[1];
	can->sbeg = pcan->item[2];
	can->send = pcan->item[3];
	can->qid = pcan->item[4];
	can->qbeg = pcan->item[5];
	can->qend = pcan->item[6];
	
	if (pcan_off_flag(*pcan)) {
		can->qoff = can->qbeg;
		can->soff = can->sbeg;
	} else {
		can->qoff = can->qend;
		can->soff = can->send;
	}
}

void
change_pcan_roles(PackedGappedCandidate* src, PackedGappedCandidate* dst)
{
	u32 item = 0;
	if (pcan_sdir(*src)) item |= PcanQdirMask;
	if (pcan_qdir(*src)) item |= PcanSdirMask;
	if (pcan_off_flag(*src)) item |= PcanOffMask;
	item |= pcan_score(*src);
	dst->item[0] = item;
	
	pcan_sid(*dst) 	= pcan_qid(*src);
	pcan_sbeg(*dst) = pcan_qbeg(*src);
	pcan_send(*dst) = pcan_qend(*src);
	pcan_qid(*dst) 	= pcan_sid(*src);
	pcan_qbeg(*dst) = pcan_sbeg(*src);
	pcan_qend(*dst) = pcan_send(*src);
}

void
normalise_pcan_sdir(PackedGappedCandidate* pcan, const u32 qsize, const u32 ssize)
{
	if (pcan_sdir(*pcan) == FWD) return;
	oc_assert(pcan_sdir(*pcan) == REV);
	
	u32 item = 0;
	if (pcan_qdir(*pcan) == FWD) item |= PcanQdirMask;
	if (!pcan_off_flag(*pcan)) item |= PcanOffMask;
	item |= pcan_score(*pcan);
	pcan->item[0] = item;
	
	u32 beg, end;
	beg = qsize - pcan_qend(*pcan);
	end = qsize - pcan_qbeg(*pcan);
	pcan_qbeg(*pcan) = beg;
	pcan_qend(*pcan) = end;
	
	beg = ssize - pcan_send(*pcan);
	end = ssize - pcan_sbeg(*pcan);
	pcan_sbeg(*pcan) = beg;
	pcan_send(*pcan) = end;
}

static int
PackedGappedCandidate_CnsScoreGT(PackedGappedCandidate a, PackedGappedCandidate b)
{
	int a_item, b_item;

	a_item = pcan_score(a);
	b_item = pcan_score(b);
	if (a_item != b_item) return a_item > b_item;
	
	a_item = pcan_qid(a);
	b_item = pcan_qid(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_qdir(a);
	b_item = pcan_qdir(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_qbeg(a);
	b_item = pcan_qbeg(b);
	if (a_item != b_item) return a_item < b_item;
	
	a_item = pcan_sbeg(a);
	b_item = pcan_sbeg(b);
	if (a_item != b_item) return a_item < b_item;
	
	return 0;
}

KSORT_INIT(PackedGappedCandidate_CnsScoreGT, PackedGappedCandidate, PackedGappedCandidate_CnsScoreGT);
