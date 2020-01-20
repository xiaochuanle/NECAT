#include "hbn_align.h"

HbnMapData*
hbn_map_data_new()
{
    HbnMapData* data = (HbnMapData*)malloc(sizeof(HbnMapData));
    data->edlib_data = new_BlockwiseEdlibData(0.5, kOcaBlockSize2048);
    data->dalign_data = new_OcDalignData(0.35);
    kv_init(data->qstr);
    kv_init(data->tstr);
    return data;
}

HbnMapData*
hbn_map_data_free(HbnMapData* data)
{
    data->edlib_data = free_BlockwiseEdlibData(data->edlib_data);
    data->dalign_data = free_OcDalignData(data->dalign_data);
    kv_destroy(data->qstr);
    kv_destroy(data->tstr);
    free(data);
    return NULL;
}

static void
validate_aligned_string(const u8* query,
    const u8* target,
    int qoff,
    int qend,
    int toff,
    int tend,
    kstring_t* qaln,
    kstring_t* taln)
{
    hbn_assert(ks_size(*qaln) == ks_size(*taln));
    int q_idx = qoff;
    int t_idx = toff;
    for (size_t i = 0; i < ks_size(*qaln); ++i) {
        char c = ks_A(*qaln, i);
        if (c != '-') {
            int c1 = query[q_idx];
            c1 = "ACGT"[c1];
            hbn_assert(c == c1);
            ++q_idx;
        }
        c = ks_A(*taln, i);
        if (c != '-') {
            int c1 = target[t_idx];
            c1 = "ACGT"[c1];
            hbn_assert(c == c1);
            ++t_idx;
        }
    }
    hbn_assert(q_idx == qend);
    hbn_assert(t_idx == tend);
}

BOOL 
hbn_cns_extend(HbnMapData* hbn_map_data,
	const u8* query,
	const u8* target,
    GappedCandidate* can,
    const int min_align_size,
    const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align,
	kstring_t* target_align)
{
    int r = 0;
    r = blockwise_edlib_align(hbn_map_data->edlib_data,
            query,
            target,
            can,
            min_align_size,
            min_ident_perc,
            qbeg,
            qend,
            tbeg,
            tend,
            ident_perc,
            query_align,
            target_align);
    if (r) validate_aligned_string(query, target, *qbeg, *qend, *tbeg, *tend, query_align, target_align);
    return r;
}

static int
left_extend(HbnMapData* map_data,
    const u8* query,
    const u8* target,
    int* qbeg,
    int* tbeg,
    kstring_t* query_align,
    kstring_t* target_align,
    size_t* a_idx_from,
    int* dist)
{
    const int kMaxHang = 300;
    const int kMatchSize = 8;
    *dist = 0;
    *a_idx_from = 0;
    int qls = *qbeg;
    int tls = *tbeg;
    int ls = hbn_min(qls, tls);
    if (ls > kMaxHang) return 0;
    if (ls == 0) return 0;
    
    ls = 0;
    hbn_assert(ks_size(*query_align) == ks_size(*target_align));
    int qi = 0;
    int ti = 0;
    size_t i = 0;
    for (; ls < kMatchSize && i < ks_size(*query_align); ++i) {
        int qc = ks_A(*query_align, i);
        int tc = ks_A(*target_align, i);
        if (qc != '-') ++qi;
        if (tc != '-') ++ti;
        if (qc == tc) {
            ++ls;
        } else {
            ls = 0;
        }
    }
    if (ls < kMatchSize) return 0;
    *a_idx_from = i + 1;

//HBN_LOG("lqi = %d, lti = %d", qi, ti);
    qls += qi;
    tls += ti;
    const u8* qs = query;
    const u8* ts = target;
    int r = 0;
    GappedCandidate can; memset(&can, 0, sizeof(GappedCandidate));
    can.qbeg = 0;
    can.qend = qls;
    can.qsize = qls;
    can.qoff = qls;
    can.sbeg = 0;
    can.send = tls;
    can.ssize = tls;
    can.soff = tls;
    //DUMP_GAPPED_CANDIDATE(fprintf, stderr, can);
    int qb, qe, tb, te;
    double ident_perc;
    r = ocda_go(map_data->dalign_data,
            qs,
            ts,
            &can,
            1,
            0.0,
            &qb,
            &qe,
            &tb,
            &te,
            &ident_perc,
            NULL,
            NULL);
    if (!r) return 0;
    if (qe != qls || te != tls) return 0;

//HBN_LOG("left diff = %d", ocda_distance(*map_data->dalign_data));
    *qbeg = qb;
    *tbeg = tb;
    *dist = ocda_distance(*map_data->dalign_data);
    return 1;
}

static int
right_extend(HbnMapData* map_data,
    const u8* query,
    const u8* target,
    int* qend,
    int qsize,
    int* tend,
    int tsize,
    kstring_t* query_align,
    kstring_t* target_align,
    size_t* a_idx_to,
    int* dist)
{
    const int kMaxHang = 300;
    const int kMatchSize = 8;
    *dist = 0;
    *a_idx_to = ks_size(*query_align);
    int qrs = qsize - (*qend);
    int trs = tsize - (*tend);
    int rs = hbn_min(qrs, trs);
    if (rs > kMaxHang) return 0;
    if (rs == 0) return 0;
    
    rs = 0;
    hbn_assert(ks_size(*query_align) == ks_size(*target_align));
    size_t i = ks_size(*query_align);
    int qi = 0;
    int ti = 0;
    while (i && rs < kMatchSize) {
        --i;
        int qc = ks_A(*query_align, i);
        int tc = ks_A(*target_align, i);
        if (qc != '-') ++qi;
        if (tc != '-') ++ti;
        if (qc == tc) {
            ++rs;
        } else {
            rs = 0;
        }
    }
    if (rs < kMatchSize) return 0;
    *a_idx_to = i;

//HBN_LOG("qi = %d, ti = %d", qi, ti);
    qrs = (*qend) - qi;
    trs = (*tend) - ti;
    const u8* qs = (query + qrs);
    const u8* ts = (target + trs);
    int r = 0;
    GappedCandidate can;
    can.qbeg = 0;
    can.qend = qsize - qrs;
    can.qsize = qsize - qrs;
    can.qoff = 0;
    can.sbeg = 0;
    can.send = tsize - trs;
    can.ssize = tsize - trs;
    can.soff = 0;
    int qb, qe, tb, te;
    double ident_perc;
    r = ocda_go(map_data->dalign_data,
            qs,
            ts,
            &can,
            1,
            0.0,
            &qb,
            &qe,
            &tb,
            &te,
            &ident_perc,
            NULL,
            NULL);
    if (!r) return 0;
    if (qb != 0 || tb != 0) return 0;

//HBN_LOG("right diff = %d", ocda_distance(*map_data->dalign_data));
    *qend = qrs + qe;
    *tend = trs + te;
    *dist = ocda_distance(*map_data->dalign_data);
    return 1;
}

static double 
fix_ident_perc(int left_diff,
    int right_diff,
    kstring_t* query_align,
    kstring_t* target_align,
    size_t align_from,
    size_t align_to,
    int qbeg,
    int qend,
    int tbeg,
    int tend)
{
    int diff = left_diff + right_diff;
    hbn_assert(ks_size(*query_align) == ks_size(*target_align));
    size_t align_size = ks_size(*query_align);
    hbn_assert(align_from < align_to, "align_from = %zu, align_to = %zu, align_size = %zu",
        align_from, align_to, align_size);
    hbn_assert(align_to <= align_size);
    for (size_t i = align_from; i < align_to; ++i) {
        int qc = ks_A(*query_align, i);
        int tc = ks_A(*target_align, i);
        if (qc != tc) ++diff;
    }
    return 100.0 - 200.0 * diff / (qend + tend - qbeg - tbeg);
}

BOOL 
hbn_map_extend(HbnMapData* hbn_map_data,
	const u8* query,
	const u8* target,
    GappedCandidate* can,
    const int min_align_size,
    const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align,
	kstring_t* target_align)
{
    int r = 0;
    r = hbn_cns_extend(hbn_map_data,
            query,
            target,
            can,
            min_align_size,
            min_ident_perc,
            qbeg,
            qend,
            tbeg,
            tend,
            ident_perc, 
            query_align,
            target_align);  
    if (r == 0) return r;
    //HBN_LOG("before0:\t[%d, %d, %d] x [%d, %d, %d], %g", *qbeg, *qend, can->qsize, *tbeg, *tend, can->ssize, *ident_perc);

    int qb = *qbeg;
    int qe = *qend;
    int tb = *tbeg;
    int te = *tend;
    size_t a_idx_from, a_idx_to;
    int left_dist, right_dist;
    int lext = left_extend(hbn_map_data, query, target, qbeg, tbeg, query_align, target_align, &a_idx_from, &left_dist);
    int rext = right_extend(hbn_map_data, query, target, qend, can->qsize, tend, can->ssize, query_align, target_align, &a_idx_to, &right_dist);
    if (lext || rext) {
        //HBN_LOG("before:\t[%d, %d, %d] x [%d, %d, %d], %g", qb, qe, can->qsize, tb, te, can->ssize, *ident_perc);
        *ident_perc = fix_ident_perc(left_dist, right_dist, query_align, target_align, a_idx_from, a_idx_to, *qbeg, *qend, *tbeg, *tend);
        //HBN_LOG("after:\t[%d, %d, %d] x [%d, %d, %d], %g", *qbeg, *qend, can->qsize, *tbeg, *tend, can->ssize, *ident_perc);
    }
    return 1;
}