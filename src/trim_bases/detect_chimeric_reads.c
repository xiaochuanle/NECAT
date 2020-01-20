#include "detect_chimeric_reads.h"

#include "largest_cover_range.h"

static int
remove_low_quality_m4(M4Record* m4v,
                    int nm4,
                    const double min_ident_perc)
{
    int k = 0;
    for (int i = 0; i < nm4; ++i) {
        if (m4v[i].ident_perc >= min_ident_perc) m4v[k++] = m4v[i];
    }
    return k;
}

int is_complete_read(M4Record* m4v,
                    int nm4,
                    const double min_ident_perc,
                    int* fbgn,
                    int* fend)
{
    nm4 = remove_low_quality_m4(m4v, nm4, min_ident_perc);
    LargestCoverRange lcr;
    lcr.size = m4v[0].ssize;
    *fbgn = 0;
    *fend = lcr.size;
    for (int i = 0; i < nm4; ++i) {
        lcr.left = m4v[i].soff;
        lcr.right = m4v[i].send;
        if (lcr_is_complete(&lcr)) return 1;
    }
    return 0;
}

/// case I
/// uncomplete read, complete target
static int
is_chimeric_read_case_i(int qb1, int qe1,
    int qb2, int qe2,
    int tb1, int te1,
    int tb2, int te2,
    int qsize, int tsize)
{
    int lqb, lqe, rqb, rqe;
    int ltb, lte, rtb, rte;
    if (qb1 < qb2) {
        lqb = qb1;
        lqe = qe1;
        rqb = qb2;
        rqe = qe2;
    } else {
        lqb = qb2;
        lqe = qe2;
        rqb = qb1;
        rqe = qe1;
    }
    if (tb1 < tb2) {
        ltb = tb1;
        lte = te1;
        rtb = tb2;
        rte = te2;
    } else {
        ltb = tb2;
        lte = te2;
        rtb = tb1;
        rte = te1;
    }

    /// read overlap size do not vary too much
    int ov1 = lqe - lqb;
    int ov2 = rqe - rqb;
    int max_ov = OC_MAX(ov1, ov2);
    int min_ov = OC_MIN(ov1, ov2);
    if (min_ov < max_ov * 0.9) return 0;

    /// read share 90% of bps
    int common_read_bps = 0;
    if (lqe > rqb) common_read_bps = (lqe - rqb);
    int r = (common_read_bps >= (lqe - lqb) * 0.9) && (common_read_bps >= (rqe - rqb) * 0.9);
    if (!r) return r;

    /// complete target
    int mapped_target_bps = rte - ltb;
    if (rtb > lte) mapped_target_bps -= (rtb - lte);
    r = mapped_target_bps >= tsize * 0.9;
    if (!r) return r;

    /// target share not too much
    if (lte > rtb) {
        int ov = lte - rtb;
        r = (ov < (lte - ltb) * 0.4) && (ov < (rte - rtb) * 0.4);
    }
    
    return r ? 1 : 0;
}

/// case II
/// complete read, uncomplete target
static int
is_chimeric_read_case_ii(int qb1, int qe1,
    int qb2, int qe2,
    int tb1, int te1,
    int tb2, int te2,
    int qsize, int tsize)
{
    int lqb, lqe, rqb, rqe;
    int ltb, lte, rtb, rte;
    if (qb1 < qb2) {
        lqb = qb1;
        lqe = qe1;
        rqb = qb2;
        rqe = qe2;
    } else {
        lqb = qb2;
        lqe = qe2;
        rqb = qb1;
        rqe = qe1;
    }
    if (tb1 < tb2) {
        ltb = tb1;
        lte = te1;
        rtb = tb2;
        rte = te2;
    } else {
        ltb = tb2;
        lte = te2;
        rtb = tb1;
        rte = te1;
    }

    /// read overlap size do not vary too much
    int ov1 = lqe - lqb;
    int ov2 = rqe - rqb;
    int max_ov = OC_MAX(ov1, ov2);
    int min_ov = OC_MIN(ov1, ov2);
    if (min_ov < max_ov * 0.9) return 0;

    /// read share 90% of bps
    int common_read_bps = 0;
    if (lqe > rqb) common_read_bps = (lqe - rqb);
    int r = (common_read_bps >= (lqe - lqb) * 0.9) && (common_read_bps >= (rqe - rqb) * 0.9);
    if (!r) return r;  


    /// complete read
    r = (lqe - lqb >= qsize * 0.9) && (rqe - rqb >= qsize * 0.9);
    if (!r) return r;

    /// target should close to each other
    if (lte > rtb) {
        r = lte - rtb <= 1000;
    } else {
        r = rtb - lte <= 1000;
    }
    return r ? 2 : 0;  
}

#define m4_record_chimeric_cmp(a, b) ( \
    ((a).qid < (b).qid) \
    || \
    ((a).qid == (b).qid && (a).qdir < (b).qdir) \
    || \
    ((a).qid == (b).qid && (a).qdir == (b).qdir && (a).vscore > (b).vscore) \
)

KSORT_INIT(m4_record_chimeric_cmp, M4Record, m4_record_chimeric_cmp)

int is_chimeric_read(M4Record* m4v,
                    int nm4,
                    const double min_ident_perc,
                    int* fbgn,
                    int* fend)
{
    nm4 = remove_low_quality_m4(m4v, nm4, min_ident_perc);
    ks_introsort_m4_record_chimeric_cmp(nm4, m4v);
    int i = 0;
    int max_size = 0;
    int num_chimeric_reads = 0;
    while (i < nm4) {
        int j = i+1;
        while (j < nm4 && m4v[j].qid == m4v[i].qid) ++j;
        int k = i + 1;
        while (k < j && m4v[k].qdir == m4v[i].qdir) ++k;
        if (k < j) {
            int r = is_chimeric_read_case_i(m4v[i].qoff, m4v[i].qend,
                        m4v[k].qoff, m4v[k].qend,
                        m4v[i].soff, m4v[i].send,
                        m4v[k].soff, m4v[k].send,
                        m4v[i].qsize, m4v[i].ssize)
                    ||
                    is_chimeric_read_case_ii(m4v[i].qoff, m4v[i].qend,
                        m4v[k].qoff, m4v[k].qend,
                        m4v[i].soff, m4v[i].send,
                        m4v[k].soff, m4v[k].send,
                        m4v[i].qsize, m4v[i].ssize);
            if (r) {
                ++num_chimeric_reads;
                if (m4v[i].send - m4v[i].soff > max_size) {
                    max_size = m4v[i].send - m4v[i].soff;
                    *fbgn = m4v[i].soff;
                    *fend = m4v[i].send;
                }
                if (m4v[k].send - m4v[k].soff > max_size) {
                    max_size = m4v[k].send - m4v[k].soff;
                    *fbgn = m4v[k].soff;
                    *fend = m4v[k].send;
                }    
                //HBN_LOG("find chimeric read case %d", r);
                //DUMP_M4_RECORD(fprintf, stderr, m4v[i]);
                //DUMP_M4_RECORD(fprintf, stderr, m4v[k]);
            }
        }
        i = j;
    }
    return (max_size > 0) && (num_chimeric_reads > 1);
}
