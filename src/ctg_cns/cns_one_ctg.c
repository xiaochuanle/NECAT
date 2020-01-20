#include "cns_one_ctg.h"

#include "chain_dp.h"
#include "mem_finder.h"
#include "cns_ctg_subseq.h"
#include "../common/makedb_aux.h"
#include "../common/packed_db.h"
#include "../common/m4_record.h"
#include "../gapped_align/oc_daligner.h"
#include "../gapped_align/oc_aligner.h"
#include "../edlib/edlib_wrapper.h"
#include "../klib/ksort.h"

static const int kCtgSegmentSize = 1000000;

#define m4_soff_lt(a, b) ((a).soff < (b).soff)
KSORT_INIT(m4_soff_lt, M4Record, m4_soff_lt);


static M4Record*
load_m4v(const char* wrk_dir, const int ctg_id, int * m4_count)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/p%d", wrk_dir, ctg_id);
    OC_LOG("loading %s", path);
    size_t s = FILE_SIZE(path);
    oc_assert(s % sizeof(M4Record) == 0);
    s /= sizeof(M4Record);
    if (s == 0) return NULL;

    M4Record* a = (M4Record*)calloc(s, sizeof(M4Record));
    DFOPEN(in, path, "rb");
    FREAD(a, sizeof(M4Record), s, in);
    FCLOSE(in);
    *m4_count = s;
    return a;
}

static int
m4_is_contained_in_ctg(const M4Record* m4, const idx ctg_from, const idx ctg_to)
{
    if (m4->soff >= ctg_to) return -1;
    if (m4->send <= ctg_from) return 0;
    if (m4->soff >= ctg_from && m4->send <= ctg_to) return 1;
    if (m4->send > ctg_from && m4->send - ctg_from >= 1000) return 1;
    if (ctg_to > m4->soff && ctg_to - m4->soff >= 1000) return 1;
    return 0;
}

static M4Record*
extract_cov_m4s(const M4Record* m4_array,
    const int m4_count,
    const idx ctg_from,
    const idx ctg_to,
    int* cov_m4_count)
{
    int n = 0;
    for (int i = 0; i < m4_count; ++i) {
        int r = m4_is_contained_in_ctg(m4_array + i, ctg_from, ctg_to);
        if (r == -1) break;
        if (r) ++n;
    }
    M4Record* a = (M4Record*)calloc(n, sizeof(M4Record));
    int m = 0;
    for (int i = 0; i < m4_count; ++i) {
        int r = m4_is_contained_in_ctg(m4_array + i, ctg_from, ctg_to);
        if (r == -1) break;
        if (r) a[m++] = m4_array[i];
    }
    oc_assert(m == n);
    *cov_m4_count = m;
    return a;
}

void
cns_one_ctg(const char* mkdb_dir,
    PackedDB* reads,
    const char* contig,
    const idx ctg_size,
    const int ctg_id,
    kstring_t* cns_ctg)
{
    int m4_count = 0;
    M4Record* m4_array = load_m4v(mkdb_dir, ctg_id, &m4_count);
    ks_introsort_m4_soff_lt(m4_count, m4_array);
    idx ctg_subseq_from = 0;
    idx ctg_subseq_to = 0;
    idx left = ctg_size;
    while (left) {
        ctg_subseq_to = ctg_subseq_from + kCtgSegmentSize;
        ctg_subseq_to = OC_MIN(ctg_size, ctg_subseq_to);
        left = ctg_size - ctg_subseq_to;
        if (left <= kCtgSegmentSize/2) {
            ctg_subseq_to = ctg_size;
            left = 0;
        }
        idx ctg_subseq_size = ctg_subseq_to - ctg_subseq_from;

        u8* ctg_subseq = (u8*)(contig + ctg_subseq_from);
        int cov_m4_count = 0;
        M4Record* cov_m4_array = extract_cov_m4s(m4_array, m4_count, ctg_subseq_from, ctg_subseq_to, &cov_m4_count);
        OC_LOG("correct %zu --- %zu, number of m4: %d", ctg_subseq_from, ctg_subseq_to, cov_m4_count);
        cns_ctg_subseq(reads, ctg_subseq, ctg_subseq_from, ctg_subseq_size, cov_m4_array, cov_m4_count, cns_ctg);
        free(cov_m4_array);
        ctg_subseq_from = ctg_subseq_to;
    }
    free(m4_array);
}