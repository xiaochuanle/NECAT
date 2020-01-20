#ifndef FC_CORRECT_ONE_READ_H
#define FC_CORRECT_ONE_READ_H

#include "../common/ontcns_aux.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"
#include "small_object_alloc.h"

#include <pthread.h>
#include <stdint.h>
#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int             cns_pos_t;
typedef uint16_t        cns_delta_t;
#define MAX_CNS_DELTA   UINT16_MAX

typedef struct {
    double          weight;
    cns_pos_t       t_pos;
    cns_pos_t       p_t_pos;
    cns_delta_t     delta;
    cns_delta_t     p_delta;
    char            q_base;
    char            p_q_base;
} AlignTag;

#define DUMP_ALIGN_TAG_STD(tag) do { \
    fprintf(stdout, "t_pos = %d, p_t_pos = %d, delta = %d, p_delta = %d, q_base = %c, p_q_base = %c", \
        (tag).t_pos, \
        (tag).p_t_pos, \
        (tag).delta, \
        (tag).p_delta, \
        (tag).q_base, \
        (tag).p_q_base); \
} while(0)

typedef kvec_t(AlignTag) vec_align_tag;

void ks_introsort_align_tag(size_t n, AlignTag* tags);

void add_align_tags(const char* qaln,
        const char* taln,
        const cns_pos_t aln_size,
        const cns_pos_t qoff,
        const cns_pos_t qend,
        const cns_pos_t toff,
        const cns_pos_t tend,
        double weight,
        vec_align_tag* tags);

typedef struct {
    double weight;
    cns_pos_t p_t_pos;
    cns_delta_t p_delta;
    char p_q_base;
    int link_cout;
} LinkInfo;

typedef struct {
    int n_link;
    int coverage;
    LinkInfo* plinks;
    cns_pos_t best_p_t_pos;
    cns_delta_t best_p_delta;
    char best_p_q_base;
    double score;
} BaseLinks;

typedef struct {
    BaseLinks links[5];
} DeltaCovInfo;

typedef struct {
    int n_delta;
    DeltaCovInfo* delta;
} BackboneItem;

typedef kvec_t(BackboneItem) vec_backbone_item;

#define backbone_item_clear(item) ((item).n_delta = 0, (item).delta = 0)

typedef kvec_t(void*) vec_voidp;

typedef struct {
    vec_align_tag tags;
    vec_backbone_item backbone;
    vec_int coverage;
    vec_voidp p_alloc;
    //void*   km;
    kstring_t cns_seq;
    vec_int t_pos;
    CnsMemoryAlloc* li_alloc;
    CnsMemoryAlloc* dci_alloc;
} CnsData;

CnsData*
cns_data_new();

CnsData* 
cns_data_free(CnsData* data);

void cns_data_clear(CnsData* data);

void
get_cns_from_align_tags(CnsData* data,
    cns_pos_t target_size,
    cns_pos_t min_cov);

void
dump_cns_seq(CnsData* data,
    const char* hdr,
    u8* template,
    const int template_size,
    const int min_size,
    const int split_reads,
    FILE* cns_out,
    FILE* raw_out,
    pthread_mutex_t* out_lock);

#ifdef __cplusplus
}
#endif
#endif // FC_CORRECT_ONE_READ_H