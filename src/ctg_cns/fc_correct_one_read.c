#include "fc_correct_one_read.h"

#include "kalloc.h"
#include "../common/oc_assert.h"

#define align_tag_plink_eq(a, b) ( \
    (a).p_t_pos == (b).p_t_pos \
    && \
    (a).p_delta == (b).p_delta \
    && \
    (a).p_q_base == (b).p_q_base \
)

CnsData*
cns_data_new()
{
    CnsData* data = (CnsData*)malloc(sizeof(CnsData));
    kv_init(data->tags);
    kv_init(data->backbone);
    kv_init(data->coverage);
    kv_init(data->p_alloc);
    //data->km = km_init();
    kstr_init(data->cns_seq);
    kv_init(data->t_pos);
    data->li_alloc = cns_mem_alloc_new(sizeof(LinkInfo));
    data->dci_alloc = cns_mem_alloc_new(sizeof(DeltaCovInfo));
    return data;
}

CnsData* 
cns_data_free(CnsData* data)
{
    kv_destroy(data->tags);
    kv_destroy(data->backbone);
    kv_destroy(data->coverage);
    kv_destroy(data->p_alloc);
    //km_destroy(data->km);
    free_kstring(data->cns_seq);
    kv_destroy(data->t_pos);
    data->li_alloc = cns_mem_alloc_free(data->li_alloc);
    data->dci_alloc = cns_mem_alloc_free(data->dci_alloc);
    free(data);
    return NULL;
}

void cns_data_clear(CnsData* data)
{
    //HBN_LOG("number of alloc: %lu", kv_size(data->p_alloc));
#if 0
    for (size_t i = 0; i < kv_size(data->p_alloc); ++i) {
        void* p = kv_A(data->p_alloc, i);
        kfree(data->km, p);
    }
#endif
    cns_mem_alloc_clear(data->li_alloc);
    cns_mem_alloc_clear(data->dci_alloc);
    kv_clear(data->tags);
    kv_clear(data->backbone);
    kv_clear(data->coverage);
    kv_clear(data->p_alloc);
    kstr_clear(data->cns_seq);
    kv_clear(data->t_pos);
}

#define AlignTag_LT(a, b) ( \
	((a).t_pos < (b).t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta < (b).delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base < (b).q_base)\
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos < (b).p_t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta < (b).p_delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base < (b).p_q_base) \
	)

KSORT_INIT(align_tag, AlignTag, AlignTag_LT);

void add_align_tags(const char* qaln,
        const char* taln,
        const cns_pos_t aln_size,
        const cns_pos_t qoff,
        const cns_pos_t qend,
        const cns_pos_t toff,
        const cns_pos_t tend,
        double weight,
        vec_align_tag* tags)
{
    char p_q_base = '.';
    int i = qoff - 1; // position in query
    int j = toff - 1; // postion in template
    int jj = 0; // number of non-gap bases in query to a gap in template
    int p_j = -1;
    int p_jj = 0;
    cns_pos_t k;
    AlignTag tag;
    tag.weight = weight;

    for (k = 0; k < aln_size; ++k) {
        if (qaln[k] != '-') {
            ++i;
            ++jj;
        }
        if (taln[k] != '-') {
            ++j;
            jj = 0;
        }
        oc_assert(i >= 0);
        oc_assert(i < qend);
        oc_assert(j >= 0);
        oc_assert(j < tend);

        if ( (jj >= MAX_CNS_DELTA) || (p_jj >= MAX_CNS_DELTA) ) break;
        
        tag.t_pos = j;
        tag.delta = jj;
        tag.p_t_pos = p_j;
        tag.p_delta = p_jj;
        tag.p_q_base = p_q_base;
        tag.q_base = qaln[k];
        
        p_j = j;
        p_jj = jj;
        p_q_base = qaln[k];

        oc_assert(tag.q_base != '.');
        kv_push(AlignTag, *tags, tag);
    }  
}

static u8 
encode_dna_base(const char c)
{
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
        default: return 4;
    }
    OC_ERROR("invalid dna base: %c", c);
    return 4;
}

static void
build_base_links(AlignTag* tags, 
    const int ntag, 
    BaseLinks* link, 
    CnsMemoryAlloc* li_alloc)
{
    int n_link = 0;
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && align_tag_plink_eq(tags[i], tags[j])) ++j;
        ++n_link;
        i = j;
    }

    link->plinks = (LinkInfo*)cns_mem_alloc_alloc(li_alloc, n_link);
    for (i = 0; i < n_link; ++i) {
        link->plinks[i].link_cout = 0;
        link->plinks[i].weight = 0.0;
    }
    link->n_link = n_link;
    link->coverage = ntag;

    n_link = 0;
    i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && align_tag_plink_eq(tags[i], tags[j])) ++j;
        LinkInfo* linfo = link->plinks + n_link;
        linfo->p_t_pos = tags[i].p_t_pos;
        linfo->p_delta = tags[i].p_delta;
        linfo->p_q_base = tags[i].p_q_base;
        linfo->link_cout = j - i;
        linfo->weight = 0;
        for (int k = i; k < j; ++k) linfo->weight += tags[k].weight;
        ++n_link;
        i = j;
    }

    oc_assert(n_link == link->n_link);
}

static void
build_delta_links(AlignTag* tags, 
    const int ntag, 
    DeltaCovInfo* dci, 
    CnsMemoryAlloc* li_alloc)
{
    for (int i = 0; i < 5; ++i) {
        dci->links[i].n_link = 0;
        dci->links[i].coverage = 0;
        dci->links[i].best_p_t_pos = -1;
        dci->links[i].best_p_delta = MAX_CNS_DELTA;
        dci->links[i].best_p_q_base = '.';
        dci->links[i].score = 0;
    }

    int i = 0;
    while (i < ntag) {
        int j = i;
        while (j < ntag && tags[i].q_base == tags[j].q_base) ++j;
        if (tags[i].q_base == '.') DUMP_ALIGN_TAG_STD(tags[i]);
        u8 c = encode_dna_base(tags[i].q_base);
        build_base_links(tags + i, j - i, dci->links + c, li_alloc);
        i = j;
    }
}

static void
build_backbone_item(AlignTag* tags,
    const int ntag,
    BackboneItem* item,
    int* coverage,
    CnsMemoryAlloc* li_alloc,
    CnsMemoryAlloc* dci_alloc)
{
    item->n_delta = tags[ntag - 1].delta + 1;
    item->delta = (DeltaCovInfo*)cns_mem_alloc_alloc(dci_alloc, item->n_delta);
    for (int i = 0; i < item->n_delta; ++i) {
        for (int j = 0; j < 5; ++j) {
            item->delta[i].links[j].n_link = 0;
            item->delta[i].links[j].coverage = 0;
            item->delta[i].links[j].score = 0.0;
        }
    }
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tags[i].delta == tags[j].delta) ++j;
        build_delta_links(tags + i, j - i, item->delta + tags[i].delta, li_alloc);
        if (tags[i].delta == 0) coverage[ tags[i].t_pos ] = j - i;
        i = j;
    }
}

void
build_backbone(AlignTag* tags,
    const int ntag,
    const int target_size,
    CnsMemoryAlloc* li_alloc,
    CnsMemoryAlloc* dci_alloc,
    vec_backbone_item* backbone,
    vec_int* coverage)
{
    kv_resize(BackboneItem, *backbone, target_size);
    for (int i = 0; i < target_size; ++i) backbone_item_clear(kv_A(*backbone, i));
    kv_resize(int, *coverage, target_size);
    kv_zero(int, *coverage);

    ks_introsort_align_tag(ntag, tags);
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tags[i].t_pos == tags[j].t_pos) ++j;
        oc_assert(tags[i].t_pos < target_size);
        build_backbone_item(tags + i, j - i, kv_data(*backbone) +  tags[i].t_pos, kv_data(*coverage), li_alloc, dci_alloc);
        i = j;
    }
}

#define fc_reverse_array(type, a, n) do { \
    size_t __l = 0, __r = (n); \
    type* array = (a); \
    while (__l < __r) { \
        type tmp = array[__l]; \
        array[__l] = array[__r-1]; \
        array[__r-1] = tmp; \
        ++__l; \
        --__r; \
    } \
} while(0)

void
get_cns_from_align_tags(CnsData* data,
    cns_pos_t target_size,
    cns_pos_t min_cov)
{
    build_backbone(kv_data(data->tags),
        kv_size(data->tags),
        target_size,
        data->li_alloc,
        data->dci_alloc,
        &data->backbone,
        &data->coverage);

    BaseLinks* g_best_aln_col = 0;
    int g_best_t_pos = -1;
    double g_best_score = -1.0;
    int kk;
    int ck;
    int best_i;
    int best_j;
    int best_b;
    double score;
    double best_score;
    BaseLinks* aln_col;  
    BackboneItem* backbone = kv_data(data->backbone);
    int* coverage = kv_data(data->coverage);

    for (cns_pos_t i = 0; i < target_size; ++i) {
        for (cns_pos_t j = 0; j < backbone[i].n_delta; ++j) {
            for (kk = 0; kk < 5; ++kk) {
                aln_col = backbone[i].delta[j].links + kk;
                best_score = -1;
                best_i = -1;
                best_j = -1;
                best_b = -1;
                for (ck = 0; ck < aln_col->n_link; ++ck) {
                    cns_pos_t pi = aln_col->plinks[ck].p_t_pos;
                    cns_pos_t pj = aln_col->plinks[ck].p_delta;
                    int pkk = encode_dna_base(aln_col->plinks[ck].p_q_base);
                    score = aln_col->plinks[ck].weight - 0.1 * coverage[i];
                    if (pi != -1) score += backbone[pi].delta[pj].links[pkk].score;

                    if (score > best_score) {
                        best_score = score;
                        aln_col->best_p_t_pos = best_i = pi;
                        aln_col->best_p_delta = best_j = pj;
                        aln_col->best_p_q_base = best_b = pkk;
                    }
                }
                aln_col->score = best_score;
                if (best_score > g_best_score) {
                    g_best_score = best_score;
                    g_best_aln_col = aln_col;
                    g_best_t_pos = i;
                }
            }
        }
    } 

    int i = g_best_t_pos;
    int j = 0;
    kk = g_best_aln_col ? g_best_aln_col->best_p_q_base : 0;
    while (i != -1) {
        char bb = '-';
        switch (kk)
        {
            case 0: bb = (coverage[i] < min_cov) ? 'a' : 'A'; break;
            case 1: bb = (coverage[i] < min_cov) ? 'c' : 'C'; break;
            case 2: bb = (coverage[i] < min_cov) ? 'g' : 'G'; break;
            case 3: bb = (coverage[i] < min_cov) ? 't' : 'T'; break;
            case 4: bb = '-'; break;
        }
        //HBN_LOG("i = %d, bb = %c", i, bb);
        if (bb != '-') {
            kputc(bb, &data->cns_seq);
            kv_push(cns_pos_t, data->t_pos, i);
        }

        i = g_best_aln_col->best_p_t_pos;
        j = g_best_aln_col->best_p_delta;
        kk = g_best_aln_col->best_p_q_base;
        if (i != -1) g_best_aln_col = backbone[i].delta[j].links + kk;
    }

    fc_reverse_array(char, kstr_data(data->cns_seq), kstr_size(data->cns_seq));
    fc_reverse_array(cns_pos_t, kv_data(data->t_pos), kv_size(data->t_pos));     
}

static void
dump_one_result(FILE* out, 
    pthread_mutex_t* out_lock, 
    const char* hdr, const char* seq,
    const int seq_size)
{
    pthread_mutex_lock(out_lock);
    FWRITE(hdr, 1, strlen(hdr), out);
    fprintf(out, "\n");
    FWRITE(seq, 1, seq_size, out);
    fprintf(out, "\n");
    pthread_mutex_unlock(out_lock);
}

#define IS_LOWER(c) ((c) >= 'a' && (c) <= 'z')
#define IS_UPPER(c) ((c) >= 'A' && (c) <= 'Z')

void
dump_cns_seq(CnsData* data,
    const char* hdr,
    u8* template,
    const int template_size,
    const int min_size,
    const int split_reads,
    FILE* cns_out,
    FILE* raw_out,
    pthread_mutex_t* out_lock)
{
    const char* cns_seq = kstr_data(data->cns_seq);
    const int* t_pos = kv_data(data->t_pos);
    const size_t n = kstr_size(data->cns_seq);
    size_t i = 0;
    for (i = 0; i < template_size; ++i) template[i] = "acgt"[template[i]];
    i = 0;
    char header[1000];

    if (split_reads) {
        int last_raw_idx = 0;
        while (i < n) {
            while (i < n && IS_LOWER(cns_seq[i])) ++i;
            size_t j = i + 1;
            while (j < n && IS_UPPER(cns_seq[j])) ++j;
            if (j - i >= min_size) {
                int curr_raw_idx = t_pos[i];
                if (curr_raw_idx - last_raw_idx >= min_size) {
                    sprintf(header, ">%s_(%d_%d_%d_%d)", 
                        hdr, last_raw_idx, curr_raw_idx, curr_raw_idx - last_raw_idx, template_size);
                    const char* tmp = (const char*)(template + last_raw_idx);
                    dump_one_result(raw_out, out_lock, header, tmp, curr_raw_idx - last_raw_idx);
                }
                sprintf(header, ">%s_(%d_%d_%d_%d)", 
                    hdr, curr_raw_idx, t_pos[j-1], (int)(j - i), template_size);
                dump_one_result(cns_out, out_lock, header, cns_seq + i, (int)(j - i));
                last_raw_idx = t_pos[j - 1] + 1;
            }
            i = j;
        }
        if (template_size - last_raw_idx >= min_size) {
            sprintf(header, ">%s_(%d_%d_%d_%d)", 
                hdr, last_raw_idx, template_size, template_size - last_raw_idx, template_size);
            const char* tmp = (const char*)(template + last_raw_idx);
            dump_one_result(raw_out, out_lock, header, tmp, template_size - last_raw_idx);
        }
    } else {
        new_kstring(out_cns_seq);
        int last_raw_idx = 0;
        while (i < n) {
            while (i < n && IS_LOWER(cns_seq[i])) ++i;
            if (i >= n) break;
            size_t j = i + 1;
            while (j < n && IS_UPPER(cns_seq[i])) ++j;
            if (j - i >= min_size) {
                int curr_raw_idx = t_pos[i];
                const char* tmp = (const char*)(template + last_raw_idx);
                if (curr_raw_idx > last_raw_idx) kputsn(tmp, curr_raw_idx - last_raw_idx, &out_cns_seq);
                kputsn(cns_seq + i, j - i, &out_cns_seq);
                last_raw_idx = t_pos[j - 1] + 1;
            }
            i = j;
        }
        const char* tmp = (const char*)(template + last_raw_idx);
        if (template_size > last_raw_idx) kputsn(tmp, template_size - last_raw_idx, &out_cns_seq);
        sprintf(header, ">%s_(%d_%d_%d_%d)", hdr, 0, template_size, (int)kstr_size(out_cns_seq), template_size);
        dump_one_result(cns_out, out_lock, header, kstr_data(out_cns_seq), kstr_size(out_cns_seq));
        free_kstring(out_cns_seq);
    }
}
