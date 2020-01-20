#include "cns_aux.h"

#include "../common/ontcns_aux.h"
#include "../common/oc_assert.h"
#include "../klib/ksort.h"

static inline u8
encode_dna_base(const char c) 
{
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
    }
    OC_ERROR("invalid dna base: %c", c);
    return 4;
}

static void
build_base_links(AlignTag* tags, const int ntag, BaseLinks* link, LinkInfoAllocator* alloc)
{
    int n_link = 0;
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && AlignTag_PLinkEq(tags[i], tags[j])) ++j;
        ++n_link;
        i = j;
    }

    link->plinks = (LinkInfo*)alloc_OcObjectAllocator(alloc, (u32)n_link);
    link->n_link = n_link;
    link->coverage = ntag;

    n_link = 0;
    i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && AlignTag_PLinkEq(tags[i], tags[j])) ++j;
        LinkInfo* linfo = link->plinks + n_link;
        linfo->p_t_pos = tags[i].p_t_pos;
        linfo->p_delta = tags[i].p_delta;
        linfo->p_q_base = tags[i].p_q_base;
        linfo->link_count = j - i;
        linfo->weight = 0;
        for (int k = i; k < j; ++k) linfo->weight += tags[k].weight;

        ++n_link;
        i = j;
    }

    oc_assert(n_link == link->n_link, "n_link = %d, link->n_link = %d", n_link, link->n_link);
}

static void
build_delta_links(AlignTag* tags, const int ntag, DeltaCovInfo* dci, LinkInfoAllocator* alloc)
{
    for (int i = 0; i < 5; ++i) {
        dci->links[i].n_link = 0;
        dci->links[i].coverage = 0;
        dci->links[i].best_p_t_pos = -1;
        dci->links[i].best_p_delta = U8_MAX;
        dci->links[i].best_p_q_base = '.';
        dci->links[i].score = 0;
    }

    int i = 0;
    while (i < ntag) {
        int j = i;
        while (j < ntag && tags[i].q_base == tags[j].q_base) ++j;
        u8 c = encode_dna_base(tags[i].q_base);
        build_base_links(tags + i, j - i, dci->links + c, alloc);
        i = j;
    }
}

static void
build_backbone_item(AlignTag* tags,
        const int ntag,
        BackboneItem* item,
        DeltaCovInfoAllocator* dci_alloc,
        LinkInfoAllocator* li_alloc,
        int* coverage)
{
    item->n_delta = tags[ntag - 1].delta + 1;
    item->delta = (DeltaCovInfo*)alloc_OcObjectAllocator(dci_alloc, (u32)item->n_delta);
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tags[i].delta == tags[j].delta) ++j;
        build_delta_links(tags + i, j - i, item->delta + tags[i].delta, li_alloc);
        if (tags[i].delta == 0) {
            coverage[ tags[i].t_pos ] = j - i;
        }
        i = j;
    }
}

void
build_backbone(AlignTag* tags,
        const int ntag,
        const int template_size,
        DeltaCovInfoAllocator* dci_alloc,
        LinkInfoAllocator* li_alloc,
        vec_backbone_item* backbone,
        vec_int* coverage)
{
    kv_resize(BackboneItem, *backbone, (size_t)template_size);
    for (int i = 0; i < template_size; ++i) clear_BackboneItem(kv_A(*backbone, i));
    kv_resize(int, *coverage, (size_t)template_size);
    kv_zero(int, *coverage);

    ks_introsort(AlignTag, (size_t)ntag, tags);
    int i = 0;
    while (i < ntag) {
        int j = i + 1;
        while (j < ntag && tags[i].t_pos == tags[j].t_pos) ++j;
        oc_assert(tags[i].t_pos < template_size, "i = %d, j = %d, t_pos = %d, template_size = %d", i, j, tags[i].t_pos, template_size);
        build_backbone_item(tags + i, j - i, kv_data(*backbone) + tags[i].t_pos, dci_alloc, li_alloc, kv_data(*coverage));
        i = j;
    }
}

void
consensus_backbone_segment(BackboneItem* backbone,
        int from,
        int to,
        int* coverage,
        kstring_t* cns_seq,
        int* cns_from,
        int* cns_to)
{
    int g_best_ck = 0;
    BaseLinks* g_best_aln_col = 0;
    int g_best_t_pos = 0;
    double g_best_score = -1.0;
    int g_best_q_base = -1;
    int kk;
    int ck;
    int best_i;
    int best_j;
    int best_b;
    int best_ck = -1;
    double score;
    double best_score;
    BaseLinks* aln_col;

    for (int i = from; i < to; ++i) {
        for (int j = 0; j < backbone[i].n_delta; ++j) {
            for (kk = 0; kk < 5; ++kk) {
                aln_col = backbone[i].delta[j].links + kk;
                if (aln_col->coverage) {
                    best_score = -1;
                    best_i = -1;
                    best_j = -1;
                    best_b = -1;

                    for (ck = 0; ck < aln_col->n_link; ++ck) {
                        int pi = aln_col->plinks[ck].p_t_pos;
                        int pj = aln_col->plinks[ck].p_delta;
                        int pkk = encode_dna_base(aln_col->plinks[ck].p_q_base);
                        if (j) score = aln_col->plinks[ck].weight - 0.4 * 0.5 * coverage[i];
                        else score = aln_col->plinks[ck].weight - 0.4 * 0.5 * coverage[i];
                        if (pi != -1) {
                            score += backbone[pi].delta[pj].links[pkk].score;
                        }

                        if (score > best_score) {
                            best_score = score;
                            aln_col->best_p_t_pos = best_i = pi;
                            aln_col->best_p_delta = best_j = pj;
                            aln_col->best_p_q_base = best_b = pkk;
                            best_ck = ck;
                        }
                    }

                    aln_col->score = best_score;
                    if (best_score > g_best_score) {
                        g_best_score = best_score;
                        g_best_aln_col = aln_col;
                        g_best_ck = best_ck;
                        g_best_t_pos = i;
                        g_best_q_base = kk;
                    }
                }
            }
        }
    }
    oc_assert(g_best_score != -1, "g_best_score = %d", g_best_score);

    char bb = '$';
    ck = g_best_q_base;
    int i = g_best_t_pos;
    int _cns_to = i + 1, _cns_from = 0;
    int j;
    kstr_clear(*cns_seq);
    while (1) {
        bb = ck;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = backbone[i].delta[j].links + ck;
        _cns_from = i;

        if (bb != 4) {
            oc_assert(bb >= 0 && bb < 4);
            kputc(bb, cns_seq);
        }
    }
	
	reverse_kstring(cns_seq);
	if (cns_from) *cns_from = _cns_from;
	if (cns_to) *cns_to = _cns_to;
}
