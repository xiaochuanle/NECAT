#ifndef CNS_AUX_H
#define CNS_AUX_H

#include "align_tags.h"
#include "../common/soa.h"
#include "../klib/kvec.h"

typedef struct {
    double weight;
    int p_t_pos;
    u8 p_delta;
    char p_q_base;
    int link_count;
} LinkInfo;

#define LinkInfo_EQ(a, b) ((a).p_t_pos == (b).p_t_pos && (a).delta == (b).delta && (a).p_q_base == (b).p_q_base)

typedef OcObjectAllocator LinkInfoAllocator;

typedef struct {
    int n_link;
    int coverage;
    LinkInfo* plinks;
    int best_p_t_pos;
    u8 best_p_delta;
    u8 best_p_q_base;
    double score;
} BaseLinks;

typedef struct {
    BaseLinks links[5];
} DeltaCovInfo;

typedef OcObjectAllocator DeltaCovInfoAllocator;

typedef struct {
    int n_delta;
    DeltaCovInfo* delta;
} BackboneItem;

typedef kvec_t(BackboneItem) vec_backbone_item;

#define clear_BackboneItem(item) ((item).n_delta = 0, (item).delta = 0)

#define AlignTag_PLinkEq(a, b) ((a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base == (b).p_q_base)

void
build_backbone(AlignTag* tags,
			   const int ntag,
			   const int template_size,
			   DeltaCovInfoAllocator* dci_alloc,
			   LinkInfoAllocator* li_alloc,
			   vec_backbone_item* backbone,
			   vec_int* coverage);

void
consensus_backbone_segment(BackboneItem* backbone,
						   int from,
						   int to,
						   int* coverage,
						   kstring_t* cns_seq,
						   int* cns_from,
						   int* cns_to);

#endif // CNS_AUX_H
