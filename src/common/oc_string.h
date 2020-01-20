#ifndef OC_STRING_H
#define OC_STRING_H

#include "ontcns_aux.h"

typedef char* OcString;

typedef struct {
    size_t len;
    size_t free;
    char buf[];
}OcsHdr;

#define ocs_hdr(ocs) ((OcsHdr*)(void*)((ocs) - sizeof(OcsHdr)))

static inline size_t
ocs_len(const OcString s)
{
    OcsHdr* sh = ocs_hdr(s);
    return sh->len;
}

static inline size_t
ocs_avail(const OcString s)
{
    OcsHdr* sh = ocs_hdr(s);
    return sh->free;
}

#endif // OC_STRING_H
