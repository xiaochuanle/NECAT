#include "oc_string.h"

#define OCS_MAX_PREALLOC (1024 * 1024)

OcString
ocs_new_len(const void* init, size_t init_len)
{
    OcsHdr* sh;
    if (init) {
        sh = (OcsHdr*)malloc(sizeof(OcsHdr) + init_len + 1);
    } else {
        sh = (OcsHdr*)calloc(1, sizeof(OcsHdr) + init_len + 1);
    }

    if (sh == NULL) return NULL;
    sh->len = init_len;
    sh->free = 0;
    if (init_len && init) memcpy(sh->buf, init, init_len);
    sh->buf[init_len] = '\0';
    return sh->buf;
}

OcString
ocs_empty(void)
{
    return ocs_new_len("", 0);
}

OcString
ocs_new(const char* init)
{
    size_t init_len = (init == NULL) ? 0 : strlen(init);
    return ocs_new_len(init, init_len);
}

OcString
ocs_free(OcString s)
{
    if (s) {
        free(ocs_hdr(s));
    }
    return 0;
}

void
ocs_clear(OcString s)
{
    OcsHdr* sh = ocs_hdr(s);
    sh->free += sh->len;
    sh->len = 0;
    sh->buf[0] = '\0';
}

OcString
ocs_make_room_for(OcString s, size_t add_len)
{
    OcsHdr *sh, *newsh;
    size_t free = ocs_avail(s);
    size_t len, newlen;

    if (free >= add_len) return s;
    len = ocs_len(s);
    sh = ocs_hdr(s);
    newlen = len + add_len;
    if (newlen < OCS_MAX_PREALLOC) {
        newlen *= 2;
    } else {
        newlen += OCS_MAX_PREALLOC;
    }    
    newsh = (OcsHdr*)realloc(sh, sizeof(OcsHdr) + newlen + 1);
    if (newsh == NULL) return NULL;

    newsh->free = newlen - len;
    return newsh->buf;
}

OcString
ocs_puts_len(OcString s, const void* t, size_t len)
{
    OcsHdr* sh;
    size_t curlen = ocs_len(s);

    s = ocs_make_room_for(s, len);
    if (s == NULL) return NULL;
    sh = ocs_hdr(s);
    memcpy(s + curlen, t, len);
    sh->len = curlen + len;
    sh->free = sh->free - len;
    s[curlen + len] = '\0';
    return s;
}

OcString
ocs_puts(OcString s, const char* t)
{
    return ocs_puts_len(s, t, strlen(t));
}
