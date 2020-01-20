#ifndef SEQ_FLAG_AUX_H
#define SEQ_FLAG_AUX_H

#include "../common/ontcns_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline void
set_seq_flag(u8* flags, const int id)
{
    const u8 flag_table[8] = {1,2,4,8,16,32,64,128};
    flags[id>>3] |= flag_table[id&7];
}

static inline int
seq_flag_is_set(const u8* flags, const int id)
{
    const u8 flag_table[8] = {1,2,4,8,16,32,64,128};
    int r = flags[id>>3] & flag_table[id&7];
    //HBN_LOG("id =  %d, r = %d", id, r);
    return r;
}

static inline u8*
load_reads_flags(const char* mkdb_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/seq_status.bin", mkdb_dir);
    size_t size = FILE_SIZE(path);
    if (!size) return NULL;
    u8* flags = (u8*)calloc(size, 1);
    DFOPEN(in, path, "rb");
    FREAD(flags, 1, size, in);
    FCLOSE(in);
    return flags;
}

#ifdef __cplusplus
}
#endif
#endif // SEQ_FLAG_AUX_H
