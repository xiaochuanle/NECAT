#ifndef CNS_ONE_CTG_H
#define CNS_ONE_CTG_H

#include "../common/packed_db.h"

#ifdef __cplusplus
extern "C" {
#endif

void
cns_one_ctg(const char* mkdb_dir,
    PackedDB* reads,
    const char* contig,
    const idx ctg_size,
    const int ctg_id,
    kstring_t* cns_ctg);

#ifdef __cplusplus
}
#endif
#endif // CNS_ONE_CTG_H