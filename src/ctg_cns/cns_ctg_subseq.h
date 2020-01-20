#ifndef CNS_CTG_SUBSEQ_H
#define CNS_CTG_SUBSEQ_H

#include "../common/packed_db.h"
#include "../common/m4_record.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_CNS_COV 10
#define MIN_CNS_COV 4
#define CNS_MEM_KMER_SIZE 15
#define CNS_MEM_WINDOW_SIZE 6
#define CNS_MEM_MEM_SIZE 20

void
cns_ctg_subseq(PackedDB* reads,
    u8* ctg_subseq,
    const idx ctg_subseq_offset,
    const int ctg_subseq_size,
    M4Record* m4_array,
    const int m4_count,
    kstring_t* cns_ctg);

#ifdef __cplusplus
}
#endif
#endif // CNS_CTG_SUBSEQ_H