#ifndef DETECT_CHIMERIC_READS_H
#define DETECT_CHIMERIC_READS_H

#include "../common/ontcns_aux.h"
#include "../common/m4_record.h"

#ifdef __cplusplus
extern "C" {
#endif

int is_complete_read(M4Record* m4v,
                    int nm4,
                    const double min_ident_perc,
                    int* fbgn,
                    int* fend);

int is_chimeric_read(M4Record* m4v,
                    int nm4,
                    const double min_ident_perc,
                    int* fbgn,
                    int* fend);

#ifdef __cplusplus
}
#endif
#endif // DETECT_CHIMERIC_READS_H
