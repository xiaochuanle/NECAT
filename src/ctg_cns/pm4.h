#ifndef NECAT_PM4_H
#define NECAT_PM4_H

#include "../common/m4_record.h"

void
part_m4(const char* wrk_dir,
    const char* m4_path,
    const int n_ctg,
    const int num_dumpped_files,
    const int num_threads);

#endif // NECAT_PM4_H