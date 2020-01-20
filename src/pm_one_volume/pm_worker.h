#ifndef PM_WORKER_H
#define PM_WORKER_H

#include "../common/map_options.h"

#define GroupThreadSize 8

void
pm_main(MapOptions* options, const int vid, const char* wrk_dir, const char* output);

void
pm_multi_group(MapOptions* options, const int vid, const char* wrk_dir, const char* output);

#endif // PM_WORKER_H
