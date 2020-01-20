#ifndef PCAN_H
#define PCAN_H

#include "pcan_options.h"

void
pcan_main(PcanOptions* options,
			const char* wrk_dir,
			const char* can_path);

#endif // PCAN_H
