#ifndef CONSENSUS_ONE_PARTITION_H
#define CONSENSUS_ONE_PARTITION_H

#include "../common/packed_db.h"
#include "cns_options.h"

#include <stdio.h>

void
consensus_one_partition(const char* pac_reads_dir,
						PackedDB* reads,
						const char* can_path,
						CnsOptions* options,
						FILE* cns_out,
						FILE* raw_out,
						const int pid);

#endif // CONSENSUS_ONE_PARTITION_H
