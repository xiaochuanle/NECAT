#ifndef MAP_OPTIONS_H
#define MAP_OPTIONS_H

#include "ontcns_defs.h"
#include "../klib/kstring.h"

#define MAP_JOB_CAN		0
#define MAP_JOB_ALN		1

typedef struct {
    int     kmer_size;
    int     scan_window;
    int     kmer_cnt_cutoff;
    int     block_size;
    int     block_score_cutoff;
    int     num_candidates;
    int     align_size_cutoff;
	double 	ddfs_cutoff;
    double  error;
    int     num_output;
    int     num_threads;
	int		job;
	int		binary_output;
	int		use_hdr_as_id;
} MapOptions;

extern const MapOptions sDefaultPairwiseMapingOptions;

extern const MapOptions sDefaultReferenceMapingOptions;

void
print_MapOptions(const MapOptions* p);

BOOL
parse_MapOptions(int argc, char* argv[], MapOptions* options);

void
describe_MapOptions(const MapOptions* defaultOptions);

void
MapOptions2String(const MapOptions* options, kstring_t* s);

#endif // MAP_OPTIONS_H
