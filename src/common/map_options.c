#include "map_options.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

static const char*
argn_list = "k:z:q:b:s:n:a:d:e:m:t:j:u:i:";

const MapOptions
sDefaultPairwiseMapingOptions = {
	15, 	// kmer_size
	10, 	// scan_window
	500, 	// kmer count cutoff
	2000, 	// block size
	3, 		// block score cutoff
	500, 	// number of candidate
	500, 	// align size cutoff
	0.25, 	// ddf score cutoff
	0.5, 	// sequencing error
	500, 	// number of output
	1, 		// number of threads
	MAP_JOB_ALN,  // job
	0,		// binary output
	1		// hdr as id
};

const MapOptions
sDefaultReferenceMapingOptions = {
	15, 	// kmer_size
	5, 	// scan_window
	500, 	// kmer count cutoff
	1000, 	// block size
	3, 		// block score cutoff
	20, 	// number of candidate
	400, 	// align size cutoff
	0.25, 	// ddf score cutoff
	0.5, 	// sequencing error
	20, 	// number of output
	1, 		// number of threads
	MAP_JOB_ALN,  // job
	0,		// binary output
	1		// hdr as id
};

void
print_MapOptions(const MapOptions* p)
{
	FILE* out = stderr;
	fprintf(out, "-k %d ", p->kmer_size);
	fprintf(out, "-z %d ", p->scan_window);
	fprintf(out, "-q %d ", p->kmer_cnt_cutoff);
	fprintf(out, "-b %d ", p->block_size);
	fprintf(out, "-s %d ", p->block_score_cutoff);
	fprintf(out, "-n %d ", p->num_candidates);
	fprintf(out, "-a %d ", p->align_size_cutoff);
	fprintf(out, "-d %f ", p->ddfs_cutoff);
	fprintf(out, "-e %f ", p->error);
	fprintf(out, "-m %d ", p->num_output);
	fprintf(out, "-t %d ", p->num_threads);
	fprintf(out, "-j %d ", p->job);
	fprintf(out, "-u %d ", p->binary_output);
	fprintf(out, "-i %d ", p->use_hdr_as_id);
	fprintf(out, "\n");
}

void
MapOptions2String(const MapOptions* p, kstring_t* s)
{
	kstring_t* out = s;
	ksprintf(out, "-k %d ", p->kmer_size);
	ksprintf(out, "-z %d ", p->scan_window);
	ksprintf(out, "-q %d ", p->kmer_cnt_cutoff);
	ksprintf(out, "-b %d ", p->block_size);
	ksprintf(out, "-s %d ", p->block_score_cutoff);
	ksprintf(out, "-n %d ", p->num_candidates);
	ksprintf(out, "-a %d ", p->align_size_cutoff);
	ksprintf(out, "-d %f ", p->ddfs_cutoff);
	ksprintf(out, "-e %f ", p->error);
	ksprintf(out, "-m %d ", p->num_output);
	ksprintf(out, "-t %d ", p->num_threads);
	ksprintf(out, "-j %d ", p->job);
	ksprintf(out, "-u %d ", p->binary_output);
	ksprintf(out, "-i %d ", p->use_hdr_as_id);
}

BOOL
parse_MapOptions(int argc, char* argv[], MapOptions* options)
{
	int c;
	while ((c = getopt(argc, argv, argn_list)) != -1) {
		switch (c) {
			case 'k':
				options->kmer_size = atoi(optarg);
				break;
			case 'z':
				options->scan_window = atoi(optarg);
				break;
			case 'q':
				options->kmer_cnt_cutoff = atoi(optarg);
				break;
			case 'b':
				options->block_size = atoi(optarg);
				break;
			case 's':
				options->block_score_cutoff = atoi(optarg);
				break;
			case 'n':
				options->num_candidates = atoi(optarg);
				break;
			case 'a':
				options->align_size_cutoff = atoi(optarg);
				break;
			case 'd':
				options->ddfs_cutoff = atof(optarg);
				break;
			case 'e':
				options->error = atof(optarg);
				break;
			case 'm':
				options->num_output = atoi(optarg);
				break;
			case 't':
				options->num_threads = atoi(optarg);
				break;
			case 'j':
				options->job = atoi(optarg);
				break;
			case 'u':
				options->binary_output = atoi(optarg);
				break;
			case 'i':
				options->use_hdr_as_id = atoi(optarg);
				break;
			case '?':
				fprintf(stderr, "invalid option '%c'\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			case ':':
				fprintf(stderr, "argument to option '%c' is not provided\n", (char)c);
				return ARG_PARSE_FAIL;
				break;
			default:
				break;
		}
	}
	return ARG_PARSE_SUCCESS;
}

void
describe_MapOptions(const MapOptions* defaultOptions)
{
	FILE* out = stderr;
	
	fprintf(out, "-k <Integer>\tkmer size\n");
	fprintf(out, "-z <Integer>\tscan window size\n");
	fprintf(out, "-q <Integer>\tkmer occurs > q times will be ignored\n");
	fprintf(out, "-b <Integer>\tblock size\n");
	fprintf(out, "-n <Integer>\tnumber of candidates\n");
	fprintf(out, "-a <Integer>\tmin align length\n");
	fprintf(out, "-d <Real>\tddf score cutoff\n");
	fprintf(out, "-e <Real>\tsequencing error\n");
	fprintf(out, "-m <Integer>\tnumber of output\n");
	fprintf(out, "-t <Integer>\tnumber of cpu threads\n");
	fprintf(out, "-j <0 or 1>\tjob: 0 = find candidates only, 1 = perform alignemnt\n");
	fprintf(out, "-u <0 or 1>\toutput binary results: 0 = no, 1 = yes\n");
	fprintf(out, "-i <0 or 1>\tuse header as sequence id: 0 = no, 1 = yes\n");
	
	fprintf(out, "\n");
	fprintf(out, "DEFAULT OPTIONS:\n");
	print_MapOptions(defaultOptions);
}
