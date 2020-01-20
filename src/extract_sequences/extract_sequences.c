#include "../common/ontcns_aux.h"
#include "../klib/kseq.h"
#include "../klib/ksort.h"

#include <assert.h>
#include <stdio.h>

KSEQ_DECLARE(gzFile)

#define INT_GT(a, b) ((a) > (b))
KSORT_INIT(INT_GT, int, INT_GT)

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s reads genome_size coverage output\n", prog);
}

int
calc_length_cutoff(const char* path, const idx genome_size, const int coverage)
{
	new_kvec(vec_int, length);
	DGZ_OPEN(in, path, "r");
	kseq_t* read = kseq_init(in);
	while (kseq_read(read) >= 0) {
		int s = kstr_size(read->seq);
		kv_push(int, length, s);
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
	
	ks_introsort_INT_GT(kv_size(length), kv_data(length));
	idx target_size = genome_size * coverage;
	idx curr = 0;
	size_t i = 0;
	for (; i < kv_size(length); ++i) {
		curr += kv_A(length, i);
		if (curr >= target_size) break;
	}
	
	int ret = (i < kv_size(length)) ? kv_A(length, i) : 0;
	free_kvec(length);
	return ret;
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* reads_path = argv[1];
	const idx genome_size = atoll(argv[2]);
	const int coverage = atoi(argv[3]);
	const char* output = argv[4];
	
	int length_cutoff = calc_length_cutoff(reads_path, genome_size, coverage);
	new_kvec(vec_int, length);
	DGZ_OPEN(in, reads_path, "r");
	kseq_t* read = kseq_init(in);
	DFOPEN(out, output, "w");
	int id = 1;
	while (kseq_read(read) >= 0) {
		int s = kstr_size(read->seq);
		if (s >= length_cutoff) {
			fprintf(out, ">%d\n", id);
			++id;
			for (size_t i = 0; i < kstr_size(read->seq); ++i) {
				char c = kstr_A(read->seq, i);
				fprintf(out, "%c", c);
			}
			fprintf(out, "\n");
		}
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
	FCLOSE(out);
}
