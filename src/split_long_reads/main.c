#include "../klib/kseq.h"
#include "../klib/kstring.h"
#include "../common/ontcns_aux.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

KSEQ_DECLARE(gzFile)

static const int sOvlpSize = 10000;

int
calc_next_subseq_size(int left, int subseq_size)
{
	int L = OC_MIN(left, subseq_size);
	left -= L;
	if (left <= sOvlpSize) {
		L += left;
		left = 0;
	}	
	return L;
}

void
mutate_subsequence(FILE* out, kstring_t* read, int from, int to)
{
	for (int i = from; i < to; ++i) {
		int N = rand() % 101;
		if (N < 2) { // mutate
			int T = rand() % 3;
			if (T == 0) { // mis
				T = rand() % 4;
				char c = "ACGT"[T];
				fprintf(out, "%c", c);
			} else if (T == 1) { // ins
				char c = kstr_A(*read, i);
				fprintf(out, "%c", c);
				int L = rand() % 4;
				for (int j = 0; j < L; ++j) {
					T = rand() % 4;
					c = "ACGT"[T];
					fprintf(out, "%c", c);
				}
			}	
		}
	}
}

static void
output_subsequence(FILE* out, kstring_t* read, int* read_id, int from, int to, int mutate)
{
	fprintf(out, ">%d\n", *read_id);
	++(*read_id);
	if (mutate) {
		mutate_subsequence(out, read, from, from + sOvlpSize);
		from += sOvlpSize;
	}
	for (int i = from; i < to; ++i) fprintf(out, "%c", kstr_A(*read, i));
	fprintf(out, "\n");
}

static void
split_one_read(FILE* out, kstring_t* read, int* read_id, int max_length)
{
	int read_size = kstr_size(*read);

	int left = read_size, from = 0, to = 0;
	while (left) {
		int L = calc_next_subseq_size(left, max_length);
		to = from + L;
		output_subsequence(out, read, read_id, from, to, from > 0);
		left -= L;
		from = to - sOvlpSize;
	}
}

int main(int argc, char* argv[])
{
 	srand((unsigned)time( NULL ) );
	if (argc != 5) {
		fprintf(stderr, "USAGE:\n");
		fprintf(stderr, "%s min-length max-length input output\n", argv[0]);
		return 1;
	}
	size_t min_length = atoll(argv[1]);
	size_t max_length = atoll(argv[2]);
	const char* input = argv[3];
	const char* output = argv[4];
	assert(max_length > sOvlpSize);
	int id = 1;
	DFOPEN(out, output, "w");
	
	DGZ_OPEN(in, input, "r");
    kseq_t* read = kseq_init(in);
    while (kseq_read(read) >= 0) {
		split_one_read(out, &read->seq, &id, max_length);
    }
    GZ_CLOSE(in);
    kseq_destroy(read);
	
	FCLOSE(out);
}
