#include "../common/ontcns_aux.h"
#include "../klib/kseq.h"
#include <stdio.h>

KSEQ_DECLARE(gzFile)

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s input output\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
		print_usage(argv[0]);
		return 1;
	}

	const char* input = argv[1];
	const char* output = argv[2];
	DGZ_OPEN(in, input, "r");
	kseq_t* read = kseq_init(in);
	DFOPEN(out, output, "w");
	int id = 1;
	while (kseq_read(read) >= 0) {
		fprintf(out, ">%d\n", id);
		++id;
		for (size_t i = 0; i < kstr_size(read->seq); ++i) {
			char c = kstr_A(read->seq, i);
			fprintf(out, "%c", c);
		}
		fprintf(out, "\n");
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
	FCLOSE(out);
	return 0;
}
