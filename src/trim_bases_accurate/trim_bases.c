#include "../common/ontcns_aux.h"
#include "range_list.h"
#include "../klib/kseq.h"

#include <assert.h>
#include <stdio.h>

KSEQ_DECLARE(gzFile)

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s reads clipped_ranges output\n", prog);
}

void
load_clipped_ranges(const char* path, vec_ClippedRange* clipped_ranges)
{
	gzFile in = safe_gzopen(__FILE__, __LINE__, path, "r");
	
	char line[1024];
	ClippedRange r;
	int id;
	while (gzgets(in, line, 1024)) {
		SAFE_SCANF(sscanf, line, 4, "%d%d%d%d", &id, &r.left, &r.right, &r.size);
		kv_push(ClippedRange, *clipped_ranges, r);
	}

	safe_gzclose(__FILE__, __LINE__, in);
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* reads_path = argv[1];
	const char* clipped_ranges_path = argv[2];
	const char* output = argv[3];
	
	new_kvec(vec_ClippedRange, clipped_ranges);
	load_clipped_ranges(clipped_ranges_path, &clipped_ranges);
	
	DFOPEN(out, output, "w");
	DGZ_OPEN(in, reads_path, "r");
	kseq_t* read = kseq_init(in);
	int id = -1;
	int out_id = 0;
	while (kseq_read(read) >= 0) {
		++id;
		ClippedRange range = kv_A(clipped_ranges, id);
		if (range.size == 0) continue;
		int s = kstr_size(read->seq);
		assert(s == range.size);
		fprintf(out, ">%d\n", out_id + 1);
		++out_id;
		for (int i = range.left; i < range.right; ++i) {
			char c = kstr_A(read->seq, i);
			fprintf(out, "%c", c);
		}
		fprintf(out, "\n");
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
	FCLOSE(out);
}
