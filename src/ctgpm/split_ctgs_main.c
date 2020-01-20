#include "split_ctgs.h"
#include <stdio.h>
#include "../klib/kseq.h"

KSEQ_DECLARE(gzFile)

void
print_usage(const char* prog)
{
	FILE* out = stderr;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s contigs contig_seqs\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
		print_usage(argv[0]);
		return 1;
	}
	
	const char* ctg_path = argv[1];
	const char* ctg_reads_path = argv[2];
	DFOPEN(out, ctg_reads_path, "w");
	DGZ_OPEN(in, ctg_path, "r");
	kseq_t* contig = kseq_init(in);
	int ctg_id = 0;
	while (kseq_read(contig) >= 0) {
		split_contig(ctg_id, &contig->seq, out);
		++ctg_id;
	}
	
	GZ_CLOSE(in);
	kseq_destroy(contig);
	FCLOSE(out);
}
