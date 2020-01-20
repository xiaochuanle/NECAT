#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"
#include "../klib/ksort.h"
#include "../klib/kvec.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>

static void
print_usage(const char* prog)
{
	FILE* out = stderr;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s cns_reads_in raw_reads_in cns_reads_out raw_reads_out\n", prog);
}

typedef struct {
	int read_id;
	int order_in_pdb;
} CnsReadInfo;
	
typedef kvec_t(CnsReadInfo) vec_CnsReadInfo;
#define CnsReadInfoLt(a, b) ((a).read_id < (b).read_id)
KSORT_INIT(CnsReadInfo, CnsReadInfo, CnsReadInfoLt)

static int
extract_cns_info_from_hdr(const char* hdr)
{
	int id = 0;
	const char* p = hdr;
	while ( (p != NULL) && (!isdigit(*p)) ) ++p;
	assert(p != NULL);
	while ( (p != NULL) && isdigit(*p) ) {
		id = id * 10 + (*p - '0');
		++p;
	}
	return id;
}

void
load_reads(const char* path, PackedDB* pdb, vec_CnsReadInfo* cns_info)
{
	DGZ_OPEN(in, path, "r");
	kseq_t* read = kseq_init(in);
	int order_in_pdb = 0;
	while (kseq_read(read) >= 0) {
		int id = extract_cns_info_from_hdr(kstr_str(read->name));
		CnsReadInfo crin;
		crin.read_id = id;
		crin.order_in_pdb = order_in_pdb;
		++order_in_pdb;
		kv_push(CnsReadInfo, *cns_info, crin);
		pdb_add_one_seq(pdb, read, TECH_NANOPORE);
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
}

static void
dump_one_read(FILE* out, const char* hdr, kstring_t* read)
{
	fprintf(out, ">%s\n", hdr);
	for (size_t i = 0; i < kstr_size(*read); ++i) {
		int c = kstr_A(*read, i);
		fprintf(out, "%c", "ACGT"[c]);
	}
	fprintf(out, "\n");
}

static void
dump_reordered_reads(PackedDB* pdb, vec_CnsReadInfo* cns_info, const char* out_path)
{
	size_t n = kv_size(*cns_info);
	CnsReadInfo* infov = kv_data(*cns_info);
	ks_introsort_CnsReadInfo(n, infov);
	DFOPEN(out, out_path, "w");
	new_kstring(read);
	for (size_t i = 0; i < n; ++i) {
		pdb_extract_sequence(pdb, infov[i].order_in_pdb, FWD, &read);
		const char* hdr = PDB_SEQ_NAME(pdb, infov[i].order_in_pdb);
		dump_one_read(out, hdr, &read);
	}
	free_kstring(read);
	FCLOSE(out);
}

static void
reorder_reads(const char* in_path, const char* out_path)
{
	PackedDB* reads = new_PackedDB();
	new_kvec(vec_CnsReadInfo, cns_info);
	load_reads(in_path, reads, &cns_info);
	dump_reordered_reads(reads, &cns_info, out_path);
	free_PackedDB(reads);
	free_kvec(cns_info);
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		return 1;
	}
	const char* cns_reads_in = argv[1];
	const char* raw_reads_in = argv[2];
	const char* cns_reads_out = argv[3];
	const char* raw_reads_out = argv[4];

	reorder_reads(cns_reads_in, cns_reads_out);
	reorder_reads(raw_reads_in, raw_reads_out);
}
