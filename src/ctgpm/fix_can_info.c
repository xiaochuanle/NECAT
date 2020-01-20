#include "split_ctgs.h"

#include "../common/gapped_candidate.h"
#include "../klib/kseq.h"

#include <assert.h>

KSEQ_DECLARE(gzFile)

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE\n");
	fprintf(out, "%s reads candidates output\n", prog);
}

void
build_hdr_info(const char* reads_path, vec_hdr_info* hdr_info)
{
	DGZ_OPEN(in, reads_path, "r");
	kseq_t* read = kseq_init(in);
	CtgSeqHdrInfo hdr;
	while (kseq_read(read) >= 0) {
		extract_ctg_seq_hdr(kstr_str(read->name), &hdr.ctg_id, &hdr.seq_size, &hdr.ctg_offset, &hdr.ctg_size);
		kv_push(CtgSeqHdrInfo, *hdr_info, hdr);
	}
	GZ_CLOSE(in);
	kseq_destroy(read);
}

static void
fix_offsets(int dir, idx beg_in, idx end_in, idx off_in, 
			int seq_size, idx ctg_offset, idx ctg_size,
			idx* beg_out, idx* end_out, idx* off_out)
{
	if (dir == FWD) {
		*beg_out = beg_in + ctg_offset;
		*end_out = end_in + ctg_offset;
		*off_out = off_in + ctg_offset;
	} else {
		assert(dir == REV);
		idx beg = seq_size - end_in;
		idx end = seq_size - beg_in;
		idx off = seq_size - 1 - off_in;
		beg += ctg_offset;
		end += ctg_offset;
		off += ctg_offset;
		
		*beg_out = ctg_size - end;
		*end_out = ctg_size - beg;
		*off_out = ctg_size - 1 - off;
	}
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		print_usage(argv[0]);
		return 1;
	}
	const char* reads_path = argv[1];
	const char* can_path = argv[2];
	const char* output = argv[3];
	new_kvec(vec_hdr_info, hdr_info);
	build_hdr_info(reads_path, &hdr_info);
	
	char line[1024];
	GappedCandidate can;
	DFOPEN(out, output, "w");
	DGZ_OPEN(in, can_path, "r");
	while (gzgets(in, line, 1024)) {
		LOAD_GAPPED_CANDIDATE(sscanf, line, can);
		
		int qctg_id = kv_A(hdr_info, can.qid).ctg_id;
		int tctg_id = kv_A(hdr_info, can.sid).ctg_id;
		if (qctg_id == tctg_id) continue;
		
		fix_offsets(can.qdir, can.qbeg, can.qend, can.qoff, 
					kv_A(hdr_info, can.qid).seq_size, 
					kv_A(hdr_info, can.qid).ctg_offset, 
					kv_A(hdr_info, can.qid).ctg_size,
					&can.qbeg, &can.qend, &can.qoff);
		
		fix_offsets(can.sdir, can.sbeg, can.send, can.soff, 
					kv_A(hdr_info, can.sid).seq_size, 
					kv_A(hdr_info, can.sid).ctg_offset, 
					kv_A(hdr_info, can.sid).ctg_size,
					&can.sbeg, &can.send, &can.soff);
		
		can.qsize = kv_A(hdr_info, can.qid).ctg_size;
		can.ssize = kv_A(hdr_info, can.sid).ctg_size;
		can.qid = qctg_id;
		can.sid = tctg_id;
		
		assert(can.qoff <= can.qsize);
		assert(can.soff <= can.ssize);
		DUMP_GAPPED_CANDIDATE(fprintf, out, can);
	}
	
	GZ_CLOSE(in);
	FCLOSE(out);
	free_kvec(hdr_info);
}
