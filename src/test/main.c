#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"
#include "../common/makedb_aux.h"
#include "../common/oc_assert.h"
#include "../common/gapped_candidate.h"
#include "../common/m4_record.h"

#include "../klib/khash.h"
#include "../klib/kvec.h"

#include "../edlib/edlib_wrapper.h"
#include "../gapped_align/oc_daligner.h"
#include "../gapped_align/oc_aligner.h"

#include <assert.h>

int main(int argc, char* argv[])
{
	assert(argc == 2);
	while (1) {
        int id = 0;
		DGZ_OPEN(in, argv[1], "r");
		kseq_t* read = kseq_init(in);
		while (kseq_read(read) >= 0) {
			idx bps = kstr_size(read->seq);
            if (bps == 37117) { OC_LOG("%d\t%d", id, bps); break; }
            ++id;
		}
		GZ_CLOSE(in);
		kseq_destroy(read);
        break;
	}
}
