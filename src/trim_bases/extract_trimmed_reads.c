#include "largest_cover_range.h"
#include "../common/ontcns_aux.h"
#include "../common/packed_db.h"

static void
print_usage(const char* prog)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s lcrv_path input_reads all_ovlps complete_reads trimmed_reads complete_ovlps\n", prog);
}

static void
add_one_header(kstring_t* headers, kstring_t* line, vec_size_type* offsets)
{
    size_t offset = kstr_size(*headers);
    kv_push(size_t, *offsets, offset);
    size_t n = 0;
    for (; n < kstr_size(*line); ++n) {
        if (isspace(kstr_A(*line, n))) break;
    }
    kputsn(kstr_str(*line), n, headers);
    kputc('\0', headers);
}

static const char*
extract_header(kstring_t* headers, vec_size_type* offsets, int id)
{
    oc_assert(id < kv_size(*offsets));
    const char* name = kstr_str(*headers) + kv_A(*offsets, id);
    return name;
}

static void
dump_complete_m4s(LargestCoverRange* lcrv,
    kstring_t* headers,
    vec_size_type* hdr_offsets,
    const char* input,
    const char* output)
{
    M4Record m4;
    DFOPEN(in, input, "rb");
    DFOPEN(out, output, "w");
    char line[2048];
    while (fread(&m4, sizeof(M4Record), 1, in)) {
        int r = lcr_is_complete(&lcrv[m4.qid]) && lcr_is_complete(&lcrv[m4.sid]);
        if (!r) continue;
oc_assert(m4.qsize == lcrv[m4.qid].size);
oc_assert(m4.ssize == lcrv[m4.sid].size);
        const char* qhdr = extract_header(headers, hdr_offsets, m4.qid);
        const char* shdr = extract_header(headers, hdr_offsets, m4.sid);
int qid = atoi(qhdr);
int sid = atoi(shdr);
if (qid != m4.qid) {
    DUMP_M4_RECORD(fprintf, stdout, m4);
    printf("%s\n", qhdr);
    exit(0);
}
if (sid != m4.sid) {
    DUMP_M4_RECORD(fprintf, stdout, m4);
    printf("%s\n", shdr);
    exit(0);
}
        DUMP_ASM_M4_HDR_ID(sprintf, line, m4);
        FWRITE(line, 1, strlen(line), out);
    }
    FCLOSE(in);
    FCLOSE(out);
}

int main(int argc, char* argv[])
{
    if (argc != 7) {
        print_usage(argv[0]);
        return 1;
    }
    const char* lcrv_path = argv[1];
    const char* input_reads = argv[2];
    const char* all_ovlps = argv[3];
    const char* complete_reads = argv[4];
    const char* trimmed_reads = argv[5];
    const char* complete_ovlps = argv[6];
    new_kvec(vec_lcr, v_lcr);
    load_lcrs(lcrv_path, &v_lcr);
    size_t num_reads = kv_size(v_lcr);
    LargestCoverRange* lcrv = kv_data(v_lcr);
    int read_id = 0;
    DFOPEN(complete_out, complete_reads, "w");
    DFOPEN(trimmed_out, trimmed_reads, "w");
    new_kstring(headers);
    new_kvec(vec_size_type, hdr_offsets);
    kv_push(size_t, hdr_offsets, 0);
    DGZ_OPEN(in, input_reads, "r");
    kseq_t* read = kseq_init(in);
    while (1) {
        int r = kseq_read(read);
        if (r == -1) break;
        ++read_id;
        add_one_header(&headers, &read->name, &hdr_offsets);
        oc_assert(read_id <= num_reads);
        int from = 0, to = 0;
        if (!lcr_is_valid(lcrv[read_id])) continue;
        oc_assert(kstr_size(read->seq) == lcrv[read_id].size, "id: %d, %d [%d, %d, %d]", 
            read_id, kstr_size(read->seq), lcrv[read_id].left, lcrv[read_id].right, lcrv[read_id].size);;
        FILE* out = NULL;
        if (lcr_is_complete(&lcrv[read_id])) {
            out = complete_out;
            from = 0;
            to = lcrv[read_id].size;
        } else {
            out = trimmed_out;
            from = lcrv[read_id].left;
            to = lcrv[read_id].right;
        }
        oc_assert(from >= 0 && from < to && to <= lcrv[read_id].size);
        int n = to - from;
        oc_assert(out);
        fprintf(out, ">");
        FWRITE(kstr_str(read->name), 1, kstr_size(read->name), out);
        fprintf(out, "\n");
        FWRITE(kstr_str(read->seq) + from, 1, n, out);
        fprintf(out, "\n");
    }

    dump_complete_m4s(lcrv, &headers, &hdr_offsets, all_ovlps, complete_ovlps);
    FCLOSE(complete_out);
    FCLOSE(trimmed_out);
    kseq_destroy(read);
    free_kstring(headers);
    kv_destroy(hdr_offsets);
    kv_destroy(v_lcr);
}
