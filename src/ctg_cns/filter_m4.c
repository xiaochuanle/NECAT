#include "../common/m4_record.h"
#include "seq_flag_aux.h"
#include "load_ctg_read_ids.h"
#include "../common/makedb_aux.h"

#include <stdio.h>

static const int kMaxLoadedM4 = 10000;

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s reads_dir ctg_read_ids m4_in m4_out\n", pn);
}

static size_t
load_records(void* array, const size_t size, const size_t count, FILE* in)
{
    size_t n = 0;
    n = fread(array, size, count, in);
    return n;
}

static int
load_next_m4(M4Record* a, int* next_i, int *max_i, FILE* in, M4Record* m4)
{
    int i = *next_i;
    int m = *max_i;
    if (i >= m) {
        m = load_records(a, sizeof(M4Record), kMaxLoadedM4, in);
        if (m == 0) return 0;
        i = 0;
        *max_i = m;
    }
    *m4 = a[i];
    //DUMP_ASM_M4(fprintf, stderr, a[i]);
    ++i;
    *next_i = i;
    return 1;
}

static int
load_next_m4v(M4Record* a, int* next_i, int *max_i, FILE* in, vec_m4* m4_list)
{
    kv_clear(*m4_list);
    M4Record m4;
    if (!load_next_m4(a, next_i, max_i, in, &m4)) return 0;
    int qid = m4.qid;
    kv_push(M4Record, *m4_list, m4);
    while (1) {
        if (!load_next_m4(a, next_i, max_i, in, &m4)) break;
        if (m4.qid != qid) {
            --(*next_i);
            break;
        }
        kv_push(M4Record, *m4_list, m4);
    }
    return 1;
}

static void
process_m4v(vec_m4* m4_list, const u8* ctg_read_ids, FILE* out)
{
    M4Record* m4_array = kv_data(*m4_list);
    const int n = kv_size(*m4_list);
    const int E = 100;
    int num_full_m4 = 0;
    int full_m4_idx = -1;
    for (int i = 0; i < n; ++i) {
        M4Record* m4 = m4_array + i;
        if (seq_flag_is_set(ctg_read_ids, m4->qid)) continue;
        if (m4->qoff <= E && m4->qsize - m4->qend <= E) {
            ++num_full_m4;
            full_m4_idx = i;
            continue;
        }
        if (m4->soff <= E && m4->ssize - m4->send <= E) {
            ++num_full_m4;
            full_m4_idx = i;
            continue;
        }
        idx qbeg = m4->qoff, qend = m4->qend;
        if (m4->qdir == REV) {
            qbeg = m4->qsize - m4->qend;
            qend = m4->qsize - m4->qoff;
        }
        idx sbeg = m4->soff, send = m4->send;
        if (m4->sdir == REV) {
            sbeg = m4->ssize - m4->send;
            send = m4->ssize - m4->soff;
        }
        if (m4->qsize - qend <= E && sbeg <= E) {
            if (qend - qbeg >= m4->qsize * 0.6 || send - sbeg >= m4->ssize * 0.6) {
                ++num_full_m4;
                full_m4_idx = i;
                continue;
            }
        }
        if (m4->ssize - send <= E && qbeg <= E) {
            if (qend - qbeg >= m4->qsize * 0.6 || send - sbeg >= m4->ssize * 0.6) {
                ++num_full_m4;
                full_m4_idx = i;
                continue;
            }            
        }
    }
    if (num_full_m4 == 1) {
        M4Record* m4 = m4_array + full_m4_idx;
        if (m4->ident_perc >= 90.0) {
            if (m4->sdir == REV) {
                m4->sdir = FWD;
                m4->qdir = 1 - m4->qdir;
            }
            FWRITE(m4, sizeof(M4Record), 1, out);
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return 1;
    }
    const char* reads_dir = argv[1];
    const char* ctg_read_ids_path = argv[2];
    const char* m4_input = argv[3];
    const char* m4_output = argv[4];

    M4Record* a = (M4Record*)calloc(kMaxLoadedM4, sizeof(M4Record));
    int a_i = 0;
    int max_a_i = 0;
    const int num_reads = load_num_reads(reads_dir);
    u8* ctg_read_ids = load_ctg_read_ids(ctg_read_ids_path, num_reads);
    DFOPEN(in, m4_input, "rb");
    DFOPEN(out, m4_output, "wb");
    new_kvec(vec_m4, m4v);
    while (load_next_m4v(a, &a_i, &max_a_i, in, &m4v)) {
        process_m4v(&m4v, ctg_read_ids, out);
    }

    kv_destroy(m4v);
    FCLOSE(out);
    FCLOSE(in);
    free(a);
    free(ctg_read_ids);
}
