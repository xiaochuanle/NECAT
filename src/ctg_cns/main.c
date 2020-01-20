#include "cns_one_ctg.h"

#include "../common/oc_assert.h"
#include "../common/makedb_aux.h"
#include "../common/nst_nt4_table.h"

#include <stdio.h>
#include <stdlib.h>

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s reads_dir reference_dir num_threads output\n", pn);
}

typedef struct {
    const char* wrk_dir;
    PackedDB* reads;
    PackedDB* contigs;
    FILE* out;
    pthread_mutex_t out_lock;
    int ctg_id;
    pthread_mutex_t ctg_id_lock;
} CtgCnsData;

CtgCnsData*
CtgCnsDataNew(const char* wrk_dir, PackedDB* reads, PackedDB* contigs, FILE* out)
{
    CtgCnsData* data = (CtgCnsData*)calloc(1, sizeof(CtgCnsData));
    data->wrk_dir = wrk_dir;
    data->reads = reads;
    data->contigs = contigs;
    data->out = out;
    pthread_mutex_init(&data->out_lock, NULL);
    data->ctg_id = 0;
    pthread_mutex_init(&data->ctg_id_lock, NULL);
    return data;
}

CtgCnsData*
CtgCnsDataFree(CtgCnsData* data)
{
    free(data);
    return NULL;
}

void*
ctg_cns_worker(void* params)
{
    CtgCnsData* data = (CtgCnsData*)(params);
    const int n_ctg = pdb_num_seqs(data->contigs);
    int ctg_id = 0;
    new_kstring(contig);
    new_kstring(cns_ctg);
    char job_name[256];
    while (1) {
        pthread_mutex_lock(&data->ctg_id_lock);
        ctg_id = data->ctg_id;
        ++data->ctg_id;
        pthread_mutex_unlock(&data->ctg_id_lock);
        if (ctg_id >= n_ctg) break;
        pdb_extract_sequence(data->contigs, ctg_id, FWD, &contig);
        sprintf(job_name, "correct contig %d, length = %zu", ctg_id, kstr_size(contig));
        TIMING_START(job_name);
        kstr_clear(cns_ctg);
        cns_one_ctg(data->wrk_dir, data->reads, kstr_data(contig), kstr_size(contig), ctg_id, &cns_ctg);
        for (size_t p = 0; p < kstr_size(cns_ctg); ++p) {
            int c = kstr_A(cns_ctg, p);
            oc_assert(isupper(c), "p = %zu, c = %c", p, c);
            c = nst_nt4_table[c];
            oc_assert(c >= 0 && c < 4);
        }
        const char* hdr = PDB_SEQ_NAME(data->contigs, ctg_id);
        pthread_mutex_lock(&data->out_lock);
        fprintf(data->out, ">%s\n", hdr);
        FWRITE(kstr_data(cns_ctg), 1, kstr_size(cns_ctg), data->out);
        fprintf(data->out, "\n");
        pthread_mutex_unlock(&data->out_lock);
        TIMING_END(job_name);
    }
    free_kstring(cns_ctg);
    free_kstring(contig);
    return NULL;
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    const char* reads_dir = argv[1];
    const char* reference_dir = argv[2];
    const int num_threads = atoi(argv[3]);
    const char* output = argv[4];
    PackedDB* reads = merge_volumes(reads_dir);
    PackedDB* reference = merge_volumes(reference_dir);
    DFOPEN(out, output, "w");
    pthread_t job_ids[num_threads];
    CtgCnsData* data = CtgCnsDataNew(reference_dir, reads, reference, out);
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(job_ids + i, NULL, ctg_cns_worker, data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(job_ids[i], NULL);
    }
    CtgCnsDataFree(data);
    FCLOSE(out);
    free_PackedDB(reference);
    free_PackedDB(reads);
}