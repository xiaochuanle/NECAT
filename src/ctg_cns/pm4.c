#include "pm4.h"

#include "../common/ontcns_aux.h"
#include "../common/oc_assert.h"

#include <pthread.h>

void make_partition_name(const char* wrk_dir, const int pid, char path[])
{
    sprintf(path, "%s/p%d", wrk_dir, pid);
}

void dump_num_partitions(const char* wrk_dir, const int np)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/np.txt", wrk_dir);
    DFOPEN(out, path, "w");
    fprintf(out, "%d\n", np);
    FCLOSE(out);
}

int load_num_partitions(const char* wrk_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/np.txt", wrk_dir);
    DFOPEN(in, path, "r");
    int np;
    SAFE_SCANF(fscanf, in, 1, "%d", &np);
    FCLOSE(in);
    return np;
}

typedef struct {
    FILE** out_list;
    int n;
} RecordWriter;

RecordWriter*
can_writer_new(const char* wrk_dir, const int pid_from, const int pid_to)
{
    RecordWriter* w = (RecordWriter*)malloc( sizeof(RecordWriter) );
    int n = pid_to - pid_from;
    w->out_list = (FILE**)malloc( sizeof(FILE*) * n );
    w->n = n;
    char path[HBN_MAX_PATH_LEN];
    for (int i = pid_from; i < pid_to; ++i) {
        make_partition_name(wrk_dir, i, path);
        int fid = i - pid_from;
        FOPEN(w->out_list[fid], path, "wb");
    }
    return w;
}

RecordWriter*
can_writer_free(RecordWriter* w)
{
    for (int i = 0; i < w->n; ++i) FCLOSE(w->out_list[i]);
    free(w->out_list);
    free(w);
    return NULL;
}

typedef struct {
    FILE*                       in;
    pthread_mutex_t*            in_lock;
    RecordWriter*               w;
    pthread_mutex_t*            w_lock;
    int                         min_ctg_id;
    int                         max_ctg_id;
} PartRecordData;

static size_t
load_records(void* array, const size_t size, const size_t count, FILE* in, pthread_mutex_t* in_lock)
{
    size_t n = 0;
    pthread_mutex_lock(in_lock);
    n = fread(array, size, count, in);
    pthread_mutex_unlock(in_lock);
    return n;
}

static void*
pcan_worker(void* param)
{
    PartRecordData* data = (PartRecordData*)(param);
    const size_t S = U64_ONE * 256 * 1024 * 1024;
    const size_t N = S / sizeof(M4Record);
    M4Record* a = (M4Record*)malloc( N * sizeof(M4Record) );
    new_kvec(vec_size_type, idx_range);

    size_t n, m;
    while ((n = load_records(a, sizeof(M4Record), N, data->in, data->in_lock))) {
        m = 0;
        for (size_t i = 0; i < n; ++i) {
            M4Record* e = a + i;
            int sid = e->sid;
            if (sid >= data->min_ctg_id && sid < data->max_ctg_id) {
                //DUMP_ASM_M4(fprintf, stderr, *e);
                if (e->sdir == REV) {
                    e->sdir = FWD;
                    e->qdir = 1 - e->qdir;
                } else {
                    oc_assert(e->sdir == FWD);
                }
                a[m++] = a[i];
            }
        }
        if (m == 0) continue;

        ks_introsort_M4Record_SidLT(m, a);
        kv_clear(idx_range);
        kv_push(size_t, idx_range, 0);
        size_t i = 0;
        n = m;
        while (i < n) {
            M4Record* e = a + i;
            const int sid = e->sid;
            size_t j = i + 1;
            while (j < n) {
                e = a + j;
                int id = e->sid;
                if (id >= sid) break;
                ++j;
            }
            kv_push(size_t, idx_range, j);
            i = j;
        }
        oc_assert(kv_back(idx_range) == n);

        pthread_mutex_lock(data->w_lock);
        for (i = 0; i < kv_size(idx_range) - 1; ++i) {
            size_t from = kv_A(idx_range, i);
            size_t to = kv_A(idx_range, i + 1);
            m = to - from;
            M4Record* e = a + from;
            int sid =e->sid;
            int fid = sid - data->min_ctg_id;
            oc_assert(fid < data->w->n);
            FWRITE(e, sizeof(M4Record), m, data->w->out_list[fid]);
        }
        pthread_mutex_unlock(data->w_lock);
    }

    free(a);
    kv_destroy(idx_range);
    return NULL;
}

void
part_m4(const char* wrk_dir,
    const char* m4_path,
    const int n_ctg,
    const int num_dumpped_files,
    const int num_threads)
{
    pthread_mutex_t in_lock;
    pthread_mutex_init(&in_lock, NULL);
    pthread_mutex_t out_lock;
    pthread_mutex_init(&out_lock, NULL);
    pthread_t job_ids[num_threads];

    for (int fid = 0; fid < n_ctg; fid += num_dumpped_files) {
        int sfid = fid;
        int efid = OC_MIN(sfid + num_dumpped_files, n_ctg);
        RecordWriter* w = can_writer_new(wrk_dir, sfid, efid);
        DFOPEN(record_in, m4_path, "rb");
        PartRecordData can_data = {
            record_in,
            &in_lock,
            w,
            &out_lock,
            sfid,
            efid
        };
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(job_ids + i, NULL, pcan_worker, &can_data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(job_ids[i], NULL);
        }
        FCLOSE(record_in);
        can_writer_free(w);
    }
}
