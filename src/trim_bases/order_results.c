#include <stdio.h>

#include "../common/ontcns_aux.h"
#include "../common/m4_record.h"
#include "../klib/kseq.h"
#include "../common/packed_db.h"
#include "../common/oc_assert.h"

static void 
print_usage(const char* prog)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s input_reads input_m4 output_reads output_m4\n", prog);
}

static int
get_number_of_reads(const char* path)
{
    DGZ_OPEN(in, path, "r");
    kseq_t* read = kseq_init(in);
    int max_id = -1;
    while (1) {
        int r = kseq_read(read);
        if (r == -1) break;
        kputc('\0', &read->name);
        int oid = atoi(kstr_str(read->name));
        max_id = OC_MAX(max_id, oid);
    }
    kseq_destroy(read);
    return max_id + 1;
}

static int*
order_reads(const char* input_reads, const char* output_reads)
{
    int num_reads = get_number_of_reads(input_reads);
    int* id_maps = (int*)malloc(sizeof(int) * (num_reads+2));
    for (int i = 0; i <= num_reads; ++i) {
        id_maps[i] = -1;
    }

    DGZ_OPEN(in, input_reads, "r");
    kseq_t* read = kseq_init(in);
    DFOPEN(out, output_reads, "w");
    int id = 1;
    while (1) {
        int r = kseq_read(read);
        if (r == -1) break;
        kputc('\0', &read->name);
        int oid = atoi(kstr_str(read->name));
        oc_assert(oid < num_reads);
        id_maps[oid] = id;
        fprintf(out, ">%d\n", id);
        ++id;
        FWRITE(kstr_str(read->seq), 1, kstr_size(read->seq), out);
        fprintf(out, "\n");
    }
    kseq_destroy(read);
    FCLOSE(out);
    return id_maps;    
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return 1;
    }

    const char* input_reads = argv[1];
    const char* input_m4 = argv[2];
    const char* output_reads = argv[3];
    const char* output_m4 = argv[4];

    int* id_maps = order_reads(input_reads, output_reads);
    new_kstring(line);
    DFOPEN(in, input_m4, "r");
    M4Record m4;
    DFOPEN(out, output_m4, "w");
    while (1) {
        kstr_clear(line);
        if (kgetline(&line, fgets, in) == EOF) break;
        kputc('\0', &line);
        LOAD_ASM_M4(sscanf, kstr_str(line), m4);
        oc_assert(id_maps[m4.qid] != -1);
if (id_maps[m4.sid] == -1) DUMP_ASM_M4(fprintf, stdout, m4);
        oc_assert(id_maps[m4.sid] != -1, "%d", m4.sid);
        m4.qid = id_maps[m4.qid];
        m4.sid = id_maps[m4.sid];
        DUMP_ASM_M4(fprintf, out, m4);        
    }
    free_kstring(line);
    free(id_maps);
    FCLOSE(out);
}
