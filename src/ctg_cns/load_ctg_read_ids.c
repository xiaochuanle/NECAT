#include "load_ctg_read_ids.h"

#include "../common/oc_assert.h"
#include "seq_flag_aux.h"

#include <stdlib.h>
#include <ctype.h>

static void
extract_read_ids(const char* line, int* id1, int* id2)
{
    const int n = strlen(line);
    int i = 0;
    while (i < n && (!isspace(line[i]))) ++i;
    ++i;
    oc_assert(i < n);
    while (i < n && (!isdigit(line[i]))) ++i;
    oc_assert(i < n);
    *id1 = atoi(line + i);
    while (i < n && line[i] != '~') ++i;
    ++i;
    oc_assert(i < n);
    oc_assert(isdigit(line[i]));
    *id2 = atoi(line + i);
}

u8* load_ctg_read_ids(const char* path, const int num_reads)
{
    int num_bytes = (num_reads + 7) / 8;
    u8* a = (u8*)calloc(num_bytes, 1);
    char line[HBN_MAX_PATH_LEN];
    DFOPEN(in, path, "r");
    int id1, id2;
    int cnt = 0;
    while (fgets(line, HBN_MAX_PATH_LEN, in)) {
        extract_read_ids(line, &id1, &id2);
        set_seq_flag(a, id1);
        set_seq_flag(a, id2);
    }
    FCLOSE(in);
    return a;
}