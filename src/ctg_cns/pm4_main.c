#include "pm4.h"
#include "../common/makedb_aux.h"

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s mkdb_dir input_m4\n", pn);
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    const char* mkdb_dir = argv[1];
    const char* m4_path = argv[2];
    int n_ctg = load_num_reads(mkdb_dir);
    part_m4(mkdb_dir, m4_path, n_ctg, 100, 1);
}