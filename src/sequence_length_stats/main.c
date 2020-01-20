#include "../klib/kseq.h"
#include "../klib/ksort.h"
#include "../klib/kvec.h"
#include "../common/ontcns_aux.h"

#include <stdint.h>
#include <stdio.h>

KSEQ_DECLARE(gzFile)

typedef int64_t idx_t;
typedef kvec_t(idx_t) vec_i64;
#define idxt_gt(a, b) ((a) > (b))
KSORT_INIT(idx_gt, idx_t, idxt_gt)

typedef struct {
    idx_t num_sequences;
    idx_t max_length;
    idx_t min_length;
    idx_t total_length;
    idx_t avg_length;
    idx_t median_length;
    idx_t n25_length;
    idx_t n25_cnt;
    idx_t n50_length;
    idx_t n50_cnt;
    idx_t n75_length;
    idx_t n75_cnt;
} LengthStats;

static void
print_digit_with_comma(FILE* out, idx_t num)
{
    const char* digits = "0123456789";
    char buf[256];
    char* p = buf;
    int c = 0;
    do {
        *p = digits[num % 10];
        ++p;
        ++c;
        num /= 10;
        if (c == 3 && num) {
            *p++ = ',';
            c = 0;
        }
    } while(num);

    while (p > buf) fprintf(out, "%c", *--p);
}

static void 
read_sequence_length(const char* path, vec_i64* length_list, const idx_t length_cutoff)
{
    DGZ_OPEN(in, path, "r");
    kseq_t* read = kseq_init(in);
    while (kseq_read(read) >= 0) {
        idx_t s = kseq_size(*read);
        if (s >= length_cutoff) {
            kv_push(idx_t, *length_list, s);
        }
    }
    GZ_CLOSE(in);
    kseq_destroy(read);

    ks_introsort(idx_gt, kv_size(*length_list), kv_data(*length_list));
}

void
calc_nxx_stats(vec_i64* length_list, const idx_t nxx_sum, idx_t* nxx_length, idx_t* nxx_cnt)
{
    idx_t sum = 0;
    idx_t cnt = 0;
    idx_t length = 0;
    for (size_t i = 0; i != kv_size(*length_list); ++i) {
        ++cnt;
        sum += kv_A(*length_list, i);
        if (sum >= nxx_sum) {
            length = kv_A(*length_list, i);
            break;
        }
    }

    *nxx_length = length;
    *nxx_cnt = cnt;
}

static void
calc_length_stats(vec_i64* length_list, LengthStats* stats)
{
    if (kv_size(*length_list) == 0) return;
    size_t num_sequences = kv_size(*length_list);
    stats->num_sequences = num_sequences;
    stats->max_length = kv_A(*length_list, 0);
    stats->min_length = kv_A(*length_list, num_sequences - 1);
    stats->median_length = kv_A(*length_list, num_sequences/2);

    idx_t sum = 0;
    for (size_t i = 0; i != kv_size(*length_list); ++i) sum += kv_A(*length_list, i);
    idx_t n25_sum = sum / 4;
    idx_t n50_sum = sum / 2;
    idx_t n75_sum = n25_sum * 3;
    stats->total_length = sum;
    stats->avg_length = sum / num_sequences;

    calc_nxx_stats(length_list, n25_sum, &stats->n25_length, &stats->n25_cnt);
    calc_nxx_stats(length_list, n50_sum, &stats->n50_length, &stats->n50_cnt);
    calc_nxx_stats(length_list, n75_sum, &stats->n75_length, &stats->n75_cnt);
}

static void
print_nxx_info(FILE* out, int xx, idx_t nxx_cnt, idx_t nxx_length)
{
    fprintf(out, "N%d stats:\t\t\t%d%% of total sequences are contained in the ", xx, xx);
    print_digit_with_comma(out, nxx_cnt);
    fprintf(out, " sequences having >= ");
    print_digit_with_comma(out, nxx_length);
    fprintf(out, " bps\n");
}

static void
print_length_stats(const LengthStats* stats)
{
    fprintf(stdout, "Number of sequences:\t\t");
    print_digit_with_comma(stdout, stats->num_sequences);
    fprintf(stdout, "\n");
    fprintf(stdout, "Number of bps:\t\t\t");
    print_digit_with_comma(stdout, stats->total_length);
    fprintf(stdout, "\n");
    fprintf(stdout, "Max length:\t\t\t");
    print_digit_with_comma(stdout, stats->max_length);
    fprintf(stdout, "\n");
    fprintf(stdout, "Min length:\t\t\t");
    print_digit_with_comma(stdout, stats->min_length);
    fprintf(stdout, "\n");
    fprintf(stdout, "Average length:\t\t\t");
    print_digit_with_comma(stdout, stats->avg_length);
    fprintf(stdout, "\n");
    fprintf(stdout, "Median length:\t\t\t");
    print_digit_with_comma(stdout, stats->median_length);
    fprintf(stdout, "\n");
    print_nxx_info(stdout, 25, stats->n25_cnt, stats->n25_length);
    print_nxx_info(stdout, 50, stats->n50_cnt, stats->n50_length);
    print_nxx_info(stdout, 75, stats->n75_cnt, stats->n75_length);
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s sequence-name sequence-length-cutoff\n", argv[0]);
        return 1;
    }

    new_kvec(vec_i64, length_list);
    LengthStats stats; memset(&stats, 0, sizeof(LengthStats));
    const char* sequence_path = argv[1];
    idx_t sequence_length_cutoff = atoll(argv[2]);
    read_sequence_length(sequence_path, &length_list, sequence_length_cutoff);
    calc_length_stats(&length_list, &stats);
    print_length_stats(&stats);
    free_kvec(length_list);
    return 0;
}
