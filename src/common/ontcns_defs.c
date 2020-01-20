#include "ontcns_defs.h"

#include "../klib/ksort.h"

KSORT_INIT(IntPair, IntPair, IntPair_LT)
#define generic_gt(a, b) ((a) > (b))
KSORT_INIT(double_gt, double, generic_gt)
