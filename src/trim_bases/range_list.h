#ifndef RANGE_LIST_H
#define RANGE_LIST_H

#include "../common/ontcns_defs.h"
#include "../klib/ksort.h"
#include "../klib/kvec.h"

typedef struct {
	int left, right, size;
} ClippedRange;

typedef kvec_t(ClippedRange) vec_ClippedRange;

typedef struct {
  int lo;
  int hi;
  int ct; // number of source ranges
  int va; // value
} CovRange;

#define CovRange_LT(lhs, rhs) (\
  ((lhs).lo < (rhs).lo) \
  || \
  ((lhs).lo == (rhs).lo && (lhs).hi < (rhs).hi) \
  ) 

typedef kvec_t(CovRange)  vec_CovRange;

#define CovRangeIsInvalid(r) ((r).lo == 0 && (r).hi == 0)

#define InvalidCovRange(r) \
  do { \
    (r).lo = 0; \
    (r).hi = 0; \
  } while(0)

typedef struct {
  BOOL          is_sorted;
  BOOL          is_merged;
  vec_CovRange  list;
} CovRangeList;

#define CovRangeListIsMerged(list) ((list).is_merged)
#define CovRangeListIsSorted(list) ((list).is_sorted)
#define CovRangeListLo(crl, i) (kv_A((crl).list, i).lo)
#define CovRangeListHi(crl, i) (kv_A((crl).list, i).hi)
#define CovRangeListCount(crl, i) (kv_A((crl).list, i).ct)
#define CovRangeListDepth(crl, i) (kv_A((crl).list, i).ct)


void
init_CovRangeList(CovRangeList* list);

void
destroy_CovRangeList(CovRangeList* list);

void
clear_CovRangeList(CovRangeList* list);

void
copy_CovRangeList(CovRangeList* dst, CovRangeList* src);

void
add_CovRangeList(CovRangeList* list, int position, int length, int val);

void
sort_CovRangeList(CovRangeList* list);

void
merge_CovRangeList(CovRangeList* list, int min_ovlp);

void
depth_from_CovRangeList(CovRangeList* dst, CovRangeList* src);

#endif // RANGE_LIST_H
