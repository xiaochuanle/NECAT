#include "range_list.h"
#include <assert.h>

KSORT_INIT(CovRange_LT, CovRange, CovRange_LT)

void
init_CovRangeList(CovRangeList* list)
{
  list->is_sorted = FALSE;
  list->is_merged = FALSE;
  kv_init(list->list);
}

void
destroy_CovRangeList(CovRangeList* list)
{
  kv_destroy(list->list);
}

void
clear_CovRangeList(CovRangeList* list)
{
  list->is_sorted = FALSE;
  list->is_merged = FALSE;
  kv_clear(list->list);
}

void
copy_CovRangeList(CovRangeList* dst, CovRangeList* src)
{
  dst->is_sorted = src->is_sorted;
  dst->is_merged = src->is_merged;
  kv_resize(CovRange, dst->list, kv_size(src->list));
  size_t cp_size = sizeof(CovRange) * kv_size(src->list);
  memcpy(kv_data(dst->list), kv_data(src->list), cp_size);
}

void
add_CovRangeList(CovRangeList* list, int position, int length, int val)
{
  CovRange range;
  range.lo = position;
  range.hi = position + length;
  range.ct = 1;
  range.va = val;
  kv_push(CovRange, list->list, range);
  list->is_sorted = 0;
  list->is_merged = 0;
}

void
sort_CovRangeList(CovRangeList* list)
{
  if (CovRangeListIsSorted(*list)) return;
  ks_introsort_CovRange_LT(kv_size(list->list), kv_data(list->list));
  list->is_sorted = TRUE;
}

void
merge_CovRangeList(CovRangeList* list, int min_ovlp)
{
  if (CovRangeListIsMerged(*list)) return;
  int curr = 0, next = 1;
  sort_CovRangeList(list);
  int nrange = kv_size(list->list);
  CovRange* crl = kv_data(list->list);

  while (next < nrange) {
    if (CovRangeIsInvalid(crl[curr])) {
      crl[curr] = crl[next];
      InvalidCovRange(crl[next]);
      ++next;
    } else {
      BOOL intersect = FALSE;
      /// contained
      if (crl[curr].lo <= crl[next].lo && crl[next].hi <= crl[curr].hi) intersect = TRUE;
      /// overlap
      if (crl[curr].hi - min_ovlp >= crl[next].lo) intersect = TRUE;

      if (intersect) {
        if (crl[curr].hi < crl[next].hi) {
          crl[curr].hi = crl[next].hi;
	}
          crl[curr].ct += crl[next].ct;
          crl[curr].va += crl[next].va;

          InvalidCovRange(crl[next]);
          ++next;
      } else {
        ++curr;
        if (curr != next) crl[curr] = crl[next];
        ++next;
      }
    }
  }

  if (curr + 1 < nrange) kv_resize(CovRange, list->list, curr + 1);
  list->is_merged = 1;
}

typedef struct {
  int pos; // position of the change in depth
  int change; // the value associated with this object; added or subtracted from va
  BOOL open; // if true, the start of a new interval
} IntervalDepthRegion;

#define IntervalDepthRegion_LT(lhs, rhs) ( \
    ((lhs).open > (rhs).open) \
    || \
    ((lhs).open == (rhs).open && (lhs).pos < (rhs).pos) \
    )

KSORT_INIT(IntervalDepthRegion_LT, IntervalDepthRegion, IntervalDepthRegion_LT)

static void
compute_depth_for_CovRangeList(CovRangeList* list, IntervalDepthRegion* id, size_t id_len)
{
  clear_CovRangeList(list);
  if (id_len == 0) return;

  ks_introsort_IntervalDepthRegion_LT(id_len, id);
  kv_reserve(CovRange, list->list, id_len);
  CovRange* crl = kv_data(list->list);
  size_t list_len = 0;
  /// init the first interval
  assert(id[0].open == TRUE);
  crl[list_len].lo = id[0].pos;
  crl[list_len].hi = id[0].pos;
  crl[list_len].ct = 1;
  crl[list_len].va = id[0].change;

  int nct;
  int nva;

  for (size_t i = 1; i < id_len; ++i) {
    crl[list_len].hi = id[i].pos;

    if (id[i].open == TRUE) {
      nct = crl[list_len].ct + 1;
      nva = crl[list_len].va + id[i].change;
    } else {
      nct = crl[list_len].ct - 1;
      nva = crl[list_len].va - id[i].change;
    }

    if (
        ((id[i - 1].pos != id[i].pos) || (crl[list_len].va != nva))
        &&
        (crl[list_len].lo != crl[list_len].hi)
       ) {
      ++list_len;
      crl[list_len].lo = id[i].pos;
      crl[list_len].ct = crl[list_len - 1].ct;
      crl[list_len].va = crl[list_len - 1].va;
    }

    crl[list_len].hi = id[i].pos;
    crl[list_len].ct = nct;
    crl[list_len].va = nva;

    if (
        (list_len > 1)
        &&
        (crl[list_len - 1].hi == crl[list_len].lo)
        &&
        (crl[list_len - 1].ct == crl[list_len].ct)
        &&
        (crl[list_len - 1].va == crl[list_len].va)
       ) {
      crl[list_len - 1].hi = crl[list_len].hi;
      --list_len;
    }
  }

  kv_resize(CovRange, list->list, list_len);
}

void
depth_from_CovRangeList(CovRangeList* dst, CovRangeList* src)
{
  size_t id_len = 2 * kv_size(src->list);
  IntervalDepthRegion* id = (IntervalDepthRegion*)malloc( sizeof(IntervalDepthRegion) * id_len );
  for (size_t i = 0; i < kv_size(src->list); ++i) {
    id[2 * i    ].pos       = CovRangeListLo(*src, i);
    id[2 * i    ].change    = kv_A(src->list, i).va;
    id[2 * i    ].open      = TRUE;

    id[2 * i + 1].pos     = CovRangeListHi(*src, i);
    id[2 * i + 1].change  = kv_A(src->list, i).va;
    id[2 * i + 1].open    = FALSE;
  }

  compute_depth_for_CovRangeList(dst, id, id_len);

  free(id);
}
