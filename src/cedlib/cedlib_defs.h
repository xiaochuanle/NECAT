#ifndef CEDLIB_DEFS_H
#define CEDLIB_DEFS_H

#include <stdint.h>
#include "../common/ontcns_aux.h"
#include "../klib/kalloc.h"

typedef uint64_t        Word;
typedef kvec_t(Word)    vec_word;
#define WORD_SIZE       64
#define WORD_1          ((Word)1)
#define HIGH_BIT_MASK   ((WORD_1) << (WORD_SIZE -1))
#define MAX_UCHAR       255
#define AlphabetSize    4

#define CALC_NUM_BLOCKS(s) (((s) + WORD_SIZE - 1) / WORD_SIZE)

#define EDLIB_STATUS_OK     0
#define EDLIB_STATUS_ERROR  1

#define EDLIB_EDOP_MATCH    0
#define EDLIB_EDOP_INSERT   1
#define EDLIB_EDOP_DELETE   2
#define EDLIB_EDOP_MISMATCH 3

typedef kvec_t(char) vec_char;

typedef enum {
    EDLIB_MODE_NW,
    EDLIB_MODE_SHW,
    EDLIB_MODE_HW
} EdlibAlignMode;

typedef enum {
    EDLIB_TASK_DISTANCE,
    EDLIB_TASK_LOC,
    EDLIB_TASK_PATH
} EdlibAlignTask;

typedef enum {
    EDLIB_CIGAR_STANDARD,
    EDLIB_CIGAR_EXTENDED
} EdlibCigarFormat;

typedef struct {
    int k;
    EdlibAlignMode mode;
    EdlibAlignTask   task;
} EdlibAlignConfig;

EdlibAlignConfig
edlibNewAlignConfig(int k, EdlibAlignMode mode, EdlibAlignTask task);

typedef struct {
    Word P;
    Word M;
    int score;
} EdlibBlock;

typedef kvec_t(EdlibBlock) vec_EdlibBlock;

typedef struct {
	/// matrix 
    vec_word               Ps;
    vec_word               Ms;
    vec_int                scores;
    vec_int                firstBlocks;
    vec_int                lastBlocks;
	// scan data
    vec_EdlibBlock         blocks;
	// peq
    vec_word                peq;
    // result
	int						status;
	int						editDistance;
	vec_int					startLocations;
	vec_int					endLocations;
	vec_uchar				alignment;
} EdlibData;

EdlibData*
new_EdlibData(const int maxQueryLength, const int maxTargetLength);

EdlibData*
free_EdlibData(EdlibData* data);

void
reset_EdlibData(EdlibData* data, const int queryLength, const int targetLength);

#endif // CEDLIB_DEFS_H
