#ifndef BLOCKWISE_EDLIB_H
#define BLOCKWISE_EDLIB_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "../common/ontcns_aux.h"
#include "../klib/kstring.h"
#include "../klib/kvec.h"
#include "../common/gapped_candidate.h"
#include "align_defs.h"

typedef uint64_t        Word;
#define WORD_SIZE       64
#define WORD_ONE        ((Word)1)
#define HIGH_BIT_MASK   (WORD_ONE << (WORD_SIZE - 1))

#define MaxSeqSize		4096
#define MaxNumBlocks	64
#define AlphabetSize    4

#define kOcaBlockSize512	512
#define kOcaBlockSize1024	1024
#define kOcaBlockSize2048	2048

#define EDLIB_OP_SUB 	('M')
#define EDLIB_OP_INS 	('I')
#define EDLIB_OP_DEL 	('D')

/**
 * Describes cigar format.
 * @see http://samtools.github.io/hts-specs/SAMv1.pdf
 * @see http://drive5.com/usearch/manual/cigar.html
 */
typedef enum{
	EEDLIB_CIGAR_STANDARD,  //!< Match: 'M', Insertion: 'I', Deletion: 'D', Mismatch: 'M'.
	EEDLIB_CIGAR_EXTENDED   //!< Match: '=', Insertion: 'I', Deletion: 'D', Mismatch: 'X'.
}SmallEdlibCigarFormat;

typedef enum 
{
    EEDLIB_MODE_NW,
    EEDLIB_MODE_SHW,
    EEDLIB_MODE_HW
}SmallEdlibAlignMode;

typedef struct {
	Word* Ps;
	Word* Ms;
	int* scores;
	int* first_blocks;
	int* last_blocks;
}EdlibAlignMatrix;

EdlibAlignMatrix*
new_EdlibAlignMatrix(int maxNumBlocks, int targetLength);

EdlibAlignMatrix*
free_EdlibAlignMatrix(EdlibAlignMatrix* data);

typedef struct {
	Word P;
	Word M;
	int score;
}EdlibBlock;

typedef struct {
	int 		edit_distance;
	vec_int		end_locations;
	vec_int 	start_locations;
	vec_u8	alignment;
}SmallEdlibAlignResult;

SmallEdlibAlignResult*
new_SmallEdlibAlignResult();

SmallEdlibAlignResult*
free_SmallEdlibAlignResult(SmallEdlibAlignResult* result);

void
clear_SmallEdlibAlignResult(SmallEdlibAlignResult* result);

typedef struct {
	int num;
	int op;
} EdlibGapAlignOp;

typedef kvec_t(EdlibGapAlignOp)	vec_EdlibGapOp;

typedef struct {
    EdlibAlignMatrix*   align_matrix;
    SmallEdlibAlignResult*   result;
    EdlibBlock*         blocks;
    Word*               peq;
    vec_EdlibGapOp      cigar;
	int					maxNumBlocks;
	int					maxTargetLength;
} EdlibAlignData;

EdlibAlignData*
new_EdlibAlignData(int maxNumBlocks, int targetLength);

EdlibAlignData*
free_EdlibAlignData(EdlibAlignData* data);

void
alloc_space_EdlibAlignData(EdlibAlignData* data, int query_length, int target_length);

typedef struct {
	EdlibAlignData*		edlib;
	int					qoff;
	int					qend;
	int					toff;
	int					tend;
	double				ident_perc;
	kstring_t			query_align;
	kstring_t			target_align;
	kstring_t			fqaln;
	kstring_t			rqaln;
	kstring_t			ftaln;
	kstring_t			rtaln;
	kstring_t			qfrag;
	kstring_t			tfrag;
	char*				qabuf;
	char*				tabuf;
	double				error;
	int					qid;
	int					tid;
	int 				block_size;
} BlockwiseEdlibData;

BlockwiseEdlibData*
new_BlockwiseEdlibData(double error, int block_size);

BlockwiseEdlibData*
free_BlockwiseEdlibData(BlockwiseEdlibData* data);

int
blockwise_edlib_align(BlockwiseEdlibData* oca_data,
	const u8* query,
	const u8* target,
	GappedCandidate* can,
	const int min_align_size,
	const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align,
	kstring_t* target_align);

#ifdef __cplusplus
}
#endif
#endif // BLOCKWISE_EDLIB_H