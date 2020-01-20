#include "blockwise_edlib.h"

#include <assert.h>

#define EDLIB_STATUS_OK     0
#define EDLIB_STATUS_ERROR  0

// Edit operations.
#define EDLIB_EDOP_MATCH 0    //!< Match.
#define EDLIB_EDOP_INSERT 1   //!< Insertion to target = deletion from query.
#define EDLIB_EDOP_DELETE 2   //!< Deletion from target = insertion to query.
#define EDLIB_EDOP_MISMATCH 3 //!< Mismatch.

#ifndef GAP_CHAR
#define GAP_CHAR '-'
#endif

EdlibAlignMatrix*
new_EdlibAlignMatrix(int maxNumBlocks, int targetLength)
{
    EdlibAlignMatrix* data = (EdlibAlignMatrix*)malloc(sizeof(EdlibAlignMatrix));
    data->Ps = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->Ms = (Word*)malloc(sizeof(Word) * maxNumBlocks * targetLength);
    data->scores = (int*)malloc(sizeof(int) * maxNumBlocks * targetLength);
    data->first_blocks = (int*)malloc(sizeof(int) * targetLength);
    data->last_blocks = (int*)malloc(sizeof(int) * targetLength);

    return data;
}

EdlibAlignMatrix*
free_EdlibAlignMatrix(EdlibAlignMatrix* data)
{
    free(data->Ps);
    free(data->Ms);
    free(data->scores);
    free(data->first_blocks);
    free(data->last_blocks);
    free(data);
    return 0;
}

SmallEdlibAlignResult*
new_SmallEdlibAlignResult()
{
    SmallEdlibAlignResult* result = (SmallEdlibAlignResult*)malloc(sizeof(SmallEdlibAlignResult));
    kv_init(result->end_locations);
    kv_init(result->start_locations);
    kv_init(result->alignment);
    return result;
}

SmallEdlibAlignResult*
free_SmallEdlibAlignResult(SmallEdlibAlignResult* result)
{
    kv_destroy(result->end_locations);
    kv_destroy(result->start_locations);
    kv_destroy(result->alignment);
    free(result);
    return 0;
}

void
clear_SmallEdlibAlignResult(SmallEdlibAlignResult* result)
{
    result->edit_distance = -1;
    kv_clear(result->end_locations);
    kv_clear(result->start_locations);
    kv_clear(result->alignment);
}

EdlibAlignData*
new_EdlibAlignData(int maxNumBlocks, int targetLength)
{
    EdlibAlignData* data = (EdlibAlignData*)malloc(sizeof(EdlibAlignData));
    data->align_matrix = new_EdlibAlignMatrix(maxNumBlocks, targetLength);
    data->result = new_SmallEdlibAlignResult();
    data->blocks = (EdlibBlock*)malloc(sizeof(EdlibBlock) * maxNumBlocks);
	data->peq = (Word*)malloc(sizeof(Word) * (AlphabetSize + 1) * maxNumBlocks);
    kv_init(data->cigar);
	data->maxNumBlocks = maxNumBlocks;
	data->maxTargetLength = targetLength;
    return data;
}

EdlibAlignData*
free_EdlibAlignData(EdlibAlignData* data)
{
    data->align_matrix = free_EdlibAlignMatrix(data->align_matrix);
    data->result = free_SmallEdlibAlignResult(data->result);
    free(data->blocks);
	free(data->peq);
    kv_destroy(data->cigar);
    free(data);
    return 0;
}

void
alloc_space_EdlibAlignData(EdlibAlignData* data, int query_length, int target_length)
{
	int newMaxNumBlocks = (query_length + WORD_SIZE - 1) / WORD_SIZE;
	int newMaxTargetLength = target_length;
	int newMaxSpace = sizeof(Word) * newMaxNumBlocks * newMaxTargetLength;
	
	if (newMaxNumBlocks > data->maxNumBlocks || newMaxTargetLength > data->maxTargetLength) {
		data->align_matrix->Ps = (Word*)realloc(data->align_matrix->Ps, newMaxSpace);
		data->align_matrix->Ms = (Word*)realloc(data->align_matrix->Ms, newMaxSpace);
		data->align_matrix->scores = (int*)realloc(data->align_matrix->scores, sizeof(int) * newMaxNumBlocks * newMaxTargetLength);
		data->align_matrix->first_blocks = (int*)realloc(data->align_matrix->first_blocks, sizeof(int) * newMaxTargetLength);
		data->align_matrix->last_blocks = (int*)realloc(data->align_matrix->last_blocks, sizeof(int) * newMaxTargetLength);
		
		data->blocks = (EdlibBlock*)realloc(data->blocks, sizeof(EdlibBlock) * newMaxNumBlocks);
		data->peq = (Word*)realloc(data->peq, sizeof(Word) * (AlphabetSize + 1) * newMaxNumBlocks);
		
		data->maxNumBlocks = newMaxNumBlocks;
		data->maxTargetLength = newMaxTargetLength;
	}
}

static inline int
calc_num_blocks(const int n)
{
    return (n + WORD_SIZE - 1) / WORD_SIZE;
}

static inline void
calc_block_cell_scores(const EdlibBlock block, int scores[])
{
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; ++i) {
        scores[i] = score;
        if (block.P & mask) --score;
        if (block.M & mask) ++score;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
}

static void
build_peq(const char* query, const int query_size, Word* peq)
{
    const int nblk = calc_num_blocks(query_size);
    for (int s = 0; s <= AlphabetSize; ++s) {
        for (int b = 0; b < nblk; ++b) {
            const int bid = s * MaxNumBlocks + b;
            if (s < AlphabetSize) {
                peq[bid] = 0;
                for (int r = (b + 1) * WORD_SIZE - 1; r >= b * WORD_SIZE; --r) {
                    peq[bid] <<= 1;
                    if (r >= query_size || query[r] == s) peq[bid] += 1;
                }
            } else {
                peq[bid] = (Word)-1;
            }
        }
    }
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word* PvOut, Word* MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = (Word)(hin >> 2) & WORD_ONE; // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= (Word)((hin + 1) >> 1);

    (*PvOut) = Mh | ~(Xv | Ph);
    (*MvOut) = Ph & Xv;

    return hout;
}

static void
calc_edit_distance_semi_global(const char* query, const int query_size,
							   const char* target, const int target_size,
							   const Word* peq,
							   int k,
							   const SmallEdlibAlignMode mode,
							   EdlibBlock* blocks,
							   int* edit_distance,
							   vec_int* end_positions)
{
    *edit_distance = -1;
    kv_clear(*end_positions);
    const int nblk = calc_num_blocks(query_size);
    const int W = nblk * WORD_SIZE - query_size;
    int fblk = 0;
    int lblk; // = min(calc_num_blocks(k + 1), nblk) - 1;
	{
		int X = calc_num_blocks(k + 1);
		lblk = hbn_min(X, nblk);
		--lblk;
	}
    EdlibBlock* blk = blocks;

    if (mode == EEDLIB_MODE_HW) {
        //k = min(query_size, k);
		k = hbn_min(query_size, k);
    }

    for (int b = 0; b <= lblk; ++b) {
        blk->score = (b + 1) * WORD_SIZE;
        blk->P = (Word)-1;
        blk->M = (Word)0;
        ++blk;
    }

    int best_d = -1;
    const int start_hout = mode == EEDLIB_MODE_HW ? 0 : 1;
    for (int c = 0; c < target_size; ++c) {
        const int tc = target[c];
        const Word* cpeq = peq + tc * MaxNumBlocks;
        int hout = start_hout;
        blk = blocks + fblk;
        cpeq += fblk;
        for (int b = fblk; b <= lblk; ++b) {
            hout = calculateBlock(blk->P, blk->M, *cpeq, hout, &blk->P, &blk->M);
            blk->score += hout;
            ++blk;
            ++cpeq;
        }
        --blk;
        --cpeq;

        if ( (lblk < nblk - 1) && (blk->score - hout <= k) && ((*(cpeq + 1) & WORD_ONE) || hout < 0) ) {
            ++lblk;
            ++blk;
            ++cpeq;
            blk->P = (Word)-1;
            blk->M = (Word)0;
            blk->score = (blk - 1)->score - hout + WORD_SIZE + calculateBlock(blk->P, blk->M, *cpeq, hout, &blk->P, &blk->M);
        } else {
            while (lblk >= fblk && blk->score >= k + WORD_SIZE) {
                --lblk;
                --blk;
                --cpeq;
            }
        }

        if (mode != EEDLIB_MODE_HW) {
            while (fblk <= lblk && blocks[fblk].score >= k + WORD_SIZE) {
                ++fblk;
            }
        }

        if (mode == EEDLIB_MODE_HW) {
            //lblk = max(0, lblk);
			lblk = hbn_max(0, lblk);
        }

        if (lblk < fblk) {
            *edit_distance = best_d;
            return;
        }

        if (lblk == nblk - 1) {
            int col_score = blk->score;
            if (col_score <= k) {
                if (best_d == -1 || col_score <= best_d) {
                    if (col_score != best_d) {
                        //end_positions.clear();
						kv_clear(*end_positions);
                        best_d = col_score;
                        k = best_d;
                    }
                    //end_positions.push_back(c - W);
					kv_push(int, *end_positions, c - W);
                }
            }
        }
    }

    if (lblk == nblk - 1) {
        int scores[WORD_SIZE];
        calc_block_cell_scores(*blk, scores);
        for (int i = 0; i < W; ++i) {
            int col_score = scores[i + 1];
            if (col_score <= k && (best_d == -1 || col_score <= best_d)) {
                if (col_score != best_d) {
                    //end_positions.clear();
					kv_clear(*end_positions);
                    k = best_d = col_score;
                }
                //end_positions.push_back(target_size - W + i);
				kv_push(int, *end_positions, target_size - W + i);
            }
        }
    }

    *edit_distance = best_d;
}

static void
calc_edit_distance_nw(const char* query, const int query_size,
					  const char* target, const int target_size,
					  const Word* peq,
					  int k,
					  const BOOL traceback,
					  EdlibAlignMatrix* align_matrix,
					  EdlibBlock* blocks,
					  const int target_stop_position,
					  int* edit_distance,
					  int* end_position)
{
    *edit_distance = -1;
    *end_position = -1;

    if (target_stop_position > -1 && traceback) {
		HBN_ERR("invalid parameters: target_stop_position = %d, traceback = %d", target_stop_position, traceback);
    }

    if (k < abs(target_size - query_size)) {
        return;
    }

    //k = min(k, max(query_size, target_size));
	{
		int X = hbn_max(query_size, target_size);
		k = hbn_min(k, X);
	}

    const int nblk = calc_num_blocks(query_size);
    const int W = nblk * WORD_SIZE - query_size;
    int fblk = 0;
    int lblk; // = min(nblk, calc_num_blocks(min(k, (k + query_size - target_size) / 2) + 1)) - 1;
	{
		int X = (k + query_size - target_size) / 2;
		int Y = hbn_min(k, X);
		int Z = calc_num_blocks(Y + 1);
		lblk = hbn_min(nblk, Z) - 1;
	}
    EdlibBlock* blk = blocks;

    for (int b = 0; b <= lblk; ++b) {
        blk->score = (b + 1) * WORD_SIZE;
        blk->P = (Word)-1;
        blk->M = (Word)0;
        ++blk;
    }

    for (int c = 0; c < target_size; ++c) {
        const int tc = target[c];
        const Word* cpeq = peq + tc * MaxNumBlocks;
        int hout = 1;
        blk = blocks + fblk;
        for (int b = fblk; b <= lblk; ++b) {
            hout = calculateBlock(blk->P, blk->M, cpeq[b], hout, &blk->P, &blk->M);
            blk->score += hout;
            ++blk;
        }
        --blk;

        //k = min(k, blk->score
        //        + max(target_size -c - 1, query_size - ((1 + lblk) * WORD_SIZE - 1) - 1)
        //        + (lblk == nblk - 1 ? W : 0));
		
		{
			int X1 = target_size - c - 1, X2 = query_size - ((1 + lblk) * WORD_SIZE - 1) - 1;
			int X = hbn_max(X1, X2);
			int Y = (lblk == nblk - 1) ? W : 0;
			int Z = X + Y + blk->score;
			k = hbn_min(k, Z);
		}

        if (lblk + 1 < nblk) {
            BOOL r = (lblk + 1) * WORD_SIZE - 1 > k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size;
            if (!r) {
                ++lblk;
                ++blk;
                blk->P = (Word)-1;
                blk->M = (Word)0;
                int newHout = calculateBlock(blk->P, blk->M, cpeq[lblk], hout, &blk->P, &blk->M);
                blk->score = (blk - 1)->score - hout + WORD_SIZE + newHout;
                hout = newHout;
            }
        }

        while (lblk >= fblk
				&&
				(blk->score >= k + WORD_SIZE
				 || ((lblk + 1) * WORD_SIZE -1 > 
					 k - blk->score + 2 * WORD_SIZE - 2 - target_size + c + query_size + 1))) {
            --lblk;
            --blk;
        }

        while (fblk <= lblk
				&&
				(blocks[fblk].score >= k + WORD_SIZE
				 || ((fblk + 1) * WORD_SIZE - 1 <
					 blocks[fblk].score - k - target_size + query_size + c))) {
            ++fblk;
        }

        if (lblk < fblk) {
            *edit_distance = *end_position = -1;
            return;
        }

        if (traceback) {
            blk = blocks + fblk;
            for (int b = fblk; b <= lblk; ++b) {
                align_matrix->Ps[MaxNumBlocks * c + b] = blk->P;
                align_matrix->Ms[MaxNumBlocks * c + b] = blk->M;
                align_matrix->scores[MaxNumBlocks * c + b] = blk->score;
                align_matrix->first_blocks[c] = fblk;
                align_matrix->last_blocks[c] = lblk;
                ++blk;
            }
        }

        if (c == target_stop_position) {
            for (int b = fblk; b <= lblk; ++b) {
                align_matrix->Ps[b] = blocks[b].P;
                align_matrix->Ms[b] = blocks[b].M;
                align_matrix->scores[b] = blocks[b].score;
                align_matrix->first_blocks[0] = fblk;
                align_matrix->last_blocks[0] = lblk;
            }
            *edit_distance = -1;
            *end_position = target_stop_position;
            return;
        }
    }

    if (lblk == nblk - 1) {
        int scores[WORD_SIZE];
        calc_block_cell_scores(blocks[lblk], scores);
        int col_score = scores[W];
        if (col_score <= k) {
            *edit_distance = col_score;
            *end_position = target_size - 1;
            return;
        }
    }

    *edit_distance = *end_position = -1;
}

/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in alignData. Consumes large amount of memory: O(queryLength * targetLength).
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] bestScore  Best score.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentTraceback(const int queryLength, 
									const int targetLength,
                                    const int bestScore, //const AlignmentData* const alignData,
									const EdlibAlignMatrix* align_matrix,
									vec_u8* alignment) {
                                    //vector<unsigned char>& alignment) {
	const int nblk = calc_num_blocks(queryLength);
	const int W = nblk * WORD_SIZE - queryLength;
	// alignment.clear();
	kv_clear(*alignment);
    int c = targetLength - 1; // index of column
    int b = nblk - 1; // index of block in column
    int currScore = bestScore; // Score of current cell
    int lScore  = -1; // Score of left cell
    int uScore  = -1; // Score of upper cell
    int ulScore = -1; // Score of upper left cell
    Word currP = align_matrix->Ps[c * MaxNumBlocks + b]; // P of current block
    Word currM = align_matrix->Ms[c * MaxNumBlocks + b]; // M of current block
    // True if block to left exists and is in band
    BOOL thereIsLeftBlock = c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1];
    // We set initial values of lP and lM to 0 only to avoid compiler warnings, they should not affect the
    // calculation as both lP and lM should be initialized at some moment later (but compiler can not
    // detect it since this initialization is guaranteed by "business" logic).
    Word lP = 0, lM = 0;
    if (thereIsLeftBlock) {
        lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b]; // P of block to the left
        lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b]; // M of block to the left
    }
    currP <<= W;
    currM <<= W;
    int blockPos = WORD_SIZE - W - 1; // 0 based index of current cell in blockPos

    // TODO(martin): refactor this whole piece of code. There are too many if-else statements,
    // it is too easy for a bug to hide and to hard to effectively cover all the edge-cases.
    // We need better separation of logic and responsibilities.
    while (1) {
        if (c == 0) {
            thereIsLeftBlock = TRUE;
            lScore = b * WORD_SIZE + blockPos + 1;
            ulScore = lScore - 1;
        }

        // TODO: improvement: calculate only those cells that are needed,
        //       for example if I calculate upper cell and can move up,
        //       there is no need to calculate left and upper left cell
        //---------- Calculate scores ---------//
        if (lScore == -1 && thereIsLeftBlock) {
            lScore = align_matrix->scores[(c - 1) * MaxNumBlocks + b]; // score of block to the left
            for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
                if (lP & HIGH_BIT_MASK) lScore--;
                if (lM & HIGH_BIT_MASK) lScore++;
                lP <<= 1;
                lM <<= 1;
            }
        }
        if (ulScore == -1) {
            if (lScore != -1) {
                ulScore = lScore;
                if (lP & HIGH_BIT_MASK) ulScore--;
                if (lM & HIGH_BIT_MASK) ulScore++;
            }
            else if (c > 0 && b-1 >= align_matrix->first_blocks[c-1] && b-1 <= align_matrix->last_blocks[c-1]) {
                // This is the case when upper left cell is last cell in block,
                // and block to left is not in band so lScore is -1.
                ulScore = align_matrix->scores[(c - 1) * MaxNumBlocks + b - 1];
            }
        }
        if (uScore == -1) {
            uScore = currScore;
            if (currP & HIGH_BIT_MASK) uScore--;
            if (currM & HIGH_BIT_MASK) uScore++;
            currP <<= 1;
            currM <<= 1;
        }
        //-------------------------------------//

        // TODO: should I check if there is upper block?

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query
        if (uScore != -1 && uScore + 1 == currScore) {
            currScore = uScore;
            lScore = ulScore;
            uScore = ulScore = -1;
            if (blockPos == 0) { // If entering new (upper) block
                if (b == 0) { // If there are no cells above (only boundary cells)
					//alignment.push_back(EDLIB_EDOP_INSERT);
					//for (int i = 0; i < c + 1; ++i) alignment.push_back(EDLIB_EDOP_DELETE);
					kv_push(unsigned char, *alignment, EDLIB_EDOP_INSERT);
					for (int i = 0; i < c + 1; ++i) kv_push(unsigned char, *alignment, EDLIB_EDOP_DELETE);
					break;
                } else {
                    blockPos = WORD_SIZE - 1;
                    b--;
                    currP = align_matrix->Ps[c * MaxNumBlocks + b];
                    currM = align_matrix->Ms[c * MaxNumBlocks + b];
                    if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                        thereIsLeftBlock = TRUE;
                        lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b]; // TODO: improve this, too many operations
                        lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
                    } else {
                        thereIsLeftBlock = FALSE;
                        // TODO(martin): There may not be left block, but there can be left boundary - do we
                        // handle this correctly then? Are l and ul score set correctly? I should check that / refactor this.
                    }
                }
            } else {
                blockPos--;
                lP <<= 1;
                lM <<= 1;
            }
            // Mark move
            ///-(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
			//alignment.push_back(EDLIB_EDOP_INSERT);
			kv_push(unsigned char, *alignment, EDLIB_EDOP_INSERT);
        }
        // Move left - deletion from target - insertion to query
        else if (lScore != -1 && lScore + 1 == currScore) {
            currScore = lScore;
            uScore = ulScore;
            lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE; // Move left
                ///int numUp = b * WORD_SIZE + blockPos + 1;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				//alignment.push_back(EDLIB_EDOP_DELETE);
				kv_push(unsigned char, *alignment, EDLIB_EDOP_DELETE);
				int numUp = b * WORD_SIZE + blockPos + 1;
				//for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				for (int i = 0; i < numUp; ++i) kv_push(unsigned char, *alignment, EDLIB_EDOP_INSERT);
				break;
            }
            currP = lP;
            currM = lM;
            if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                thereIsLeftBlock = TRUE;
                lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b];
                lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = TRUE;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = FALSE;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
			//alignment.push_back(EDLIB_EDOP_DELETE);
			kv_push(unsigned char, *alignment, EDLIB_EDOP_DELETE);
        }
        // Move up left - (mis)match
        else if (ulScore != -1) {
            unsigned char moveCode = ulScore == currScore ? EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH;
            currScore = ulScore;
            uScore = lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                ///(*alignment)[(*alignmentLength)++] = moveCode; // Move left
                ///int numUp = b * WORD_SIZE + blockPos;
                ///for (int i = 0; i < numUp; i++) // Move up until end
                ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                ///break;
				
				//alignment.push_back(moveCode);
				kv_push(unsigned char, *alignment, moveCode);
				int numUp = b * WORD_SIZE + blockPos;
				//for (int i = 0; i < numUp; i++) alignment.push_back(EDLIB_EDOP_INSERT);
				for (int i = 0; i < numUp; ++i) kv_push(unsigned char, *alignment, EDLIB_EDOP_INSERT);
				break;
            }
            if (blockPos == 0) { // If entering upper left block
                if (b == 0) { // If there are no more cells above (only boundary cells)
                    ///(*alignment)[(*alignmentLength)++] = moveCode; // Move up left
                    ///for (int i = 0; i < c + 1; i++) // Move left until end
                    ///    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
                    ///break;
					
					//alignment.push_back(moveCode);
					//for (int i = 0; i < c + 1; i++) alignment.push_back(EDLIB_EDOP_DELETE);
					kv_push(unsigned char, *alignment, moveCode);
					for (int i = 0; i < c + 1; ++i) kv_push(unsigned char, *alignment, EDLIB_EDOP_DELETE);
					break;
                }
                blockPos = WORD_SIZE - 1;
                b--;
                currP = align_matrix->Ps[c * MaxNumBlocks + b];
                currM = align_matrix->Ms[c * MaxNumBlocks + b];
            } else { // If entering left block
                blockPos--;
                currP = lP;
                currM = lM;
                currP <<= 1;
                currM <<= 1;
            }
            // Set new left block
            if (c > 0 && b >= align_matrix->first_blocks[c-1] && b <= align_matrix->last_blocks[c-1]) {
                thereIsLeftBlock = TRUE;
                lP = align_matrix->Ps[(c - 1) * MaxNumBlocks + b];
                lM = align_matrix->Ms[(c - 1) * MaxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = TRUE;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = FALSE;
                }
            }
            // Mark move
            ///(*alignment)[(*alignmentLength)++] = moveCode;
			//alignment.push_back(moveCode);
			kv_push(unsigned char, *alignment, moveCode);
        } else {
            // Reached end - finished!
            break;
        }
        //----------------------------------//
    }

    //reverse(alignment.begin(), alignment.end());
	{
		int align_size = kv_size(*alignment);
		int mid = align_size / 2;
		unsigned char* align_str = kv_data(*alignment);
		for (int i = 0; i < mid; ++i) {
			int j = align_size - 1 - i;
			char tmpc = align_str[i];
			align_str[i] = align_str[j];
			align_str[j] = tmpc;
		}
	}
    return EDLIB_STATUS_OK;
}

static void
edlibAlignmentToCigar(const unsigned char* alignment, 
					  int alignmentLength, 
					  SmallEdlibCigarFormat cigarFormat,
					  vec_EdlibGapOp* cigar)
{
	kv_clear(*cigar);
	if (cigarFormat != EEDLIB_CIGAR_EXTENDED && cigarFormat != EEDLIB_CIGAR_STANDARD) {
        return;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EEDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
	EdlibGapAlignOp op;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
			op.num = numOfSameMoves;
			op.op = lastMove;
			//cigar.push_back(op);
			kv_push(EdlibGapAlignOp, *cigar, op);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    //cigar.clear();
					kv_clear(*cigar);
                    return;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
}

static void
cigar2AlignedString(vec_EdlibGapOp* cigar,
					const char* query,
					const char* target,
					char* query_align,
					char* target_align,
					int* _qend,
					int* _tend)
{
	int aidx = 0;
	int qend = 0;
	int tend = 0;
	int n_cigar = kv_size(*cigar);
	EdlibGapAlignOp* oplist = kv_data(*cigar);
	
	for (int i = 0; i < n_cigar; ++i) {
		int n = oplist[i].num;
		char op = oplist[i].op;
		switch (op) {
			case 'M':
				for (int c = 0; c < n; ++c) {
					int qc = *query++;
					int tc = *target++;
					assert(qc >= 0 && qc < 4);
					assert(qc >= 0 && qc < 4);
					query_align[aidx] = "ACGT"[qc]; //DecodeDNA(qc);
					target_align[aidx] = "ACGT"[tc]; //DecodeDNA(tc);
					++aidx;
					++qend;
					++tend;
				}
				break;
			case 'I':
				for (int c = 0; c < n; ++c) {
					int qc = *query++;
					assert(qc >= 0 && qc < 4);
					query_align[aidx] = "ACGT"[qc]; //DecodeDNA(qc);
					target_align[aidx] = GAP_CHAR;
					++aidx;
					++qend;
				}
				break;
			case 'D':
				for (int c = 0; c < n; ++c) {
					int tc = *target++;
					assert(tc >= 0 && tc < 4);
					query_align[aidx] = GAP_CHAR;
					target_align[aidx] = "ACGT"[tc]; //DecodeDNA(tc);
					++aidx;
					++tend;
				}
				break;
			default:
				HBN_ERR("invalid operation: %c", op);
				break;
		}
	}
	query_align[aidx] = '\0';
	target_align[aidx] = '\0';
	*_qend = qend;
	*_tend = tend;
}

int
Edlib_align(const char* query,
			const int query_size,
			const char* target,
			const int target_size,
			EdlibAlignData* align_data,
			const double error,
			char* query_align,
			char* target_align,
			int* qend,
			int* tend)
{
	query_align[0] = '\0';
	target_align[0] = '\0';
	*qend = 0;
	*tend = 0;
	clear_SmallEdlibAlignResult(align_data->result);
	
	build_peq(query, query_size, align_data->peq);
	int k = hbn_min(query_size, target_size) * error * 1.1;
	calc_edit_distance_semi_global(query,
								   query_size,
								   target,
								   target_size,
								   align_data->peq,
								   k,
								   EEDLIB_MODE_SHW,
								   align_data->blocks,
								   &align_data->result->edit_distance,
								   &align_data->result->end_locations);
	if (align_data->result->edit_distance == -1) return 0;
	assert(kv_size(align_data->result->end_locations) > 0);
	
	int edit_distance;
	int end_position;
	calc_edit_distance_nw(query,
						  query_size,
						  target,
						  kv_A(align_data->result->end_locations, 0) + 1,
						  align_data->peq,
						  align_data->result->edit_distance,
						  TRUE,
						  align_data->align_matrix,
						  align_data->blocks,
						  -1,
						  &edit_distance,
						  &end_position);
	assert(edit_distance == align_data->result->edit_distance);
	assert(end_position == kv_A(align_data->result->end_locations, 0));
	
	obtainAlignmentTraceback(query_size,
							 kv_A(align_data->result->end_locations, 0) + 1,
							 edit_distance,
							 align_data->align_matrix,
							 &align_data->result->alignment);
	edlibAlignmentToCigar(kv_data(align_data->result->alignment),
						  kv_size(align_data->result->alignment),
						  EEDLIB_CIGAR_STANDARD,
						  &align_data->cigar);
	cigar2AlignedString(&align_data->cigar,
						query,
						target,
						query_align,
						target_align,
						qend,
						tend);
	
	return 1;
}

static const int kOcaMatCnt = 8;
//static const int kOcaBlockSize = 512;

#define kMatchCnt1 8
#define kMatchCnt2	8

BlockwiseEdlibData*
new_BlockwiseEdlibData(double error, int block_size)
{
	BlockwiseEdlibData* data = (BlockwiseEdlibData*)malloc(sizeof(BlockwiseEdlibData));
	data->edlib = new_EdlibAlignData(MaxNumBlocks, MaxSeqSize);
	ks_init(data->query_align);
	ks_init(data->target_align);
	ks_init(data->fqaln);
	ks_init(data->rqaln);
	ks_init(data->ftaln);
	ks_init(data->rtaln);
	ks_init(data->qfrag);
	ks_init(data->tfrag);
	data->qabuf = (char*)malloc(100000);
	data->tabuf = (char*)malloc(100000);
	data->error = error;
	data->block_size = block_size;
	
	return data;
}

BlockwiseEdlibData*
free_BlockwiseEdlibData(BlockwiseEdlibData* data)
{
	if (data) {
		if (data->edlib) data->edlib = free_EdlibAlignData(data->edlib);
		ks_destroy(data->query_align);
		ks_destroy(data->target_align);
		ks_destroy(data->fqaln);
		ks_destroy(data->rqaln);
		ks_destroy(data->ftaln);
		ks_destroy(data->rtaln);
		ks_destroy(data->qfrag);
		ks_destroy(data->tfrag);
		free(data->qabuf);
		free(data->tabuf);
		free(data);
	}
	return 0;
}

static void
validate_aligned_string(int qid,
						const char* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const char* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	//return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
			const char qc1 = DecodeDNA(right_extend ? query[x] : query[-x]);
			hbn_assert(qc == qc1, "qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu",
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
					  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
			const char tc1 = DecodeDNA(right_extend ? target[y] : target[-y]);
			hbn_assert(tc == tc1, "qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d",
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);

			++y;
		}
	}
}

static BOOL
get_next_sequence_block(const char* query,
						int qidx,
						const int qsize,
						const char* target,
						int tidx,
						const int tsize,
						const int desired_block_size,
						const BOOL right_extend,
						kstring_t* qfrag,
						kstring_t* tfrag)
{
	BOOL last_block = FALSE;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	int qblk;
	int tblk;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = tleft * 1.3;
		qblk = hbn_min(qblk, qleft);
		tblk = qleft * 1.3;
		tblk = hbn_min(tblk, tleft);
		last_block = TRUE;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = FALSE;
	}
	
	ks_clear(*qfrag);
	ks_clear(*tfrag);
	if (right_extend) {
		const char* Q = query + qidx;
		for (int i = 0; i < qblk; ++i) kputc(Q[i], qfrag);
		const char* R = target + tidx;
		for (int i = 0; i < tblk; ++i) kputc(R[i], tfrag);
	} else {
		const char* Q = query - qidx;
		for (int i = 0; i < qblk; ++i) kputc(Q[-i], qfrag);
		const char* R = target - tidx;
		for (int i = 0; i < tblk; ++i) kputc(R[-i], tfrag);
	}
	
	return last_block;
}

static void
oca_extend(const char* query,
		   const int query_size,
		   const char* target,
		   const int target_size,
		   EdlibAlignData* align_data,
		   const int block_size,
		   const double error,
		   const BOOL right_extend,
		   kstring_t* qfrag,
		   kstring_t* tfrag,
		   char* qabuf,
		   char* tabuf,
		   kstring_t* qaln,
		   kstring_t* taln)
{
	ks_clear(*qaln);
	ks_clear(*taln);
	int qidx = 0, tidx = 0;
	
	//OC_LOG("query_size = %d, target_size = %d", query_size, target_size);
	
	while (1) {
		int qfae, tfae, qfrag_size, tfrag_size;
		BOOL last_block = get_next_sequence_block(query,
												  qidx,
												  query_size,
												  target,
												  tidx,
												  target_size,
												  block_size,
												  right_extend,
												  qfrag,
												  tfrag);
		qfrag_size = ks_size(*qfrag);
		tfrag_size = ks_size(*tfrag);
		if (qfrag_size == 0 || tfrag_size == 0) break;
		
		Edlib_align(ks_s(*qfrag),
					qfrag_size,
					ks_s(*tfrag),
					tfrag_size,
					align_data,
					error,
					qabuf,
					tabuf,
					&qfae,
					&tfae);
		
		//OC_LOG("qsize = %d, tsize = %d, qfae = %d, tfae = %d", qfrag_size, tfrag_size, qfae, tfae);
		
		//OC_LOG("validating align string");
		validate_aligned_string(0,
								ks_s(*qfrag),
								0,
								qfae,
								qabuf,
								0,
								ks_s(*tfrag),
								0,
								tfae,
								tabuf,
								strlen(qabuf),
								TRUE);
		
		BOOL done = last_block;
		int acnt = 0, qcnt = 0, tcnt = 0;
		if (qfrag_size - qfae > 30 && tfrag_size - tfae > 30) done = 1;
		const int M = done ? kMatchCnt2 : kMatchCnt1;
		int align_size = strlen(qabuf);
		int k = align_size - 1, m = 0;
		while (k >= 0) {
			const char qc = qabuf[k];
			const char tc = tabuf[k];
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == M) break;
			--k;
		}
		
		if (m != M || k < 1) {
			align_size = 0;
			for (int i = 0; i < qfrag_size && i < tfrag_size; ++i) {
				const char qc = ks_A(*qfrag, i);
				const char tc = ks_A(*tfrag, i);
				if (qc != tc) break;
				qabuf[align_size] = DecodeDNA(qc);
				tabuf[align_size] = DecodeDNA(tc);
				++align_size;
			}
			done = TRUE;
		} else {
			align_size -= acnt;
			qidx += (qfae - qcnt);
			tidx += (tfae - tcnt);
			if (done) align_size += M;
		}
	
		kputsn(qabuf, align_size, qaln);
		kputsn(tabuf, align_size, taln);
		if (done) break;
	}
	
	if (1) {
		int qend = 0, tend = 0;
		for (size_t i = 0; i != ks_size(*qaln); ++i) {
			if (ks_A(*qaln, i) != GAP_CHAR) ++qend;
			if (ks_A(*taln, i) != GAP_CHAR) ++tend;
		}
		//OC_LOG("validating align string");
		validate_aligned_string(0,
							query,
							0,
							qend,
							ks_s(*qaln),
							0,
							target,
							0,
							tend,
							ks_s(*taln),
							ks_size(*qaln),
							right_extend);
	}
}

static double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
	return 100.0 * n / align_size;
}

int
blockwise_edlib_align(BlockwiseEdlibData* oca_data,
	const u8* u8_query,
	const u8* u8_target,
	GappedCandidate* can,
	const int min_align_size,
	const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align_,
	kstring_t* target_align_)
{
	const char* query = (const char*)(u8_query);
	const char* target = (const char*)(u8_target);
	kstring_t* query_align = &oca_data->query_align;
	kstring_t* target_align = &oca_data->target_align;
	ks_clear(*query_align);
	ks_clear(*target_align);
	const int query_start = can->qoff;
	const int target_start = can->soff;
	const int query_size = can->qsize;
	const int target_size = can->ssize;
	int QS = query_start;
	int TS = target_start;
	BOOL right_extend;
	
	right_extend = FALSE;
	oca_extend(query + QS - 1,
			   QS,
			   target + TS -1,
			   TS,
			   oca_data->edlib,
			   oca_data->block_size,
			   oca_data->error,
			   right_extend,
			   &oca_data->qfrag,
			   &oca_data->tfrag,
			   oca_data->qabuf,
			   oca_data->tabuf,
			   &oca_data->rqaln,
			   &oca_data->rtaln);
	int rqcnt = 0, rtcnt = 0;
	{
		int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
		int raln_size = ks_size(oca_data->rqaln);
		for (int i = 0; i < raln_size; ++i) {
			const char qc = ks_A(oca_data->rqaln, i);
			const char tc = ks_A(oca_data->rtaln, i);
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == kOcaMatCnt) break;
		}
		if (m == kOcaMatCnt) {
			QS -= qcnt;
			TS -= tcnt;
			for (int i = raln_size; i > acnt; --i) {
				const char qc = ks_A(oca_data->rqaln, i - 1);
				kputc(qc, query_align);
				if (qc != GAP_CHAR) ++rqcnt;
				const char tc = ks_A(oca_data->rtaln, i - 1);
				kputc(tc, target_align);
				if (tc != GAP_CHAR) ++rtcnt;
			}
		}
	}
	
	right_extend = TRUE;
	oca_extend(query + QS,
			   query_size - QS,
			   target + TS,
			   target_size - TS,
			   oca_data->edlib,
			   oca_data->block_size,
			   oca_data->error,
			   right_extend,
			   &oca_data->qfrag,
			   &oca_data->tfrag,
			   oca_data->qabuf,
			   oca_data->tabuf,
			   &oca_data->fqaln,
			   &oca_data->ftaln);
	int fqcnt = 0, ftcnt = 0;
	if (ks_size(*query_align) == 0) {
		int qcnt = 0, tcnt = 0, acnt = 0, m = 0;
		int faln_size = ks_size(oca_data->fqaln);
		for (int i = 0; i < faln_size; ++i) {
			const char qc = ks_A(oca_data->fqaln, i);
			const char tc = ks_A(oca_data->ftaln, i);
			if (qc != GAP_CHAR) ++qcnt;
			if (tc != GAP_CHAR) ++tcnt;
			if (qc == tc) {
				++m;
			} else {
				m = 0;
			}
			++acnt;
			if (m == kOcaMatCnt) break;
		}
		if (m == kOcaMatCnt) {
			acnt -= kOcaMatCnt;
			qcnt -= kOcaMatCnt;
			tcnt -= kOcaMatCnt;
			QS += qcnt;
			TS += tcnt;
			for (int i = acnt; i < faln_size; ++i) {
				const char qc = ks_A(oca_data->fqaln, i);
				if (qc != GAP_CHAR) ++fqcnt;
				kputc(qc, query_align);
				const char tc = ks_A(oca_data->ftaln, i);
				if (tc != GAP_CHAR) ++ftcnt;
				kputc(tc, target_align);
			}
		}
	} else {
		int faln_size = ks_size(oca_data->fqaln);
		for (int i = 0; i < faln_size; ++i) {
			const char qc = ks_A(oca_data->fqaln, i);
			if (qc != GAP_CHAR) ++fqcnt;
			kputc(qc, query_align);
			const char tc = ks_A(oca_data->ftaln, i);
			if (tc != GAP_CHAR) ++ftcnt;
			kputc(tc, target_align);
		}
	}
	
	oca_data->qoff = QS - rqcnt;
	oca_data->qend = QS + fqcnt;
	oca_data->toff = TS - rtcnt;
	oca_data->tend = TS + ftcnt;
	
	//OC_LOG("validating align string");
	validate_aligned_string(oca_data->qid,
							query,
							oca_data->qoff,
							oca_data->qend,
							ks_s(*query_align),
							oca_data->tid,
							target,
							oca_data->toff,
							oca_data->tend,
							ks_s(*target_align),
							ks_size(*target_align),
						    TRUE);
	int align_size = ks_size(*target_align);
	oca_data->ident_perc = calc_ident_perc(ks_s(*query_align), ks_s(*target_align), align_size);
	int r = (align_size >= min_align_size) && (oca_data->ident_perc >= min_ident_perc);
	ks_clear(*query_align_);
	ks_clear(*target_align_);
	if (r) {
		*qbeg = oca_data->qoff;
		*qend = oca_data->qend;
		*tbeg = oca_data->toff;
		*tend = oca_data->tend;
		*ident_perc = oca_data->ident_perc;
		kputsn(ks_s(*query_align), ks_size(*query_align), query_align_);
		kputsn(ks_s(*target_align), ks_size(*target_align), target_align_);
	}
	return r;
}
