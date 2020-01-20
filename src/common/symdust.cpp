/*  $Id: symdust.cpp 464807 2015-04-14 16:29:40Z vakatov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Aleksandr Morgulis
 *
 * File Description:
 *   Implementation file for CSymDustMasker class.
 *
 */

//#include <ncbi_pch.hpp>

//#include <algo/dustmask/symdust.hpp>

#include "symdust.hpp"

BEGIN_NCBI_SCOPE

//------------------------------------------------------------------------------
CSymDustMasker::triplets::triplets( 
    size_type window, Uint1 low_k,
    perfect_list_type & perfect_list, thres_table_type & thresholds )
    : start_( 0 ), stop_( 0 ), max_size_( window - 2 ), low_k_( low_k ),
      L( 0 ), P( perfect_list ), thresholds_( thresholds ),
      r_w( 0 ), r_v( 0 ), num_diff( 0 )
{
    std::fill( c_w, c_w + 64, 0 );
    std::fill( c_v, c_v + 64, 0 );
}

//------------------------------------------------------------------------------
bool CSymDustMasker::triplets::shift_high( triplet_type t )
{
    triplet_type s = triplet_list_.back();
    triplet_list_.pop_back();
    rem_triplet_info( r_w, c_w, s );
    if( c_w[s] == 0 ) --num_diff;
    ++start_;

    triplet_list_.push_front( t );
    if( c_w[t] == 0 ) ++num_diff;
    add_triplet_info( r_w, c_w, t );
    ++stop_;

    if( num_diff <= 1 ) {
        P.insert( P.begin(), perfect( start_, stop_ + 1, 0, 0 ) );
        return false;
    }
    else return true;
}

//------------------------------------------------------------------------------
bool CSymDustMasker::triplets::shift_window( triplet_type t )
{
    // shift the left end of the window, if necessary
    if( triplet_list_.size() >= max_size_ ) {
        if( num_diff <= 1 ) {
            return shift_high( t );
        }

        triplet_type s = triplet_list_.back();
        triplet_list_.pop_back();
        rem_triplet_info( r_w, c_w, s );
        if( c_w[s] == 0 ) --num_diff;

        if( L == start_ ) {
            ++L;
            rem_triplet_info( r_v, c_v, s );
        }

        ++start_;
    }

    triplet_list_.push_front( t );
    if( c_w[t] == 0 ) ++num_diff;
    add_triplet_info( r_w, c_w, t );
    add_triplet_info( r_v, c_v, t );

    if( c_v[t] > low_k_ ) {
        Uint4 off = triplet_list_.size() - (L - start_) - 1;

        do {
            rem_triplet_info( r_v, c_v, triplet_list_[off] );
            ++L;
        }while( triplet_list_[off--] != t );
    }

    ++stop_;

    if( triplet_list_.size() >= max_size_ && num_diff <= 1 ) {
        P.clear();
        P.insert( P.begin(), perfect( start_, stop_ + 1, 0, 0 ) );
        return false;
    }else return true;
}

//------------------------------------------------------------------------------
inline void CSymDustMasker::triplets::find_perfect()
{
    typedef perfect_list_type::iterator perfect_iter_type;
    counts_type counts;

    Uint4 count = stop_ - L; // count is the suffix length

    // we need a local copy of triplet counts
    std::copy( c_v, c_v + 64, counts );

    Uint4 score = r_v; // and of the partial sum
    perfect_iter_type perfect_iter = P.begin();
    Uint4 max_perfect_score = 0;
    size_type max_len = 0;
    size_type pos = L - 1; // skipping the suffix
    impl_citer_type it = triplet_list_.begin() + count; // skipping the suffix
    impl_citer_type iend = triplet_list_.end();

    for( ; it != iend; ++it, ++count, --pos ) {
        Uint1 cnt = counts[*it];
        add_triplet_info( score, counts, *it );

        if( cnt > 0 && score*10 > thresholds_[count] ) {
            // found the candidate for the perfect interval
            // get the max score for the existing perfect intervals within
            //      current suffix
            while(    perfect_iter != P.end() 
                   && pos <= perfect_iter->bounds_.first ) {
                if(    max_perfect_score == 0 
                    || max_len*perfect_iter->score_ 
                       > max_perfect_score*perfect_iter->len_ ) {
                    max_perfect_score = perfect_iter->score_;
                    max_len = perfect_iter->len_;
                }

                ++perfect_iter;
            }

            // check if the current suffix score is at least as large
            if(    max_perfect_score == 0 
                || score*max_len >= max_perfect_score*count ) {
                max_perfect_score = score;
                max_len = count;
                perfect_iter = P.insert( 
                        perfect_iter, perfect( pos, stop_ + 1, 
                        max_perfect_score, count ) );
            }
        }
    }
}
    
//------------------------------------------------------------------------------
CSymDustMasker::CSymDustMasker( 
    Uint4 level, size_type window, size_type linker )
    : level_( (level >= 2 && level <= 64) ? level : DEFAULT_LEVEL ), 
      window_( (window >= 8 && window <= 64) ? window : DEFAULT_WINDOW ), 
      linker_( (linker >= 1 && linker <= 32) ? linker : DEFAULT_LINKER ),
      low_k_( level_/5 )
{
    thresholds_.reserve( window_ - 2 );
    thresholds_.push_back( 1 );

    for( size_type i = 1; i < window_ - 2; ++i )
        thresholds_.push_back( i*level_ );
}

//------------------------------------------------------------------------------
inline void CSymDustMasker::save_masked_regions( 
        TMaskList & res, size_type wstart, size_type start )
{
    if( !P.empty() ) {
        TMaskedInterval b = P.back().bounds_;
        
        if( b.first < wstart ) {
            TMaskedInterval b1( b.first + start, b.second + start );

            if( !res.empty() ) {
                size_type s = res.back().second;

                if( s + linker_ >= b1.first ) {
                    res.back().second = std::max( s, b1.second );
                }else {
                    res.push_back( b1 );
                }
            }else {
                res.push_back( b1 );
            }

            while( !P.empty() && P.back().bounds_.first < wstart ) {
                P.pop_back();
            }
        }
    }
}

//------------------------------------------------------------------------------
void
CSymDustMasker::operator()( const sequence_type & seq, 
							size_type seq_len,
                            size_type start, size_type stop,
						    TMaskList& res)
{	
	res.clear();
	
	if (seq_len == 0) return;
	
	if (stop >= seq_len) stop = seq_len - 1;

    if( start > stop )
        start = stop;
	
#define _GET_POS(_start, _stop) ((_stop) - (_start))

    while( stop > 2 + start )    // there must be at least one triplet
    {
        // initializations
        P.clear();
        triplets w( window_, low_k_, P, thresholds_ );

        ///seq_citer_type it(seq, start);
		seq_citer_type it = seq + start;
		seq_citer_type it_start = seq;

        char c1 = *it, c2 = *++it;
        triplet_type t = (converter_( c1 )<<2) + converter_( c2 );

        //it.SetPos(start + w.stop() + 2);
		it = it_start + (start + w.stop() + 2);

        bool done = false;
        while( !done && _GET_POS(it_start, it) <= stop )
        {
            save_masked_regions( res, w.start(), start );

            // shift the window
            t = ((t<<2)&TRIPLET_MASK) + (converter_( *it )&0x3);
            ++it;

            if( w.shift_window( t ) ) {
                if( w.needs_processing() ) {
                    w.find_perfect();
                }
            }else {
                while( _GET_POS(it_start, it) <= stop ) {
                    save_masked_regions( res, w.start(), start );
                    t = ((t<<2)&TRIPLET_MASK) + (converter_( *it )&0x3);

                    if( w.shift_window( t ) ) {
                        done = true;
                        break;
                    }

                    ++it;
                }
            }
        }

        // append the rest of the perfect intervals to the result
        {
            size_type wstart = w.start();

            while( !P.empty() ) {
                save_masked_regions( res, wstart, start );
                ++wstart;
            }
        }

        if( w.start() > 0 ) start += w.start();
        else break;
    }
}

//------------------------------------------------------------------------------
void
CSymDustMasker::operator()( const sequence_type & seq, size_type seq_len, TMaskList& res)
{ (*this)( seq, seq_len, 0, seq_len - 1, res ); }

//------------------------------------------------------------------------------
void CSymDustMasker::GetMaskedLocs( 
    const sequence_type & seq, 
	size_type seq_len,
    std::vector< std::pair<size_type, size_type > >& locs )
{
    // typedef std::vector< CConstRef< objects::CSeq_loc > > locs_type;
	TMaskList res;
    (*this)( seq, seq_len, res );
    locs.clear();
    locs.reserve( res.size() );

    for( TMaskList::const_iterator it = res.begin(); it != res.end(); ++it )
		locs.push_back( std::make_pair(it->first, it->second) );
}

END_NCBI_SCOPE
