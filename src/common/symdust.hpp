/*  $Id: symdust.hpp 464803 2015-04-14 16:29:37Z vakatov $
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
 *   Header file for CSymDustMasker class.
 *
 */

#ifndef C_SYM_DUST_MASKER_HPP
#define C_SYM_DUST_MASKER_HPP

//#include <corelib/ncbitype.h>
//#include <corelib/ncbistr.hpp>
//#include <corelib/ncbiobj.hpp>

//#include <objects/seqloc/Seq_loc.hpp>
//#include <objects/seqloc/Packed_seqint.hpp>
//#include <objmgr/seq_vector.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <memory>
#include <deque>
#include <list>

typedef uint8_t Uint1;
typedef uint32_t Uint4;

//#include "mask_macros.h"

#define BEGIN_NCBI_SCOPE
#define END_NCBI_SCOPE
#define NCBI_XALGODUSTMASK_EXPORT


BEGIN_NCBI_SCOPE

/**
    \brief Looks for low complexity parts of sequences according to the
           symmetric version of DUST algorithm.
 */
class NCBI_XALGODUSTMASK_EXPORT CSymDustMasker
{
    private:

        /** \internal
            \brief Function class responsible for conversion from IUPACNA
                   to NCBI2NA coding.
         */
        struct CIupac2Ncbi2na_converter
        {
            /** \internal
                \brief Operator performing the actual conversion.
                \param r base letter in IOPACNA encoding
                \return the same letter in NCBI2NA encoding
             */
            Uint1 operator()( Uint1 r ) const
            {
                switch( r )
                {
                    case 67: return 1;
                    case 71: return 2;
                    case 84: return 3;
                    default: return 0;
                }
            }
        };

        //typedef objects::CSeqVector seq_t;          /**<\internal Sequence type. */
		typedef char* seq_t;
        typedef CIupac2Ncbi2na_converter convert_t; /**<\internal Converter type. */

    public:

        /**\brief Public sequence type. */
        typedef seq_t sequence_type;    
        /**\brief Integer size type corresponding to sequence_type. */
        //typedef sequence_type::size_type size_type; 
		typedef Uint4 size_type;
        /**\brief Type respresenting an interval selected for masking. */
        typedef std::pair< size_type, size_type > TMaskedInterval;
        /**\brief Type representing a list of masked intervals. */
        typedef std::vector< TMaskedInterval > TMaskList;

        static const Uint4 DEFAULT_LEVEL  = 20; /**< Default value of score threshold. */
        static const Uint4 DEFAULT_WINDOW = 64; /**< Default window size. */
        static const Uint4 DEFAULT_LINKER = 1;  /**< Default value of the longest distance between
                                                     consequtive masked intervals at which they
                                                     should be merged. */

        // These (up to constructor) are public to work around the bug in SUN C++ compiler.

        /** \internal
            \brief Type representing a perfect interval.
         */
        struct perfect
        {
            TMaskedInterval bounds_;    /**<\internal The actual interval. */
            Uint4 score_;               /**<\internal The score of the interval. */
            size_type len_;             /**<\internal The length of the interval. */

            /** \internal   
                \brief Object constructor.
                \param start position of the left end
                \param stop position of the right end
                \param score the score
                \param len the length
             */
            perfect( size_type start, size_type stop, Uint4 score, size_type len )
                : bounds_( start, stop ), score_( score ), len_( len )
            {}
        };

        /**\brief Type representing a list of perfect intervals. */
        typedef std::list< perfect > perfect_list_type;
        /**\brief Table type to store score sum thresholds for each window length. */
        typedef std::vector< Uint4 > thres_table_type;
        /**\brief Type representing a triplet value. */
        typedef Uint1 triplet_type;

        /**\brief Selects the significant bits in triplet_type. */
        static const triplet_type TRIPLET_MASK = 0x3F;

        /**
            \brief Object constructor.
            \param level score threshold
            \param window max window size
            \param linker max distance at which to merge consequtive masked intervals
         */
        CSymDustMasker( Uint4 level = DEFAULT_LEVEL, 
                        size_type window 
                            = static_cast< size_type >( DEFAULT_WINDOW ),
                        size_type linker 
                            = static_cast< size_type >( DEFAULT_LINKER ) );

        /**
            \brief Mask a sequence.
            \param seq a sequence to mask
            \return list of masked intervals
         */
        void operator()( const sequence_type & seq, size_type seq_len, TMaskList& res);

        /**
            \brief Mask a part of the sequence.
            \param seq the sequence to mask
            \param start beginning position of the subsequence to mask
            \param stop ending position of the subsequence to mask
            \return list of masked intervals
         */
        void operator()( const sequence_type & seq,
											   size_type seq_len,
                                              size_type start, size_type stop,
											  TMaskList& res);

        /**
            \brief Mask a sequence and return result as a sequence of CSeq_loc
                   objects.
            \param seq_id sequence id
            \param seq the sequence
            \param [out] vector of const (smart) references to CSeq_loc
         */
        void GetMaskedLocs( 
            const sequence_type & seq,
			size_type seq_len,
            std::vector< std::pair<size_type, size_type> >& locs );

        /**\brief Mask a sequence and return result as a CPacked_seqint
                  instance.
           \param seq_id sequence id
           \param seq the sequence
          */
        ///CRef< objects::CPacked_seqint > GetMaskedInts( 
        ///    objects::CSeq_id & seq_id, const sequence_type & seq );

    private:

        /**\internal Sequence iterator type. */
        //typedef sequence_type::const_iterator seq_citer_type;
		typedef const char* seq_citer_type;

        /** \internal
            \brief Class representing the set of triplets in a window.
         */
        class triplets
        {
            public:
                
                /** \internal
                    \brief Object constructor.
                    \param window max window size
                    \param low_k max triplet multiplicity that guarantees that
                                 the window score is not above the threshold
                    \param perfect_list [in/out] current list of perfect intervals
                    \param thresholds table of threshold values for each window size
                 */
                triplets( size_type window, 
                          Uint1 low_k,
                          perfect_list_type & perfect_list,
                          thres_table_type & thresholds );

                size_type start() const { return start_; }  /**<\internal Get position of the first triplet. */
                size_type stop() const { return stop_; }    /**<\internal Get position of the last triplet. */
                size_type size() const { return triplet_list_.size(); } /**<\internal Get the number of triplets. */

                /** \internal
                    \brief Update the list of perfect intervals with with suffixes
                           of the current window.
                 */
                void find_perfect(); 

                /** \internal
                    \brief Shift the window one base to the right using triplet
                           value t.
                    \param t the triplet value to add to the right end of the
                             triplet list
                    \return false, if the new window contains a single triplet value;
                            true otherwise
                */
                bool shift_window( triplet_type t );

                /** \internal
                    \brief Shift a single triplet window using a triplet value t.
                    \param t the triplet value to add to the right end of the
                             triplet list
                    \return false, if the new window contains a single triplet value;
                            true otherwise
                */
                bool shift_high( triplet_type t );

                /** \internal
                    \brief Check the condition of Proposition 2 allowing
                           to skip the window processing.
                    \return true if the window requires suffix processing,
                            false otherwise
                */
                bool needs_processing() const
                {
                  Uint4 count = stop_ - L; 
                  return count < triplet_list_.size() && 
                         10*r_w > thresholds_[count];
                }

            private:
                
                /**\internal Implementation type for triplets list. */
                typedef std::deque< triplet_type > impl_type;
                /**\internal Triplets list iterator type. */
                typedef impl_type::const_iterator impl_citer_type;
                /**\internal Type for triplet counts tables. */
                typedef Uint1 counts_type[64];

                /** \internal
                    \brief Recompute the value of the running sum
                           and the triplet counts when a new triplet
                           is added.
                    \param r the running sum
                    \param c the triplet counts
                    \param t the new triplet value
                */
                void add_triplet_info( 
                        Uint4 & r, counts_type & c, triplet_type t )
                { r += c[t]; ++c[t]; }

                /** \internal
                    \brief Recompute the value of the running sum
                           and the triplet counts when a triplet
                           is removed.
                    \param r the running sum
                    \param c the triplet counts
                    \param t the triplet value being removed
                */
                void rem_triplet_info( 
                        Uint4 & r, counts_type & c, triplet_type t )
                { --c[t]; r -= c[t]; }

                impl_type triplet_list_;            /**<\internal The triplet list. */

                size_type start_;                   /**<\internal Position of the first triplet in the window. */
                size_type stop_;                    /**<\internal Position of the last triplet in the window. */
                size_type max_size_;                /**<\internal Maximum window size. */

                Uint1 low_k_;                       /**<\internal Max triplet multiplicity that guarantees that
                                                                  that the window score is not above the threshold. */
                Uint4 L;                            /**<\internal Position of the start of the window suffix
                                                                  corresponding to low_k_. */

                perfect_list_type & P;              /**<\internal Current list of perfect subintervals. */
                thres_table_type & thresholds_;     /**<\internal Table containing thresholds for each 
                                                                  value of window length. */

                counts_type c_w;             /**<\internal Table of triplet counts for the whole window. */
                counts_type c_v;             /**<\internal Table of triplet counts for the window suffix. */
                Uint4 r_w;                   /**<\internal Running sum for the whole window. */
                Uint4 r_v;                   /**<\internal Running sum for the window suffix. */
                Uint4 num_diff;              /**<\internal Number of different triplets values. */
        };

        /** \internal
            \brief Merge perfect intervals into the result list.
            \param res the result list
            \param w the list of perfect intervals
            \param start the start position of the subsequence
        */
        void save_masked_regions( 
                TMaskList & res, size_type w, size_type start );

        Uint4 level_;       /**<\internal Score threshold. */
        size_type window_;  /**<\internal Max window size. */
        size_type linker_;  /**<\internal Max distance at which consequtive masked intervals should be merged. */

        Uint1 low_k_;   /**<\internal max triplet multiplicity guaranteeing not to exceed score threshold. */

        perfect_list_type P;            /**<\internal List of perfect intervals within current window. */
        thres_table_type thresholds_;   /**<\internal Table containing score thresholds for each window size. */

        convert_t converter_;   /**\internal IUPACNA to NCBI2NA converter object. */
};

END_NCBI_SCOPE

#endif
