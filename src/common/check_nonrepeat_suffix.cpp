#include "check_nonrepeat_suffix.h"

#include "symdust.hpp"

#include <vector>
using namespace std;

extern "C" int
is_nonrepeat_sequence(char* seq, const idx seq_size)
{
    CSymDustMasker masker;
    typedef CSymDustMasker::size_type loc_type;
    typedef vector< pair<loc_type, loc_type> >::iterator loc_iter;
    vector< pair<loc_type, loc_type> > masked_locs;
    masker.GetMaskedLocs(seq, seq_size, masked_locs);
    idx msize = 0;
    for (loc_iter iter = masked_locs.begin(); iter != masked_locs.end(); ++iter) {
        msize += iter->second - iter->first;
    }
    return msize + 200 < seq_size;
}
