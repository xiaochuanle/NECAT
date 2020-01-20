#include "sequence.hpp"

namespace fsa {

DnaSerialTable DnaSeq::s_SerialTable;
DnaComplementTable DnaSeq::s_ComplTable;

std::string Seq::ReverseComplement(const std::string &seq) {
    std::string r(seq.rbegin(), seq.rend());

    for (size_t i = 0; i < r.length(); ++i) {
        switch (r[i]) {
        case 'A': r[i] = 'T'; break;
        case 'C': r[i] = 'G'; break;
        case 'G': r[i] = 'C'; break;
        case 'T': r[i] = 'A'; break;
        case 'a': r[i] = 't'; break;
        case 'c': r[i] = 'g'; break;
        case 'g': r[i] = 'c'; break;
        case 't': r[i] = 'a'; break;
        default: break;
        }
    }
    return r;
}

} // namespace fsa {
    