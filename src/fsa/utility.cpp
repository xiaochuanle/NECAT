#include "utility.hpp"

#include <algorithm>

namespace fsa {

std::vector<std::string> SplitStringBySpace(const std::string &str) {
    std::vector<std::string> substrs;

    auto IsSpace = [](char c) { return ::isspace(c); };
    auto IsNotSpace = [](char c) { return !::isspace(c); };

    auto begin = std::find_if(str.begin(), str.end(), IsNotSpace);

    while (begin != str.end()) {
        auto end = std::find_if(begin, str.end(), IsSpace);
        substrs.push_back(std::string(begin, end));
        begin = std::find_if(end, str.end(), IsNotSpace);
    }

    return substrs;
}



size_t FindLongestXHeap(std::vector<int> &lengths, long long base_size) {

    assert(base_size > 0);

    long long accu = 0;
    size_t index = 0;
    auto cmp = [](size_t a, size_t b) { return a > b; };
    for (size_t i = 0; i< lengths.size(); ++i) {
        if (accu < base_size || lengths[0] < lengths[i]) {
            accu += lengths[i];
            std::swap(lengths[index], lengths[i]);
            index++;
            std::push_heap(lengths.begin(), lengths.begin()+index, cmp);
        }

        while (accu - lengths[0] >= base_size) {
            accu -= lengths[0];
            std::pop_heap(lengths.begin(), lengths.begin()+index, cmp);
            index--;
        }
    }
    assert(index > 0);
    
    return index;
}

size_t FindLongestXSort(std::vector<int> &lengths, long long base_size) {

    long long accu = 0;

    std::sort(lengths.begin(), lengths.end(), [](int a, int b) { return a > b; });

    for (size_t i = 0; i< lengths.size(); ++i) {
        if (accu < base_size) {
            accu += lengths[i];
        } else {
            return i;
        }
    }

    return lengths.size();
}

} // namespace fsa {
    