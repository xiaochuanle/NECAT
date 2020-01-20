#include "truncate_end_spaces.h"

#include <cassert>
#include <cctype>
#include <string>

using namespace std;

extern "C"
void truncate_end_spaces(char* line)
{
    string tmp_line(line);
    string::size_type from = 0, size = tmp_line.size(), to = size;
    if (tmp_line.empty()) {
        line[0] = '\0';
        return;
    }

    while (from < size) {
        if (!isspace(tmp_line[from])) break;
        ++from;
    }
    if (from == size) {
        line[0] = '\0';
        return;
    }

    while (to > from) {
        if (!isspace(tmp_line[to-1])) break;
        --to;
    }
    assert(to > from);

    size_t i = 0;
    for (string::size_type k = from; k < to; ++k) {
        line[i++] = tmp_line[k];
    }
    assert(i > 0);
    line[i] = '\0';
}