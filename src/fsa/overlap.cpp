#include "overlap.hpp"

#include <sstream>
#include <vector>

#include "utility.hpp"

namespace fsa {

std::string Overlap::ToM4Line() const {
    std::ostringstream oss;

    oss << a_.id + 1 << " " << b_.id + 1 << " " << identity_ << " " << score_ << " "
        << a_.strand << " " << a_.start << " " << a_.end << " " << a_.len << " "
        << b_.strand << " " << b_.start << " " << b_.end << " " << b_.len << "\n";

    return oss.str();
}

bool Overlap::FromM4Line(const std::string &line) {

    std::vector<std::string> items = SplitStringBySpace(line);

    if (items.size() >= 12) {

        // M4文件的Id就是read在fasta文件的序号。
        a_.id = atoi(items[0].c_str()) - 1;
        b_.id = atoi(items[1].c_str()) - 1;

        identity_ = atof(items[2].c_str());
        score_ = atoi(items[3].c_str());

        a_.strand = atoi(items[4].c_str());
        a_.start = atoi(items[5].c_str());
        a_.end = atoi(items[6].c_str());
        a_.len = atoi(items[7].c_str());

        b_.strand = atoi(items[8].c_str());
        b_.start = atoi(items[9].c_str());
        b_.end = atoi(items[10].c_str());
        b_.len = atoi(items[11].c_str());

        // 调整strand，保证a.strand = 0, 先设置b，再设置a
        b_.strand = a_.strand == b_.strand ? 0 : 1;
        a_.strand = 0;

        score_ = -((a_.end - a_.start) + (b_.end - b_.start)) / 2;

        return true;
    }
    else {
        return false;
    }
}


bool Overlap::CheckEnd(int error) const {
    int a_start = a_.strand == 0 ? a_.start : a_.len - a_.end;
    int a_end = a_.strand == 0 ? a_.end : a_.len - a_.start;

    int b_start = b_.strand == 0 ? b_.start : b_.len - b_.end;
    int b_end = b_.strand == 0 ? b_.end : b_.len - b_.start;

    return !((a_start > error && b_start > error) ||
        (a_end < a_.len - error && b_end < b_.len - error));
}

Overlap::Loc Overlap::ReverseLocation(Loc loc, bool direct) {
    switch (loc) {
        case Loc::Left: return direct ? Loc::Right : Loc::Left;
        case Loc::Right: return direct ? Loc::Left : Loc::Right;
        case Loc::Equal: return Loc::Equal;
        case Loc::Containing: return Loc::Contained;
        case Loc::Contained: return Loc::Containing;
        case Loc::Abnormal: return Loc::Abnormal;
        default:    { assert(!"never come here"); return Loc::Abnormal; }
    }
}   

Overlap::Loc Overlap::Location(int err) const {
    auto compare = [err](int a, int b) { 
        if (a < b - err) return -1;
        else if (a > b + err) return 1;
        else return 0;
    };

    auto a_start_loc = compare(a_.start, 0);
    auto a_end_loc = compare(a_.end, a_.len);
    auto b_start_loc = compare(b_.start, 0);
    auto b_end_loc = compare(b_.end, b_.len);
    if (SameDirect()) {

        // a: --------->
        // b:       -------->
        if (a_start_loc > 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc < 0) {
                return Loc::Left;
        }

        // a:       -------->
        // b:  ---------->
        if (a_start_loc == 0 && a_end_loc < 0 && b_start_loc > 0 && b_end_loc == 0) {
                return Loc::Right;
        }

        // a: -------->
        // b: -------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Equal;
        }

        // a: --------->
        // b:    ---->
        if (a_start_loc >= 0 && a_end_loc <= 0 &&(a_start_loc!=0 || a_end_loc !=0) && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Containing;
        }

        // a:   ----->
        // b: ---------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc >= 0 && b_end_loc <= 0 && (b_start_loc != 0 || b_end_loc != 0)) {
            return Loc::Contained;
        }
    } else {

        // query :  <---------
        // target:       -------->
        if (a_start_loc == 0 && a_end_loc < 0 && b_start_loc == 0 && b_end_loc < 0) {
                return Loc::Left;
        }

        // query :       <---------
        // target:  ---------->
        if (a_start_loc > 0 && a_end_loc == 0 && b_start_loc > 0 && b_end_loc == 0) {
                return Loc::Right;
        }

        // query : <--------
        // target: -------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Equal;
        }

        // query : <---------
        // target:    ---->
        if (a_start_loc >= 0 && a_end_loc <= 0 &&(a_start_loc != 0 || a_end_loc != 0)&& b_start_loc == 0 && b_end_loc == 0) {
            return Loc::Containing;
        }

        // query :   <----
        // target: ---------->
        if (a_start_loc == 0 && a_end_loc == 0 && b_start_loc >= 0 && b_end_loc <= 0 && (b_start_loc != 0 || b_end_loc != 0)) {
            return Loc::Contained;
        }
    }

    return Loc::Abnormal;
}


Overlap::Loc Overlap::Location(Seq::Id id, int err) const {
    assert(id == b_.id || id == a_.id);
    Loc loc = Location(err);
    return id == a_.id ? loc : ReverseLocation(loc, SameDirect());
}

bool Overlap::IsContaining(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Equal || loc == Loc::Containing;
}

bool Overlap::IsContained(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Equal || loc == Loc::Contained;
}


bool Overlap::IsProper(int err) const {
    Loc loc = Location(err);
    return loc == Loc::Left || loc == Loc::Right;
}


std::array<int, 2> Overlap::Overhang(const Overlap &o, Loc loc) {
    std::array<int,2> overhang({-1, -1}); // -1 mean no overhang
    
    if (o.SameDirect()) {
        if (loc == Overlap::Loc::Left) {
            //int oh = std::max(o.b_.start-0, o.a_.len-o.a_.end);
            overhang[0] = o.a_.len - o.a_.end;
            overhang[1] = o.b_.start;
        } else if (loc == Overlap::Loc::Right) {
            //int oh = std::max(o.a_.start-0, o.b_.len-o.b_.end);
            overhang[0] = o.a_.start;
            overhang[1] = o.b_.len - o.b_.end;
        } else if (loc == Overlap::Loc::Contained) {
            //int oh = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
        } else if (loc == Overlap::Loc::Containing) {
            //int oh = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        } else if (loc == Overlap::Loc::Equal) {
            //int oh1 = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            //int oh2 = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            //int oh = std::max(oh1, oh2);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        }
    } else {
        if (loc == Overlap::Loc::Left) {
            //int oh = std::max(o.b_.start-0, o.a_.start-0);
            overhang[0] = o.a_.start;
            overhang[1] = o.b_.start;
        } else if (loc == Overlap::Loc::Right) {
            //int oh = std::max(o.a_.len-o.a_.end, o.b_.len-o.b_.end);
            overhang[0] = o.a_.len - o.a_.end;
            overhang[1] = o.b_.len - o.b_.end;
        } else if (loc == Overlap::Loc::Contained) {
            //int oh = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
        } else if (loc == Overlap::Loc::Containing) {
            //int oh = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        } else if (loc == Overlap::Loc::Equal) {
            //int oh1 = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            //int oh2 = std::max(o.b_.start-0, o.b_.len-o.b_.end);
            //int oh = std::max(oh1, oh2);
            overhang[0] = std::max(o.a_.start-0, o.a_.len-o.a_.end);
            overhang[1] = std::max(o.b_.start-0, o.b_.len-o.b_.end);
        }
    }
    return overhang;   
}


std::array<int, 2> Overlap::Overhang() const {
    std::array<int,2> result;
    int a0 = a_.start - 0;
    int a1 = a_.len - a_.end;
    int b0 = b_.start - 0;
    int b1 = b_.len - b_.end;

    if (SameDirect()) {
        std::vector<int> oh(4);
        oh[0] = std::max(a0, b1);
        oh[1] = std::max(a1, b0);
        oh[2] = std::max(a0, a1);
        oh[3] = std::max(b0, b1);

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{a0, b1};
        } else if (min_oh == 1) {
            result = std::array<int,2>{a1, b0};
        } else if (min_oh == 2) {
            result = std::array<int,2> {oh[2], -1};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {-1, oh[3]};
        }
    } else {
        std::vector<int> oh(4);
        oh[0] = std::max(a0, b0);
        oh[1] = std::max(a1, b1);
        oh[2] = std::max(a0, a1);
        oh[3] = std::max(b0, b1);

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{a0, b0};
        } else if (min_oh == 1) {
            result = std::array<int,2>{a1, b1};
        } else if (min_oh == 2) {
            result = std::array<int,2> {oh[2], -1};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {-1, oh[3]};
        }
    }
    //printf("-- %s", ToM4Line().c_str());
    //printf("   %d, %d\n", result[0], result[1]);
    return result;
}

std::array<int, 2> Overlap::Overhang2() const {
    std::array<int,2> result;   // 0表示没有判断为overhang，1表示左端为overhnag 2表示右端为overhang 3表示两端为overhang
    int a0 = a_.start - 0;
    int a1 = a_.len - a_.end;
    int b0 = b_.start - 0;
    int b1 = b_.len - b_.end;

    if (SameDirect()) {
        std::vector<int> oh(4);
        oh[0] = std::max(a0, b1);
        oh[1] = std::max(a1, b0);
        oh[2] = std::max(a0, a1);
        oh[3] = std::max(b0, b1);

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{1, 2};
        } else if (min_oh == 1) {
            result = std::array<int,2>{2, 1};
        } else if (min_oh == 2) {
            result = std::array<int,2> {3, 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, 3};
        }
    } else {
        std::vector<int> oh(4);
        oh[0] = std::max(a0, b0);
        oh[1] = std::max(a1, b1);
        oh[2] = std::max(a0, a1);
        oh[3] = std::max(b0, b1);

        auto min_oh = std::min_element(oh.begin(), oh.end()) - oh.begin();
        if (min_oh == 0) {
            result = std::array<int,2>{1, 1};
        } else if (min_oh == 1) {
            result = std::array<int,2>{2, 2};
        } else if (min_oh == 2) {
            result = std::array<int,2> {3, 0};
        } else {
            assert(min_oh == 3);
            result = std::array<int,2> {0, 3};
        }
    }
    return result;
}


} // namespace fsa {