#ifndef FSA_OVERLAP_SHOW_HPP
#define FSA_OVERLAP_SHOW_HPP


#include<climits>
#include <sstream>
#include <array>

#include "overlap_store.hpp"
#include "argument_parser.hpp"

class OverlapShow {
public:
    OverlapShow();

    bool ParseArgument(int argc, const char *const argv[]);

    void Run();
    void Usage();

protected:
    ArgumentParser GetArgumentParser();



    void LoadOverlaps(const std::string &fname);
    void SaveOverlaps(const std::string &fname);

    std::string ToString(const Overlap &o) const;

protected:
    std::string ifname_;                        //!< 输入文件
    std::string type_;
    std::string read_name_;
    OverlapStore ol_store_;
    

};

#endif // FSA_OVERLAP_SHOW_HPP
