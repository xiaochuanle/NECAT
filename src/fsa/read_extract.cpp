#include "read_extract.hpp"

#include <numeric> 
#include <algorithm>
#include <list>
#include <iostream>
#include <cstdio>

#include "logger.hpp"
#include "read_store.hpp"


namespace fsa {
bool ReadExtract::ParseArgument(int argc, const char* const argv[]) {
    return GetArgumentParser().ParseArgument(argc, argv);
}


void ReadExtract::Usage() {
    std::cout << GetArgumentParser().Usage();
}

ArgumentParser ReadExtract::GetArgumentParser() {
    ArgumentParser ap;
    //ap.AddPositionOption(action_, "action", "");
    ap.AddPositionOption(ifname_, "ifname", "input read file name");
    ap.AddPositionOption(ofname_, "ofname", "output read file name");
    ap.AddNamedOption(action_, "action", "");
    ap.AddNamedOption(thread_size_, "thread_size", "Number of threads");
    ap.AddNamedOption(base_size_, "base_size", "");
    return ap;
}

void ReadExtract::Run() {
    PrintArguments();

    LOG(INFO)("Start");

    if (action_ == "longest") {
        ExtractLongest(ifname_, ofname_, base_size_);
    } else if (action_ == "split") {
        Split(ifname_, ofname_, base_size_);
    } else if (action_ == "copy") {
        Copy(ifname_, ofname_);
    } else {

    }

    LOG(INFO)("END");
}


void ReadExtract::ExtractLongest(const std::string &ifname, const std::string &ofname, long long goal) {

    if (goal == 0) LOG(FATAL)("base_size should be greater than 0");

    std::vector<int> lengths;
    LoadReadFile(ifname, "", [&lengths](const SeqReader::Item& item) {
        lengths.push_back((int)item.seq.size());
    });

    LOG(INFO)("length size = %zd", lengths.size());
    size_t min_length = 0;
    
    size_t index = FindLongestXHeap(lengths, goal);
    min_length = lengths[0];
    

    LOG(INFO)("min_length = %zd", min_length);

    FastaWriter writer(ofname);
    
    LoadReadFile(ifname, "", [min_length, &writer](const SeqReader::Item& item) {
        if (item.seq.size() >= min_length) {
            writer.Write(item);
        }
    });
    

}


size_t ReadExtract::FindLongestXHeap(std::vector<int> &lengths, long long goal) {

    assert(goal > 0);

    long long accu = 0;
    size_t index = 0;
    auto cmp = [](size_t a, size_t b) {
        return a > b;
    };
    for (size_t i = 0; i< lengths.size(); ++i) {
        if (accu < goal || lengths[0] < lengths[i]) {
            accu += lengths[i];
            std::swap(lengths[index], lengths[i]);
            index++;
            std::push_heap(lengths.begin(), lengths.begin()+index, cmp);
        }

        while (accu - lengths[0] >= goal) {
            accu -= lengths[0];
            std::pop_heap(lengths.begin(), lengths.begin()+index, cmp);
            index--;
        
        }
    }
    assert(index > 0);
    

    return index;
}

void ReadExtract::Split(const std::string &ifname, const std::string &opattern, long long size) {
    
    if (size == 0) LOG(FATAL)("base_size should be greater than 0");

    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        printf("%s\n", result.c_str());
        return result;
    };

    int index = 0;
    FastaWriter *writer = nullptr;
    long long accu = 0;
    
    std::vector<int> lengths;
    LoadReadFile(ifname, "", [&](const SeqReader::Item& item) {
        if (writer == nullptr) {
            writer = new FastaWriter(format(opattern, index++));
        }

        if (accu < size)  {
            writer->Write(item);
            accu += item.seq.size();

        }

        if (accu >= size) {
            delete writer;
            writer = nullptr;
            accu = 0;
        }
    });

    if (writer != nullptr) delete writer;
    
}

void ReadExtract::Copy(const std::string &ifname, const std::string &ofname) {
    FastaWriter writer(ofname);
    
    LoadReadFile(ifname, "", [&](const SeqReader::Item& item) {
        writer.Write(item);
    });

    
}

void ReadExtract::PrintArguments() {
    LOG(INFO)("Arguments: \n%s", GetArgumentParser().PrintOptions().c_str());
}

} // namespace fsa {
