#include "read_tools.hpp"

#include "overlap_store.hpp"
#include "file_io.hpp"

namespace fsa {

ArgumentParser ReadTools::GetArgumentParser() {
    ArgumentParser ap("fsa_rd_tools", "tools about reads", "1.0");

    ap.AddPositionOption(cmd_, "cmd", "command");
    ap.AddNamedOption(ifname_, "ifname", "input file name");
    ap.AddNamedOption(size_, "size", "");
    ap.AddNamedOption(block_size_, "block_size", "");
    ap.AddNamedOption(base_size_, "base_size", "");
    ap.AddNamedOption(ofname_, "ofname", "output file name");
    ap.AddNamedOption(thread_size_, "thread_size", "number of threads");
    ap.AddNamedOption(names_, "names", "names of reads");
    ap.AddNamedOption(min_length_, "min_length", "mininum length of reads");
    ap.AddNamedOption(genome_size_, "genome_size", "genome size");
    ap.AddNamedOption(id2name_, "id2name", "id2name");
    ap.AddNamedOption(discard_illegal_read_, "discard_illegal_read", "the read is discarded if it contains a illegal base");
    return ap;
}

void ReadTools::Running() {
    if (cmd_ == "n50") {
        StatN50(ifname_);
    } else if (cmd_ == "split") {
        Split(ifname_, ofname_, block_size_, base_size_);
    } else if (cmd_ == "split_name") {
        SplitName(ifname_, ofname_, block_size_, base_size_);
    } else if (cmd_ == "longest") {
        Longest(ifname_, ofname_, base_size_, min_length_);
    } else if (cmd_ == "check") {
        Check(ifname_);
    } else {
        LOG(ERROR)("Not recognizing the command: %s", cmd_.c_str());
    }
    
}

void ReadTools::StatN50(const std::string& ifname) {
    LOG(INFO)("Load read file %s", ifname_.c_str());
    ReadStore rs;
    rs.Load(ifname_, "", 4);

    std::vector<int> read_lens(rs.GetIdRange()[1]);
    for (size_t i=0; i<read_lens.size(); ++i) {
        read_lens[i] = rs.GetSeqLength(i);
    }

    std::sort(read_lens.begin(), read_lens.end(), [](int a, int b) { return a > b; });

    long long total_length = genome_size_ == 0 ? std::accumulate(read_lens.begin(), read_lens.end(), (long long)0) : genome_size_;

    std::cout << "Count: " << read_lens.size() << "\n";
    std::cout << "Total: " << total_length << "\n";
    std::cout << "Max: " << read_lens.front() << "\n";
    std::cout << "Min: " << read_lens.back() << "\n";
    
    long long accu = 0;
    int ns[] = { 25, 50, 75};
    size_t ins = 0;
    for (size_t i=0; i<read_lens.size(); ++i) {
        accu += read_lens[i];
        for (; ins < sizeof(ns) / sizeof(ns[0]) && accu > (long long)ns[ins]*1.0/100 * total_length; ++ins) {
            std::cout << "N" << ns[ins] << ": " << read_lens[i] << "\n";
            std::cout << "L" << ns[ins] << ": " << i+1 << "\n";
        }
        if (ins >= sizeof(ns) / sizeof(ns[0])) {
            break;
        }
    }
}

long long FileLength(const std::string &fname) {
    std::ifstream in(fname);
	in.seekg(0, std::ios::end);
	return in.tellg(); 
}

void ReadTools::Split(const std::string &ifname, const std::string &opattern, long long block_size, long long base_size) {
    
    if (block_size == 0) LOG(FATAL)("block_size should be greater than 0");


    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        return result;
    };

    size_t min_length = 0;
    if (base_size > 0) {

        std::vector<int> lengths;
        LoadReadFile(ifname, "", [&lengths](const SeqReader::Item& item) {
            lengths.push_back((int)item.seq.size());
        });

        LOG(INFO)("length size = %zd", lengths.size());
    
        FindLongestXHeap(lengths, base_size);
        min_length = lengths[0];
    }

    int index = 0;
    FastaWriter *writer = nullptr;
    long long accu = 0;
    
    std::vector<int> lengths;
    LoadReadFile(ifname, "", [&](const SeqReader::Item& item) {
        if (writer == nullptr) {
            writer = new FastaWriter(format(opattern, index++));
        }

        assert(accu < block_size);
        
        if (item.seq.size() >= min_length)  {
            writer->Write(item);
            accu += item.seq.size();
        } 
        
        if (accu >= block_size) {
            delete writer;
            writer = nullptr;
            accu = 0;
        }
    });

    if (writer != nullptr) {
        delete writer;
    }

}

void ReadTools::SplitName(const std::string &ifname, const std::string &opattern, long long block_size, long long base_size) {
    
    if (block_size == 0) LOG(FATAL)("block_size should be greater than 0");

    auto format = [](const std::string& pattern, int d) {
        std::string result = pattern;
        std::string s = std::to_string(d);
        result.replace(pattern.find("{}"), 2, s);
        return result;
    };

    size_t min_length = 0;
    if (base_size > 0) {

        std::vector<int> lengths;
        LoadReadFile(ifname, "", [&lengths](const SeqReader::Item& item) {
            lengths.push_back((int)item.seq.size());
        });

        LOG(INFO)("length size = %zd", lengths.size());
    
        FindLongestXHeap(lengths, base_size);
        min_length = lengths[0];
    }


    int index = 0;
    GzFileWriter *writer = nullptr;
    long long accu = 0;
    
    std::vector<int> lengths;
    LoadReadFile(ifname, "", [&](const SeqReader::Item& item) {
        if (writer == nullptr) {
            writer = new GzFileWriter(format(opattern, index++));
        }

        assert(accu < block_size);

        if (item.seq.size() >= min_length)  {
            writer->Write(item.head);
            writer->Write("\n");
            accu += item.seq.size();
        } 

        if (accu >= block_size) {
            delete writer;
            writer = nullptr;
            accu = 0;

        }
    });

    if (writer != nullptr) {
        delete writer;
    }   
}

void ReadTools::Longest(const std::string &ifname, const std::string &ofname, long long base_size, int min_length) {
   
    if (base_size > 0) {

        std::vector<int> lengths;
        LoadReadFile(ifname, "", [&lengths, this](const SeqReader::Item& item) {
            if (!discard_illegal_read_ || DnaSeq::Check(item.seq)) {
                lengths.push_back((int)item.seq.size());
            }
        });

        LOG(INFO)("length size = %zd", lengths.size());
        
        FindLongestXHeap(lengths, base_size);
        if (lengths[0] > min_length)  min_length = lengths[0];
    }

    LOG(INFO)("min_length = %zd", min_length);

    FastaWriter writer(ofname);

    if (id2name_.empty()) {
        LoadReadFile(ifname, "", [min_length, this, &writer](const SeqReader::Item& item) {
            if ((int)item.seq.size() >= min_length && (!discard_illegal_read_ || DnaSeq::Check(item.seq))) {
                writer.Write(item);
            }
        });
    } else {
        GzFileWriter i2n(id2name_);
        size_t index = 0;
        LoadReadFile(ifname, "", [&index, min_length, this, &writer, &i2n](SeqReader::Item& item) {        
            char buff[24];
            sprintf(buff, "%zd ", index);
            i2n.Write(buff);
            i2n.Write(item.head);
            i2n.Write("\n");

            if ((int)item.seq.size() >= min_length && (!discard_illegal_read_ || DnaSeq::Check(item.seq))) {    
                item.head = buff;
                writer.Write(item);
            }
            index++;
        });
    }
    
}



void ReadTools::Check(const std::string &ifname) {

    DnaSerialTable table;
    LoadReadFile(ifname, "", [&table](const SeqReader::Item& item) {
        for (auto c : item.seq) {
            if (table[c] == (decltype(table[c]))-1) {
                LOG(ERROR)("Found bad base %c(%d) in %s", c, (int)c, item.head.c_str());
            }
        }
    });
}

} // namespace fsa
