#ifndef FSA_FILE_IO_HPP
#define FSA_FILE_IO_HPP

#include <fstream>
#include <zlib.h>
#include <vector>


namespace fsa {

class GzFileReader {
public:
    GzFileReader(const std::string &fname);
    ~GzFileReader() { if (in_ != nullptr) gzclose(in_); }
    bool Valid() const { return in_ != nullptr; }
    bool IsEnd() const { return gzeof(in_); }
    std::string GetStrippedLine();
    const std::string &QueryStrippedLine() {
        if (!HasStrippedLine()) {
            int64_t pos = Tell();   // TODO record pos
            next_line_ = GetStrippedLine();
            next_line_pos_ = pos;
        }
        return next_line_;
    }
    void ConsumeStrippedLine() { next_line_pos_ = -1; }     
    bool HasStrippedLine() const { return next_line_pos_ >= 0; }

    std::string GetLine();
    bool GetLine(std::string &s);
    size_t GetLines(std::vector<std::string> &lines);

    int64_t Tell() { return !HasStrippedLine() ? gztell(in_) : next_line_pos_; }
    void Seek(int64_t pos) { 
        next_line_.clear(); 
        next_line_pos_ = -1;   
        gzseek(in_, pos, SEEK_SET);
    }
    
protected:
    gzFile in_ { nullptr };

    std::string next_line_ { "" };
    int64_t next_line_pos_ {-1} ;
};

class GzFileWriter {
public:
    GzFileWriter(const std::string &fname, bool plain) : plain_(plain) {
        /*

        */
    }

    GzFileWriter(const std::string &fname) {
        plain_ = !(fname.size() >= 3 && fname.substr(fname.size()-3) == ".gz");
        if (plain_) {
            out_plain_ = fopen(fname.c_str(), "wb");
        } else {
            out_compress_ = gzopen(fname.c_str(), "wb");
        }
    }

    void Write(const std::string &str) {
        if (plain_) {
            fputs(str.c_str(), out_plain_);
        } else {
            gzputs(out_compress_, str.c_str());
        }
    }

    ~GzFileWriter() { 
        if (out_compress_ != nullptr) gzclose(out_compress_); 
        if (out_plain_ != nullptr) fclose(out_plain_); 
    }
    
protected:
    gzFile out_compress_ {nullptr};
    FILE* out_plain_ {nullptr};
    bool plain_ ;

};


} // namespace fsa {

#endif // FSA_FILE_IO_HPP

