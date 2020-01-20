#ifndef FSA_SEQUENCE_IO_HPP
#define FSA_SEQUENCE_IO_HPP

#include "file_io.hpp"

namespace fsa {

class SeqReader {
public:
    using ItemId = int64_t;
    struct Item {
        std::string head;
        std::string sub_head;
        std::string seq;
        std::string quality;
        ItemId id;            // location in file
    };
public:
    virtual bool Next(Item &item) = 0;
    virtual bool Get(ItemId id, Item &item) = 0;

    static std::string NextNonEmptyLine(std::ifstream &in);
};

class SeqWriter {
public:
    virtual void Write(const SeqReader::Item& item)=0;
};


class FastaReader : public SeqReader {
public:
    FastaReader(const std::string &fname) : in_(fname) {}
    virtual ~FastaReader() {}

    bool IsValid() const { return in_.Valid(); }
    bool IsFileEnd() { return in_.IsEnd(); }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item);}

protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq);

    ItemId Tell() { return in_.Tell(); }
    void Seek(ItemId id) { in_.Seek(id); }

protected:
    GzFileReader in_;
};



class FastqReader : public SeqReader {
public:
    FastqReader(const std::string &fname) : in_(fname) { }
    virtual ~FastqReader() {}

    bool IsValid() const { return in_.Valid(); }
    bool IsFileEnd() { return in_.IsEnd(); }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item); }

    
protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq) { seq = in_.GetStrippedLine(); return true; }
    bool GetHead1() { std::string line = in_.GetStrippedLine(); return line[0] == '+'; }
    bool GetQuality(std::string &qua) { qua = in_.GetStrippedLine(); return true;  }

    ItemId Tell() { return in_.Tell(); }
    void Seek(ItemId id) { in_.Seek(id); }
protected:
    GzFileReader in_;
};


class FastaWriter : public SeqWriter {
public:
    FastaWriter(const std::string &fname) : out_(fname) {
    }
    virtual ~FastaWriter() {}
    
    void Write(const SeqReader::Item& item) {
        out_.Write(">");
        out_.Write(item.head);
        out_.Write(" ");
        out_.Write(item.sub_head);
        out_.Write("\n");
        out_.Write(item.seq);
        out_.Write("\n");
    }

protected:
    GzFileWriter out_;
};

} // namespace fsa {

#endif // FSA_SEQUENCE_IO_HPP
