#ifndef FSA_LOGGER_HPP
#define FSA_LOGGER_HPP


#include <stdio.h>
#include <map>
#include <string>
#include <memory>
#include <vector>
#include <mutex>

namespace fsa {

class Dumper {
public:

    class Stream {
    public:
        Stream(const std::string &fname="");
        ~Stream();
        void operator() (const char* const format, ...);
    protected:
        FILE * file_{ nullptr };
    };

    Stream& operator [](const std::string &name);

    void SetDirectory(const std::string &dir) { directory_ = dir; }
    void SetLevel(int level) { level_ = level; }
protected:
    std::string directory_;
    int level_{ 0 };
    std::map<std::string, std::unique_ptr<Stream>> streams_;
    Stream empty_;
};

extern Dumper DUMPER;

class Logger {
public:
    enum Level {
        L_DEBUG=0,
        L_INFO,
        L_WARNING,
        L_ERROR,
        L_FATAL
    };
    Logger();

    class Stream {
    public:
        Stream(Logger& logger, Level level) :logger_(logger), level_(level){}
        void operator() (const char* const format, ...);
    protected:
        Logger & logger_;
        Level level_;
    };

    void SetFileName(const std::string &fname);
    void SetLevel(Level level) { level_ = level; }
    void Log(Level level, const char* format, va_list arglist);
protected:
    std::mutex mutex_;
    std::string filename_;
    Level level_{ L_INFO };
    std::vector<std::string> levelname_{ "DEBUG", "INFO", "WARNING", "ERROR", "FATAL" };
    FILE *file_;
};

#define LOG(s) Logger::Stream(LOGGER, Logger::L_##s)

extern Logger LOGGER;

} // namespace fsa {

#endif // FSA_LOGGER_HPP
