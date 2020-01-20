#include "logger.hpp"

#include<stdarg.h>  
#include <time.h>

namespace fsa {

Dumper DUMPER;

Dumper::Stream::Stream(const std::string &fname) {
    if (fname != "") {
        file_ = fopen(fname.c_str(), "w");
    }
}

Dumper::Stream::~Stream() {
    if (file_ != nullptr && (file_ != stdout || file_ != stderr)) {
        fclose(file_);
    }
}

void Dumper::Stream::operator()(const char *const format, ...) {
    if (file_ != nullptr) {
        va_list arglist;
        va_start(arglist, format);
        vfprintf(file_, format, arglist);
        va_end(arglist);
        fflush(file_);
    }
}

Dumper::Stream& Dumper::operator [](const std::string &name) {

    if (level_ > 0) {
        auto it = streams_.find(name);
        if (it != streams_.end()) {
            return *it->second;
        }
        else {
            Stream *s = new Stream(name);
            streams_[name] = std::unique_ptr<Stream>(s);
            return *s;
        }
    }
    else {
        return empty_;
    }
}



Logger LOGGER;

void Logger::Stream::operator()(const char *const format, ...) {

    va_list arglist;
    va_start(arglist, format);
    logger_.Log(level_, format, arglist);
    va_end(arglist);

}


Logger::Logger() {
    file_ = stderr;
    SetLevel(level_);
}

void Logger::SetFileName(const std::string &fname) {
    if (file_ != nullptr && file_ != stdout && file_ != stderr) {
        fclose(file_);
    }
    filename_ = fname;

    if (fname != "") {
        file_ = fopen(filename_.c_str(), "w");
    }

    if (file_ == nullptr) {
        file_ = stderr;
    }
}

void Logger::Log(Level level, const char* format, va_list arglist) {

    std::lock_guard<std::mutex> lock(mutex_);
    if (level >= level_) {

        time_t timep;
        time(&timep);
        char tmp[64];
        strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S", localtime(&timep));
        fprintf(file_, "%s", tmp);

        fprintf(file_, " [%s] ", levelname_[level].c_str());

        vfprintf(file_, format, arglist);

        fprintf(file_, "\n");
        fflush(file_);
    }

    if (level == Logger::L_FATAL) {
        abort();
    }
}

} // namespace fsa {