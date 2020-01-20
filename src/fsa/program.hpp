#ifndef FSA_PROGRAM_HPP
#define FSA_PROGRAM_HPP

#include <iostream>

#include "argument_parser.hpp"
#include "logger.hpp"

namespace fsa {
    
class Program {
public:
    virtual bool ParseArgument(int argc, char *const argv[]) {
        return  GetArgumentParser().ParseArgument(argc, argv);
    }

    virtual void Run() {
        LOG(INFO)("Start");
        PrintArguments();
        Running();
        LOG(INFO)("End");
    };
    
    virtual void Usage() {
        std::cout << GetArgumentParser().Usage();
    }

    virtual void PrintArguments() {
        LOG(INFO)("Arguments: \n%s", GetArgumentParser().PrintOptions().c_str());
    }

    virtual ArgumentParser GetArgumentParser() = 0;
    virtual void Running() = 0;
};

#define PROGRAM_INSTANCE(cls) \
    int main(int argc, char *argv[]) { \
        fsa::cls runner;   \
        if (runner.ParseArgument(argc, argv)) {    \
            runner.Run();   \
        } else {    \
            runner.Usage(); \
        }   \
        return 0;   \
    }

} // namespace fsa {

#endif // FSA_PROGRAM_HPP
