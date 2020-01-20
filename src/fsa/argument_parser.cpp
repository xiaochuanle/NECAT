#include "argument_parser.hpp"

#include <cassert>
#include <sstream>
#include <getopt.h>
#include "logger.hpp"

namespace fsa {

std::string ArgumentParser::Usage() {
    const std::string indent="  ";
    std::ostringstream oss;
    oss << description_ << "\n\n";
    oss << "Usage: " << program_ <<" [options] ";
    for (const auto &o : positional_options_) {
        oss << o.name << " ";
    }
    oss << "\n";

    if (positional_options_.size() > 0) {
        //oss << "\n";
        for (const auto &o : positional_options_) {
            oss << indent << o.name  << "\n"
                << indent << "    " << o.desc << "\n";
        }
    }

    if (named_options_.size() > 0) {
        oss << "\nOptions:\n";
        for (const auto &i : named_options_) {
            if (i.type == ValueType::BOOL) {
                oss << indent  << "--" << i.name  <<"\n";

            } else {
            oss << indent  << "--" << i.name  << "=" << (i.DesiredValue().empty() ? i.TypeString() : i.DesiredValue()) <<"\n";

            }
            oss  << indent << "    " << i.desc << "\n";
            switch (i.type) {
            case ValueType::INT:        oss << indent << "    default: "<< *i.value.i << "\n"; break;
            case ValueType::DOUBLE:     oss << indent << "    default: "<< *i.value.d << "\n"; break;
            case ValueType::STRING:     oss << indent << "    default: \""<< *i.value.s << "\"\n"; break;
            case ValueType::BOOL:       break;
            case ValueType::LONGLONG:   oss << indent << "    default: " << *i.value.lli << "\n"; break;
            default:                    assert(!"never come here");
            }

        }
    }

    return oss.str();
}




bool ArgumentParser::ParseArgument(int argc,  const char *const argv[]) {

    std::unordered_map<std::string, std::string> named_arguments_;

    std::vector<std::string> position_arguments_;

    int opt = 0;
    int option_index = 0;
    const char *short_options = "";

    const size_t long_options_capability = 100;
    struct option long_options[long_options_capability];
    assert(long_options_capability >= named_options_.size() + 1);

    for (size_t i = 0; i < named_options_.size(); ++i) {
        long_options[i].name = named_options_[i].name.c_str();
        long_options[i].flag = NULL;
        long_options[i].has_arg = named_options_[i].type != ValueType::BOOL ? required_argument : no_argument;
        long_options[i].val = (int)i;
    }
    long_options[named_options_.size()] = { 0, 0, 0, 0 };

    int oldopterr = opterr;
    opterr = 0; // Don't terminate the program
    while ((opt = getopt_long(argc, (char *const*)argv, short_options, long_options, &option_index)) != -1)
    {
        if (opt != '?') {
            if (opt <(int)(sizeof(long_options) / sizeof(long_options[0]))) {
                if (named_options_[opt].type != ValueType::BOOL) {
                    named_arguments_[named_options_[opt].name] = optarg;
                } else {
                    named_arguments_[named_options_[opt].name] = "";
                }
            }
        } else {
            LOG(WARNING)("unrecognized option %s", argv[optind-1]);   
        }
    }
    opterr = oldopterr;
    for (; optind < argc; optind++) {
        position_arguments_.push_back(argv[optind]);
    }

    return ParsePositionalOptions(position_arguments_) && ParseNamedOptions(named_arguments_);
}



bool ArgumentParser::ParsePositionalOptions(const std::vector<std::string> &positions) {
    if (positions.size() == positional_options_.size()) {
        
        bool r = true;
        for (size_t i=0; r && i<positional_options_.size(); ++i) {
            r = SetValue(positional_options_[i], positions[i]);
        }

        return r;
    }
    else {
        return false;
    }
}

bool ArgumentParser::ParseNamedOptions(const std::unordered_map<std::string, std::string> &named) {
    for (const auto& i : named) {
        for (auto& opt : named_options_) {

            if (opt.name == i.first) {
                if (!SetValue(opt, i.second))
                    return false;
            }
        }
    }
    return true;
}

std::string ArgumentParser::PrintOptions(const std::string& indent) const {
    std::ostringstream oss;
    auto print_opt = [&indent](std::ostringstream &oss, const Option &opt) {
        oss << indent << opt.name << " = ";
        switch (opt.type) {
        case ValueType::INT:        oss << *opt.value.i; break;
        case ValueType::DOUBLE:     oss << *opt.value.d; break;
        case ValueType::STRING:     oss << *opt.value.s; break;
        case ValueType::BOOL:       oss << *opt.value.b; break;
        case ValueType::LONGLONG:   oss << *opt.value.lli; break;
        default:                    assert(!"never come here");
        }

        oss << "\n";
    };

    for (const auto &i : positional_options_) {
        print_opt(oss, i);
    }

    for (const auto &i : named_options_) {
        print_opt(oss, i);
    }

    return oss.str();
}

bool ArgumentParser::SetValue(ArgumentParser::Option &opt, const std::string &value) {
    switch (opt.type) {
    case ValueType::INT:
        return opt.toFunc.toInt(value, opt.value.i);
    case ValueType::STRING:
        return opt.toFunc.toString(value, opt.value.s);
    case ValueType::BOOL:
        return opt.toFunc.toBool(value, opt.value.b);
    case ValueType::DOUBLE:
        return opt.toFunc.toDouble(value, opt.value.d);
    case ValueType::LONGLONG:
        return opt.toFunc.toLonglong(value, opt.value.lli);
    default:
        assert(!"never come here");
        return false;
    }
}

bool ArgumentParser::ToInt(const std::string &str, int *v) {
    *v = std::atoi(str.c_str());
    return true;
}

bool ArgumentParser::ToString(const std::string &str, std::string *v) {
    *v = str;
    return true;
}

bool ArgumentParser::ToBool(const std::string &str, bool *v) {
    *v = true;
    return true;
}

bool ArgumentParser::ToDouble(const std::string &str, double *v) {
    *v = std::atof(str.c_str());
    return true;
}

bool ArgumentParser::ToLonglong(const std::string &str, long long int *v) {
    *v = atoll(str.c_str());
    return true;
}

} // namespace fsa {