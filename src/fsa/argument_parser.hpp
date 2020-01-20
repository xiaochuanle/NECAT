#ifndef ARGUMENT_PARSER_HPP
#define ARGUMENT_PARSER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace fsa {

class ArgumentParser {
public:
    ArgumentParser(const std::string &prog="", const std::string &desc="", const std::string &version="") 
        : program_(prog), description_(desc), version_(version) {}

    bool ParseArgument(int argc, const char* const argv[]);
    std::string Usage();
    std::string PrintOptions(const std::string& indent="  ") const;


    template<typename T>
    void AddPositionOption(T &value, const std::string &name, const std::string &desc, const std::string &desired="") {
        positional_options_.push_back(Option(value, name, desc, desired));
    }

    template<typename T>
    void AddNamedOption(T &value, const std::string &name, const std::string &desc, const std::string &desired="", bool (*to)(const std::string&str, T*v) = nullptr) {
        named_options_.push_back(Option(value, name, desc, desired, to));
    }

    bool ParsePositionalOptions(const std::vector<std::string> &positional);
    bool ParseNamedOptions(const std::unordered_map<std::string, std::string> &named);


protected:
    enum class ValueType {
        INT, STRING, DOUBLE, BOOL, LONGLONG,
    };

    struct Option {
        Option(int &v, const std::string &n, const std::string &d, const std::string &de, bool (*toInt)(const std::string &str, int* v) = nullptr) 
            : name(n), desc(d), desired(de), type(ValueType::INT){
                
            value.i = &v;

            toFunc.toInt = toInt == nullptr ? ToInt : toInt;
        }
        Option(std::string &v, const std::string &n, const std::string &d, const std::string &de, bool (*toString)(const std::string &str, std::string* v) = nullptr) 
            : name(n), desc(d), desired(de), type(ValueType::STRING){
            value.s = &v;
            toFunc.toString = toString == nullptr ? ToString : toString;
        }
        Option(double &v, const std::string &n, const std::string &d, const std::string &de, bool (*toDouble)(const std::string &str, double* v) = nullptr) 
            : name(n), desc(d), desired(de), type(ValueType::DOUBLE){
            value.d = &v;
            toFunc.toDouble = toDouble == nullptr ? ToDouble : toDouble;
        }
        Option(bool &v, const std::string &n, const std::string &d, const std::string &de, bool (*toBool)(const std::string &str, bool* v) = nullptr) 
            : name(n), desc(d), desired(de), type(ValueType::BOOL){
            value.b = &v;
            toFunc.toBool = toBool == nullptr ? ToBool : toBool;
        }

        Option(long long &v, const std::string &n, const std::string &d, const std::string &de, bool (*toLonglong)(const std::string &str, long long* v) = nullptr) 
            : name(n), desc(d), desired(de), type(ValueType::LONGLONG){
            value.lli = &v;
            toFunc.toLonglong = toLonglong == nullptr ? ToLonglong : toLonglong;
        }

        std::string name;
        std::string desc;
        std::string desired;
        ValueType type;
        char short_name;
        union {
            long long int *lli;
            int *i;
            std::string *s;
            bool *b;
            double *d;
        } value;

        union {
            bool (*toInt)(const std::string &str, int* v);
            bool (*toString)(const std::string &str, std::string* v);
            bool (*toDouble)(const std::string &str, double* v);
            bool (*toBool)(const std::string &str, bool* v);
            bool (*toLonglong)(const std::string &str, long long* v);
        } toFunc;

        std::string DesiredValue() const {
            return !desired.empty() ? desired : TypeString();
        }

        std::string TypeString() const {
            switch (type) {
            case ValueType::INT: return "INT";
            case ValueType::LONGLONG: return "INT";
            case ValueType::DOUBLE: return "DOUBLE";
            case ValueType::BOOL: return "BOOL";
            case ValueType::STRING: return "STRING";
            default:
                return "";
            }
        }
    };

protected:    
    static bool ToInt(const std::string &str, int *v);
    static bool ToString(const std::string &str, std::string *v);
    static bool ToBool(const std::string &str, bool *v);
    static bool ToDouble(const std::string &str, double *v);
    static bool ToLonglong(const std::string &str, long long int *v);
    bool SetValue(Option& opt, const std::string &value);
protected:
    std::string program_;
    std::string description_;
    std::string version_;


    std::vector<Option> positional_options_;
    std::vector<Option> named_options_;
};

} // namespace fsa {
#endif // ARGUMENT_PARSER_HPP  
