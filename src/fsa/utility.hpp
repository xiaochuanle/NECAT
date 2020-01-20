#ifndef FSA_UTILITY_HPP
#define FSA_UTILITY_HPP

#include <string>
#include <array>
#include <vector>
#include <thread>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <numeric>

namespace fsa {

template<typename T>
auto SplitConstIterater(size_t sz, const T& container) -> std::vector<std::array<typename T::const_iterator, 2>> {
    assert(sz >= 1);

    std::vector<std::array<typename T::const_iterator, 2>> result(sz);

    size_t sub_size = (container.size() + sz - 1) / sz;
    size_t index = 0;
    for (typename T::const_iterator it = container.begin(); it != container.end(); ++it, ++index) {
        if (index % sub_size == 0) {
            result[index/sub_size][0] = it;
        }
    }

    for (size_t i=0; i<result.size()-1; ++i) {
        result[i][1] = result[i+1][0];
    }
    result.back()[1] = container.end();

    return result;
}

template<typename M>
auto SplitMapKeys(size_t sz, const M& map) -> std::vector<std::vector<typename M::key_type>> {
    std::vector<std::vector<typename M::key_type>> results(sz);
    size_t sub_size = (map.size() + sz - 1) / sz;
    typename M::size_type index = 0;
    for (const auto &i : map) {
        results[index++ / sub_size].push_back(i.first);
    }
    return results;
}

template<typename S>
auto SplitSet(size_t sz, const S& set) -> std::vector<S> {
    std::vector<S> results(sz);

    size_t sub_size = (set.size() + sz - 1) / sz;
    typename S::size_type index = 0;
    for (const auto &i : set) {
        results[index++ / sub_size].insert(i);
    }
    return results;
}

template<typename V>
auto SplitVectorKeys(size_t sz, const V& container) ->std::vector<std::array<size_t, 2>> {
    std::vector<std::array<size_t, 2>> result;
    
    size_t bsize = container.size() / sz;
    for (size_t i = 1; i < sz; ++i) {
        result.push_back(std::array<size_t, 2>({ (i - 1)*bsize, i*bsize }));
    }
    result.push_back(std::array<size_t, 2>({ (sz-1)*bsize, container.size() }));

    return result;
}


template<typename T=int>
auto SplitRange(size_t sz, T low, T high) -> std::vector<std::array<T,2>>{
    assert(sz >= 1 && high >= low);
    std::vector<std::array<T,2>> result(sz);
    T inv = (high - low) / sz;

    for (size_t i=0; i<sz-1; ++i) {
        result[i][0] = low + inv*i;
        result[i][1] = low + inv*(i+1);
    }
    result[sz-1][0] = low + inv*(sz-1);
    result[sz-1][1] = high;
    return result;
}

template<typename T>
auto MoveCombineVector(const std::vector<T> &sub_output) -> T {
    T output;
    for (auto &so : sub_output) {
        // TODO output.insert(output.end(), std::make_move_iterator(so.begin()), std::make_move_iterator(so.end()));
        output.insert(output.end(), so.begin(), so.end());
    }
    return output;
}
template<typename T>
auto MoveCombineMapOrSet(const std::vector<T> &sub_output) -> T {
    T output;
    for (auto &so : sub_output) {
        // TODO output.insert(output.end(), std::make_move_iterator(so.begin()), std::make_move_iterator(so.end()));
        output.insert(so.begin(), so.end());
    }
    return output;
}

template<typename W>
void MultiThreadRun(size_t thread_size, W work_func) {
    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(work_func, i));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}

template<typename S, typename W>
void MultiThreadRun(size_t thread_size, S split_func, W work_func) {
    auto sub_inputs = split_func();
    assert(sub_inputs.size() == thread_size);

    std::vector<std::thread> workers;
    for (size_t i = 0; i < sub_inputs.size(); ++i) {
        workers.push_back(std::thread(work_func, std::ref(sub_inputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }

}

// TODO simplify return type
template<typename I, typename S, typename W, typename C>
auto MultiThreadRun(size_t thread_size, const I &inputs, S split_func, W work_func, C combine_func)
-> decltype(combine_func(std::vector<decltype(work_func(split_func(thread_size, inputs)[0]))>(1, work_func(split_func(thread_size, inputs)[0])))) {
    auto sub_inputs = split_func(thread_size, inputs);
    assert(sub_inputs.size() == thread_size);
    std::vector<decltype(work_func(sub_inputs[0]))> sub_outputs(thread_size);

    // TODO "const decltype(sub_inputs[0])& in" get compiler error
    auto thread_work = [&work_func](decltype(sub_inputs[0])& in, decltype(sub_outputs[0])& out) {
        out = work_func(in);
    };

    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(thread_work, std::ref(sub_inputs[i]), std::ref(sub_outputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }

    return combine_func(sub_outputs);
}

template<typename I, typename S, typename W>
void MultiThreadRun(size_t thread_size, I &inputs, S split_func, W work_func) {
    auto sub_inputs = split_func(thread_size, inputs);
    assert(sub_inputs.size() == thread_size);

    std::vector<std::thread> workers;
    for (size_t i = 0; i < thread_size; ++i) {
        workers.push_back(std::thread(work_func, std::ref(sub_inputs[i])));
    }

    for (size_t i = 0; i < workers.size(); ++i) {
        workers[i].join();
    }
}



template<typename T, size_t N>
struct ArrayHash {
    size_t operator()(const std::array<T, N> & r) const {
        size_t h = 17;
        for (size_t i = 0; i < N; ++i) {
            //h ^= std::hash<T>()(r[i]*(T)(i+1)* (T)122777);   // TODO test the efficiency of number
            h = h * 31 + r[i];
        }
        return h;
    }
};


template<typename T, size_t N>
struct ArrayCompare {
    size_t operator()(const std::array<T, N> &a, const std::array<T, N> & b) const {
        for (size_t i = 0; i < N; ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }
};

template<typename M>
auto FindIntersectKeys(const M &a, const M &b) -> std::unordered_set<typename M::key_type> {
    std::unordered_set<typename M::key_type> result;

    for (const auto &i : a) {
        if (b.find(i.first) != b.end()) result.insert(i.first);
    }

    return result;
}


std::vector<std::string> SplitStringBySpace(const std::string &str);

template<typename T>
void DeletePtrContainer(T & c) {
    for (auto e : c) {
        delete e;
    }
    c.clear();
}

size_t FindLongestXHeap(std::vector<int> &lengths, long long base_size);
size_t FindLongestXSort(std::vector<int> &lengths, long long base_size);

template<typename T>
void ComputeMedianAbsoluteDeviation(const std::vector<T>& data_, T &median, T &mad) {
    std::vector<T> data = data_;
    assert(data.size() > 0);
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    median = data[data.size()/2];

    for (size_t i=0; i<data.size(); ++i) {
        data[i] = std::abs(data[i]-median);
    }
    std::nth_element(data.begin(), data.begin() + data.size()/2, data.end());
    mad = data[data.size()/2];
}

template<typename T>
void ComputeMedianAbsoluteDeviation(const std::vector<std::array<T,2>>& data_, T &median, T &mad) {
    std::vector<std::array<T,2>> data = data_;
    assert(data.size() > 0);

    auto find_median = [](std::vector<std::array<T,2>> &data) {
        std::sort(data.begin(), data.end(), [](const std::array<T,2> &a, const std::array<T,2> &b) {
            return a[0] < b[0];
        });

        T total = std::accumulate(data.begin(), data.end(), 0, [](T a, const std::array<T,2> &b) {
            return a + b[1];
        });
        
        T accu = 0;
        size_t index = 0;
        for (; index < data.size(); ++index) {
             accu += data[index][1];
             if (accu >= total / 2) break;
        }

        return data[index][0];
    };

    median = find_median(data);
    for (size_t i=0; i<data.size(); ++i) {
        data[i][0] = std::abs(data[i][0]-median);
    }

    mad = find_median(data);
}



template<typename T>
void ComputeMeanAbsoluteDeviation(std::vector<T>& data, T &mean, T &mad) {
    assert(data.size() > 0);
    mean = std::accumulate(data.begin(), data.end()) / data.size();

    for (size_t i=0; i<data.size(); ++i) {
        data[i] = std::abs(data[i]-mean);
    }
    
    mad = std::accumulate(data.begin(), data.end()) / data.size();;
}


template<typename T>
void ComputeMeanAbsoluteDeviation(std::vector<std::array<T,2>>& data, T &mean, T &mad) {
    assert(data.size() > 0);

    auto find_mean = [](std::vector<std::array<T,2>> &data) {
        std::sort(data.begin(), data.end(), [](const std::array<T,2> &a, const std::array<T,2> &b) {
            return a[0] < b[0];
        });

        T total = std::accumulate(data.begin(), data.end(), (T)0, [](T a, const std::array<T,2> &b) {
            return a + b[0]*b[1];
        });
        
        T weight = std::accumulate(data.begin(), data.end(), (T)0, [](T a, const std::array<T,2> &b) {
            return a + b[1];
        }); 

        return total / weight;
    };

    mean = find_mean(data);
    for (size_t i=0; i<data.size(); ++i) {
        data[i][0] = std::abs(data[i][0]-mean);
    }

    mad = find_mean(data);
}


} // namespace fsa {

#endif // FSA_UTILITY_HPP  
