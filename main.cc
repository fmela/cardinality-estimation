#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>

#include "murmurhash3.h"

using namespace std;

class CardinalityEstimator
{
public:
    CardinalityEstimator() = default;
    virtual ~CardinalityEstimator() = default;

    virtual void item(const string& s) = 0;
    virtual size_t count() = 0;
};

class ExactCardinalityEstimator : public CardinalityEstimator
{
protected:
    unordered_set<string> items_;

public:
    ExactCardinalityEstimator() = default;

    void item(const string& s) { items_.insert(s); }
    size_t count() { return items_.size(); }
};

namespace PCSA {
    const double kPhi = .77351;

    // Number of ones in the binary representation of x
    uint32_t p(uint32_t x) {
        return __builtin_popcount(x);
    }

    // 2^r(x)
    uint32_t R(uint32_t x) {
        return ~x & (x+1);
    }

    // Number of trailing ones in the binary representation of x
    uint32_t r(uint32_t x) {
        return p(R(x) - 1);
    }
};

class PCSACardinalityEstimator : public CardinalityEstimator
{
protected:
    hash<string> hash_;
    uint32_t sketch_;

public:
    PCSACardinalityEstimator() : CardinalityEstimator(), sketch_(0) {}

    void item(const string& s) {
        sketch_ |= PCSA::R(hash_(s));
    }
    size_t count() { return PCSA::R(sketch_) / PCSA::kPhi; }
};

template<const unsigned M>
class StochasticAveragingCardinalityEstimator : public CardinalityEstimator
{
protected:
    hash<string> hash_;
    uint32_t sketch_[M];

public:
    StochasticAveragingCardinalityEstimator() : CardinalityEstimator() {
        for (unsigned k = 0; k < M; ++k)
            sketch_[k] = 0;
    }

    void item(const string& s) {
        const unsigned k = murmurhash3_32((const uint8_t*)s.c_str(), s.size()) % M;
        sketch_[k] |= PCSA::R(hash_(s));
    }

    size_t count() {
        uint32_t sum = 0;
        for (unsigned k = 0; k < M; ++k)
            sum += PCSA::r(sketch_[k]);
        double mean = sum / (double) M;
        return M * pow(2, mean) / PCSA::kPhi;
    }
};

template<const unsigned M>
class LogLogCardinalityEstimator : public CardinalityEstimator
{
protected:
    hash<string> hash_;
    uint32_t sketch_[M];

public:
    LogLogCardinalityEstimator() : CardinalityEstimator() {
        for (unsigned k = 0; k < M; ++k)
            sketch_[k] = 0;
    }

    void item(const string& s) {
        const unsigned k = murmurhash3_32((const uint8_t*)s.c_str(), s.size()) % M;
        const uint32_t r = PCSA::r(hash_(s));
        if (sketch_[k] < r)
            sketch_[k] = r;
    }

    size_t count() {
        uint32_t sum = 0;
        for (unsigned k = 0; k < M; ++k)
            sum += PCSA::R(sketch_[k]);
        double mean = sum / (double) M;
        return M * pow(2, mean + 1.0) * PCSA::kPhi;
    }
};

template<const unsigned M>
class HyperLogLogCardinalityEstimator : public CardinalityEstimator
{
protected:
    hash<string> hash_;
    uint32_t sketch_[M];

public:
    HyperLogLogCardinalityEstimator() : CardinalityEstimator() {
        for (unsigned k = 0; k < M; ++k)
            sketch_[k] = 0;
    }

    void item(const string& s) {
        const unsigned k = murmurhash3_32((const uint8_t*)s.c_str(), s.size()) % M;
        const uint32_t r = PCSA::r(hash_(s));
        if (sketch_[k] < r)
            sketch_[k] = r;
    }

    size_t count() {
        double sum = 0;
        for (unsigned k = 0; k < M; ++k)
            sum += pow(2, -1.0 - sketch_[k]);
        return M * M * PCSA::kPhi / sum;
    }
};

template<typename generator, typename distribution>
void generate_random_string(string* dest,
                            size_t length, generator& gen, distribution& dis)
{
    dest->resize(length);
    for (size_t i = 0; i < length; ++i) {
        dest->push_back(dis(gen));
    }
}

int main()
{
    const size_t kStringLength = 6;

    // mt19937 generator((random_device())());
    mt19937 generator(0xdeadbeef);
    uniform_int_distribution<> distribution('A', 'z');

    struct {
        string name;
        CardinalityEstimator* estimator;
    } estimators[] = {
        { "   exact", new ExactCardinalityEstimator },
        { "    pcsa", new PCSACardinalityEstimator },

        { "    sa_5", new StochasticAveragingCardinalityEstimator<5> },
        { "   sa_29", new StochasticAveragingCardinalityEstimator<29> },
        { "   sa_73", new StochasticAveragingCardinalityEstimator<73> },
        { "  sa_257", new StochasticAveragingCardinalityEstimator<257> },
        { " sa_1531", new StochasticAveragingCardinalityEstimator<1531> },

        { "    ll_5", new LogLogCardinalityEstimator<5> },
        { "   ll_29", new LogLogCardinalityEstimator<29> },
        { "   ll_73", new LogLogCardinalityEstimator<73> },
        { "  ll_257", new LogLogCardinalityEstimator<257> },
        { " ll_1531", new LogLogCardinalityEstimator<1531> },

        { "   hll_5", new HyperLogLogCardinalityEstimator<5> },
        { "  hll_29", new HyperLogLogCardinalityEstimator<29> },
        { "  hll_73", new HyperLogLogCardinalityEstimator<73> },
        { " hll_257", new HyperLogLogCardinalityEstimator<257> },
        { "hll_1531", new HyperLogLogCardinalityEstimator<1531> },
    };

    string s;
    for (size_t i = 0; i < 5*1000*1000; ++i) {
        generate_random_string(&s, kStringLength, generator, distribution);

        for (auto& e : estimators)
            e.estimator->item(s);
    }

    for (auto& e : estimators)
        cout << e.name << " count: " << e.estimator->count() << endl;

    return 0;
}
