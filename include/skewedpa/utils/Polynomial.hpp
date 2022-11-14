#pragma once
#include <random>

struct Polynomial {
    double alpha = 1;

public:
    using value_type = double;
    using distribution_type = std::uniform_real_distribution<double>;

    Polynomial() : alpha(1) { }

    Polynomial(double alpha_) : alpha(alpha_) { }

    value_type operator() (uint64_t v) const {
        return std::pow(v, alpha);
    }

    static value_type dist_max(value_type weight) {
        return weight;
    }
};