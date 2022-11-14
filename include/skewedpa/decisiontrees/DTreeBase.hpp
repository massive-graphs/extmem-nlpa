#pragma once

namespace skewedpa {

template <typename IndexType, typename FunctionType>
class DTreeBase {
public:
    using index_type = IndexType;
    using weight_type = typename FunctionType::value_type;
    using function_type = FunctionType;
    using distribution_type = typename FunctionType::distribution_type;

    virtual size_t add(index_type value, int64_t n) = 0;

    [[nodiscard]] virtual size_t number_of_entries() const = 0;

    virtual weight_type weight() const = 0;
};

}