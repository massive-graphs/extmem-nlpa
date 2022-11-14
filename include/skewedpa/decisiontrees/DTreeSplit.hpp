#include <cmath>
#include <tlx/define/likely.hpp>
#include <skewedpa/decisiontrees/DTreeBase.hpp>

namespace skewedpa {

template <typename IndexType, typename FunctionType, typename DTreeSmall, typename DTreeLarge>
class DTreeSplit : DTreeBase<IndexType, FunctionType> {
public:
    static constexpr bool bounded = false;

    using base_type         = DTreeBase<IndexType, FunctionType>;
    using index_type        = typename base_type::index_type;
    using weight_type       = typename base_type::weight_type;
    using function_type     = typename base_type::function_type;
    using distribution_type = typename base_type::distribution_type;

    explicit DTreeSplit(FunctionType f, index_type sup_index) :
        threshold_(static_cast<index_type>(std::sqrt(sup_index))),
        dt_small_(f, threshold_),
        dt_large_(f, threshold_)
    {
        static_assert(!DTreeLarge::bounded);
    }

    explicit DTreeSplit(index_type sup_index) :
        threshold_(static_cast<index_type>(std::sqrt(sup_index))),
        dt_small_(threshold_),
        dt_large_(threshold_)
    {
        static_assert(!DTreeLarge::bounded);
    }

    size_t add(index_type value, int64_t n = 1) {
        if (value < threshold_) {
            dt_small_.add(value, n);
        } else {
            dt_large_.add(value, n);
        }

        return 0;
    }

    size_t increment(index_type value, size_t n = 1) {
        assert(value > 0);
        if (TLX_LIKELY(value < threshold_ - 1)) {
            dt_small_.increment(value, n);
            return 0;
        }

        if (TLX_LIKELY(value >= threshold_)) {
            dt_large_.increment(value, n);
            return 0;
        }

        dt_small_.add(value, -n);
        dt_large_.add(threshold_, n);

        return 0;
    }

    template <bool sequential_attachment, typename Generator, typename SingleArgCallback>
    size_t sample(Generator &&gen, SingleArgCallback &&cb)  {
        auto dist = distribution_type{0, function_type::dist_max(weight())};
        auto weighted_coin = dist(gen);

        if (weighted_coin < dt_small_.weight() || dt_large_.fully_sampled()) {
            auto overflow = dt_small_.template sample<sequential_attachment>(std::forward<Generator>(gen), std::forward<SingleArgCallback>(cb));
            if (TLX_UNLIKELY(overflow)) {
                dt_large_.add(threshold_, overflow);
            }
        } else {
            dt_large_.template sample<sequential_attachment>(std::forward<Generator>(gen), std::forward<SingleArgCallback>(cb));
        }

        return 0;
    }

    template <bool sequential_attachment>
    size_t update() {
        // if not sequential attachment, correct decision tree after all entries are sampled
        const auto overflows = dt_small_.template update<sequential_attachment>();
        dt_large_.template update<sequential_attachment>();
        if (TLX_UNLIKELY(overflows))
            dt_large_.add(threshold_, overflows);

        return 0;
    }

    [[nodiscard]] size_t number_of_entries() const {
        return dt_small_.number_of_entries() + dt_large_.number_of_entries();
    }

    weight_type weight() const {
        return dt_small_.weight() + dt_large_.weight();
    }

    void print(const std::string& label) const {
        std::cout << label << std::endl;
        dt_small_.print("");
        dt_large_.print("");
    }

    static std::string tree_name() {
        return "DTreeSplit<" + DTreeSmall::tree_name() + ", " + DTreeLarge::tree_name() + ">";
    }

protected:
    const index_type threshold_;
    DTreeSmall dt_small_;
    DTreeLarge dt_large_;
};

}
