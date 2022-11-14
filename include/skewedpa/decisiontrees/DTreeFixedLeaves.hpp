#include <iostream>
#include <memory>
#include <vector>
#include <tlx/math/round_to_power_of_two.hpp>
#include <tlx/container/btree_map.hpp>
#include <skewedpa/decisiontrees/DTreeBase.hpp>

namespace skewedpa {

template <typename IndexType, typename FunctionType>
class DTreeFixedLeaves final : DTreeBase<IndexType, FunctionType> {
public:
    static constexpr bool bounded = true;

    using base_type         = DTreeBase<IndexType, FunctionType>;
    using index_type        = typename base_type::index_type;
    using weight_type       = typename base_type::weight_type;
    using function_type     = typename base_type::function_type;
    using distribution_type = typename base_type::distribution_type;

    using value_type        = std::pair<index_type, weight_type>;

    explicit DTreeFixedLeaves(index_type sup_index) :
        f_(),
        first_leaf_(tlx::round_up_to_power_of_two(sup_index)),
        guard_leaf_(sup_index),
        tree_(first_leaf_),
        counts_(sup_index)
    {
        tree_one_indexed_ = std::addressof(tree_.front()) - 1;
    }

    explicit DTreeFixedLeaves(FunctionType f, index_type sup_index) :
        f_(f),
        first_leaf_(tlx::round_up_to_power_of_two(sup_index)),
        guard_leaf_(sup_index),
        tree_(first_leaf_),
        counts_(sup_index)
    {
        tree_one_indexed_ = std::addressof(tree_.front()) - 1;
    }

    size_t add(index_type value, int64_t n = 1) {
        assert(value <= guard_leaf_); // out-of-bounds otherwise
        if (value >= guard_leaf_)
            return n;

        assert(value > 0); // groups with value 0 are by design prohibited, they indicate non-existence
        assert(n != 0); // adding nothing is useless

        // update global params
        const auto delta = n * f_(value);
        total_entries_ += n;
        total_weight_ += delta;

        // update counts
        counts_[value] += n;

        auto i = first_leaf_ + value;
        do {
            auto step = [&] {
                const auto parent = i/2;
                const auto is_left = !(i & 1);
                tree_one_indexed_[parent].first += is_left * n;
                tree_one_indexed_[parent].second += is_left * delta;
                i = parent;
            };

            if (i > 8) {
                step(); step(); step();
            }

            // TODO @Manuel you had 3x steps in one?
            step();
        } while (i);
        assert(tree_one_indexed_[first_leaf_ / 2].first == 0);

        return 0;
    }

    size_t increment(index_type value, size_t n = 1) {
        assert(value > 0);
        const auto to_value = value + 1;

        add(value, -n);
        return add(to_value, n);
    }

    template <bool sequential_attachment, typename Generator, typename SingleArgCallback>
    size_t sample(Generator &&gen, SingleArgCallback &&cb)  {
        constexpr bool with_rejection = !sequential_attachment;
        const auto i = sample_group<with_rejection>(std::forward<Generator>(gen));
        cb(i);
        if constexpr(sequential_attachment) {
            return increment(i);
        } else {
            return 0;
        }
    }

    template <bool sequential_attachment>
    size_t update() {
        // if not sequential attachment, correct decision tree after all entries are sampled
        if constexpr(!sequential_attachment) {
            // update decision tree for the neighbors
            size_t num_overflows = 0;
            for (const auto[g, cg] : sampled_indices_) {
                num_overflows += increment(g, cg);
            }

            // clear data structure
            sampled_indices_.clear();

            return num_overflows;
        } else {
            return 0;
        }
    }

    [[nodiscard]] size_t number_of_entries() const {
        return total_entries_;
    }

    [[nodiscard]] weight_type weight() const {
        return total_weight_;
    }

    [[nodiscard]] const std::vector<value_type>& tree() const {
        return tree_;
    }

    [[nodiscard]] const std::vector<size_t>& counts() const {
        return counts_;
    }

    static std::string tree_name() {
        return "DTreeFixedLeaves";
    }

    void print(const std::string& label) const {
        std::cout << label << std::endl;
        for (size_t i = 0; i < counts_.size(); ++i) {
            if (!counts_[i]) continue;
            std::cout << "[" << i << "," << counts_[i] << "," << counts_[i] * f_(i) << "] ";
        }
        std::cout << std::endl;
    }

protected:
    const function_type f_;
    const index_type first_leaf_;
    const index_type guard_leaf_;

    std::vector<value_type> tree_;
    std::vector<size_t> counts_;
    value_type *tree_one_indexed_;
    size_t total_entries_ = 0;
    weight_type total_weight_ = 0;

    tlx::btree_map<index_type, size_t> sampled_indices_;

    template <bool with_rejection, typename Generator>
    index_type sample_group(Generator &&gen) {
        if constexpr(with_rejection) {
            index_type i = 0;
            do {
                i = sample_group<false>(std::forward<Generator>(gen));
                const auto ci = counts_[i];

                if (!sampled_indices_.exists(i)) {
                    sampled_indices_[i]++;
                    break;
                } else {
                    auto accept = std::uniform_int_distribution<size_t>{1, ci};
                    if (accept(gen) > sampled_indices_[i]) {
                        sampled_indices_[i]++;
                        break;
                    }
                }
            } while (true);

            return i;
        } else {
            auto dist = distribution_type{0, FunctionType::dist_max(total_weight_)};
            auto value = dist(gen);

            assert(!counts_[0]);
            index_type i = 1;
            do {
                auto step = [&] {
                    const auto[left_count, left_weight] = tree_one_indexed_[i];
                    tlx::unused(left_count);

                    auto to_right = (value >= left_weight);
                    value -= to_right * left_weight;
                    i = 2*i + to_right;
                };

                if (8*i + 7 < first_leaf_) {
                    step(); step(); step();
                }

                step();
            } while (i < first_leaf_);
            assert(i > first_leaf_);

            return i - first_leaf_;
        }
    }
};

}