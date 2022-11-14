#include <iostream>
#include <deque>
#include <tuple>
#include <tlx/unused.hpp>
#include <tlx/container/btree_map.hpp>
#include <tlx/container/btree_set.hpp>
#include <skewedpa/decisiontrees/DTreeBase.hpp>

namespace skewedpa {

template <typename IndexType, typename FunctionType>
class DTreeDeque final : DTreeBase<IndexType, FunctionType> {
public:
    static constexpr bool bounded = false;
    static constexpr bool assert_thoroughly = true;

    using base_type         = DTreeBase<IndexType, FunctionType>;
    using index_type        = typename base_type::index_type;
    using weight_type       = typename base_type::weight_type;
    using function_type     = typename base_type::function_type;
    using distribution_type = typename base_type::distribution_type;
    struct value_type {
        index_type g;
        int64_t n;
        weight_type w;
        value_type(index_type g_, int64_t n_, weight_type w_) :  g(g_), n(n_), w(w_) { }
        void add(int64_t n_, weight_type w_) { n += n_; w += w_; }
        bool operator== (const value_type &o) const { return g == o.g && n == o.n && w == o.w; }
    };

    template <typename Type>
    explicit DTreeDeque(const Type &o) : f_() { tlx::unused(o); }

    template <typename Type>
    explicit DTreeDeque(FunctionType f, const Type &o) : f_(f) { tlx::unused(o); }

    size_t add(index_type value, int64_t n = 1) {
        const auto it = std::lower_bound(deque_.begin(), deque_.end(), value, [](const auto &a, const auto &b) { return a.g < b; });
        const auto delta = n * f_(value);
        total_entries_ += n;
        total_weight_ += delta;

        // new largest group incoming
        if (it == deque_.end()) {
            deque_.emplace_back(value, n, delta);
        } else {
            auto &group = *it;
            assert(group.n > 0);
            if (group.g == value)
                group.add(n, delta);
            else
                deque_.emplace(it, value, n, delta);
        }

        return 0;
    }

    // TODO optimize, let the decision tree keep track of the sequential attachment process
    void increment(index_type i, size_t n = 1) {
        // must be present
        auto it = std::find_if(deque_.begin(), deque_.end(), [&](const auto &t) { return t.g == i; });
        assert(it != deque_.end());
        
        const auto ix = std::distance(deque_.begin(), it);
        increment_index(ix, n);
    }

    template <bool sequential_attachment, typename Generator, typename SingleArgCallback>
    size_t sample(Generator &&gen, SingleArgCallback &&cb)  {
        if constexpr(assert_thoroughly) assert(std::is_sorted(deque_.begin(), deque_.end(), [](const auto &a, const auto &b) { return a.g < b.g; }));

        constexpr bool with_rejection = !sequential_attachment;
        const auto i = sample_group_index<with_rejection>(std::forward<Generator>(gen));
        const auto group = deque_[i];
        cb(group.g);
        assert(i < deque_.size());
        assert(group.n > 0);

        if constexpr(sequential_attachment) {
            increment_index(i);
        }

        return false;
    }

    template <bool sequential_attachment>
    size_t update() {
        if constexpr(assert_thoroughly) assert(std::is_sorted(deque_.begin(), deque_.end(), [](const auto &a, const auto &b) { return a.g < b.g; }));

        // if not sequential attachment, correct decision tree after all entries are sampled
        if constexpr(!sequential_attachment) {
            // update decision tree for the neighbors
            size_t erased_offset = 0;
            size_t emplaced_offset = 0;

            assert(std::is_sorted(sampled_indices_.begin(), sampled_indices_.end()));
            for (auto [i, ci] : sampled_indices_) {
                assert(static_cast<size_t>(deque_[i + emplaced_offset - erased_offset].n) >= ci);
                const auto [erased, emplaced] = increment_index(i + emplaced_offset - erased_offset, ci);
                erased_offset += erased;
                emplaced_offset += emplaced;
            }

            // clear data structure
            sampled_indices_.clear();
            num_sampled_ = 0;
        }

        return 0;
    }

    [[nodiscard]] size_t number_of_entries() const {
        return total_entries_;
    }

    [[nodiscard]] bool fully_sampled() const {
        return num_sampled_ == number_of_entries();
    }

    [[nodiscard]] weight_type weight() const {
        return total_weight_;
    }

    [[nodiscard]] const std::deque<value_type>& deque() const {
        return deque_;
    };

    static std::string tree_name() {
        return "DTreeDeque";
    }

    void print(const std::string& label) const {
        std::cout << label << std::endl;
        for (const auto group : deque_) {
            std::cout << "[" << group.g << "," << group.n << "," << group.w << "] ";
        }
        std::cout << std::endl;
    }

protected:
    const function_type f_;
    std::deque<value_type> deque_;

    tlx::btree_map<index_type, size_t> sampled_indices_;
    size_t num_sampled_ = 0;
    size_t total_entries_ = 0;
    weight_type total_weight_ = 0;

    template <bool with_rejection, typename Generator>
    size_t sample_group_index(Generator &&gen) {
        if constexpr(with_rejection) {
            size_t i = 0;
            do {
                i = sample_group_index<false>(std::forward<Generator>(gen));
                const auto group = deque_[i];

                if (!sampled_indices_.exists(i)) {
                    sampled_indices_[i]++;
                    ++num_sampled_;
                    break;
                } else {
                    auto accept = std::uniform_int_distribution<size_t>{1, static_cast<size_t>(group.n)};
                    if (accept(gen) > sampled_indices_[i]) {
                        sampled_indices_[i]++;
                        ++num_sampled_;
                        break;
                    }
                }
            } while (true);

            return i;
        } else {
            auto dist = distribution_type{0, FunctionType::dist_max(total_weight_)};
            auto value = dist(gen);

            for (size_t i = 0; i < deque_.size(); ++i) {
                const auto group = deque_[i];
                if (value < group.w)
                    return i;
                else
                    value -= group.w;
            }
            return deque_.size() - 1;
        }
    }

public:
    std::pair<bool, bool> increment_index(index_type ix, int64_t n = 1) {
        const auto group = deque_[ix];
        assert(group.n > 0);

        bool erased = false;
        bool emplaced = false;
        // increment immediate adjacent group
        if (ix + 1ul == deque_.size()) {
            deque_.emplace_back(group.g + 1, n, n * f_(group.g + 1));
        } else {
            const auto [h, ch, wh] = deque_[ix + 1];
            tlx::unused(ch, wh);
            if (h == group.g + 1) {
                deque_[ix + 1].add(n, n * f_(h));
            } else {
                deque_.insert(deque_.begin() + ix + 1, {group.g + 1, n, static_cast<weight_type>(n) * f_(group.g + 1)});
                emplaced = true;
            }
        }

        if (group.n == n) {
            deque_.erase(deque_.begin() + ix);
            erased = true;
        } else {
            deque_[ix].add(-n, -n * f_(group.g));
        }

        total_weight_ = total_weight_ - n * f_(group.g) + n * f_(group.g + 1);

        return std::make_pair(erased, emplaced);
    }
};

}
