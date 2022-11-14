#pragma once

#include <stxxl/comparator>
#include <stxxl/sequence>
#include <stxxl/sorter>
#include <stxxl/priority_queue>

using node_t = uint64_t;
using edge_t = std::pair<node_t, node_t>;

template <typename T>
struct Scale {
    static constexpr T K = static_cast<T>(1000LLU);   ///< Kilo (base 10)
    static constexpr T M = K * K;                     ///< Mega (base 10)
    static constexpr T G = K * K * K;                 ///< Giga (base 10)
    static constexpr T P = K * K * K * K;             ///< Peta (base 10)

    static constexpr T Ki = static_cast<T>(1024LLU);  ///< Kilo (base 2)
    static constexpr T Mi = Ki* Ki;                   ///< Mega (base 2)
    static constexpr T Gi = Ki* Ki* Ki;               ///< Giga (base 2)
    static constexpr T Pi = Ki* Ki* Ki* Ki;           ///< Peta (base 2)
};
using IntScale = Scale<int64_t>;
using UIntScale = Scale<uint64_t>;

constexpr size_t SORTER_MEM = 2 * UIntScale::Gi;
constexpr size_t PQ_POOL_MEM = 128 * UIntScale::Mi;
constexpr size_t INTERNAL_PQ_MEM = 128 * UIntScale::Mi;
constexpr size_t MAX_PQ_SIZE = UIntScale::Gi; // is multiplied by 1024 according to docs

namespace skewedpa {

namespace DegreeNodeQuery {

template <bool sequential_attachment>
struct msg {
    size_t time;
    node_t node;
    node_t deg;

    msg() = default;
    msg(size_t time_, node_t node_, node_t deg_, bool pot_selfloop = false) :
        time(time_),
        node(node_ | ~(~static_cast<node_t>(0) >> (pot_selfloop && sequential_attachment))),
        deg(deg_)
    { }
    msg& operator= (const msg&) = default;

    [[nodiscard]] bool possibly_selfloop() const {
        return sequential_attachment && (node & ~(~static_cast<node_t>(0) >> 1));
    }
};

struct fixed_msg {
    size_t time;
    node_t deg;
    fixed_msg() = default;
    fixed_msg(size_t time_, node_t deg_) : time(time_), deg(deg_) { }
    fixed_msg& operator= (const fixed_msg&) = default;
};

struct key_extractor {
    template <bool sequential_attachment>
    auto operator() (msg<sequential_attachment> &a) const {
        return std::tie(a.deg, a.time, a.node);
    }

    template <bool sequential_attachment>
    auto operator() (const msg<sequential_attachment> &a) const {
        return std::tie(a.deg, a.time, a.node);
    }
};

struct fixed_key_extractor {
    auto operator() (fixed_msg &a) const {
        return std::tie(a.deg, a.time);
    }

    auto operator() (const fixed_msg &a) const {
        return std::tie(a.deg, a.time);
    }
};

template <bool sequential_attachment>
using cmp = stxxl::struct_comparator<msg<sequential_attachment>, key_extractor, stxxl::direction::Less, stxxl::direction::Less, stxxl::direction::Less>;
template <bool sequential_attachment>
using sorter_type = stxxl::sorter<msg<sequential_attachment>, cmp<sequential_attachment>>;

using fixed_cmp = stxxl::struct_comparator<fixed_msg, fixed_key_extractor, stxxl::direction::Less, stxxl::direction::Less>;
using fixed_sorter_type = stxxl::sorter<fixed_msg, fixed_cmp>;

}

template <bool sequential_attachment>
inline std::ostream& operator<< (std::ostream &os, const DegreeNodeQuery::msg<sequential_attachment> &m) {
    const bool pot_selfloop = m.possibly_selfloop();
    os << "qry[time: " << m.time << ",\tnode: " << (pot_selfloop ? m.node & (~static_cast<node_t>(0) >> 1) : m.node) << ",\tdeg: " << m.deg << "]";
    if constexpr(sequential_attachment)
        return os;
    else {
        os << "\t<sl?: " << pot_selfloop << ">";
        return os;
    }
}

namespace DegreeNodeInfo {

struct msg {
    size_t time;
    node_t node;
    node_t deg;

    msg() = default;
    msg(size_t time_, node_t node_, node_t deg_) : time(time_), node(node_), deg(deg_) { }
    msg& operator= (const msg&) = default;
};

struct key_extractor {
    auto operator() (msg &a) const {
        return std::tie(a.deg, a.time, a.node);
    }

    auto operator() (const msg &a) const {
        return std::tie(a.deg, a.time, a.node);
    }
};

// pq related definitions
using gt_cmp = stxxl::struct_comparator<
    msg,
    key_extractor,
    stxxl::direction::Greater,
    stxxl::direction::Greater,
    stxxl::direction::Greater
>;

using pq_type = typename stxxl::PRIORITY_QUEUE_GENERATOR<msg, gt_cmp, INTERNAL_PQ_MEM, MAX_PQ_SIZE>::result;

using block_type = typename pq_type::block_type;

using pq_pool_type = foxxll::read_write_pool<block_type>;

enum { pq_pool_mem = (PQ_POOL_MEM / 2) / block_type::raw_size };

// sorter related definitions
using lt_cmp = stxxl::struct_comparator<
    msg,
    key_extractor,
    stxxl::direction::Less,
    stxxl::direction::Less,
    stxxl::direction::Less
>;

using sorter_type = typename stxxl::sorter<msg, lt_cmp>;

// degree ommitted data types
struct deg_omit_msg {
    size_t time;
    node_t node;

    deg_omit_msg() = default;
    deg_omit_msg(size_t time_, node_t node_) : time(time_), node(node_) { }
    deg_omit_msg& operator= (const deg_omit_msg&) = default;
};

struct deg_omit_key_extracter {
    auto operator() (deg_omit_msg &a) const {
        return std::tie(a.time, a.node);
    }

    auto operator() (const deg_omit_msg &a) const {
        return std::tie(a.time, a.node);
    }
};

using deg_omit_sequence_type = stxxl::sequence<deg_omit_msg>;

}

inline std::ostream& operator<< (std::ostream& os, const DegreeNodeInfo::msg &m) {
    os << "inf[time: " << m.time << ",\tnode: " << m.node << ",\tdeg: " << m.deg << "]";
    return os;
}

inline std::ostream& operator<< (std::ostream& os, const DegreeNodeInfo::deg_omit_msg &m) {
    os << "inf[time: " << m.time << ",\tnode: " << m.node << ",\tdeg: -]";
    return os;
}

namespace DegreeNodeSample {

struct msg {
    size_t rand;
    node_t node;

    msg() = default;
    msg(size_t rand_, node_t node_) : rand(rand_), node(node_) { }
    msg& operator= (const msg&) = default;
};
struct key_extractor {
    auto operator() (msg &a) const {
        return std::tie(a.rand, a.node);
    }

    auto operator() (const msg &a) const {
        return std::tie(a.rand, a.node);
    }
};

using cmp = stxxl::struct_comparator<msg, key_extractor, stxxl::direction::Greater, stxxl::direction::Greater>;
using pq_type = typename stxxl::PRIORITY_QUEUE_GENERATOR<msg, cmp, INTERNAL_PQ_MEM, MAX_PQ_SIZE>::result;
using block_type = typename pq_type::block_type;
using pq_pool_type = foxxll::read_write_pool<block_type>;
enum { pq_pool_mem = (PQ_POOL_MEM / 2) / block_type::raw_size };

}

}