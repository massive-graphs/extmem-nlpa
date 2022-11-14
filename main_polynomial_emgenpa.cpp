#include <tlx/cmdline_parser.hpp>

#include <skewedpa/utils/Polynomial.hpp>
#include <skewedpa/decisiontrees/DTreeDeque.hpp>
#include <skewedpa/decisiontrees/DTreeFixedLeaves.hpp>
#include <skewedpa/decisiontrees/DTreeSplit.hpp>

#include <skewedpa/algorithms/externalmemory/GeneralizedPATFP.hpp>
#include <skewedpa/ScopedTimer.hpp>

template <typename ValueType>
struct Identity {
    using value_type = ValueType;
    using distribution_type = std::uniform_int_distribution<value_type>;

    value_type operator() (value_type v) const {
        return v;
    }

    static value_type dist_max(value_type weight) {
        return weight - 1;
    }
};

struct EdgeOutMock {
    size_t count = 0;
public:
    EdgeOutMock() {}

    void push(edge_t e) {
        ++count;
        tlx::unused(e);
    }

    size_t size() const {
        return count;
    }
};

int main(int argc, const char** argv) {
    tlx::CmdlineParser cp;
    cp.set_description("Benchmark for EM-GenPA");

    size_t ell = 10;
    cp.add_size_t("ell", ell, "Number of incoming edges per node");

    size_t ldexp = 36;
    cp.add_size_t("ldexp", ldexp, "Double of lower exponent");

    size_t udexp = 56;
    cp.add_size_t("udexp", udexp, "Double of upper exponent");

    size_t repeats = 3;
    cp.add_size_t("repeats", repeats, "Repeats");

    double alpha = 1.;
    cp.add_double("alpha", alpha, "Polynomial exponent");

    if (!cp.process(argc, argv)) {
        return -1;
    }

    NetworKit::Graph g(2 * ell);
    g.addEdge(0, 2 * ell - 1);
    for (size_t c = 1; c < 2 * ell; ++c) {
        g.addEdge(c - 1, c);
    }

    using DTreeFixedLeavesType = skewedpa::DTreeFixedLeaves<node_t, Polynomial>;
    using DTreeDequeType = skewedpa::DTreeDeque<node_t, Polynomial>;
    using DTreeSplitType = skewedpa::DTreeSplit<node_t, Polynomial, DTreeFixedLeavesType, DTreeDequeType>;
    using AlgoType = skewedpa::GeneralizedPATFP<DTreeSplitType, EdgeOutMock, 16 * sizeof(node_t)>;
    std::random_device rd;
    Polynomial p(alpha);
    for (size_t exp = ldexp; exp <= udexp; ++exp) {
        for (size_t j = 0; j < repeats; ++j) {
            double halfexp = exp / 2.;
            const auto n = static_cast<node_t>(std::pow(2, halfexp));
            EdgeOutMock out_mock;
            AlgoType skewedpatfp(g, n, ell, p, out_mock);

            size_t algoseed = rd();
            std::cout << algoseed << std::endl;
            incpwl::ScopedTimer timer("[seq?: " + std::to_string(false) + ", n: " + std::to_string(n) + ", d: " + std::to_string(ell) + "]\n");
            std::mt19937_64 gen(algoseed);
            skewedpatfp.run(gen);
            std::cout << "Edges produced: " << out_mock.size() << std::endl;
            std::cout << "EM-GenPA, " << alpha << ", " << n << ", " << (n * ell) << ", " << timer.elapsedSeconds() << std::endl;
        }
    }

    return 0;
}
