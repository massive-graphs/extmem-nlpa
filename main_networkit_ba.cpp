#include <tlx/cmdline_parser.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/generators/BarabasiAlbertGenerator.hpp>
#include <skewedpa/ScopedTimer.hpp>
#include <skewedpa/defs.hpp>

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

    if (!cp.process(argc, argv)) {
        return -1;
    }

    NetworKit::Graph g(2 * ell);
    g.addEdge(0, 2 * ell - 1);
    for (size_t c = 1; c < 2 * ell; ++c) {
        g.addEdge(c - 1, c);
    }

    for (size_t exp = ldexp; exp <= udexp; ++exp) {
        for (size_t j = 0; j < repeats; ++j) {
            double halfexp = exp / 2.;
            const auto n = static_cast<node_t>(std::pow(2, halfexp));
            NetworKit::BarabasiAlbertGenerator nk_ba(ell, n, g);
            incpwl::ScopedTimer timer(
                "[seq?: " + std::to_string(false) + ", n: " + std::to_string(n) + ", d: " + std::to_string(ell) +
                "]\n");
            nk_ba.generate();
            std::cout << "NK-BA, " << n << ", " << (n * ell) << ", " << timer.elapsedSeconds() << std::endl;
        }
    }

    return 0;
}
