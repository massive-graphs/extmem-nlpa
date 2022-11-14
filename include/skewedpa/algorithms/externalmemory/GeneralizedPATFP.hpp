#include <queue>
#include "../../ScopedTimer.hpp"
#include "../../defs.hpp"
#include <networkit/graph/Graph.hpp>
#include <foxxll/mng/read_write_pool.hpp>
#include <stxxl/sequence>
#include <stxxl/sorter>
#include <stxxl/priority_queue>

namespace skewedpa {

template <typename IMDecisionTree, typename Out, size_t semi_external_raw_size = 0>
class GeneralizedPATFP {
public:
    static constexpr bool verbose = false;

    // data structures and types for externally handled requests
    using DegreeNodeQueryMsg                = DegreeNodeQuery::msg<false>;
    using DegreeNodeQuerySorter             = DegreeNodeQuery::sorter_type<false>;

    // data structures and types for externally handled sampling
    using DegreeNodeSamplePQPool            = DegreeNodeSample::pq_pool_type;
    using DegreeNodeSamplePQ                = DegreeNodeSample::pq_type;

    // data structures and types for externally handled existence infos
    using DegreeNodeInfoSorter              = DegreeNodeInfo::sorter_type;
    using DegreeOmmittedNodeInfoMsg         = DegreeNodeInfo::deg_omit_msg;
    using DegreeOmmittedNodeInfoSequence    = DegreeNodeInfo::deg_omit_sequence_type;

    GeneralizedPATFP(const NetworKit::Graph &seed_graph,
                     node_t num_new_nodes,
                     node_t out_deg,
                     IMDecisionTree::function_type f,
                     Out &out) :
    // seed graph related members
        n0_(seed_graph.numberOfNodes()),
        m0_(seed_graph.numberOfEdges()),

        // generation related members
        n_(num_new_nodes),
        d_(out_deg),
        out_(out),

        // decision tree
        dt_(f, seed_graph.numberOfNodes() + n_),

        // data structure for handling externally held requests
        deg_node_qrys_(DegreeNodeQuery::cmp<false>(), SORTER_MEM),

        // data structure for forwarding seed degrees and incoming degrees
        deg_node_infos_(DegreeNodeInfo::lt_cmp(), SORTER_MEM),

        // data structure for forwarding to degree added by one
        deg_node_forwarded_infos_ptr_(std::make_unique<DegreeOmmittedNodeInfoSequence>()),
        deg_node_forwarded_infos_next_ptr_(std::make_unique<DegreeOmmittedNodeInfoSequence>()),

        // data structures for sampling random vertices with same degrees
        deg_node_samples_pool_(DegreeNodeSample::pq_pool_mem, DegreeNodeSample::pq_pool_mem),
        deg_node_samples_ptr_(std::make_unique<DegreeNodeSamplePQ>(deg_node_samples_pool_))
    {
        if constexpr(verbose) std::cout << "Initialize" << std::endl;

        seed_graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
            out.push( { std::min(u, v), std::max(u, v) } );
            if constexpr(verbose) std::cout << "    @out[u: " << std::min(u, v) << ",\tv: " << std::max(u, v) << "]" << std::endl;
        });
        seed_graph.forNodes([&](NetworKit::node u) {
            assert(seed_graph.degree(u) > 0);
            dt_.add(seed_graph.degree(u));
            deg_node_infos_.push( { 0, u, seed_graph.degree(u) } );
        });
    }

    GeneralizedPATFP(const NetworKit::Graph &seed_graph,
                     node_t num_new_nodes,
                     node_t out_deg,
                     Out &out) :
        // seed graph related members
        n0_(seed_graph.numberOfNodes()),
        m0_(seed_graph.numberOfEdges()),

        // generation related members
        n_(num_new_nodes),
        d_(out_deg),
        out_(out),

        // decision tree
        dt_(seed_graph.numberOfNodes() + n_),

        // data structure for handling externally held requests
        deg_node_qrys_(DegreeNodeQuery::cmp<false>(), SORTER_MEM),

        // data structure for forwarding seed degrees and incoming degrees
        deg_node_infos_(DegreeNodeInfo::lt_cmp(), SORTER_MEM),

        // data structure for forwarding to degree added by one
        deg_node_forwarded_infos_ptr_(std::make_unique<DegreeOmmittedNodeInfoSequence>()),
        deg_node_forwarded_infos_next_ptr_(std::make_unique<DegreeOmmittedNodeInfoSequence>()),

        // data structures for sampling random vertices with same degrees
        deg_node_samples_pool_(DegreeNodeSample::pq_pool_mem, DegreeNodeSample::pq_pool_mem),
        deg_node_samples_ptr_(std::make_unique<DegreeNodeSamplePQ>(deg_node_samples_pool_))
    {
        if constexpr(verbose) std::cout << "Initialize" << std::endl;

        seed_graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
            out.push( { std::min(u, v), std::max(u, v) } );
            if constexpr(verbose) std::cout << "    @out[u: " << std::min(u, v) << ",\tv: " << std::max(u, v) << "]" << std::endl;
        });
        seed_graph.forNodes([&](NetworKit::node u) {
            assert(seed_graph.degree(u) > 0);
            dt_.add(seed_graph.degree(u));
            deg_node_infos_.push( { 0, u, seed_graph.degree(u) } );
        });
    }

    template <typename Generator>
    void run(Generator &&gen) {
        std::cout << "Run algorithm" << std::endl;

        // phase 1
        for (node_t new_node = 0; new_node < n_; ++new_node) {
            for (node_t kth_edge = 0; kth_edge < d_; ++kth_edge) {
                // sample and push query
                const auto cb = [&](auto deg) {
                    deg_node_qrys_.push( { new_node * d_ + kth_edge + 1, n0_ + new_node, deg } );
                    if constexpr(verbose) std::cout << "+" << DegreeNodeQueryMsg{ new_node * d_ + kth_edge + 1, n0_ + new_node, deg } << std::endl;
                };
                dt_.template sample<false>(std::forward<Generator>(gen), cb);
            }

            // update the potentially delayed updates
            dt_.template update<false>();

            // add new node to decision tree if non-sequential attachment
            dt_.add(d_);

            // generate existence infos
            deg_node_infos_.push( { new_node * d_ + 1, n0_ + new_node, d_ } );
        }
        deg_node_qrys_.sort();
        deg_node_infos_.sort();
        assert(!deg_node_qrys_.empty());
        assert(!deg_node_infos_.empty());

        //! PHASE 2

        std::cout << "  - fully-external" << std::endl;

        // - fully-external

        // track the number of incoming deg info messages, and switch to semi-external once smaller than a threshold
        size_t deg_node_info_counts = std::numeric_limits<size_t>::max();

        // save last request
        auto prev_qry = DegreeNodeQuery::cmp<false>{}.max_value();

        // distribution related
        size_t max_U = std::numeric_limits<size_t>::max();
        auto dist_U = std::uniform_int_distribution<size_t>{ 0, max_U };

        do {
            // retrieve newest query in degree then in time,
            // if new degree sampling priority queue must be emptied
            const auto qry = *deg_node_qrys_;
            assert(qry.deg > 0);

            if constexpr(verbose) std::cout << "@" << qry << std::endl;

            // processing new degree group
            if (qry.deg != prev_qry.deg) {
                // previous samples no longer needed
                deg_node_samples_ptr_.reset();
                deg_node_samples_ptr_ = std::make_unique<DegreeNodeSamplePQ>(deg_node_samples_pool_);
                assert(deg_node_samples_ptr_->empty());

                // previous existence infos to smaller degree group no longer needed
                while (!deg_node_infos_.empty()) {
                    const auto elem = *deg_node_infos_;
                    if (elem.deg >= qry.deg)
                        break;
                    else {
                        ++deg_node_info_counts;
                        ++deg_node_infos_;
                    }
                }
                // previous forwarded existence infos no longer needed and swap with newly forwarded ones
                deg_node_forwarded_infos_stream_ptr_.reset();
                deg_node_forwarded_infos_ptr_.reset();
                deg_node_forwarded_infos_ptr_ = std::make_unique<DegreeOmmittedNodeInfoSequence>();
                std::swap(deg_node_forwarded_infos_ptr_, deg_node_forwarded_infos_next_ptr_);

                // if new degree group is larger by at least 2 forfeit tracked infos
                if (qry.deg > prev_qry.deg + 1) {
                    deg_node_forwarded_infos_ptr_.reset();
                    deg_node_forwarded_infos_ptr_ = std::make_unique<DegreeOmmittedNodeInfoSequence>();
                }

                deg_node_forwarded_infos_stream_ptr_ = std::make_unique<DegreeOmmittedNodeInfoSequence::stream>(deg_node_forwarded_infos_ptr_->get_stream());

                // break and proceed semi-externally if sufficiently small
                if (qry.deg > d_ && deg_node_info_counts < semi_external_raw_size / sizeof(node_t))
                    break;

                // reset counter
                deg_node_info_counts = 0;
            }

            // retrieve possible nodes to fulfill query from PQ
            while (!deg_node_infos_.empty()) {
                const auto elem = *deg_node_infos_;
                if constexpr(verbose) std::cout << "  @top+" << elem << std::endl;

                if (elem.deg > qry.deg)
                    break;
                assert(elem.deg == qry.deg);
                if (elem.time > qry.time)
                    break;

                // push to sampling priority queue
                deg_node_samples_ptr_->push( { dist_U(gen), elem.node } );

                // pop from info priority queue
                ++deg_node_info_counts;
                ++deg_node_infos_;
            }

            // retrieve possible nodes to fulfill query from Stream
            DegreeOmmittedNodeInfoSequence::stream &deg_node_forwarded_infos_stream = *deg_node_forwarded_infos_stream_ptr_;
            while (!deg_node_forwarded_infos_stream.empty()) {
                const auto stream_front = *deg_node_forwarded_infos_stream;
                if constexpr(verbose) std::cout << "  @stream+" << stream_front << std::endl;

                if (stream_front.time > qry.time)
                    break;

                // push to sampling priority queue
                deg_node_samples_ptr_->push( { dist_U(gen), stream_front.node } );

                // move stream
                ++deg_node_forwarded_infos_stream;
                ++deg_node_info_counts;
            }
            assert(!deg_node_samples_ptr_->empty());

            // get consistent output edge
            assert(!deg_node_samples_ptr_->empty());
            const auto[sample_rand, sample_node] = deg_node_samples_ptr_->top();
            deg_node_samples_ptr_->pop();
            out_.push( {sample_node, qry.node} );
            assert(sample_node <= qry.node);

            if constexpr(verbose) std::cout << "    @out[u: " << sample_node << ",\tv: " << qry.node << "]" << std::endl;

            // forward the node with incremented degree
            deg_node_forwarded_infos_next_ptr_->push_back({qry.time, sample_node } );

            // update sampling distribution, as remaining entries are skewed
            dist_U = std::uniform_int_distribution<size_t>{ sample_rand, max_U };

            // move query entry
            ++deg_node_qrys_;

            // update last previous query
            prev_qry = qry;
        } while (!deg_node_qrys_.empty());

        // already done, skip semi-external processing
        if (deg_node_qrys_.empty()) {
            if constexpr(verbose)
                std::cout << "Expected: " << m0_ + n_ * d_ << "\n"
                          << "  Actual: " << out_.size() << std::endl;
            assert(out_.size() == m0_ + n_ * d_);

            return;
        }

        // - semi-external processing
        std::cout << "  - semi-external" << std::endl;

        std::vector<node_t> deg_nodes(deg_node_info_counts);
        std::vector<DegreeOmmittedNodeInfoMsg> deg_nodes_tfp;
        std::vector<DegreeOmmittedNodeInfoMsg> deg_nodes_tfp_inc;
        size_t deg_nodes_ptr = 0;
        size_t deg_nodes_tfp_ptr = 0;
        size_t deg_nodes_num_sampled = 0;

        // forward degree info messages to semi-external data structures
        DegreeOmmittedNodeInfoSequence::stream deg_node_forwarded_infos_stream = deg_node_forwarded_infos_ptr_->get_stream();
        for (; !deg_node_forwarded_infos_stream.empty(); ++deg_node_forwarded_infos_stream) {
            const auto elem = *deg_node_forwarded_infos_stream;
            deg_nodes_tfp_inc.emplace_back( elem.time, elem.node );
        }

        do {
            // retrieve newest query in degree then in time
            const auto qry = *deg_node_qrys_;

            if constexpr(verbose) std::cout << "@" << qry << std::endl;

            // processing new degree group
            if (qry.deg != prev_qry.deg) {
                // previous existence infos to smaller degree group no longer needed
                // - externally held
                while (!deg_node_infos_.empty()) {
                    const auto elem = *deg_node_infos_;
                    if (elem.deg >= qry.deg)
                        break;
                    else
                        ++deg_node_infos_;
                }

                // - internally held, always have the previous degree added by one
                std::swap(deg_nodes_tfp, deg_nodes_tfp_inc);
                if (qry.deg != prev_qry.deg + 1) {
                    deg_nodes_tfp.clear();
                }
                deg_nodes_tfp_inc.clear();

                // reset
                deg_nodes_ptr = 0;
                deg_nodes_tfp_ptr = 0;
                deg_nodes_num_sampled = 0;
            }

            // retrieve possible nodes to fulfill query
            // - externally held
            while (!deg_node_infos_.empty()) {
                const auto elem = *deg_node_infos_;

                if constexpr(verbose) std::cout << "  @top+" << elem << std::endl;

                if (elem.deg > qry.deg)
                    break;
                assert(elem.deg == qry.deg);
                if (elem.time > qry.time)
                    break;

                // insert to sampling vector
                deg_nodes[deg_nodes_ptr++] = elem.node;

                // pop from info priority queue
                ++deg_node_infos_;
            }

            // - internally held
            auto deg_nodes_tfp_it = deg_nodes_tfp.cbegin() + deg_nodes_tfp_ptr;
            while (deg_nodes_tfp_it != deg_nodes_tfp.cend()) {
                const auto top = *deg_nodes_tfp_it;

                if constexpr(verbose) std::cout << "  @top+" << top << std::endl;

                if (top.time > qry.time)
                    break;

                // insert to sampling vector
                deg_nodes[deg_nodes_ptr++] = top.node;
                ++deg_nodes_tfp_it;
            }
            deg_nodes_tfp_ptr = std::distance(deg_nodes_tfp.cbegin(), deg_nodes_tfp_it);

            //!! sample

            // - if no self-loop possible at least a single entry must exist that can be sampled
            assert(deg_nodes_ptr > 0);
            assert(deg_nodes_num_sampled <= deg_nodes_ptr - 1);
            assert(deg_nodes_ptr - 1 < deg_nodes.size());
            auto sample_dist = std::uniform_int_distribution<size_t>{ deg_nodes_num_sampled, deg_nodes_ptr - 1 };
            if constexpr(verbose) std::cout << "    -params(" << deg_nodes_num_sampled << ", " << (deg_nodes_ptr - 1) << ")" << std::endl;

            [&](auto qry, auto rand_index) {
                const auto sample_node = deg_nodes[rand_index];
                deg_nodes[rand_index] = deg_nodes[deg_nodes_num_sampled++];

                // output immediately
                out_.push({sample_node, qry.node});
                assert(sample_node <= qry.node); // preferential attachment does not connect to the future

                if constexpr(verbose)
                    std::cout << "    @out[u: " << sample_node << ",\tv: " << qry.node << "]" << std::endl;

                // write to vector for next iteration, update sampling structure
                deg_nodes_tfp_inc.emplace_back(qry.time, sample_node);
            }(qry, sample_dist(gen));

            // move query entry
            ++deg_node_qrys_;

            // update last previous query
            prev_qry = qry;
        } while (!deg_node_qrys_.empty());


        if constexpr(verbose)
            std::cout << "Expected: " << m0_ + n_ * d_ << "\n"
                      << "  Actual: " << out_.size() << std::endl;
        assert(out_.size() == m0_ + n_ * d_);
    }

private:
    const node_t n0_;
    const node_t m0_;
    const node_t n_;
    const node_t d_;

    Out& out_;

    IMDecisionTree dt_;

    DegreeNodeQuerySorter deg_node_qrys_;

    DegreeNodeInfoSorter deg_node_infos_;

    std::unique_ptr<DegreeOmmittedNodeInfoSequence> deg_node_forwarded_infos_ptr_;
    std::unique_ptr<DegreeOmmittedNodeInfoSequence::stream> deg_node_forwarded_infos_stream_ptr_;
    std::unique_ptr<DegreeOmmittedNodeInfoSequence> deg_node_forwarded_infos_next_ptr_;

    DegreeNodeSamplePQPool deg_node_samples_pool_;
    std::unique_ptr<DegreeNodeSamplePQ> deg_node_samples_ptr_;
};

}

