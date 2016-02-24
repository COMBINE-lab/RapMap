#pragma once

#include <cassert>
#include <cstdint>
#include <tuple>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "common.hpp"
#include "hypergraph.hpp"
#include "perfutils.hpp"
#include "packed_edge_list.hpp"
#include "bitstream.hpp"

namespace emphf {

    template <typename NodeType, typename MemoryModel>
    class hypergraph_sorter_scan {
    public:
        typedef hypergraph<NodeType> hg;
        typedef typename hg::node_t node_t;
        typedef typename hg::hyperedge hyperedge;
        typedef typename hg::xored_adj_list xored_adj_list;

        hypergraph_sorter_scan()
            : m_memory_model()
            , m_edges(m_memory_model, 0, 0)
        {}

        template <typename Range, typename EdgeGenerator>
        bool try_generate_and_sort(Range const& input_range,
                                   EdgeGenerator const& edge_gen,
                                   size_t n,
                                   size_t hash_domain)
        {
            auto m = hash_domain * 3;
            uint64_t node_bits = msb(m - 1) + 1;
            logger() << "Using " << node_bits << " bits per node" << std::endl;

            const hyperedge sentinel_edge = hg::sentinel();

            bitstream<MemoryModel> adj_lists(m_memory_model);
            typedef typename bitstream<MemoryModel>::writer bs_writer;
            typedef typename bitstream<MemoryModel>::reader bs_reader;

            auto write_adj_list = [&](bs_writer& wr,
                                      node_t& last_v0,
                                      std::pair<node_t, xored_adj_list> const& adj) {
                assert(adj.first + 1 > last_v0 + 1);
                wr.write_gamma(adj.first - last_v0 - 1);
                last_v0 = adj.first;
                assert(adj.second.degree >= 1);
                wr.write_unary(adj.second.degree - 1);
                wr.write(adj.second.v1s, node_bits);
                wr.write(adj.second.v2s, node_bits);
            };

            auto read_adj_list = [&](bs_reader& rr, node_t& last_v0) {
                node_t delta_v0 = (node_t)rr.read_gamma();
                node_t v0 = (last_v0 += delta_v0 + 1);
                node_t degree = (node_t)rr.read_unary() + 1;
                node_t v1s = (node_t)rr.read(node_bits);
                node_t v2s = (node_t)rr.read(node_bits);
                return std::make_pair(v0, xored_adj_list(degree, v1s, v2s));
            };

            // do all the allocations upfront
            adj_lists.resize(2 * m + // delta-v0
                             3 * n + // unary degrees
                             2 * node_bits * m + // v1 and v2
                             64); // padding

            // to save space, the m_edges vector is used to
            // store all the initial edges, then the round hinges and
            // the deletion lists, and finally the peeling order
            packed_edge_list<MemoryModel, hg>
                (m_memory_model, node_bits,
                 3 * n) // one for each of three orientations
                .swap(m_edges);

            m_round_boundaries.clear();

            auto edges_begin = m_edges.begin();
            auto edges_end = edges_begin;

            // generate canonical edges
            logger() << "Generating hyperedges" << std::endl;
            for (auto const& val: input_range) {
                auto edge = edge_gen(val);
                // canonical by construction
                *edges_end++ = edge;
            }

            // we temporarily the first round hinges after the edge list
            auto round_hinges_begin = edges_end;
            auto round_hinges_end = round_hinges_begin;
            size_t initial_nodes = 0;
            size_t remaining_nodes = 0;
            bs_writer wr(adj_lists);
            node_t last_v0 = -1;

            // we exploit the tripartition by populating the adjacency
            // lists of each orientation of the tripartition
            // separately
            for (unsigned o = 0; o < 3; ++o) {
                logger() << "Sorting " << o << "-orientation edges" << std::endl;
                auto sorted_edges =
                    m_memory_model.make_sorter(edges_begin, edges_end,
                                               [](hyperedge const& e) { return e.v0; },
                                               [&](hyperedge const& e, size_t k) {
                                                   uint64_t b = e.v0 - o * hash_domain;
                                                   return size_t(b * k / hash_domain);
                                               });

                logger() << "Populating " << o << "-orientation lists" << std::endl;

                xored_adj_list adj;
                auto it = sorted_edges.begin();
                hyperedge next_edge = *it;
                auto next_orientation = edges_begin;

                while (it != sorted_edges.end()) {
                    hyperedge edge = next_edge;
                    next_edge = (++it != sorted_edges.end()) ? *it : sentinel_edge;

                    node_t v0 = edge.v0;
                    assert(orientation(edge) == o);

                    adj.add_edge(edge);

                    // if last edge of run of edges starting with same v0,
                    // commit adjacency list
                    if (v0 != next_edge.v0) {
                        initial_nodes += 1;
                        if (adj.degree > 1) {
                            write_adj_list(wr, last_v0, std::make_pair(v0, adj));
                            remaining_nodes += 1;
                        } else {
                            *round_hinges_end++ = edge;
                        }

                        adj = xored_adj_list();
                    }

                    // prepare next orientation
                    if (o == 0) {
                        std::swap(edge.v0, edge.v1);
                        *next_orientation++ = edge;
                    } else if (o == 1) {
                        std::swap(edge.v0, edge.v2);
                        *next_orientation++ = edge;
                    }
                }
            }
            wr.flush();

            auto peeling_order_begin = m_edges.begin();
            auto peeling_order_end = peeling_order_begin;


            // we now move the initial hinges to the beginning of the
            // peeling_order vector, to leave room for the deletion
            // lists (can be as many as 2*n edges)
            round_hinges_end = std::copy(round_hinges_begin, round_hinges_end,
                                         peeling_order_begin);
            round_hinges_begin = peeling_order_begin;


            // iterate rounds until no more hinges are found
            for (size_t round = 0; round_hinges_begin != round_hinges_end; ++round) {
                logger() << "Round " << round << ", "
                         << 100 * double(remaining_nodes) / double(initial_nodes)
                         << "% nodes remaining" << std::endl;

                // sort hinges by their canonical orientation
                auto sorted_hinges =
                    m_memory_model.make_sorter(round_hinges_begin, round_hinges_end,
                                               [](hyperedge const& e) {
                                                   return canonicalize_edge(e);
                                               },
                                               [&](hyperedge const& e, size_t k) {
                                                   uint64_t c_v0 = std::min(e.v0, std::min(e.v1, e.v2));
                                                   return size_t(c_v0 * k / m);
                                               });

                m_round_boundaries.push_back((size_t)(std::distance(peeling_order_begin, peeling_order_end)));
                // we store the deletion list after the round hinges
                auto round_deletion_begin = round_hinges_end;
                auto round_deletion_end = round_deletion_begin;
                // prepare the edge deletion list; in round_hinges,
                // the same edge may appear with different
                // orientations, since we already removed the hinges
                // from their adjacency lists we only need to delete
                // the missing orientations
                bool orientations[] = {false, false, false};
                auto it = sorted_hinges.begin();
                hyperedge next_edge = *it;

                while (it != sorted_hinges.end()) {
                    hyperedge edge = next_edge;
                    next_edge = (++it != sorted_hinges.end()) ? *it : sentinel_edge;

                    orientations[orientation(edge)] = true;
                    auto edge_canonical = canonicalize_edge(edge);

                    if (edge_canonical != canonicalize_edge(next_edge)) {
                        // append the current hinge to the global
                        // peeling order
                        *peeling_order_end++ = edge;

                        // add to the deletion list the missing orientations
                        auto cur_edge = edge_canonical;
                        if (!orientations[0]) *round_deletion_end++ = cur_edge;
                        orientations[0] = false;

                        std::swap(cur_edge.v0, cur_edge.v1);
                        assert(orientation(cur_edge) == 1);
                        if (!orientations[1]) *round_deletion_end++ = cur_edge;
                        orientations[1] = false;

                        std::swap(cur_edge.v0, cur_edge.v2);
                        assert(orientation(cur_edge) == 2);
                        if (!orientations[2]) *round_deletion_end++ = cur_edge;
                        orientations[2] = false;
                    }
                }

                auto sorted_deletions =
                    m_memory_model.make_sorter(round_deletion_begin, round_deletion_end,
                                               [](hyperedge const& e) { return e.v0; },
                                               [&](hyperedge const& e, size_t k) {
                                            return size_t(uint64_t(e.v0) * k / m);
                                        });

                // finally delete the current hinge edges, and find
                // the new hinges
                round_hinges_begin = peeling_order_end;
                round_hinges_end = round_hinges_begin;

                size_t new_remaining_nodes = 0;

                auto cur_del = sorted_deletions.begin();
                bs_reader rr(adj_lists);
                node_t last_v0_r = -1;
                bs_writer wr(adj_lists);
                node_t last_v0_w = -1;

                for (size_t i = 0; i < remaining_nodes; ++i) {
                    auto adj = read_adj_list(rr, last_v0_r);

                    assert(cur_del == sorted_deletions.end() ||
                           hyperedge(*cur_del).v0 >= adj.first);

                    while (cur_del != sorted_deletions.end() &&
                           hyperedge(*cur_del).v0 == adj.first) {
                        adj.second.delete_edge(*cur_del);
                        ++cur_del;
                    }

                    if (adj.second.degree > 1) {
                        write_adj_list(wr, last_v0_w, adj);
                        new_remaining_nodes += 1;
                    } else if (adj.second.degree == 1) {
                        *round_hinges_end++ = adj.second.edge_from(adj.first);
                    }
                }
                remaining_nodes = new_remaining_nodes;
                wr.flush();
            }

            if (remaining_nodes > 0) {
                logger() << "Hypergraph is not peelable" << std::endl;
                return false;
            }

            assert((size_t)std::distance(peeling_order_begin, peeling_order_end) == n);
            m_round_boundaries.push_back(n);

            return true;
        }

        struct peeling_iterator
        {
            typedef hyperedge value_type;

            peeling_iterator()
                : m_sorter(nullptr)
                , m_cur_round(0)
                , m_cur_pos(0)
            {}

            peeling_iterator(hypergraph_sorter_scan const& sorter)
                : m_sorter(&sorter)
                , m_cur_round(m_sorter->m_round_boundaries.size() - 2)
                , m_cur_pos(m_sorter->m_round_boundaries[m_cur_round])
            {
                refresh();
            }

            bool operator!=(peeling_iterator const& other) const
            {
                if (!m_sorter || !other.m_sorter) {
                    return m_sorter != other.m_sorter;
                } else {
                    assert(m_sorter == other.m_sorter);
                    return
                        m_cur_round != other.m_cur_round ||
                        m_cur_pos != other.m_cur_pos;
                }
            }

            peeling_iterator& operator++()
            {
                ++m_cur_pos;
                if (m_cur_pos == m_sorter->m_round_boundaries[m_cur_round + 1]) {
                    if (m_cur_round == 0) {
                        m_sorter = nullptr;
                    } else {
                        m_cur_round -= 1;
                        m_cur_pos = m_sorter->m_round_boundaries[m_cur_round];
                    }
                }

                refresh();
                return *this;
            }

            value_type const& operator*() const
            {
                return m_cur;
            }

            value_type const* operator->() const
            {
                return &m_cur;
            }


        private:

            void refresh()
            {
                if (m_sorter) {
                    m_cur = m_sorter->m_edges[m_cur_pos];
                }
            }

            hypergraph_sorter_scan const* m_sorter;
            size_t m_cur_round;
            size_t m_cur_pos;
            hyperedge m_cur;
        };

        std::pair<peeling_iterator, peeling_iterator>
        get_peeling_order() const
        {
            return std::make_pair(peeling_iterator(*this), peeling_iterator());
        }

    private:

        MemoryModel m_memory_model;
        packed_edge_list<MemoryModel, hg> m_edges;
        std::vector<size_t> m_round_boundaries;
    };
}
