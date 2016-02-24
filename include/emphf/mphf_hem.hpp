#pragma once

#include <random>
#include <cmath>

#include "bitpair_vector.hpp"
#include "ranked_bitpair_vector.hpp"
#include "perfutils.hpp"
#include "base_hash.hpp"
#include "hypergraph_sorter_seq.hpp"

namespace emphf {

    template <typename BaseHasher, typename OffsetType = uint64_t>
    class mphf_hem {
    public:

        mphf_hem()
        {}

        template <typename MemoryModel, typename Range, typename Adaptor>
        mphf_hem(MemoryModel& mm, size_t n, Range const& input_range,
                 Adaptor adaptor, double gamma = 1.23,
                 size_t log2_expected_bucket = 8)
            : m_n(n)
        {
            using std::get;

            std::mt19937_64 rng(37); // deterministic seed

            auto hashes = mm.make_vector(hash_triple_t());
            hashes.reserve(m_n);

            auto high_bits = int(std::ceil(std::log2(n >> log2_expected_bucket)));
            m_chunk_shift = size_t(sizeof(get<0>(hash_triple_t())) * 8 - high_bits);

            while (true) {
                m_hasher = BaseHasher::generate(rng);

                hashes.clear();

                logger() << "Generating hashes" << std::endl;
                hash_triple_t max_hash;
                for (auto const& val: input_range) {
                    auto h = m_hasher(adaptor(val));
                    hashes.push_back(h);
                    max_hash = std::max(max_hash, h);
                }

                logger() << "Sorting hashes" << std::endl;
                auto sorted_hashes =
                    mm.make_sorter(hashes.begin(),
                                   hashes.end(),
                                   [](hash_triple_t const& h) {
                                       return h;
                                   },
                                   [&](hash_triple_t const& h, size_t k) {
                                       return uint64_t(this->chunk_of(h) * k) >> high_bits;
                                   });

                bitpair_vector bv;
                logger() << "Generating chunk functions" << std::endl;
                if (try_generate_outer_hashes(sorted_hashes, max_hash,
                                              bv, rng, gamma)) {
                    m_bv.build(std::move(bv));
                    break;
                }
            }
        }

        uint64_t size() const
        {
            return m_n;
        }

        BaseHasher const& base_hasher() const
        {
            return m_hasher;
        }

        template <typename T, typename Adaptor>
        uint64_t lookup(T val, Adaptor adaptor)
        {
            using std::get;
            auto hashes = m_hasher(adaptor(val));
            auto chunk = chunk_of(hashes);

            auto const& chunk_hasher = m_chunk_hashers[chunk];
            auto offset = m_offsets[chunk];
            auto nodes_domain = m_offsets[chunk + 1] - offset;

            auto ih = chunk_hasher(hashes);
            auto hd = nodes_domain / 3;

            uint64_t nodes[3] = {offset + get<0>(ih) % hd,
                                 offset + hd + (get<1>(ih) % hd),
                                 offset + 2 * hd + (get<2>(ih) % hd)};

            uint64_t hidx = (m_bv[nodes[0]] + m_bv[nodes[1]] + m_bv[nodes[2]]) % 3;
            return m_bv.rank(nodes[hidx]);
        }

        void swap(mphf_hem& other)
        {
            std::swap(m_n, other.m_n);
            std::swap(m_chunk_shift, other.m_chunk_shift);

            m_hasher.swap(other.m_hasher);
            m_bv.swap(other.m_bv);
            m_offsets.swap(other.m_offsets);
            m_chunk_hashers.swap(other.m_chunk_hashers);
        }

        void save(std::ostream& os) const
        {
            os.write(reinterpret_cast<char const*>(&m_n), sizeof(m_n));
            os.write(reinterpret_cast<char const*>(&m_chunk_shift),
                     sizeof(m_chunk_shift));
            m_hasher.save(os);
            m_bv.save(os);

            uint64_t chunks = m_chunk_hashers.size();
            os.write(reinterpret_cast<char const*>(&chunks),
                     sizeof(chunks));

            os.write(reinterpret_cast<const char*>(m_offsets.data()),
                     (std::streamsize)(sizeof(m_offsets[0]) * (chunks + 1)));

            for (size_t i = 0; i < chunks; ++i) {
                m_chunk_hashers[i].save(os);
            }


        }

        void load(std::istream& is)
        {
            is.read(reinterpret_cast<char*>(&m_n), sizeof(m_n));
            is.read(reinterpret_cast<char*>(&m_chunk_shift),
                    sizeof(m_chunk_shift));
            m_hasher.load(is);
            m_bv.load(is);

            uint64_t chunks;
            is.read(reinterpret_cast<char*>(&chunks),
                    sizeof(chunks));

            m_offsets.resize(chunks + 1);
            is.read(reinterpret_cast<char*>(m_offsets.data()),
                    (std::streamsize)(sizeof(m_offsets[0]) * (chunks + 1)));

            m_chunk_hashers.resize(chunks);
            for (size_t i = 0; i < chunks; ++i) {
                m_chunk_hashers[i].load(is);
            }
        }


    private:

        typedef typename BaseHasher::hash_triple_t hash_triple_t;
        typedef OffsetType offset_type;

        template <typename HashesVector, typename Rng>
        bool try_generate_outer_hashes(HashesVector& hashes,
                                       hash_triple_t max_hash,
                                       bitpair_vector& bv,
                                       Rng& rng, double gamma)
        {
            size_t n_chunks = chunk_of(max_hash) + 1;
            m_offsets.resize(n_chunks + 1);
            m_chunk_hashers.resize(n_chunks);

            auto hash_it = hashes.begin();
            auto last_hash = *hash_it;
            std::vector<hash_triple_t> chunk_hashes;

            while (hash_it != hashes.end()) {
                chunk_hashes.clear();
                chunk_hashes.push_back(*hash_it);
                auto chunk = chunk_of(chunk_hashes.front());

                while (++hash_it != hashes.end() &&
                       chunk_of(*hash_it) == chunk) {
                    if (*hash_it == last_hash) {
                        logger() << "Duplicate hash value found; restarting."
                                 << std::endl;
                        return false;
                    }
                    chunk_hashes.push_back(*hash_it);
                    last_hash = *hash_it;
                }


                while (true) {
                    m_chunk_hashers[chunk] = BaseHasher::generate(rng);
                    if (try_generate_inner_hashes(chunk, bv,
                                                  chunk_hashes.begin(),
                                                  chunk_hashes.end(),
                                                  gamma)) {
                        break;
                    }
                }
            }

            return true;
        }

        template <typename Iterator>
        bool try_generate_inner_hashes(size_t chunk, bitpair_vector& bv,
                                       Iterator begin, Iterator end,
                                       double gamma)
        {
            typedef uint32_t node_t;
            typedef hypergraph_sorter_seq<hypergraph<node_t>> sorter_t;
            sorter_t sorter;
            typedef sorter_t::hyperedge hyperedge;


            auto chunk_size = (size_t)std::distance(begin, end);

            auto chunk_hash_domain = (size_t(std::ceil(double(chunk_size)
                                                       * gamma)) + 2) / 3;
            auto chunk_nodes = chunk_hash_domain * 3;

            auto const& hasher = m_chunk_hashers[chunk];

            auto edge_gen = [&](hash_triple_t t) {
                using std::get;
                auto h = hasher(t);
                return hyperedge((node_t)(get<0>(h) % chunk_hash_domain),
                                 (node_t)(chunk_hash_domain +
                                          (get<1>(h) % chunk_hash_domain)),
                                 (node_t)(2 * chunk_hash_domain +
                                          (get<2>(h) % chunk_hash_domain)));
            };

            if (!sorter.try_generate_and_sort(range(begin, end),
                                              edge_gen,
                                              chunk_size, chunk_hash_domain,
                                              false)) {
                return false;
            }

            size_t offset = bv.size();
            m_offsets[chunk] = (offset_type)offset;
            m_offsets[chunk + 1] = (offset_type)(offset + chunk_nodes);
            bv.resize(m_offsets[chunk + 1]);

            auto peeling_order = sorter.get_peeling_order();

            for (auto edge = peeling_order.first;
                 edge != peeling_order.second;
                 ++edge) {

                uint64_t target = orientation(*edge);
                uint64_t assigned = bv[offset + edge->v1] + bv[offset + edge->v2];

                // "assigned values" must be nonzeros to be ranked, so
                // if the result is 0 we assign 3
                bv.set(offset + edge->v0, ((target - assigned + 9) % 3) ?: 3);
            }

            return true;
        }

        size_t chunk_of(hash_triple_t const& h) const
        {
            return std::get<0>(h) >> m_chunk_shift;
        }

        uint64_t m_n;
        uint64_t m_chunk_shift;
        BaseHasher m_hasher;
        ranked_bitpair_vector m_bv;

        std::vector<offset_type> m_offsets;
        std::vector<BaseHasher> m_chunk_hashers;
    };
}
