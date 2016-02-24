#pragma once

#include "common.hpp"

namespace emphf {

    template <typename MemoryModel>
    class packed_vector {

        template <typename T>
        using vector = typename MemoryModel::template vector<T>;

    public:

        packed_vector()
            : m_size(0)
            , m_k(0)
            , m_mask(~0ULL)
            , m_bits()
        {}

        packed_vector(MemoryModel& mm, uint64_t k, uint64_t n = 0)
            : m_size(0)
            , m_k(k)
            , m_mask(m_k == 64 ? ~0ULL : ~(~0ULL << m_k))
            , m_bits(mm.make_vector(uninitialized_uint64()))
        {
            resize(n);
        }

        void resize(uint64_t n)
        {
            // can only grow, for now
            assert(n >= size());
            m_size = n;
            m_bits.resize((m_size * m_k + 63) / 64);
        }

        size_t size() const
        {
            return m_size;
        }

        void set(uint64_t i, uint64_t val)
        {
            assert(i < size());
            assert(m_k == 64 || (val >> m_k) == 0);

            uint64_t pos = i * m_k;
            uint64_t word_pos = pos / 64;
            uint64_t offset_pos = pos % 64;

            m_bits[word_pos] &= ~(m_mask << offset_pos);
            m_bits[word_pos] |= val << offset_pos;

            uint64_t stored = 64 - offset_pos;
            if (stored < m_k) {
                m_bits[word_pos + 1] &= ~(m_mask >> stored);
                m_bits[word_pos + 1] |= val >> stored;
            }
        }

        uint64_t operator[](uint64_t i) const
        {
            assert(i < size());
            uint64_t pos = i * m_k;
            uint64_t word_pos = pos / 64;
            uint64_t offset_pos = pos % 64;

            uint64_t val = m_bits[word_pos] >> offset_pos;
            uint64_t read = 64 - offset_pos;

            if (read < m_k) {
                val |= m_bits[word_pos + 1] << read;
            }

            val &= m_mask;
            return val;
        }

        void swap(packed_vector& other)
        {
            std::swap(m_size, other.m_size);
            std::swap(m_k, other.m_k);
            std::swap(m_mask, other.m_mask);
            m_bits.swap(other.m_bits);
        }

        void save(std::ostream& os) const
        {
            os.write(reinterpret_cast<char const*>(&m_size), sizeof(m_size));
            os.write(reinterpret_cast<char const*>(&m_k), sizeof(m_k));
            os.write(reinterpret_cast<char const*>(m_bits.data()),
                     (std::streamsize)(sizeof(m_bits[0]) * m_bits.size()));
        }

        void load(std::istream& is)
        {
            is.read(reinterpret_cast<char*>(&m_size), sizeof(m_size));
            is.read(reinterpret_cast<char*>(&m_k), sizeof(m_k));
            m_bits.resize((m_size * m_k + 63) / 64);
            is.read(reinterpret_cast<char*>(m_bits.data()),
                    (std::streamsize)(sizeof(m_bits[0]) * m_bits.size()));
        }

        vector<uint64_t> const& data() const
        {
            return m_bits;
        }

    protected:
        uint64_t m_size;
        uint64_t m_k;
        uint64_t m_mask;
        vector<uninitialized_uint64> m_bits;
    };

}
