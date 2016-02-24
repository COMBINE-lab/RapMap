#pragma once

#include <cassert>
#include "common.hpp"

namespace emphf {

    template <typename MemoryModel>
    class bitstream {

    public:

        template <typename T>
        using vector = typename MemoryModel::template vector<T>;

        bitstream()
            : m_bits()
        {}

        bitstream(MemoryModel& mm)
            : m_bits(mm.make_vector(uninitialized_uint64()))
        {}

        void resize(size_t capacity_bits)
        {
            m_bits.resize((capacity_bits + 63) / 64);
        }

        void swap(bitstream& other)
        {
            m_bits.swap(other.m_bits);
        }

        size_t capacity() const
        {
            return m_bits.size() * 64;
        }

        struct reader {
            reader()
            {}

            reader(bitstream const& bs)
                : m_bs(&bs)
                , m_next(0)
                , m_buf(0)
                , m_avail(0)
            {}

            uint64_t read(uint64_t l)
            {
                uint64_t word = m_buf;
                uint64_t left = l;

                if (l > m_avail) {
                    assert(m_next < m_bs->m_bits.size());
                    left = l - m_avail;
                    m_buf = m_bs->m_bits[m_next++];
                    word |= m_buf << m_avail;
                    m_avail = 64;
                }

                m_avail -= left;
                m_buf = (left == 64) ? 0 : (m_buf >> left);
                word &= (l == 64) ? ~0ULL : ~(~0ULL << l);

                return word;
            }

            uint64_t read_unary()
            {
                uint64_t zs = 0;
                while (!m_buf) {
                    zs += m_avail;
                    m_avail = 64;
                    m_buf = m_bs->m_bits[m_next++];
                }

                uint64_t l = __builtin_ctzll(m_buf);
                m_buf >>= l;
                m_buf >>= 1;
                m_avail -= l + 1;
                return zs + l;
            }

            uint64_t read_gamma()
            {
                uint64_t l = read_unary();
                return (((uint64_t(1) << l) | read(l)) - 1);
            }

        private:

            bitstream const* m_bs;
            uint64_t m_next;
            uint64_t m_buf;
            uint64_t m_avail;
        };

        struct writer {
            writer()
            {}

            writer(bitstream& bs)
                : m_bs(&bs)
                , m_next(0)
                , m_buf(0)
                , m_pos_in_buf(0)
            {}

            void write(uint64_t word, uint64_t l)
            {
                assert(m_pos_in_buf < 64);
                m_buf |= word << m_pos_in_buf;
                uint64_t left = l;

                if (m_pos_in_buf + l >= 64) {
                    assert(m_next < m_bs->m_bits.size());
                    m_bs->m_bits[m_next++] = m_buf;
                    left = l - (64 - m_pos_in_buf);
                    if (m_pos_in_buf) {
                        m_buf = word >> (64 - m_pos_in_buf);
                    } else {
                        m_buf = 0;
                    }
                    m_pos_in_buf = 0;
                }

                m_pos_in_buf += left;
            }

            void write_zeros(uint64_t n)
            {
                assert(m_pos_in_buf < 64);
                while (n >= (64 - m_pos_in_buf)) {
                    assert(m_next < m_bs->m_bits.size());
                    m_bs->m_bits[m_next++] = m_buf;
                    m_buf = 0;
                    n -= 64 - m_pos_in_buf;
                    m_pos_in_buf = 0;
                }

                m_pos_in_buf += n;
            }

            void write_unary(uint64_t n)
            {
                write_zeros(n);
                write(1, 1);
            }

            void write_gamma(uint64_t n)
            {
                n += 1;
                uint64_t l = msb(n);
                write_zeros(l);
                uint64_t v = n ^ (1ULL << l);
                write((v << 1) | 1, l + 1);
            }

            void flush()
            {
                assert(m_next < m_bs->m_bits.size());
                m_bs->m_bits[m_next] = m_buf;
            }

        private:

            bitstream* m_bs;
            uint64_t m_next;
            uint64_t m_buf;
            uint64_t m_pos_in_buf;
        };

    private:
        vector<uninitialized_uint64> m_bits;
    };

}
