#pragma once

#include <stddef.h>

#include "packed_vector.hpp"

namespace emphf {

    template <typename MemoryModel, typename Hypergraph>
    class packed_edge_list {

    protected:
        // forward declarations
        struct value_reference;

    public:
        typedef typename Hypergraph::hyperedge value_type;
        typedef typename Hypergraph::node_t node_t;

        packed_edge_list()
        {}

        packed_edge_list(MemoryModel& mm, uint64_t bits, size_t n)
            : m_bits(bits)
            , m_n(n)
            , m_v(mm, m_bits, n * 3)
        {}

        size_t size() const
        {
            return m_n;
        }

        void set(uint64_t pos, value_type v)
        {
            assert(pos < size());
            m_v.set(pos * 3 + 0, v.v0);
            m_v.set(pos * 3 + 1, v.v1);
            m_v.set(pos * 3 + 2, v.v2);
        }

        value_type get(uint64_t pos) const
        {
            assert(pos < size());
            node_t v0, v1, v2;
            v0 = (node_t)m_v[pos * 3 + 0];
            v1 = (node_t)m_v[pos * 3 + 1];
            v2 = (node_t)m_v[pos * 3 + 2];

            return value_type(v0, v1, v2);
        }

        value_reference operator[](uint64_t pos)
        {
            return *iterator(this, pos);
        }

        value_type operator[](uint64_t pos) const
        {
            return get(pos);
        }

        struct iterator : public std::iterator<std::random_access_iterator_tag,
                                               value_type>
        {
            iterator()
            {}

            bool operator==(iterator const& rhs) const
            {
                assert(m_l == rhs.m_l);
                return m_pos == rhs.m_pos;
            }

            ptrdiff_t operator-(iterator const& rhs) const
            {
                return ptrdiff_t(m_pos) - ptrdiff_t(rhs.m_pos);
            }

            bool operator<(iterator const& rhs) const
            {
                assert(m_l == rhs.m_l);
                return m_pos < rhs.m_pos;
            }

            bool operator>(iterator const& rhs) const
            {
                assert(m_l == rhs.m_l);
                return m_pos > rhs.m_pos;
            }

            bool operator>=(iterator const& rhs) const
            {
                return !(*this < rhs);
            }

            bool operator<=(iterator const& rhs) const
            {
                return !(*this > rhs);
            }

            bool operator!=(iterator const& rhs) const
            {
                return !(*this == rhs);
            }

            value_reference operator*() const
            {
                return value_reference(m_l, m_pos);
            }

            iterator& operator+=(ptrdiff_t n)
            {
                m_pos += n;
                return *this;
            }

            iterator operator+(ptrdiff_t n) const
            {
                iterator copy(*this);
                return (copy += n);
            }

            iterator& operator-=(ptrdiff_t n)
            {
                return (*this += -n);
            }

            iterator operator-(ptrdiff_t n) const
            {
                iterator copy(*this);
                return (copy -= n);
            }

            iterator& operator++()
            {
                return (*this += 1);
            }

            iterator operator++(int)
            {
                iterator copy(*this);
                ++(*this);
                return copy;
            }

            iterator& operator--()
            {
                return (*this -= 1);
            }

            iterator operator--(int)
            {
                iterator copy(*this);
                --(*this);
                return copy;
            }

        private:
            friend class packed_edge_list;

            iterator(packed_edge_list* l, uint64_t pos)
                : m_l(l)
                , m_pos(pos)
            {}

            packed_edge_list* m_l;
            uint64_t m_pos;
        };

        iterator begin()
        {
            return iterator(this, 0);
        }

        iterator end()
        {
            return iterator(this, this->size());
        }

        void swap(packed_edge_list& other)
        {
            std::swap(m_bits, other.m_bits);
            std::swap(m_n, other.m_n);
            m_v.swap(other.m_v);
        }

    protected:

        // the idea of proxy object is similar to that of vector<bool>
        struct value_reference {

            operator value_type() const
            {
                return m_l->get(m_pos);
            }

            value_reference& operator=(value_type v)
            {
                assert(m_l);
                m_l->set(m_pos, v);
                return *this;
            }

            value_reference& operator=(value_reference const& rhs)
            {
                return *this = value_type(rhs);
            }

            bool operator==(value_reference const& rhs) const
            {
                return value_type(*this) == value_type(rhs);
            }

            bool operator<(value_reference const& rhs) const
            {
                return value_type(*this) < value_type(rhs);
            }

            friend inline
            void swap(value_reference lhs, value_reference rhs)
            {
                value_type tmp = lhs;
                lhs = value_type(rhs);
                rhs = value_type(tmp);
            }

            friend inline
            void swap(value_type& lhs, value_reference rhs)
            {
                value_type tmp = lhs;
                lhs = value_type(rhs);
                rhs = value_type(tmp);
            }

            friend inline
            void swap(value_reference lhs, value_type& rhs)
            {
                value_type tmp = lhs;
                lhs = value_type(rhs);
                rhs = value_type(tmp);
            }

        private:
            value_reference()
                : m_l(nullptr)
                , m_pos(-1)
            {}

            value_reference(packed_edge_list* l, uint64_t pos = 0)
                : m_l(l)
                , m_pos(pos)
            {}

            value_reference& operator&() = delete;
            value_reference& operator&() const = delete;

            friend struct iterator;

            packed_edge_list* m_l;
            uint64_t m_pos;
        };

        uint64_t m_bits;
        size_t m_n;
        packed_vector<MemoryModel> m_v;
    };

}
