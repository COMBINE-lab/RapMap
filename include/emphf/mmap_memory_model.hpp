#pragma once

#include <map>
#include <vector>
#include <algorithm>
#include <memory>
#include <stddef.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <unistd.h>

#include "internal_memory_model.hpp"

namespace emphf {

    typedef std::map<void*, std::pair<int, size_t>> mappings_map;

    template <typename T>
    struct mmap_allocator : std::allocator<T>
    {
        typedef T * pointer;
        typedef const T * const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T value_type;
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;

        template<typename U>
        struct rebind {
            typedef mmap_allocator<U> other;
        };

        mmap_allocator()
            : m_mappings(nullptr)
        {}

        mmap_allocator(mappings_map& mappings)
            : std::allocator<T>()
            , m_mappings(&mappings)
        {}

        pointer allocate(size_type n, const void* /* hint */)
        {
            return allocate(n);
        }

        pointer allocate(size_type n)
        {
            assert(m_mappings);
            if (!n) return nullptr;

            size_type size = n * sizeof(T);

            // create temporary file
            char tmpl[] = "mphf.temp.XXXXXX";
            int fd = mkstemp(tmpl);
            if (fd == -1) throw std::runtime_error("Impossible to create temp file");

            int ret;
            ret = unlink(tmpl);
            if (ret) {
                std::cerr << "WARNING: Error unlinking temporary file " << tmpl << std::endl;
            }

            ret = ftruncate(fd, (off_t)size);
            if (ret) throw std::runtime_error("Impossible to resize temp file");

            void* addr = mmap(nullptr, size,
                              PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
            if (!addr) throw std::runtime_error("Impossible to create mapping");

            ret = posix_madvise(addr, size, POSIX_MADV_SEQUENTIAL);
            if (ret) logger() << "Error calling madvice: " << errno << std::endl;

            (*m_mappings)[addr] = std::make_pair(fd, size);

            return static_cast<pointer>(addr);
        }

        void deallocate(pointer p, size_type /* n */)
        {
            assert(m_mappings);
            if (!p) return;
            auto mapping = m_mappings->find(static_cast<void*>(p));
            if (mapping == m_mappings->end()) {
                throw std::runtime_error("Trying to deallocate non-existent mapping");
            }
            int ret = munmap(mapping->first, mapping->second.second);
            if (ret) throw std::runtime_error("Error unmapping file");

            ret = close(mapping->second.first);
            if (ret) throw std::runtime_error("Error closing file");

            m_mappings->erase(mapping);
        }

    private:

        mappings_map* m_mappings;
    };

    struct mmap_memory_model {

        mmap_memory_model()
        {}

        ~mmap_memory_model()
        {
            if (!m_mappings.empty()) {
                std::cerr << "ERROR: leaked mappings in mmap_allocator";
            }
        }

        template <typename T>
        using vector = std::vector<T, mmap_allocator<T>>;

        template <typename T>
        vector<T> make_vector(T const&)
        {
            return vector<T>(mmap_allocator<T>(m_mappings));
        }


        template <typename Iterator,
                  typename KeyFunctor,
                  typename PartitionFunctor>
        struct sorter {

            typedef typename std::iterator_traits<Iterator>::value_type value_type;

            sorter(Iterator begin, Iterator end,
                   KeyFunctor const& kf,
                   PartitionFunctor const& pf)
                : m_cur_part(0)
                , m_buf_pos(0)
                , m_partition_sorted(false)
                , m_begin(begin)
                , m_kf(kf)
            {
                using std::swap;

                const size_t memory = 512 * 1024 * 1024;
                auto len = (size_t)std::distance(begin, end);

                if (len * sizeof(value_type) <= memory) {
                    // create a single partition
                    m_partitions.push_back(0);
                    m_partitions.push_back(len);
                    m_sort_buf.reserve(len);
                    return;
                }

                // use ceil(size / memory) * 2 buckets so that maximum
                // bucket is smaller than memory with high probability
                size_t k = (len * sizeof(value_type) + memory - 1) / memory * 2;
                m_partitions.resize(k + 1);

                // distribution count
                std::for_each(begin, end, [&](value_type const& v) {
                        m_partitions[pf(v, k) + 1] += 1;
                    });

                size_t largest_part = *std::max_element(m_partitions.begin(),
                                                        m_partitions.end() - 1);
                m_sort_buf.reserve(largest_part);

                // cumulative sum
                std::partial_sum(m_partitions.begin(), m_partitions.end(),
                                 m_partitions.begin());
                assert((size_t)m_partitions[k] == len);

                // partitioning
                size_t buffer_bytes = 1024 * 1024;
                std::vector<buffered_cursor<Iterator>> cursors;
                cursors.reserve(k);

                for (size_t p = 0; p < k; ++p) {
                    cursors.emplace_back(begin + m_partitions[p],
                                         begin + m_partitions[p + 1],
                                         buffer_bytes);
                }

                for (size_t p = k - 1; p + 1 > 0; --p) {
                    auto& cur = cursors[p];
                    while (!cur.empty()) {
                        size_t other_p;
                        while ((other_p = pf(cur.value(), k)) != p) {
                            swap(cur.value(), cursors[other_p].value());
                            cursors[other_p].advance();
                        }
                        cur.advance();
                    }
                }
            }

            struct iterator : std::iterator<std::forward_iterator_tag,
                                            value_type> {

                iterator(sorter* s = nullptr)
                    : m_s(s)
                {
                    if (m_s) {
                        advance_part();
                    }
                }

                value_type const& operator*() const
                {
                    assert(m_s);
                    assert(m_s->m_cur_part < m_s->m_partitions.size() - 1);

                    if (!m_s->m_partition_sorted) {
                        ptrdiff_t part_begin = m_s->m_partitions[m_s->m_cur_part];
                        ptrdiff_t part_end = m_s->m_partitions[m_s->m_cur_part + 1];
                        m_s->m_sort_buf.clear();
                        std::copy(m_s->m_begin + part_begin, m_s->m_begin + part_end,
                                  std::back_inserter(m_s->m_sort_buf));
                        std::sort(m_s->m_sort_buf.begin(), m_s->m_sort_buf.end(),
                                  [&](value_type const& lhs, value_type const& rhs) {
                                      return m_s->m_kf(lhs) < m_s->m_kf(rhs);
                                  });
                        m_s->m_partition_sorted = true;
                    }

                    assert(m_s->m_buf_pos < m_s->m_sort_buf.size());
                    return m_s->m_sort_buf[m_s->m_buf_pos];
                }

                iterator& operator++()
                {
                    assert(m_s);
                    ++m_s->m_buf_pos;
                    advance_part();
                    return *this;
                }

                bool operator==(iterator const& rhs) const
                {
                    return m_s == rhs.m_s;
                }

                bool operator!=(iterator const& rhs) const
                {
                    return !(*this == rhs);
                }

            private:
                void advance_part()
                {
                    while (m_s->m_partitions[m_s->m_cur_part] + m_s->m_buf_pos ==
                           m_s->m_partitions[m_s->m_cur_part + 1]) {
                        ++m_s->m_cur_part;
                        m_s->m_buf_pos = 0;
                        m_s->m_partition_sorted = false;

                        if (m_s->m_cur_part == m_s->m_partitions.size() - 1) {
                            m_s = nullptr;
                            return;
                        }
                    }
                }

                sorter* m_s;
            };

            iterator begin()
            {
                return iterator(this);
            }

            iterator end()
            {
                return iterator();
            }


        private:
            std::vector<size_t> m_partitions;
            std::vector<value_type> m_sort_buf;
            size_t m_cur_part;
            size_t m_buf_pos;
            bool m_partition_sorted;

            Iterator m_begin;
            KeyFunctor m_kf;
        };

        template <typename Iterator,
                  typename KeyFunctor,
                  typename PartitionFunctor>
        sorter<Iterator, KeyFunctor, PartitionFunctor>
        make_sorter(Iterator begin, Iterator end,
                    KeyFunctor const& kf,
                    PartitionFunctor const& pf) const
        {
            return sorter<Iterator, KeyFunctor, PartitionFunctor>(begin, end, kf, pf);
        }

        template <typename Iterator,
                  typename KeyFunctor,
                  typename PartitionFunctor>
        void sort(Iterator begin, Iterator end,
                  KeyFunctor const& kf,
                  PartitionFunctor const& pf) const
        {
            auto s = make_sorter(begin, end, kf, pf);
            std::copy(s.begin(), s.end(), begin);
        }

    private:

        template <typename Iterator>
        struct buffered_cursor {
            typedef typename Iterator::value_type value_type;

            buffered_cursor(Iterator begin, Iterator end, size_t buf_bytes)
                : m_base_cur(begin)
                , m_base_end(end)
                , m_max_bufsize((buf_bytes + sizeof(value_type) - 1) / sizeof(value_type))
                , m_buf_pos(0)
            {
                m_buffer.reserve(m_max_bufsize);
                flush();
            }

            ~buffered_cursor()
            {
                assert(empty());
                assert(m_buffer.empty());
            }

            bool empty() const
            {
                return m_base_cur == m_base_end;
            }

            void advance()
            {
                m_buf_pos += 1;
                if (m_buf_pos == m_buffer.size()) {
                    flush();
                }
            }

            value_type& value()
            {
                return m_buffer[m_buf_pos];
            }


        private:
            typedef typename std::iterator_traits<Iterator>::difference_type diff_t;

            void flush()
            {
                // flush current buffer
                assert(m_buf_pos == m_buffer.size());
                std::copy(m_buffer.begin(), m_buffer.end(), m_base_cur);
                m_base_cur += m_buffer.size();
                m_buffer.clear();
                m_buf_pos = 0;

                // fill buffer
                auto avail = std::min(m_max_bufsize,
                                      std::distance(m_base_cur, m_base_end));
                assert(avail >= 0);
                m_buffer.assign(m_base_cur, m_base_cur + avail);

            }

            Iterator m_base_cur, m_base_end;
            std::vector<value_type> m_buffer;
            diff_t m_max_bufsize;
            size_t m_buf_pos;
        };

        mappings_map m_mappings;
    };

}
