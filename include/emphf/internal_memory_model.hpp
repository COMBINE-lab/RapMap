#pragma once

#include <vector>
#include <algorithm>

namespace emphf {

    struct internal_memory_model {

        template <typename T>
        using vector = std::vector<T>;

        template <typename T>
        vector<T> make_vector(T const&) const
        {
            return vector<T>();
        }

        // for internal memory model, this is just a wrapper for an
        // already sorted range
        template <typename Iterator,
                  typename KeyFunctor,
                  typename PartitionFunctor>
        struct sorter {

            typedef typename std::iterator_traits<Iterator>::value_type value_type;
            typedef Iterator iterator;

            sorter(Iterator begin, Iterator end,
                   KeyFunctor const&,
                   PartitionFunctor const&)
                : m_begin(begin)
                , m_end(end)
            {}

            iterator begin() const
            {
                return m_begin;
            }

            iterator end() const
            {
                return m_end;
            }


        private:
            iterator m_begin;
            iterator m_end;
        };

        template <typename Iterator,
                  typename KeyFunctor,
                  typename PartitionFunctor>
        sorter<Iterator, KeyFunctor, PartitionFunctor>
        make_sorter(Iterator begin, Iterator end,
                    KeyFunctor const& kf,
                    PartitionFunctor const& pf) const
        {
            typedef typename std::iterator_traits<Iterator>::value_type value_type;
            std::sort(begin, end,
                      [&](value_type const& lhs, value_type const& rhs) {
                          return kf(lhs) < kf(rhs);
                      });
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
        }
    };

}
