#pragma once

#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <memory>
#include <ctime>
#include <cstring>
#include <cassert>

#include "emphf_config.hpp"

namespace emphf {

    std::ostream& logger();

    // XXX(ot): the following I/O code is adapted from succinct
    // library, avoiding the dependency for now
    typedef std::pair<uint8_t const*, uint8_t const*> byte_range_t;

    struct identity_adaptor
    {
        byte_range_t operator()(byte_range_t s) const
        {
            return s;
        }
    };

    struct stl_string_adaptor
    {
        byte_range_t operator()(std::string const& s) const
        {
            const uint8_t* buf = reinterpret_cast<uint8_t const*>(s.c_str());
            const uint8_t* end = buf + s.size() + 1; // add the null terminator
            return byte_range_t(buf, end);
        }
    };

    class line_iterator
        : public std::iterator<std::forward_iterator_tag, const std::string> {

    public:
        line_iterator()
            : m_is(nullptr)
            , m_buf(nullptr)
        {}

        line_iterator(FILE* is)
            : m_is(is)
            , m_pos(0)
            , m_buf(nullptr)
            , m_buf_len(0)
        {
            advance();
        }

        ~line_iterator()
        {
            free(m_buf);
        }

        value_type const& operator*() const {
            return m_line;
        }

        line_iterator& operator++() {
            advance();
            return *this;
        }

        friend bool operator==(line_iterator const& lhs, line_iterator const& rhs)
        {
            if (!lhs.m_is || !rhs.m_is) {
                if (!lhs.m_is && !rhs.m_is) {
                    return true;
                } else {
                    return false;
                }
            }

            assert(lhs.m_is == rhs.m_is);

            return rhs.m_pos == lhs.m_pos;
        }

        friend bool operator!=(line_iterator const& lhs, line_iterator const& rhs)
        {
            return !(lhs == rhs);
        }

    private:
        void advance()
        {
            assert(m_is);
            fseek(m_is, m_pos, SEEK_SET);

            // this is significantly faster than std::getline on C++
            // streams
            auto avail = getline(&m_buf, &m_buf_len, m_is);
            if (avail == -1) {
                m_is = nullptr;
                return;
            }
            m_pos = ftell(m_is);

            // trim newline character
            if (avail && m_buf[avail - 1] == '\n') {
                avail -= 1;
            }

            m_line.assign(m_buf, m_buf + avail);
        }

        FILE* m_is;
        long m_pos;
        std::string m_line;
        char* m_buf;
        size_t m_buf_len;
    };

    class file_lines
    {
    public:
        file_lines(const char* filename)
        {
            m_is = fopen(filename, "rb");
            if (!m_is) {
                throw std::invalid_argument("Error opening " + std::string(filename));
            }
        }

        ~file_lines()
        {
            fclose(m_is);
        }

        line_iterator begin() const
        {
            return line_iterator(m_is);
        }

        line_iterator end() const { return line_iterator(); }

        size_t size() const
        {
            size_t lines = 0;
            fseek(m_is, 0, SEEK_SET);
            static const size_t buf_size = 4096;
            char buf[buf_size];
            size_t avail;
            bool last_is_newline = false;
            while ((avail = fread(buf, 1, buf_size, m_is))) {
                for (size_t i = 0; i < avail; ++i) {
                    if (buf[i] == '\n') lines += 1;
                }
                last_is_newline = (buf[avail - 1] == '\n');
            }

            if (!last_is_newline) lines += 1;

            return lines;
        }

    private:
        // noncopyble
        file_lines(file_lines const&);
        file_lines& operator=(file_lines const&);

        FILE* m_is;
    };

    template <typename Iterator>
    struct iter_range
    {
        iter_range(Iterator b, Iterator e)
            : m_begin(b)
            , m_end(e)
        {}

        Iterator begin() const
        { return m_begin; }

        Iterator end() const
        { return m_end; }

        Iterator m_begin, m_end;
    };

    template <typename Iterator>
    iter_range<Iterator> range(Iterator begin, Iterator end)
    {
        return iter_range<Iterator>(begin, end);
    }

    uint64_t nonzero_pairs(uint64_t x);

    inline uint64_t msb(uint64_t x)
    {
        assert(x);
        return 63 - __builtin_clzll(x);
    }

    struct uninitialized_uint64 {
        uninitialized_uint64() {}

        uninitialized_uint64& operator=(uint64_t v)
        {
            m_val = v;
            return *this;
        }

        operator uint64_t&()
        {
            return m_val;
        }

        operator uint64_t const&() const
        {
            return m_val;
        }

    private:
        uint64_t m_val;
    };

}
