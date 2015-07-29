#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cstring>

#include "RSDic.hpp"
#include "ScopedTimer.hpp"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

//#include "HitManager.hpp"
#include "SIMDCompressionAndIntersection/intersection.h"
#include "xxhash.h"

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/hash_counter.hpp"

#include "tclap/CmdLine.h"

/*extern "C" {
#include "kseq.h"
}
*/
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#include "stringpiece.h"

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "RapMapFileSystem.hpp"
#include "RapMapConfig.hpp"
#include "ScopedTimer.hpp"
#include "SpinLock.hpp"

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using TranscriptList = std::vector<uint32_t>;
using PositionList = std::vector<uint32_t>;
using KmerIndex = std::unordered_map<uint64_t, TranscriptList, rapmap::utils::KmerKeyHasher>;
using IntervalIndex = std::unordered_map<uint64_t, rapmap::utils::KmerInterval, rapmap::utils::KmerKeyHasher>;
using OccList = std::vector<uint64_t>;
using KmerInfoList = std::vector<rapmap::utils::KmerInfo>;
using EqClassList = std::vector<rapmap::utils::EqClass>;
using EqClassLabelVec = std::vector<uint32_t>;
using PositionListHelper = rapmap::utils::PositionListHelper;
#if defined __APPLE__
using SpinLockT = SpinLock;
#else
using SpinLockT = std::mutex;
#endif

using HitCounters = rapmap::utils::HitCounters;
using MateStatus = rapmap::utils::MateStatus;
using HitInfo = rapmap::utils::HitInfo;
using ProcessedHit = rapmap::utils::ProcessedHit;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
using FixedWriter = rapmap::utils::FixedWriter;

class SASearcher {
    public:
        SASearcher(RapMapSAIndex* rmi) :
            rmi_(rmi), seq_(&rmi->seq), sa_(&rmi->SA) {}

        int cmp(std::string::iterator abeg,
                std::string::iterator aend,
                std::string::iterator bbeg,
                std::string::iterator bend) {
            auto ait = abeg;
            auto bit = bbeg;
            //size_t la = a.length();
            //size_t lb = b.length();
            while (ait < aend and bit < bend) {
                if (*ait < *bit) {
                    return -1;
                } else if (*ait > *bit) {
                    return 1;
                }
                ++ait;
                ++bit;
            }
            if (bit == bend and ait < aend) {
                return 1;
            }
            return 0;
        }

        enum class SearchDirection : uint8_t {
            UP = 0, DOWN
        };
        struct BoundSearchResult {
            int maxLen;
            int bound;
            SearchDirection dir;
        };
        // http://www.cs.jhu.edu/~langmea/resources/lecture_notes/suffix_arrays.pdf
        // templated on the iterator type so we can use a forward or revers iterator
        template <typename IteratorT>
        std::tuple<int, int, int> extendSearch(
                int lbIn, // The lower bound for the search
                int ubIn, // The upper bound for the search
                int startAt, // The offset at which to start looking
                IteratorT qb, // Iterator to the beginning of the query
                IteratorT qe, // Iterator to the end of the query
                bool complementBases=false // True if bases should be complemented
                                           // before comparison
                ) {

            std::vector<int>& SA = *sa_;
            std::string& seq = *seq_;

            int m = std::distance(qb, qe);
            size_t n = seq.length();

            auto sb = seq.begin();
            auto se = seq.end();

            // If the bounds are already trivial, just figure how long
            // of a prefix we share and return the interval.
            if (ubIn - lbIn == 2) {
                lbIn += 1;
                auto i = startAt;
                while (i < m and SA[lbIn] + i < n) {
                    char queryChar = ::toupper(*(qb + i));
                    // If we're reverse complementing
                    if (complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }
                    if ( queryChar < *(sb + SA[lbIn] + i) ) {
                        break;
                    } else if ( queryChar > *(sb + SA[lbIn] + i)) {
                        break;
                    }
                    ++i;
                }
                return std::make_tuple(lbIn, ubIn, i);
            }

            BoundSearchResult res1, res2;

            char smallest = '#';
            char largest = '}';
            char sentinel = smallest;

            int l = lbIn, r = ubIn;
            int lcpLP = startAt, lcpRP = startAt;
            int c{0};
            int i{0};
            int maxI{startAt};
            int prevI = startAt;
            int prevILow = startAt;
            int prevIHigh = startAt;
            int validBoundLow = ubIn;
            int validBoundHigh = lbIn;
            int validBound = 0;
            bool plt{true};
            bool prevPLT{true};
            //std::cerr << "lbIn = " << lbIn << ", ubIn = " << ubIn << "\n";
            // Reduce the search interval until we hit a border
            // i.e. until c == r - 1 or c == l + 1
            while (true) {
                c = (l + r) / 2;
                //std::cerr << "l = " << l << ", r = " << r << ", c = " << c << '\n';
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    char queryChar = ::toupper(*(qb + i));
                    // If we're reverse complementing
                    if (complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }

                    if ( queryChar < *(sb + SA[c] + i) ) {
                        if (i > prevIHigh) {
                            prevIHigh = i;
                            validBoundHigh = c;
                        } else if (i == prevIHigh) {
                            validBoundHigh = c < validBoundHigh ? c : validBoundHigh;
                        }
                        //std::cerr << "(l = " << l << ", r = " << r << ") pattern < SA[" << c << "]\n";
                        //std::cerr << "(i = " << i << ", m = " << m << ") " << queryChar << " < " <<  *(sb + SA[c] + i) << "\n";

                        break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        if (i > prevILow) {
                            prevILow = i;
                            validBoundLow = c;
                        } else if (i == prevILow) {
                            validBoundLow = c > validBoundLow ? c : validBoundLow;
                        }
                        //std::cerr << "(l = " << l << ", r = " << r << ") pattern > SA[" << c << "]\n";
                        //std::cerr << "(i = " << i << ", m = " << m << ") " << queryChar << " > " <<  *(sb + SA[c] + i) << "\n";
                        plt = false;
                        break;
                    }

                    ++i;
                }
                if (i > prevILow) {
                    prevILow = i;
                    validBoundLow = c;
                } else if (i == prevILow) {
                    validBoundLow = c < validBoundLow ? c : validBoundLow;
                }

                if (plt) {
                    if (c == l + 1) {
                        //std::cerr << "path 1\n";
                        auto maxI = std::max(std::max(i, prevILow), prevIHigh);
                        res1.maxLen = maxI;
                        if (maxI == m) {
                            res1.dir = SearchDirection::DOWN;
                            res1.bound = c;
                        } else {
                            validBound = (prevILow >= prevIHigh) ? validBoundLow : validBoundHigh;
                            res1.bound = validBound;
                            res1.dir = (res1.bound == validBoundLow) ? SearchDirection::DOWN : SearchDirection::UP;
                        }
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        //std::cerr << "path 2\n";
                        maxI = std::max(std::max(i, prevILow), prevIHigh);
                        res1.maxLen = maxI;
                        validBound = (prevILow >= prevIHigh) ? validBoundLow : validBoundHigh;
                        if (maxI == m) {
                            res1.bound = r;
                        } else {
                            res1.bound = validBound;
                        }
                        res1.dir = (res1.bound == validBoundLow) ? SearchDirection::DOWN : SearchDirection::UP;
                        break;
                    }
                    l = c;
                    lcpLP = i;
                }
            }

            /*
            if (ubIn - lower < 2) {
                return std::make_tuple(lower, ubIn);
            }
            */

            bool knownValid{true};
            m = res1.maxLen + 1;

            switch (res1.dir) {
                case SearchDirection::UP:
                    sentinel = '#';
                    r = res1.bound;
                    l = lbIn;
                    //std::cerr << "direction was UP; lb = " << l << ", ub = " << r << "\n";
                    break;
                case SearchDirection::DOWN:
                    sentinel = '{';
                    r = ubIn;
                    l = res1.bound;
                    //std::cerr << "direction was DOWN; lb = " << l << ", ub = " << r << "\n";
                    break;
            }

            if (r - l < 2) {
                if (r == l) { r += 1; }
                //std::cerr << "early exit!\n";
                return std::make_tuple(l, r, res1.maxLen);
            }


            lcpLP = startAt;
            lcpRP = startAt;
            c = 0;
            plt = true;
            prevPLT = true;
            prevI = 0;
            prevILow = 0;
            prevIHigh = 0;
            i = startAt;
            validBound = 0;
            validBoundLow = ubIn;
            validBoundHigh = lbIn;
            while (true) {
                c = (l + r) / 2;
                //prevI = std::max(lcpLP, lcpRP);
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    char queryChar = (i < m - 1) ? ::toupper(*(qb + i)) : sentinel;
                    // If we're reverse complementing
                    if (queryChar != sentinel and complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }

                    if ( queryChar < *(sb + SA[c] + i) ) {
                        if (i > prevILow) {
                            prevILow = i;
                            validBoundLow = c;
                        } else if (i == prevILow) {
                            validBoundLow = c < validBoundLow ? c : validBoundLow;
                        }
                        break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        if (i > prevIHigh) {
                            prevIHigh = i;
                            validBoundHigh = c;
                        } else if (i == prevIHigh) {
                            validBoundHigh = c > validBoundHigh ? c : validBoundHigh;
                        }
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) {
                        auto maxI = std::max(std::max(i, prevILow), prevIHigh);
                        res2.maxLen = maxI;
                        if (maxI == m) {
                            res2.dir = SearchDirection::DOWN;
                            res2.bound = c;
                        } else {
                            validBound = (prevILow >= prevIHigh) ? validBoundLow : validBoundHigh;
                            res2.bound = validBound;
                            res2.dir = (res2.bound == validBoundLow) ? SearchDirection::DOWN : SearchDirection::UP;
                        }
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        //std::cerr << "path 2\n";
                        maxI = std::max(i, prevI);
                        res2.maxLen = maxI;
                        validBound = (prevILow >= prevIHigh) ? validBoundLow : validBoundHigh;
                        if (maxI == m) {
                            res2.bound = r;
                        } else {
                            res2.bound = validBound;
                        }
                        res2.dir = (res2.bound == validBoundLow) ? SearchDirection::DOWN : SearchDirection::UP;
                        break;
                    }
                    l = c;
                    lcpLP = i;
                }
            }

            auto bound1 = std::min(res1.bound, res2.bound);
            auto bound2 = std::max(res1.bound, res2.bound);
            // Must occur at least once!
            if (bound1 == bound2) { bound2 += 1; }
            return std::make_tuple(bound1, bound2, res1.maxLen);
        }



        // http://www.cs.jhu.edu/~langmea/resources/lecture_notes/suffix_arrays.pdf
        std::tuple<int, int> querySimpleAccel(std::string::iterator qb,
                                              std::string::iterator qe) {
            std::vector<int>& SA = *sa_;
            std::string& seq = *seq_;
            //ForwardIt it;
            auto sb = seq.begin();
            auto se = seq.end();

            size_t n = seq.length();
            size_t m = std::distance(qb, qe);
            size_t l = 0, r = n;
            size_t lcpLP = 0, lcpRP = 0;
            size_t c{0};
            size_t i{0};
            bool plt{true};
            size_t lower{0};
            while (true) {
                c = (l + r) / 2;
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    if ( *(qb + i) < *(sb + SA[c] + i) ) {
                        break;
                    } else if ( *(qb + i) > *(sb + SA[c] + i)) {
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) { lower = c; break; }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) { lower = r; break; }
                    l = c;
                    lcpLP = i;
                }
            }

            i = 0;
            l = 0;
            r = n;
            lcpLP = 0;
            lcpRP = 0;
            size_t upper{0};
            while (true) {
                c = (l + r) / 2;
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    if ( *(qb + i) < *(sb + SA[c] + i) ) {
                        break;
                    } else if ( *(qb + i) > *(sb + SA[c] + i)) {
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) { upper = c; break; }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) { upper = r; break; }
                    l = c;
                    lcpLP = i;
                }
            }
            return std::make_tuple(lower, upper);
        }

        uint32_t lce(int p1, int p2,
                     int startAt=0,
                     int stopAt=std::numeric_limits<int>::max(),
                     bool verbose=false) {
            std::string& seq = *seq_;
            std::vector<int>& SA = *sa_;
            size_t len{0};
            auto o1 = SA[p1] + startAt;
            auto o2 = SA[p2] + startAt;
            auto maxIndex = std::max(o1, o2);
            //std::stringstream ext;
            while (maxIndex + len < textLen_ and seq[o1+len] == seq[o2+len]) {
                if (seq[o1+len] == '$') { break; }
                if (len >= stopAt) { break; }
               // ext << seq[o1+len];
                ++len;
            }
            /*
            if (verbose) {
                std::cerr << "lce is " << ext.str() << "\n";
            }
            */
            return len;
        }
    private:
        RapMapSAIndex* rmi_;
        std::string* seq_;
        std::vector<int>* sa_;
        uint32_t textLen_;
};



class SACollector {
    public:
    SACollector(RapMapSAIndex* rmi) : rmi_(rmi) {}
    bool operator()(std::string& read,
                    std::vector<QuasiAlignment>& hits,
                    SASearcher& saSearcher,
                    MateStatus ms,
                    std::atomic<uint64_t>& diffCount) {

        auto& posIDs = rmi_->positionIDs;
        auto& SA = rmi_->SA;
        auto& khash = rmi_->khash;
        auto salen = SA.size();
        auto k = rapmap::utils::my_mer::k();
        auto readStartIt = read.begin();
        auto readEndIt = read.end();

        auto readRevStartIt = read.rbegin();
        auto readRevEndIt = read.rend();

        auto rb = read.begin();
        auto re = rb + 31;
        int lbLeftFwd = 0, ubLeftFwd = 0;
        int lbLeftRC = 0, ubLeftRC = 0;
        int lbRightFwd = 0, ubRightFwd = 0;
        int lbRightRC = 0, ubRightRC = 0;
        int matchedLen;
        bool leftFwdHit = false;
        bool leftRCHit = false;
        bool leftHit = false;
        bool rightFwdHit = false;
        bool rightRCHit = false;
        bool rightHit = false;
        bool isRev = false;
        rapmap::utils::my_mer mer;
        rapmap::utils::my_mer rcMer;
        std::vector<uint32_t> leftTxps, leftTxpsRC;
        std::vector<uint32_t> rightTxps, rightTxpsRC;
        int maxInterval{1000};

        // Find a hit within the read
        // While we haven't fallen off the end
        while (re < read.end()) {
            // Get the k-mer at the current start position.

            auto pos = std::distance(readStartIt, rb);

            auto invalidPos = read.find_first_of("nN", pos);
            if (invalidPos <= pos + k) {
                rb = read.begin() + invalidPos + 1;
                re = rb + k;
                continue;
            }

            mer = rapmap::utils::my_mer(read.c_str() + pos);
            rcMer = mer.get_reverse_complement();
            // See if we can find this k-mer in the hash
            auto merIt = khash.find(mer.get_bits(0, 2*k));
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

            if (merIt != khash.end()) {
                int lb = merIt->second.begin;
                int ub = merIt->second.end;

                // lb must be 1 *less* then the current lb
                auto lbRestart = std::max(static_cast<int>(0), lb-1);
                //std::cerr << "before extend search . . . ";
                std::tie(lbLeftFwd, ubLeftFwd, matchedLen) =
                    saSearcher.extendSearch(lbRestart, ub, k, rb, readEndIt);
                //std::cerr << "after extend search\n";

                bool verbose{false};
                if (verbose) {
                    if (lb == lbLeftFwd and ub == ubLeftFwd) {
                        //std::cerr << "couldn't narrow interval at all!\n";
                    } else {
			diffCount += 1;
                        int lbNaive, ubNaive, blah;

                        //std::cerr << "before naive search . . . ";
                        std::tie(lbNaive, ubNaive, blah) =
                            saSearcher.extendSearch(0, SA.size(), 0, rb, readEndIt);
                        //std::cerr << "after naive search\n";

                        if (ubNaive - lbNaive > ubLeftFwd - lbLeftFwd and blah >= matchedLen) {//lbLeftFwd != lbNaive or ubLeftFwd != ubNaive) {
                            diffCount += 1;
                        //}
                        //if (ubNaive - lbNaive > ubLeftFwd - lbLeftFwd and !verbose) {
                            std::cerr << "narrowed interval from "
                                << "[ " << lb << ", " << ub << ") to "
                                << "[ "<< lbLeftFwd << ", " << ubLeftFwd << ")\n";

                            std::cerr << "mer is    : " << mer << "\n";
                            std::cerr << "string is : ";
                            for (auto it = rb; it != rb + matchedLen; ++it) {
                                std::cerr << *it;
                            } std::cerr << "\n";
                            std::cerr << "read is   : " << read << "\n";

                            std::cerr << "narrowed len = " << matchedLen << "\n";
                            std::cerr << "naive len = " << blah << "\n";

                            std::cerr << "entries at narrowed interval:\n";
                            for (auto x = lbLeftFwd - 1; x < ubLeftFwd + 1; ++x) {
                                std::cerr << "\tT[SA[" << x << "] = "
                                    << rmi_->seq.substr(SA[x], matchedLen) << "\n";
                            }

                            std::cerr << "naive interval is [" << lbNaive << ", " << ubNaive << ")\n";

                            std::cerr << "entries at naive interval:\n";
                            for (auto x = lbNaive - 1; x < ubNaive + 1; ++x) {
                                std::cerr << "\tT[SA[" << x << "] = "
                                    << rmi_->seq.substr(SA[x], blah) << "\n";
                            }
                            std::exit(1);
                        }
                    }

                }
                int diff = ubLeftFwd - lbLeftFwd;
                if (ubLeftFwd > lbLeftFwd and diff < maxInterval) {
                    leftFwdHit = true;
                }
            }

            if (rcMerIt != khash.end()) {
                lbLeftRC = rcMerIt->second.begin;
                ubLeftRC = rcMerIt->second.end;
                int diff = ubLeftRC - lbLeftRC;
                if (ubLeftRC > lbLeftRC and diff < maxInterval) {
                    leftRCHit = true;
                }
            }


            if (leftFwdHit or leftRCHit) {
                leftHit = true;
                break;
            }
           ++rb; ++re;
        }

        if (!leftHit) { return false; }

        if (leftFwdHit) {
            // Find the longest common extension between the suffixes
            // at lb and ub.
            auto lce = saSearcher.lce(lbLeftFwd, ubLeftFwd-1, k, read.length() - k);

            if (lce > read.length() - k) {
                size_t pos{0};
                auto invalidPos = read.length() - k - 1;
                do {
                    rb = read.begin() + invalidPos + 1;
                    re = rb + k;
                    pos = std::distance(readStartIt, rb);
                    invalidPos = read.find_first_of("nN", pos);
                } while (invalidPos <= pos + k and re <= readEndIt);

                if (re <= readEndIt) {
                    mer = rapmap::utils::my_mer(read.c_str() + pos);
                    auto merIt = khash.find(mer.get_bits(0, 2*k));

                    if (merIt != khash.end()) {
                        lbRightFwd = merIt->second.begin;
                        ubRightFwd = merIt->second.end;
                        // lb must be 1 *less* then the current lb
                        lbRightFwd = std::max(0, lbRightFwd - 1);
                        std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                            saSearcher.extendSearch(lbRightFwd, ubRightFwd,
                                                    k, rb, readEndIt);
                        int diff = ubRightFwd - lbRightFwd;
                        if (ubRightFwd > lbRightFwd and diff < maxInterval) {
                            rightFwdHit = true;
                        }
                    }

                }
            } else {
                size_t pos{0};
                auto invalidPos = std::max(static_cast<int>(0),
                                           static_cast<int>(lce) - static_cast<int>(k));
                do {
                    rb = read.begin() + invalidPos + 1;
                    re = rb + k;
                    pos = std::distance(readStartIt, rb);
                    invalidPos = read.find_first_of("nN", pos);
                } while (invalidPos <= pos + k and re <= readEndIt);

                if (re <= readEndIt) {
                    mer = rapmap::utils::my_mer(read.c_str() + pos);

                    auto merIt = khash.find(mer.get_bits(0, 2*k));

                    if (merIt != khash.end()) {
                        lbRightFwd = merIt->second.begin;
                        ubRightFwd = merIt->second.end;
                        // lb must be 1 *less* then the current lb
                        lbRightFwd = std::max(0, lbRightFwd - 1);
                        std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                            saSearcher.extendSearch(lbRightFwd, ubRightFwd, k, rb, readEndIt);
                        int diff = ubRightFwd - lbRightFwd;
                        if (ubRightFwd > lbRightFwd and diff < maxInterval) {
                            rightFwdHit = true;
                        }
                    }
                }
            }
        }

        if (leftRCHit) {
            rb = read.end() - k;
            re = rb + k;
            auto pos = std::distance(readStartIt, rb);
            mer = rapmap::utils::my_mer(read.c_str() + pos);
            rcMer = mer.get_reverse_complement();

            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

            if (rcMerIt != khash.end()) {
                lbRightRC = rcMerIt->second.begin;
                ubRightRC = rcMerIt->second.end;
                // lb must be 1 *less* then the current lb
                // We can't move any further in the reverse complement direction
                lbRightRC = std::max(0, lbRightRC - 1);
                auto queryBeg = read.rbegin();
                auto queryEnd = read.rend();

                std::tie(lbRightRC, ubRightRC, matchedLen) =
                    saSearcher.extendSearch(lbRightRC, ubRightRC, k,
                                            queryBeg, queryEnd, true);
                int diff = ubRightRC - lbRightRC;
                if (ubRightRC > lbRightRC and diff < maxInterval) {
                    rightRCHit = true;
                }
            }
        }

        std::vector<uint32_t> jointTxps;
        if (leftFwdHit and rightFwdHit) {
            for (auto i = lbLeftFwd; i < ubLeftFwd; ++i) {
                size_t txpNum = posIDs[SA[i]];
                leftTxps.push_back(txpNum);
            }
            std::sort(leftTxps.begin(), leftTxps.end());
            auto newEnd = std::unique(leftTxps.begin(), leftTxps.end());
            leftTxps.resize(std::distance(leftTxps.begin(), newEnd));

            for (auto i = lbRightFwd; i < ubRightFwd; ++i) {
                size_t txpNum = posIDs[SA[i]];
                rightTxps.push_back(txpNum);
            }
            std::sort(rightTxps.begin(), rightTxps.end());
            newEnd = std::unique(rightTxps.begin(), rightTxps.end());
            rightTxps.resize(std::distance(rightTxps.begin(), newEnd));

            std::set_intersection(leftTxps.begin(), leftTxps.end(),
                    rightTxps.begin(), rightTxps.end(),
                    std::back_inserter(jointTxps));
        } else if (leftFwdHit) {
		for (auto i = lbLeftFwd; i < ubLeftFwd; ++i) {
			size_t txpNum = posIDs[SA[i]];
			leftTxps.push_back(txpNum);
		}
		std::sort(leftTxps.begin(), leftTxps.end());
		auto newEnd = std::unique(leftTxps.begin(), leftTxps.end());
		leftTxps.resize(std::distance(leftTxps.begin(), newEnd));
		jointTxps = std::move(leftTxps);
	} else if (rightFwdHit) {
		for (auto i = lbRightFwd; i < ubRightFwd; ++i) {
			size_t txpNum = posIDs[SA[i]];
			rightTxps.push_back(txpNum);
		}
		std::sort(rightTxps.begin(), rightTxps.end());
		auto newEnd = std::unique(rightTxps.begin(), rightTxps.end());
		rightTxps.resize(std::distance(rightTxps.begin(), newEnd));
		jointTxps = std::move(rightTxps);
	}

        if (leftRCHit and rightRCHit) {
            for (auto i = lbLeftRC; i < ubLeftRC; ++i) {
                size_t txpNum = posIDs[SA[i]];
                leftTxpsRC.push_back(txpNum);
            }
            std::sort(leftTxpsRC.begin(), leftTxpsRC.end());
            auto newEnd = std::unique(leftTxpsRC.begin(), leftTxpsRC.end());
            leftTxpsRC.resize(std::distance(leftTxpsRC.begin(), newEnd));

            for (auto i = lbRightRC; i < ubRightRC; ++i) {
                size_t txpNum = posIDs[SA[i]];
                rightTxpsRC.push_back(txpNum);
            }
            std::sort(rightTxpsRC.begin(), rightTxpsRC.end());
            newEnd = std::unique(rightTxpsRC.begin(), rightTxpsRC.end());
            rightTxpsRC.resize(std::distance(rightTxpsRC.begin(), newEnd));

            std::set_intersection(leftTxpsRC.begin(), leftTxpsRC.end(),
                    rightTxpsRC.begin(), rightTxpsRC.end(),
                    std::back_inserter(jointTxps));
        } else if (leftRCHit) {
		for (auto i = lbLeftRC; i < ubLeftRC; ++i) {
			size_t txpNum = posIDs[SA[i]];
			leftTxpsRC.push_back(txpNum);
		}
		std::sort(leftTxpsRC.begin(), leftTxpsRC.end());
		auto newEnd = std::unique(leftTxpsRC.begin(), leftTxpsRC.end());
		leftTxpsRC.resize(std::distance(leftTxpsRC.begin(), newEnd));
		jointTxps .insert(jointTxps.end(), leftTxpsRC.begin(), leftTxpsRC.end());
	} else if (rightRCHit) {
		for (auto i = lbRightRC; i < ubRightRC; ++i) {
			size_t txpNum = posIDs[SA[i]];
			rightTxpsRC.push_back(txpNum);
		}
		std::sort(rightTxpsRC.begin(), rightTxpsRC.end());
		auto newEnd = std::unique(rightTxpsRC.begin(), rightTxpsRC.end());
		rightTxpsRC.resize(std::distance(rightTxpsRC.begin(), newEnd));
		jointTxps.insert(jointTxps.end(), rightTxpsRC.begin(), rightTxpsRC.end());
	}

        return jointTxps.size() > 0;
    }

    private:
        RapMapSAIndex* rmi_;
};


template <typename CollectorT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
        RapMapSAIndex& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
        std::ostream& outStream,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {
}

template <typename CollectorT, typename MutexT>
void processReadsPairSA(paired_parser* parser,
        RapMapSAIndex& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
        std::ostream& outStream,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {
    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpOffsets = rmi.txpOffsets;
    uint32_t n{0};
    //uint32_t k = rapmap::utils::my_mer::k();
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;

    size_t readLen{0};
	bool tooManyHits{false};
    uint16_t flags1, flags2;
    // 1000-bp reads are max here (get rid of hard limit later).
    std::string read1Temp(1000, 'A');
    std::string qual1Temp(1000, '~');
    std::string read2Temp(1000, 'A');
    std::string qual2Temp(1000, '~');

    char buff1[1000];
    char buff2[1000];
    FixedWriter cigarStr1(buff1, 1000);
    FixedWriter cigarStr2(buff2, 1000);

    SASearcher saSearcher(&rmi);

    // 0 means properly aligned
    // 0x1 means only alignments for left read
    // 0x2 means only alignments for right read
    // 0x3 means "orphaned" alignments for left and right
    // (currently not treated as orphan).
    std::atomic<uint64_t> diffCount{0};
    uint32_t orphanStatus{0};
    size_t lb, ub;
    while(true) {
        typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
		    tooManyHits = false;
            readLen = j->data[i].first.seq.length();
            ++hctr.numReads;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            //std::string query(j->data[i].first.seq.substr(0,31));
            //std::tie(lb, ub) = saSearcher.query(query);
            //query = j->data[i].second.seq.substr(0,31);
            //std::tie(lb, ub) = saSearcher.query(query);
    	    bool lh = hitCollector(j->data[i].first.seq,
                        leftHits, saSearcher,
                        MateStatus::PAIRED_END_LEFT, diffCount);
            bool rh = hitCollector(j->data[i].second.seq,
                        rightHits, saSearcher,
                        MateStatus::PAIRED_END_RIGHT, diffCount);
	    bool peHit = (lh and rh);
	    if (peHit) {
	    	hctr.peHits += 1;
	    } else if (lh or rh) {
		std::cerr << "WHAT?!\n";
		hctr.seHits += 1;
	    }
        /*
            hitCollector(j->data[i].second.seq,
                        rightHits, MateStatus::PAIRED_END_RIGHT);
            if (leftHits.size() > 0) {
                auto leftIt = leftHits.begin();
                auto leftEnd = leftHits.end();
                auto leftLen = std::distance(leftIt, leftEnd);
                if (rightHits.size() > 0) {
                    auto rightIt = rightHits.begin();
                    auto rightEnd = rightHits.end();
                    auto rightLen = std::distance(rightIt, rightEnd);
                    size_t numHits{0};
                    jointHits.reserve(std::min(leftLen, rightLen));
					uint32_t leftTxp, rightTxp;
                    while (leftIt != leftEnd && rightIt != rightEnd) {
                        leftTxp = leftIt->tid;
                        rightTxp = rightIt->tid;
                        if (leftTxp < rightTxp) {
                            ++leftIt;
                        } else {
                            if (!(rightTxp < leftTxp)) {
                                int32_t startRead1 = leftIt->pos;
                                int32_t startRead2 = rightIt->pos;
                                int32_t fragStartPos = std::min(leftIt->pos, rightIt->pos);
                                int32_t fragEndPos = std::max(leftIt->pos, rightIt->pos) + readLen;
                                uint32_t fragLen = fragEndPos - fragStartPos;
                                jointHits.emplace_back(leftTxp,
                                                       startRead1,
                                                       leftIt->fwd,
                                                       leftIt->readLen,
                                                       fragLen, true);
                                // Fill in the mate info
                                auto& qaln = jointHits.back();
                                qaln.mateLen = rightIt->readLen;
                                qaln.matePos = startRead2;
                                qaln.mateIsFwd = rightIt->fwd;
                                jointHits.back().mateStatus = MateStatus::PAIRED_END_PAIRED;

                                ++numHits;
                                if (numHits > maxNumHits) { tooManyHits = true; break; }
                                ++leftIt;
                            }
                            ++rightIt;
                        }
                    }
                }
				if (tooManyHits) { jointHits.clear(); ++hctr.tooManyHits; }
            }

			// If we had proper paired hits
            if (jointHits.size() > 0) {
                hctr.peHits += jointHits.size();
                orphanStatus = 0;
            } else if (leftHits.size() + rightHits.size() > 0 and !tooManyHits) {
		    // If there weren't proper paired hits, then either
			// there were too many hits, and we forcibly discarded the read
			// or we take the single end hits.
					auto numHits = leftHits.size() + rightHits.size();
					hctr.seHits += numHits;
					orphanStatus = 0;
					orphanStatus |= (leftHits.size() > 0) ? 0x1 : 0;
					orphanStatus |= (rightHits.size() > 0) ? 0x2 : 0;
					jointHits.insert(jointHits.end(),
									std::make_move_iterator(leftHits.begin()),
									std::make_move_iterator(leftHits.end()));
					jointHits.insert(jointHits.end(),
									std::make_move_iterator(rightHits.begin()),
									std::make_move_iterator(rightHits.end()));
			}

            if (jointHits.size() > 0 and !noOutput) {
                auto& readName = j->data[i].first.header;
                auto& mateName = j->data[i].second.header;
                // trim /1 and /2 from pe read names
                if (readName.length() > 2 and
                        readName[readName.length() - 2] == '/') {
                    readName[readName.length() - 2] = '\0';
                }
                if (mateName.length() > 2 and
                        mateName[mateName.length() - 2] == '/') {
                    mateName[mateName.length() - 2] = '\0';
                }


#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& trueTxpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
				uint32_t trueHitCtr{0};
				QuasiAlignment* firstTrueHit{nullptr};
                for (auto& qa : jointHits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    if (qa.isPaired) {
                        getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            flags2 |= 0x900;
                        }
                        adjustOverhang(qa, txpLens[qa.tid], cigarStr1, cigarStr2);

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            reverseRead(j->data[i].first.seq,
                                        j->data[i].first.qual,
                                        read1Temp,
                                        qual1Temp);
                        }
                        if (!qa.mateIsFwd) {
                            reverseRead(j->data[i].second.seq,
                                        j->data[i].second.qual,
                                        read2Temp,
                                        qual2Temp);
                        }

                        sstream << readName.c_str() << '\t' // QNAME
                                << flags1 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr1.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.matePos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << j->data[i].first.seq << '\t' // SEQ
                                << j->data[i].first.qual << '\n';

                        sstream << mateName.c_str() << '\t' // QNAME
                                << flags2 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.matePos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr2.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.pos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << j->data[i].second.seq << '\t' // SEQ
                                << j->data[i].first.qual << '\n';
                    } else {
                        getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            // If this is the first alignment for this read
                            // If the left end is mapped, set 0x900 on the right end
                            if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) {
                                flags2 |= 0x900;
                            } else {
                            // Otherwise, set 0x900 on the left end
                                flags1 |= 0x900;
                            }
                        }

                        std::string* readSeq{nullptr};
                        std::string* unalignedSeq{nullptr};

                        uint32_t flags, unalignedFlags;
                        std::string* qstr{nullptr};
                        std::string* unalignedQstr{nullptr};
                        std::string* unalignedName{nullptr};
                        FixedWriter* cigarStr;
                        if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) { // left read
                            readName = j->data[i].first.header;
                            unalignedName = &j->data[i].second.header;

                            readSeq = &(j->data[i].first.seq);
                            unalignedSeq = &(j->data[i].second.seq);

                            qstr = &(j->data[i].first.qual);
                            unalignedQstr = &(j->data[i].second.qual);

                            flags = flags1;
                            unalignedFlags = flags2;

                            cigarStr = &cigarStr1;
                        } else { // right read
                            readName = j->data[i].second.header;
                            unalignedName = &(j->data[i].first.header);

                            readSeq = &(j->data[i].second.seq);
                            unalignedSeq = &(j->data[i].first.seq);

                            qstr = &(j->data[i].second.qual);
                            unalignedQstr = &(j->data[i].first.qual);

                            flags = flags2;
                            unalignedFlags = flags1;

                            cigarStr = &cigarStr2;
                        }

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            reverseRead(*readSeq, *qstr,
                                        read1Temp, qual1Temp);
                        }

                        adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], *cigarStr);
                        sstream << readName.c_str() << '\t' // QNAME
                                << flags << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr->c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << 0 << '\t' // PNEXT (only 1 read in templte)
                                << 0 << '\t' // TLEN (spec says 0, not read len)
                                << *readSeq << '\t' // SEQ
                                << *qstr << '\n';

                        // Output the info for the unaligned mate.
                        sstream << unalignedName->c_str() << '\t' // QNAME
                            << unalignedFlags << '\t' // FLAGS
                            << transcriptName << '\t' // RNAME (same as mate)
                            << qa.pos + 1 << '\t' // POS (same as mate)
                            << 0 << '\t' // MAPQ
                            << readSeq->length() << 'M' << '\t' // CIGAR
                            << '=' << '\t' // RNEXT
                            << 0 << '\t' // PNEXT (only 1 read in template)
                            << 0 << '\t' // TLEN (spec says 0, not read len)
                            << *unalignedSeq << '\t' // SEQ
                            << *unalignedQstr << '\n';
                    }
                    ++alnCtr;
                    // === SAM
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (transcriptName == trueTxpName) {
							if (trueHitCtr == 0) {
									++hctr.trueHits;
									++trueHitCtr;
									firstTrueHit = &qa;
							} else {
									++trueHitCtr;
									std::cerr << "Found true hit " << trueHitCtr << " times!\n";
									std::cerr << transcriptName << '\t' << firstTrueHit->pos
											<< '\t' << firstTrueHit->fwd << '\t' << firstTrueHit->fragLen
											<< '\t' << (firstTrueHit->isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(firstTrueHit->mateStatus);
								    std::cerr << '\n';
									std::cerr << transcriptName << '\t' << qa.pos
											  << '\t' << qa.fwd << '\t' << qa.fragLen
										      << '\t' << (qa.isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(qa.mateStatus);
								    std::cerr << '\n';

							}
					}
#endif //__DEBUG__
                }
            }
*/
            if (hctr.numReads > hctr.lastPrint + 100000) {
		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                        std::cerr << "\033[F\033[F\033[F\033[F";
#else
                        std::cerr << "\033[F\033[F\033[F";
#endif // __DEBUG__
                    }
                    std::cerr << "saw " << hctr.numReads << " reads (diffCount = " << diffCount << ")\n";
                    std::cerr << "# pe hits per read = "
                        << hctr.peHits / static_cast<float>(hctr.numReads) << "\n";
                    std::cerr << "# se hits per read = "
                        << hctr.seHits / static_cast<float>(hctr.numReads) << "\n";
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    std::cerr << "The true hit was in the returned set of hits "
                        << 100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads))
                        <<  "% of the time\n";
#endif // __DEBUG__
                    iomutex->unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        if (!noOutput) {
            iomutex->lock();
            outStream << sstream.str();
            iomutex->unlock();
            sstream.clear();
        }

    } // processed all reads

}

/*
int SAQuery(int argc, char* argv[]) {
    std::string indDir(argv[1]);

    std::vector<int> SA;
    std::vector<int> LCP;
    size_t n{0};
    std::ifstream saStream(indDir + "/sa.bin");
    {
        cereal::BinaryInputArchive saArchive(saStream);
        saArchive(SA);
        saArchive(LCP);
    }
    saStream.close();

    std::ifstream seqStream(indDir + "/txpInfo.bin");
    std::vector<std::string> txpNames;
    std::vector<uint32_t> txpOffsets;
    std::string seq;
    {
        cereal::BinaryInputArchive seqArchive(seqStream);
        seqArchive(txpNames);
        seqArchive(txpOffsets);
        seqArchive(seq);
    }
    seqStream.close();

    std::ifstream rsdStream(indDir + "/rsd.bin");
    rsdic::RSDic rsd;
    {
        rsd.Load(rsdStream);
    }
    rsdStream.close();


    SASearcher saSearcher(seq, SA);

    std::cerr << "Ready to search: \n";
    std::string query;
    size_t lb, ub;
    //std::vector inds;
    while (std::getline(std::cin, query)) {
        if (query == "quit") {
            std::cerr << "exiting\n";
            std::exit(0);
        }
        auto ql = query.length();
        //inds.resize(ql);
        //inds = std::iota(inds.begin(), inds.begin()+ql, 0);
        //sacomp.setQuery(query);
        std::tie(lb, ub) = saSearcher.query(query);
        //auto lb = std::lower_bound(SA.begin(), SA.end(), query);
        //auto ub = std::upper_bound(SA.begin(), SA.end(), query);
        std::cerr << "query in = [" << lb << ", " << ub << ")\n";
        if (ub - lb > 1) {
            std::cerr << "LCE = " << saSearcher.lce(lb, ub-1, true) << "\n";
        }
        for (auto i = lb; i < ub; ++i) {
            size_t txpNum = rsd.Rank(SA[i], 1);
            if (txpNum > 0) { txpNum -= 1; }
            std::cerr << "txpNum = " << txpNum << "\n";
            std::cerr << "offset[" << txpNum << "] = " << txpOffsets[txpNum] << "\n";
            std::cerr << "SA[" << i << "] = " << SA[i] << "\n";
            size_t pos = SA[i] - txpOffsets[txpNum];
            std::cerr << "txp = " << txpNames[txpNum]
                      << ", pos = " << pos << "\n";
        }
    }
    return 0;
}
*/

int rapMapSAMap(int argc, char* argv[]) {
    std::cerr << "RapMap Mapper (SA-based)\n";

    std::string versionString = rapmap::version;
    TCLAP::CmdLine cmd(
            "RapMap Mapper",
            ' ',
            versionString);
    cmd.getProgramName() = "rapmap";

    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> unmatedReads("r", "unmatedReads", "The location of single-end reads", false, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
    TCLAP::ValueArg<uint32_t> maxNumHits("m", "maxNumHits", "Reads mapping to more than this many loci are discarded", false, 200, "positive integer");
    TCLAP::ValueArg<std::string> outname("o", "output", "The output file (default: stdout)", false, "", "path");
    TCLAP::SwitchArg endCollectorSwitch("e", "endCollector", "Use the simpler (and faster) \"end\" collector as opposed to the more sophisticated \"skipping\" collector", false);
    TCLAP::SwitchArg noout("n", "noOutput", "Don't write out any alignments (for speed testing purposes)", false);
    cmd.add(index);
    cmd.add(noout);

    cmd.add(read1);
    cmd.add(read2);
    cmd.add(unmatedReads);
    cmd.add(outname);
    cmd.add(numThreads);
    cmd.add(maxNumHits);
    cmd.add(endCollectorSwitch);

    auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
    auto consoleLog = spdlog::create("stderrLog", {consoleSink});

    try {

	cmd.parse(argc, argv);
	bool pairedEnd = (read1.isSet() or read2.isSet());
	if (pairedEnd and (read1.isSet() != read2.isSet())) {
	    consoleLog->error("You must set both the -1 and -2 arguments to align "
		    "paired end reads!");
	    std::exit(1);
	}

	if (pairedEnd and unmatedReads.isSet()) {
	    consoleLog->error("You cannot specify both paired-end and unmated "
		    "reads in the input!");
	    std::exit(1);
	}

	if (!pairedEnd and !unmatedReads.isSet()) {
	    consoleLog->error("You must specify input; either both paired-end "
			      "or unmated reads!");
	    std::exit(1);

	}

	std::string indexPrefix(index.getValue());
	if (indexPrefix.back() != '/') {
	    indexPrefix += "/";
	}

	if (!rapmap::fs::DirExists(indexPrefix.c_str())) {
	    consoleLog->error("It looks like the index you provided [{}] "
		    "doesn't exist", indexPrefix);
	    std::exit(1);
	}

	RapMapSAIndex rmi;
	rmi.load(indexPrefix);

	std::cerr << "\n\n\n\n";

	// from: http://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile
	// set either a file or cout as the output stream
	std::streambuf* outBuf;
	std::ofstream outFile;
	bool haveOutputFile{false};
	if (outname.getValue() == "") {
	    outBuf = std::cout.rdbuf();
	} else {
	    outFile.open(outname.getValue());
	    outBuf = outFile.rdbuf();
	    haveOutputFile = true;
	}
	// Now set the output stream to the buffer, which is
	// either std::cout, or a file.
	std::ostream outStream(outBuf);

	uint32_t nthread = numThreads.getValue();
	std::unique_ptr<paired_parser> pairParserPtr{nullptr};
	std::unique_ptr<single_parser> singleParserPtr{nullptr};

    /*
	if (!noout.getValue()) {
	    writeSAMHeader(rmi, outStream);
	}
    */

	SpinLockT iomutex;
	{
	    ScopedTimer timer;
	    HitCounters hctrs;
	    consoleLog->info("mapping reads . . . \n\n\n");
        if (pairedEnd) {
            std::vector<std::thread> threads;
            std::vector<std::string> read1Vec = rapmap::utils::tokenize(read1.getValue(), ',');
            std::vector<std::string> read2Vec = rapmap::utils::tokenize(read2.getValue(), ',');

            if (read1Vec.size() != read2Vec.size()) {
                consoleLog->error("The number of provided files for "
                        "-1 and -2 must be the same!");
                std::exit(1);
            }

            size_t numFiles = read1Vec.size() + read2Vec.size();
            char** pairFileList = new char*[numFiles];
            for (size_t i = 0; i < read1Vec.size(); ++i) {
                pairFileList[2*i] = const_cast<char*>(read1Vec[i].c_str());
                pairFileList[2*i+1] = const_cast<char*>(read2Vec[i].c_str());
            }
            size_t maxReadGroup{1000}; // Number of reads in each "job"
            size_t concurrentFile{2}; // Number of files to read simultaneously
            pairParserPtr.reset(new paired_parser(4 * nthread, maxReadGroup,
                        concurrentFile,
                        pairFileList, pairFileList+numFiles));

            SACollector saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPairSA<SACollector, SpinLockT>,
                        pairParserPtr.get(),
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
                        std::ref(outStream),
                        std::ref(hctrs),
                        maxNumHits.getValue(),
                        noout.getValue());
            }

            for (auto& t : threads) { t.join(); }
            delete [] pairFileList;
        } else {
            std::vector<std::thread> threads;
            std::vector<std::string> unmatedReadVec = rapmap::utils::tokenize(unmatedReads.getValue(), ',');
            size_t maxReadGroup{1000}; // Number of reads in each "job"
            size_t concurrentFile{1};
            stream_manager streams( unmatedReadVec.begin(), unmatedReadVec.end(),
                    concurrentFile);
            singleParserPtr.reset(new single_parser(4 * nthread,
                        maxReadGroup,
                        concurrentFile,
                        streams));

            /** Create the threads depending on the collector type **/
            SACollector saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingleSA<SACollector, SpinLockT>,
                        singleParserPtr.get(),
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
                        std::ref(outStream),
                        std::ref(hctrs),
                        maxNumHits.getValue(),
                        noout.getValue());
            }
            for (auto& t : threads) { t.join(); }
        }
        consoleLog->info("done mapping reads.");
	    consoleLog->info("Discarded {} reads because they had > {} alignments",
		    hctrs.tooManyHits, maxNumHits.getValue());

	}

	if (haveOutputFile) {
	    outFile.close();
	}
	return 0;
    } catch (TCLAP::ArgException& e) {
	consoleLog->error("Exception [{}] when parsing argument {}", e.error(), e.argId());
	return 1;
    }

}



