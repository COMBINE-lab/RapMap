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

#include "ScopedTimer.hpp"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "HitManager.hpp"
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
#include "stringpiece.h"

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "RapMapFileSystem.hpp"
#include "RapMapConfig.hpp"
#include "ScopedTimer.hpp"
#include "SpinLock.hpp"

//#define __TRACK_CORRECT__

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
		if (i == m or SA[c] + i == n) {
			if (i > prevIHigh) {
				prevIHigh = i;
				validBoundHigh = c;
			} else if (i == prevIHigh) {
				validBoundHigh = c < validBoundHigh ? c : validBoundHigh;
			}
		}

                if (plt) {
                    if (c == l + 1) {
                        std::cerr << "path 1\n";
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
                        std::cerr << "path 2\n";
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
                    std::cerr << "direction was UP; lb = " << l << ", ub = " << r << "\n";
                    std::cerr << "direction was UP; origLb = " << lbIn << ", origUb = " << ubIn << "\n";
                    break;
                case SearchDirection::DOWN:
                    sentinel = '{';
                    r = ubIn;
                    l = res1.bound;
                    std::cerr << "direction was DOWN; lb = " << l << ", ub = " << r << "\n";
                    std::cerr << "direction was UP; origLb = " << lbIn << ", origUb = " << ubIn << "\n";
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
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    char queryChar = (i < m - 1) ? ::toupper(*(qb + i)) : sentinel;
                    // If we're reverse complementing
                    if (queryChar != sentinel and complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }

                    if ( queryChar < *(sb + SA[c] + i) ) {
                     	break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) {
                        res2.dir = SearchDirection::DOWN;
                        res2.bound = c;
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        res2.bound = r;
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

	/**
	 * OK!  It should be (is) possible to figure out what we need with only two binary
	 * searches.  However, that seems to have some tricky corner cases and has been
	 * somewhat illusive so far.  This "naive" version performs *3* binary searches.
	 * The first determines the length of the maximum mappable prefix (MMP).  The second
	 * finds the lower bound for the query interval and the third finds the upper bound.
	 * The final binary search *is* optimized (it has a lower bound given by the value)
	 * returned by second search.  However, this method is likely a bit slower than the
	 * one above (when it can be made to work correctly at all times).
	 */
        template <typename IteratorT>
        std::tuple<int, int, int> extendSearchNaive(
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
            // Reduce the search interval until we hit a border
            // i.e. until c == r - 1 or c == l + 1
            while (true) {
                c = (l + r) / 2;
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

                        break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        if (i > prevILow) {
                            prevILow = i;
                            validBoundLow = c;
                        } else if (i == prevILow) {
                            validBoundLow = c > validBoundLow ? c : validBoundLow;
                        }
                        plt = false;
                        break;
                    }

                    ++i;
		}
		if (i == m or SA[c] + i == n) {
			if (i > prevIHigh) {
				prevIHigh = i;
				validBoundHigh = c;
			} else if (i == prevIHigh) {
				validBoundHigh = c < validBoundHigh ? c : validBoundHigh;
			}
		}

                if (plt) {
                    if (c == l + 1) {
                        auto maxI = std::max(std::max(i, prevILow), prevIHigh);
                        res1.maxLen = maxI;
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        maxI = std::max(std::max(i, prevILow), prevIHigh);
			res1.maxLen = maxI;
                        break;
                    }
                    l = c;
                    lcpLP = i;
                }
            }

            bool knownValid{true};
            m = res1.maxLen + 1;

	    // first search for the lower bound
            sentinel = '#';
	    l = lbIn;
	    r = ubIn;

            lcpLP = startAt;
            lcpRP = startAt;
            c = 0;
            plt = true;
            i = startAt;
            while (true) {
                c = (l + r) / 2;
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    char queryChar = (i < m - 1) ? ::toupper(*(qb + i)) : sentinel;
                    // If we're reverse complementing
                    if (queryChar != sentinel and complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }

                    if ( queryChar < *(sb + SA[c] + i) ) {
                     	break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) {
                        res1.bound = c;
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        res1.bound = r;
                        break;
                    }
                    l = c;
                    lcpLP = i;
                }
            }

	    // then search for the upper bound
            sentinel = '{';
	    l = res1.bound - 1;
	    r = ubIn;

            lcpLP = startAt;
            lcpRP = startAt;
            c = 0;
            plt = true;
            i = startAt;
            while (true) {
                c = (l + r) / 2;
                plt = true;
                i = std::min(lcpLP, lcpRP);
                while (i < m and SA[c] + i < n) {
                    char queryChar = (i < m - 1) ? ::toupper(*(qb + i)) : sentinel;
                    // If we're reverse complementing
                    if (queryChar != sentinel and complementBases) {
                        queryChar = rapmap::utils::my_mer::complement(queryChar);
                    }

                    if ( queryChar < *(sb + SA[c] + i) ) {
                     	break;
                    } else if ( queryChar > *(sb + SA[c] + i)) {
                        plt = false;
                        break;
                    }
                    ++i;
                }
                if (plt) {
                    if (c == l + 1) {
                        res2.bound = c;
                        break;
                    }
                    r = c;
                    lcpRP = i;
                } else {
                    if (c == r - 1) {
                        res2.bound = r;
                        break;
                    }
                    l = c;
                    lcpLP = i;
                }
            }

            // Must occur at least once!
            if (res1.bound == res2.bound) { res2.bound += 1; }
            return std::make_tuple(res1.bound, res2.bound, res1.maxLen);
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
            size_t len{startAt};
            auto o1 = SA[p1] + startAt;
            auto o2 = SA[p2] + startAt;
            auto maxIndex = std::max(o1, o2);
            while (maxIndex + len < textLen_ and seq[o1+len] == seq[o2+len]) {
                if (seq[o1+len] == '$') { break; }
                if (len >= stopAt) { break; }
                ++len;
            }
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
            MateStatus mateStatus,
            std::atomic<uint64_t>& diffCount) {

        //auto& posIDs = rmi_->positionIDs;
        auto& rankDict = rmi_->rankDict;
        auto& txpStarts = rmi_->txpOffsets;
        auto& SA = rmi_->SA;
        auto& khash = rmi_->khash;
        auto& text = rmi_->seq;
        uint32_t sampFactor{1};
        auto salen = SA.size();

        auto readLen = read.length();
        auto maxDist = 1.5 * readLen;
        auto k = rapmap::utils::my_mer::k();
        auto readStartIt = read.begin();
        auto readEndIt = read.end();

        auto readRevStartIt = read.rbegin();
        auto readRevEndIt = read.rend();

        auto rb = read.begin();
        auto re = rb + k;
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

        using rapmap::utils::SAIntervalHit;

        std::vector<SAIntervalHit> fwdSAInts;
        std::vector<SAIntervalHit> rcSAInts;

        std::vector<uint32_t> leftTxps, leftTxpsRC;
        std::vector<uint32_t> rightTxps, rightTxpsRC;
        int maxInterval{1000};

        // Find a hit within the read
        // While we haven't fallen off the end
        while (re < read.end()) {

            // Get the k-mer at the current start position.
            // And make sure that it's valid (contains no Ns).
            auto pos = std::distance(readStartIt, rb);
            auto invalidPos = read.find_first_of("nN", pos);
            if (invalidPos <= pos + k) {
                rb = read.begin() + invalidPos + 1;
                re = rb + k;
                continue;
            }

            // If the next k-bases are valid, get the k-mer and
            // reverse complement k-mer
            mer = rapmap::utils::my_mer(read.c_str() + pos);
            rcMer = mer.get_reverse_complement();

            // See if we can find this k-mer in the hash
            auto merIt = khash.find(mer.get_bits(0, 2*k));
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

            // If we can find the k-mer in the hash, get its SA interval
            if (merIt != khash.end()) {
                int lb = merIt->second.begin;
                int ub = merIt->second.end;

                // lb must be 1 *less* then the current lb
                auto lbRestart = std::max(static_cast<int>(0), lb-1);
                // Extend the SA interval using the read sequence as far as
                // possible
                std::tie(lbLeftFwd, ubLeftFwd, matchedLen) =
                    saSearcher.extendSearchNaive(lbRestart, ub, k, rb, readEndIt);
                /*
                   if (read.substr(pos, matchedLen) == text.substr(SA[lbLeftFwd - 1], matchedLen)) {
                   std::cerr << "INTERVAL LOOKS WRONG!\n";
                   std::cerr << "Start interval - 2 = " << text.substr(SA[lbLeftFwd - 2], matchedLen) << '\n';
                   std::cerr << "Start interval - 1 = " << text.substr(SA[lbLeftFwd - 1], matchedLen) << '\n';
                   std::cerr << "Interval start     = " << text.substr(SA[lbLeftFwd], matchedLen) << '\n';
                   }

                   if (read.substr(pos, matchedLen) == text.substr(SA[ubLeftFwd], matchedLen)) {
                   std::cerr << "INTERVAL LOOKS WRONG!\n";
                   std::cerr << "kmer                = " << mer << '\n';
                   std::cerr << "read substring      = " << read.substr(pos, matchedLen) << '\n';
                   std::cerr << "Before interval end = " << text.substr(SA[ubLeftFwd - 1], matchedLen) << '\n';
                   std::cerr << "Interval end        = " << text.substr(SA[ubLeftFwd], matchedLen) << '\n';
                   std::cerr << "Interval end + 1    = " << text.substr(SA[ubLeftFwd + 1], matchedLen) << '\n';
                   std::exit(1);
                   }
                //matchedLen = k;
                //lbLeftFwd = lb; ubLeftFwd = ub;
                */

                // If the SA interval is valid, and not too wide, then record
                // the hit.
                int diff = ubLeftFwd - lbLeftFwd;
                if (ubLeftFwd > lbLeftFwd and diff < maxInterval) {
                    auto queryStart = std::distance(read.begin(), rb);
                    fwdSAInts.emplace_back(lbLeftFwd, ubLeftFwd, matchedLen, queryStart, false);
                    leftFwdHit = true;
                }
            }

            // See if the reverse complement k-mer is in the hash
            if (rcMerIt != khash.end()) {
                lbLeftRC = rcMerIt->second.begin;
                ubLeftRC = rcMerIt->second.end;
                int diff = ubLeftRC - lbLeftRC;
                if (ubLeftRC > lbLeftRC and diff < maxInterval) {
                    /* Don't actually push here since we can't extend
                    auto queryStart = std::distance(read.begin(), rb);
                    rcSAInts.emplace_back(lbLeftRC, ubLeftRC, k, queryStart, true);
                    */
                    leftRCHit = true;
                }
            }

            // If we had a hit with either k-mer then we can
            // break out of this loop to look for the next informative position
            if (leftFwdHit or leftRCHit) {
                leftHit = true;
                break;
            }
            ++rb; ++re;
        }

        // If we went the entire length of the read without finding a hit
        // then we can bail.
        if (!leftHit) { return false; }

        // If we had a hit on the forward strand
        if (leftFwdHit) {

            // The length of this match
            auto matchLen = fwdSAInts.front().len;
            // The iterator to where this match began
            rb = read.begin() + fwdSAInts.front().queryPos;

            // [lb, ub) is the suffix array interval for the MMP (maximum mappable prefix)
            // of the k-mer we found.  The NIP (next informative position) in the sequence
            // is the position after the LCE (longest common extension) of
            // T[SA[lb]:] and T[SA[ub-1]:]
            auto remainingLength = std::distance(rb, readEndIt);
            auto lce = saSearcher.lce(lbLeftFwd, ubLeftFwd-1, matchLen, remainingLength);

            size_t nextInformativePosition = std::min(
                    std::max(0, static_cast<int>(readLen)- static_cast<int>(k)),
                    static_cast<int>(std::distance(readStartIt, rb) + lce)
                    );

            rb = read.begin() + nextInformativePosition;
            re = rb + k;

            while (re <= readEndIt) {
                // The offset into the string
                auto pos = std::distance(readStartIt, rb);
                //std::cerr << "readLen = " << readLen << ", lce = " << lce << "\n";
                //std::cerr << "here; pos = " << pos << "\n";
                // The position of the first N in the k-mer (if there is one)
                auto invalidPos = read.find_first_of("nN", pos);

                // If the first N is within k bases, then this k-mer is invalid
                if (invalidPos < pos + k) {
                    // A valid k-mer can't start until after the 'N'
                    nextInformativePosition = invalidPos + 1;
                    rb = read.begin() + nextInformativePosition;
                    re = rb + k;
                    // Go to the next iteration of the while loop
                    continue;
                }

                // If the current end position is valid
                if (re <= readEndIt) {

                    mer = rapmap::utils::my_mer(read.c_str() + pos);
                    auto merIt = khash.find(mer.get_bits(0, 2*k));

                    if (merIt != khash.end()) {
                        lbRightFwd = merIt->second.begin;
                        ubRightFwd = merIt->second.end;

                        // lb must be 1 *less* then the current lb
                        lbRightFwd = std::max(0, lbRightFwd - 1);
                        std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                            saSearcher.extendSearchNaive(lbRightFwd, ubRightFwd,
                                    k, rb, readEndIt);

                        int diff = ubRightFwd - lbRightFwd;
                        if (ubRightFwd > lbRightFwd and diff < maxInterval) {
                            auto queryStart = std::distance(read.begin(), rb);
                            fwdSAInts.emplace_back(lbRightFwd, ubRightFwd, matchedLen, queryStart, false);
                            rightFwdHit = true;
                            //break;
                        }

                        rb += matchedLen;
                        re = rb + k;

                        if (re <= readEndIt) {
                            rb -= matchedLen;
                            auto remainingDistance = std::distance(rb, readEndIt);
                            auto lce = saSearcher.lce(lbRightFwd, ubRightFwd-1, matchedLen, remainingDistance);
                            rb += lce;
                            re = rb + k;
                        }

                    } else {
                        rb += sampFactor;//re;
                        re = rb + k;
                    }
                }
            }
        }

        if (leftRCHit) {
            size_t pos{read.length() - k};

            auto revReadEndIt = read.rend();

            auto revRB = read.rbegin();
            auto revRE = revRB + k;

            while (revRE <= revReadEndIt){

                revRE = revRB + k;
                if (revRE >= revReadEndIt) { break; }

                // See if this k-mer would contain an N
                auto invalidPosIt = std::find_if(revRB, revRE,
                                                 [](const char c) -> bool {
                                                     return c == 'n' or c == 'N';
                                                 });
                // If we found an N before the end of the k-mer
                if (invalidPosIt < revRE) {
                    // Skip to the k-mer starting at the next position
                    // (i.e. right past the N)
                    revRB = invalidPosIt + 1;
                    continue;
                }

                // The distance from the beginning of the read to the
                // start of the k-mer
                pos = std::distance(revRE, revReadEndIt);

                // Get the k-mer and query it in the hash
                mer = rapmap::utils::my_mer(read.c_str() + pos);
                rcMer = mer.get_reverse_complement();
                auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

                // If we found the k-mer
                if (rcMerIt != khash.end()) {
                    lbRightRC = rcMerIt->second.begin;
                    ubRightRC = rcMerIt->second.end;

                    // lb must be 1 *less* then the current lb
                    // We can't move any further in the reverse complement direction
                    lbRightRC = std::max(0, lbRightRC - 1);
                    std::tie(lbRightRC, ubRightRC, matchedLen) =
                        saSearcher.extendSearchNaive(lbRightRC, ubRightRC, k,
                                revRB, revReadEndIt, true);

                    int diff = ubRightRC - lbRightRC;
                    if (ubRightRC > lbRightRC and diff < maxInterval) {
                        auto queryStart = std::distance(revRB + matchedLen, revReadEndIt);
                        rcSAInts.emplace_back(lbRightRC, ubRightRC, matchedLen, queryStart, true);
                        rightRCHit = true;
                        //break;
                    }
                    revRB += matchedLen;
                    revRE = revRB + k;

                    if (revRE <= revReadEndIt) {
                        revRB -= matchedLen;
                        auto remainingDistance = std::distance(revRB + matchedLen, revReadEndIt);
                        auto lce = saSearcher.lce(lbRightRC, ubRightRC-1, matchedLen, remainingDistance);
                        revRB += lce;
                        revRE = revRE + k;
                    }


                } else {
                    revRB += sampFactor;//= revRE;
                    revRE = revRB + k;
                }
            }
        }

        auto fwdHitsStart = hits.size();
        // If we had > 1 forward hit
        if (fwdSAInts.size() > 1) {
                auto processedHits = rapmap::hit_manager::intersectSAHits(fwdSAInts, *rmi_);
                rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, mateStatus);
        } else if (fwdSAInts.size() == 1) { // only 1 hit!
                auto& saIntervalHit = fwdSAInts.front();
                auto initialSize = hits.size();
                for (int i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
                        auto globalPos = SA[i];
                        //auto txpID = posIDs[globalPos];
			//auto txpID = rankDict.Rank(globalPos, 1);
			auto txpID = rmi_->transcriptAtPosition(globalPos);
                        // the offset into this transcript
                        auto pos = globalPos - txpStarts[txpID];
                        hits.emplace_back(txpID, pos, true, readLen);
                        hits.back().mateStatus = mateStatus;
                }
                // Now sort by transcript ID (then position) and eliminate
                // duplicates
                auto sortStartIt = hits.begin() + initialSize;
                auto sortEndIt = hits.end();
                std::sort(sortStartIt, sortEndIt,
                                [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                if (a.tid == b.tid) {
                                return a.pos < b.pos;
                                } else {
                                return a.tid < b.tid;
                                }
                                });
                auto newEnd = std::unique(hits.begin() + initialSize, hits.end(),
                                [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                return a.tid == b.tid;
                                });
                hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto fwdHitsEnd = hits.size();

        auto rcHitsStart = fwdHitsEnd;
        // If we had > 1 rc hit
        if (rcSAInts.size() > 1) {
            auto processedHits = rapmap::hit_manager::intersectSAHits(rcSAInts, *rmi_);
            rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist, hits, mateStatus);
        } else if (rcSAInts.size() == 1) { // only 1 hit!
            auto& saIntervalHit = rcSAInts.front();
            auto initialSize = hits.size();
            for (int i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
                auto globalPos = SA[i];
                //auto txpID = posIDs[globalPos];
		//auto txpID = rankDict.Rank(globalPos, 1);
		auto txpID = rmi_->transcriptAtPosition(globalPos);
                // the offset into this transcript
                auto pos = globalPos - txpStarts[txpID];
                hits.emplace_back(txpID, pos, false, readLen);
                hits.back().mateStatus = mateStatus;
            }
            // Now sort by transcript ID (then position) and eliminate
            // duplicates
            auto sortStartIt = hits.begin() + rcHitsStart;
            auto sortEndIt = hits.end();
            std::sort(sortStartIt, sortEndIt,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    if (a.tid == b.tid) {
                    return a.pos < b.pos;
                    } else {
                    return a.tid < b.tid;
                    }
                    });
            auto newEnd = std::unique(sortStartIt, sortEndIt,
                    [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid == b.tid;
                    });
            hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto rcHitsEnd = hits.size();

        // If we had both forward and RC hits, then merge them
        if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
            // Merge the forward and reverse hits
            std::inplace_merge(hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd, hits.begin() + rcHitsEnd,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid < b.tid;
                    });
            // And get rid of duplicate transcript IDs
            auto newEnd = std::unique(hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
                    [] (const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                    return a.tid == b.tid;
                    });
            hits.resize(std::distance(hits.begin(), newEnd));
        }
        // Return true if we had any valid hits and false otherwise.
        return (rcHitsEnd > fwdHitsStart);
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

    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpOffsets = rmi.txpOffsets;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};
    //uint32_t k = rapmap::utils::my_mer::k();
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> hits;

    size_t readLen{0};
	bool tooManyHits{false};
    uint16_t flags;
    // 1000-bp reads are max here (get rid of hard limit later).
    std::string readTemp(1000, 'A');
    std::string qualTemp(1000, '~');

    char buff1[1000];
    FixedWriter cigarStr(buff1, 1000);

    SASearcher saSearcher(&rmi);

    // 0 means properly aligned
    // 0x1 means only alignments for left read
    // 0x2 means only alignments for right read
    // 0x3 means "orphaned" alignments for left and right
    // (currently not treated as orphan).
    std::atomic<uint64_t> diffCount{0};
    uint32_t orphanStatus{0};
    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            readLen = j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            hitCollector(j->data[i].seq, hits, saSearcher, MateStatus::SINGLE_END, diffCount);
            auto numHits = hits.size();
            hctr.totHits += numHits;

            if (hits.size() > 0 and hits.size() < maxNumHits) {
                auto& readName = j->data[i].header;
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& txpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
                for (auto& qa : hits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    rapmap::utils::getSamFlags(qa, flags);
                    if (alnCtr != 0) {
                        flags |= 0x900;
                    }

                    std::string* readSeq = &(j->data[i].seq);
                    std::string* qstr = &(j->data[i].qual);
                    rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], cigarStr);

                    sstream << readName << '\t' // QNAME
                        << flags << '\t' // FLAGS
                        << transcriptName << '\t' // RNAME
                        << qa.pos + 1 << '\t' // POS (1-based)
                        << 255 << '\t' // MAPQ
                        << cigarStr.c_str() << '\t' // CIGAR
                        << '*' << '\t' // MATE NAME
                        << 0 << '\t' // MATE POS
                        << qa.fragLen << '\t' // TLEN
                        << *readSeq << '\t' // SEQ
                        << *qstr << '\n';
                    ++alnCtr;
                    // === SAM
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (txpNames[qa.tid] == txpName) { ++hctr.trueHits; }
#endif //__DEBUG__
                }
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()){
                    if (hctr.numReads > 0) {
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                        std::cerr << "\033[F\033[F\033[F";
#else
                        std::cerr << "\033[F\033[F";
#endif // __DEBUG__
                    }
                    std::cerr << "saw " << hctr.numReads << " reads\n";
                    std::cerr << "# hits per read = "
                        << hctr.totHits / static_cast<float>(hctr.numReads) << "\n";
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
        iomutex->lock();
        outStream << sstream.str();
        iomutex->unlock();
        sstream.clear();

    } // processed all reads


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
    auto& txpLens = rmi.txpLens;
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

            bool lh = hitCollector(j->data[i].first.seq,
                        leftHits, saSearcher,
                        MateStatus::PAIRED_END_LEFT, diffCount);
            bool rh = hitCollector(j->data[i].second.seq,
                        rightHits, saSearcher,
                        MateStatus::PAIRED_END_RIGHT, diffCount);
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
				    // The left and right transcipt ids
				    leftTxp = leftIt->tid;
				    rightTxp = rightIt->tid;

				    // They don't point to the same transcript
				    if (leftTxp < rightTxp) {
					    ++leftIt;
				    } else {

					    // The left and right iterators point to the same transcript
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
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            flags2 |= 0x900;
                        }
                        rapmap::utils::adjustOverhang(qa, txpLens[qa.tid], cigarStr1, cigarStr2);

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            rapmap::utils::reverseRead(j->data[i].first.seq,
                                        j->data[i].first.qual,
                                        read1Temp,
                                        qual1Temp);
                        }
                        if (!qa.mateIsFwd) {
                            rapmap::utils::reverseRead(j->data[i].second.seq,
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
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
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
                            rapmap::utils::reverseRead(*readSeq, *qstr,
                                        read1Temp, qual1Temp);
                        }

                        rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], *cigarStr);
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

            if (hctr.numReads > hctr.lastPrint + 1000000) {
		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads (diffCount = " << diffCount << ") : "
                              << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    std::cerr << ": true hit \% = "
                        << (100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads)));
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


	if (!noout.getValue()) {
        rapmap::utils::writeSAMHeader(rmi, outStream);
	}


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
	std::cerr << "\n\n";
        consoleLog->info("done mapping reads.");
	/*
	    consoleLog->info("Discarded {} reads because they had > {} alignments",
		    hctrs.tooManyHits, maxNumHits.getValue());
		    */

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


/*
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
*/



