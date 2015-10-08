#ifndef SA_COLLECTOR_HPP
#define SA_COLLECTOR_HPP

#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "SASearcher.hpp"

#include <algorithm>
#include <iterator>

template <typename RapMapIndexT>
class SACollector {
    public:
    using OffsetT = typename RapMapIndexT::IndexType;

    SACollector(RapMapIndexT* rmi) : rmi_(rmi) {}
    bool operator()(std::string& read,
            std::vector<rapmap::utils::QuasiAlignment>& hits,
            SASearcher<OffsetT>& saSearcher,
            rapmap::utils::MateStatus mateStatus) {

        using QuasiAlignment = rapmap::utils::QuasiAlignment;
        using MateStatus = rapmap::utils::MateStatus;

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
        OffsetT lbLeftFwd = 0, ubLeftFwd = 0;
        OffsetT lbLeftRC = 0, ubLeftRC = 0;
        OffsetT lbRightFwd = 0, ubRightFwd = 0;
        OffsetT lbRightRC = 0, ubRightRC = 0;
        OffsetT matchedLen;

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
        OffsetT maxInterval{1000};

        // The number of bases that a new query position (to which
        // we skipped) should overlap the previous extension. A
        // value of 0 means no overlap (the new search begins at the next
        // base) while a value of (k - 1) means that k-1 bases (one less than
        // the k-mer size) must overlap.
        OffsetT skipOverlap = k-1;
        // Number of nucleotides to skip when encountering a homopolymer k-mer.
        OffsetT homoPolymerSkip = k/2;

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
            if (mer.is_homopolymer()) { rb += homoPolymerSkip; re += homoPolymerSkip; continue; }
            rcMer = mer.get_reverse_complement();

            // See if we can find this k-mer in the hash
            auto merIt = khash.find(mer.get_bits(0, 2*k));
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));

            // If we can find the k-mer in the hash, get its SA interval
            if (merIt != khash.end()) {
                OffsetT lb = merIt->second.begin;
                OffsetT ub = merIt->second.end;

                // lb must be 1 *less* then the current lb
                auto lbRestart = std::max(static_cast<OffsetT>(0), lb-1);
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
                OffsetT diff = ubLeftFwd - lbLeftFwd;
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
                OffsetT diff = ubLeftRC - lbLeftRC;
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
            auto remainingLength = std::distance(rb + matchLen, readEndIt);
            auto lce = saSearcher.lce(lbLeftFwd, ubLeftFwd-1, matchLen, remainingLength);
            //auto fwdSkip = matchLen + 1 - skipOverlap;//lce - skipOverlap;
            auto fwdSkip = std::max(matchLen + 1 - skipOverlap, lce - skipOverlap);

            size_t nextInformativePosition = std::min(
                    std::max(0, static_cast<OffsetT>(readLen)- static_cast<OffsetT>(k)),
                    static_cast<OffsetT>(std::distance(readStartIt, rb) + fwdSkip)
                    );

            rb = read.begin() + nextInformativePosition;
            re = rb + k;

            size_t invalidPos{0};
            while (re <= readEndIt) {
                // The offset into the string
                auto pos = std::distance(readStartIt, rb);

                // The position of the first N in the k-mer (if there is one)
                // If we have already verified there are no Ns in the remainder
                // of the string (invalidPos is std::string::npos) then we can
                // skip this test.
                if (invalidPos != std::string::npos) {
                    invalidPos = read.find_first_of("nN", pos);
                }

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
                    if (mer.is_homopolymer()) { rb += homoPolymerSkip; re = rb + k; continue; }
                    auto merIt = khash.find(mer.get_bits(0, 2*k));

                    if (merIt != khash.end()) {
                        /*
                        if (!leftRCHit) {
                            auto rcMer = mer.get_reverse_complement();
                            auto rcMerIt = khash.find(rcMer.get_bits(0, 2*k));
                            if (rcMerIt != khash.end()) { leftRCHit = true; }
                        }
                        */

                        lbRightFwd = merIt->second.begin;
                        ubRightFwd = merIt->second.end;

                        // lb must be 1 *less* then the current lb
                        lbRightFwd = std::max(0, lbRightFwd - 1);
                        std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                            saSearcher.extendSearchNaive(lbRightFwd, ubRightFwd,
                                    k, rb, readEndIt);

                        OffsetT diff = ubRightFwd - lbRightFwd;
                        if (ubRightFwd > lbRightFwd and diff < maxInterval) {
                            auto queryStart = std::distance(read.begin(), rb);
                            fwdSAInts.emplace_back(lbRightFwd, ubRightFwd, matchedLen, queryStart, false);
                            rightFwdHit = true;
                        }

                        //fwdSkip = matchedLen - skipOverlap;

                        auto mismatchIt = rb + matchedLen;
                        if ((mismatchIt + 1 - skipOverlap) <= readEndIt) {
                            auto remainingDistance = std::distance(mismatchIt, readEndIt);
                            auto lce = saSearcher.lce(lbRightFwd, ubRightFwd-1, matchedLen, remainingDistance);
                            //rb += lce - skipOverlap;
                            //rb += std::max(static_cast<uint32_t>(fwdSkip), lce - skipOverlap);

                            // Where we would jump if we just used the MMP
                            auto skipMatch = mismatchIt + 1 - skipOverlap;
                            // Where we would jump if we used the LCE
                            auto skipLCE = rb + lce - skipOverlap;
                            // Pick the larger of the two
                            rb = std::max(skipLCE, skipMatch);
                            re = rb + k;
                        } else {
                            rb = mismatchIt + 1 - skipOverlap;
                            re = rb + k;
                        }

                    } else {
                        rb += sampFactor;
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

            auto invalidPosIt = revRB;
            while (revRE <= revReadEndIt){

                revRE = revRB + k;
                if (revRE > revReadEndIt) { break; }

                // See if this k-mer would contain an N
                // only check if we don't yet know that there are no remaining
                // Ns
                if (invalidPosIt != revReadEndIt) {
                    invalidPosIt = std::find_if(revRB, revRE,
                                                 [](const char c) -> bool {
                                                     return c == 'n' or c == 'N';
                                                 });
                }

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
                if (mer.is_homopolymer()) { revRB += homoPolymerSkip; revRE += homoPolymerSkip; continue; }
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

                    OffsetT diff = ubRightRC - lbRightRC;
                    if (ubRightRC > lbRightRC and diff < maxInterval) {
                        auto queryStart = std::distance(revRB + matchedLen, revReadEndIt);
                        rcSAInts.emplace_back(lbRightRC, ubRightRC, matchedLen, queryStart, true);
                        rightRCHit = true;
                    }

                    //auto fwdSkip = matchedLen - skipOverlap;
                    //auto fwdSkip = matchedLen + 1;
                    auto mismatchIt = revRB + matchedLen;

                    if ((mismatchIt + 1 - skipOverlap) <= revReadEndIt) {
                    //if (revRE <= revReadEndIt) {
                        auto remainingDistance = std::distance(mismatchIt, revReadEndIt);
                        auto lce = saSearcher.lce(lbRightRC, ubRightRC-1, matchedLen, remainingDistance);

                        // Where we would jump if we just used the MMP
                        auto skipMatch = mismatchIt + 1 - skipOverlap;
                        // Where we would jump if we used the lce
                        auto skipLCE = revRB + lce - skipOverlap;
                        // Choose the larger of the two
                        revRB = std::max(skipLCE, skipMatch);
                        revRE = revRB + k;
                    } else {
                        revRB = mismatchIt + 1 - skipOverlap;
                        revRE = revRB + k;
                    }

                } else {
                    revRB += sampFactor;
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
                for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
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
            for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
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

#endif // SA_COLLECTOR_HPP
