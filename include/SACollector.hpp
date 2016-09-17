//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef SA_COLLECTOR_HPP
#define SA_COLLECTOR_HPP

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>

template <typename RapMapIndexT> class SACollector {
public:
  using OffsetT = typename RapMapIndexT::IndexType;

  /** Disable NIP skipping **/
  void disableNIP() { disableNIP_ = true; }

  /** Enable NIP skipping --- default state **/
  void enableNIP() { disableNIP_ = false; }
  
  /** Require a coverage fraction of at least req for all reported mappings **/
  void setCoverageRequirement(double req) { covReq_ = req; }

  /** Get the current coverage requirement for mappings (0 means no requirement) **/
  double getCoverageRequirement() const { return covReq_; }

  /** If any hit has a suffix array interval of this length or larger, just skip it **/
  void setMaxInterval(OffsetT maxInterval) { maxInterval_ = maxInterval; }

  /** Get the maximum allowable suffix array interval **/
  OffsetT getMaxInterval(OffsetT maxInterval) const { return maxInterval_; }
  
    bool getStrictCheck() const { return strictCheck_; };
    void setStrictCheck(bool sc) { strictCheck_ = sc; }
    
  /** Construct an SACollector given an index **/
  SACollector(RapMapIndexT* rmi) : rmi_(rmi), maxInterval_(1000) {}

    enum HitStatus { ABSENT = -1, UNTESTED = 0, PRESENT = 1 };
    // Record if k-mers are hits in the
    // fwd direction, rc direction or both
    struct KmerDirScore {
      KmerDirScore(rapmap::utils::my_mer kmerIn, int32_t kposIn,
                   HitStatus fwdScoreIn, HitStatus rcScoreIn)
          : kmer(kmerIn), kpos(kposIn), fwdScore(fwdScoreIn),
            rcScore(rcScoreIn) {}
      KmerDirScore() : kpos(0), fwdScore(UNTESTED), rcScore(UNTESTED) {}
      bool operator==(const KmerDirScore& other) const {
        return kpos == other.kpos;
      }
      bool operator<(const KmerDirScore& other) const {
        return kpos < other.kpos;
      }
      void print() {
        std::cerr << "{ " << kmer.to_str() << ", " << kpos << ", "
                  << ((fwdScore) ? "PRESENT" : "ABSENT") << ", "
                  << ((rcScore) ? "PRESENT" : "ABSENT") << "}\t";
      }
      rapmap::utils::my_mer kmer;
      int32_t kpos;
      HitStatus fwdScore;
      HitStatus rcScore;
    };





    bool operator()(std::string& read,
                    std::vector<rapmap::utils::QuasiAlignment>& hits,
                    SASearcher<RapMapIndexT>& saSearcher,
                    rapmap::utils::MateStatus mateStatus,
                    bool strictCheck = false, bool consistentHits = false) {
        
    using QuasiAlignment = rapmap::utils::QuasiAlignment;
    using MateStatus = rapmap::utils::MateStatus;

    auto& rankDict = rmi_->rankDict;
    auto& txpStarts = rmi_->txpOffsets;
    auto& SA = rmi_->SA;
    auto& khash = rmi_->khash;
    auto& text = rmi_->seq;
    auto salen = SA.size();
    auto hashEnd = khash.end();
    auto readLen = read.length();

    auto maxDist = 1.5 * readLen;
    auto k = rapmap::utils::my_mer::k();
    auto readStartIt = read.begin();
    auto readEndIt = read.end();

    auto rb = read.begin();
    auto re = rb + k;
    
    uint32_t fwdHit{0};
    uint32_t rcHit{0};

    size_t fwdCov{0};
    size_t rcCov{0};

    bool foundHit = false;
    bool isRev = false;

    rapmap::utils::my_mer mer;
    rapmap::utils::my_mer rcMer;

    // This allows implementing our heurisic for comparing
    // forward and reverse-complement strand matches
    std::vector<KmerDirScore> kmerScores;

    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;

    std::vector<SAIntervalHit> fwdSAInts;
    std::vector<SAIntervalHit> rcSAInts;

    std::vector<uint32_t> leftTxps, leftTxpsRC;
    std::vector<uint32_t> rightTxps, rightTxpsRC;

    // The number of bases that a new query position (to which
    // we skipped) should overlap the previous extension. A
    // value of 0 means no overlap (the new search begins at the next
    // base) while a value of (k - 1) means that k-1 bases (one less than
    // the k-mer size) must overlap.

    OffsetT skipOverlap = k-1;
    //OffsetT skipOverlap = 0;

    // Number of nucleotides to skip when encountering a homopolymer k-mer.
    OffsetT homoPolymerSkip = k / 2;
    auto merIt = hashEnd;
    auto rcMerIt = hashEnd;

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
      if (mer.is_homopolymer()) {
        rb += homoPolymerSkip;
        re += homoPolymerSkip;
        continue;
      }
      rcMer = mer.get_reverse_complement();

      // See if we can find this k-mer in the hash
      merIt = khash.find(mer.get_bits(0, 2 * k));
      rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));

      // If we can find the k-mer in the hash, get its SA interval
      if (merIt != hashEnd) {
          if (strictCheck) {
            ++fwdHit;
            // If we also match this k-mer in the rc direction
            if (rcMerIt != hashEnd) {
              ++rcHit;
              kmerScores.emplace_back(mer, pos, PRESENT, PRESENT);
            } else { // Otherwise it doesn't match in the rc direction
              kmerScores.emplace_back(mer, pos, PRESENT, ABSENT);
            }
          } else { // no strict check
            ++fwdHit;
            if (rcMerIt != hashEnd) {
              ++rcHit;
            }
          }
        }

      // See if the reverse complement k-mer is in the hash
      if (rcMerIt != hashEnd) {
          // The original k-mer didn't match in the foward direction
          if (!fwdHit) {
              ++rcHit;
              if (strictCheck) {
                  kmerScores.emplace_back(mer, pos, ABSENT, PRESENT);
              }
          }
      }

      // If we had a hit with either k-mer then we can
      // break out of this loop to look for the next informative position
      if (fwdHit + rcHit > 0) {
        foundHit = true;
        break;
      }
      ++rb;
      ++re;
    }

    // If we went the entire length of the read without finding a hit
    // then we can bail.
    if (!foundHit) {
      return false;
    }

    // If we had a hit on the forward strand
    if (fwdHit) {
        getSAHits_(
                  saSearcher,
                  read, // the read
                  rb, // where to start the search
                  &(merIt->second), // pointer to the search interval
                  fwdCov,
                  fwdHit,
                  rcHit,
                  fwdSAInts,
                  kmerScores,
                  false);
    }

    if (rcHit >= fwdHit) {
        rapmap::utils::reverseRead(read, rcBuffer_);
        getSAHits_(
                  saSearcher,
                  rcBuffer_, // the read
                  rcBuffer_.begin(), // where to start the search
                  nullptr, // pointer to the search interval
                  rcCov,
                  rcHit,
                  fwdHit,
                  rcSAInts,
                  kmerScores,
                  true);
    }

    if (strictCheck) {
      // The first two conditions shouldn't happen
      // but I'm just being paranoid here
      if (fwdHit > 0 and rcHit == 0) {
        rcSAInts.clear();
      } else if (rcHit > 0 and fwdHit == 0) {
        fwdSAInts.clear();
      } else {
        std::sort(kmerScores.begin(), kmerScores.end());
        auto e = std::unique(kmerScores.begin(), kmerScores.end());
        // Compute the score for the k-mers we need to
        // test in both the forward and rc directions.
        int32_t fwdScore{0};
        int32_t rcScore{0};
        // For every kmer score structure
        // std::cerr << "[\n";
        for (auto kmsIt = kmerScores.begin(); kmsIt != e;
             ++kmsIt) { //: kmerScores) {
          auto& kms = *kmsIt;
          // If the forward k-mer is untested, then test it
          if (kms.fwdScore == UNTESTED) {
            auto merIt = khash.find(kms.kmer.get_bits(0, 2 * k));
            kms.fwdScore = (merIt != hashEnd) ? PRESENT : ABSENT;
          }
          // accumulate the score
          fwdScore += kms.fwdScore;

          // If the rc k-mer is untested, then test it
          if (kms.rcScore == UNTESTED) {
            rcMer = kms.kmer.get_reverse_complement();
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));
            kms.rcScore = (rcMerIt != hashEnd) ? PRESENT : ABSENT;
          }
          // accumulate the score
          rcScore += kms.rcScore;
          // kms.print();
          // std::cerr << "\n";
        }
        // std::cerr << "]\n";
        // If the forward score is strictly greater
        // then get rid of the rc hits.
        if (fwdScore > rcScore) {
          rcSAInts.clear();
        } else if (rcScore > fwdScore) {
          // If the rc score is strictly greater
          // get rid of the forward hits
          fwdSAInts.clear();
        }
      }
    }

    // Coverage requirements only make sense if
    // we have disabled NIP skipping.
    if (covReq_ > 0.0 and disableNIP_) {
      double fwdFrac{0.0};
      double rcFrac{0.0};
      if (fwdSAInts.size() > 0) {
        fwdFrac = fwdCov / static_cast<double>(readLen);
        if (fwdFrac < covReq_) {
          fwdSAInts.clear();
        }
      }
      if (rcSAInts.size() > 0) {
        rcFrac = rcCov / static_cast<double>(readLen);
        if (rcFrac < covReq_) {
          rcSAInts.clear();
        }
      }
      /*
      if (fwdFrac + rcFrac > 0.0) {
          iomutex_.lock();
          std::cerr << fwdFrac << ", " << rcFrac << "\n";
          iomutex_.unlock();
      }
      */
    }

    auto fwdHitsStart = hits.size();
    // If we had > 1 forward hit
    if (fwdSAInts.size() > 1) {
      auto processedHits = rapmap::hit_manager::intersectSAHits(
          fwdSAInts, *rmi_, readLen, consistentHits);
      rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                               hits, mateStatus);
    } else if (fwdSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = fwdSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        hits.emplace_back(txpID, hitPos, true, readLen);
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
      auto newEnd = std::unique(
          hits.begin() + initialSize, hits.end(),
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto fwdHitsEnd = hits.size();

    auto rcHitsStart = fwdHitsEnd;
    // If we had > 1 rc hit
    if (rcSAInts.size() > 1) {
      auto processedHits = rapmap::hit_manager::intersectSAHits(
          rcSAInts, *rmi_, readLen, consistentHits);
      rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                               hits, mateStatus);
    } else if (rcSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = rcSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        hits.emplace_back(txpID, hitPos, false, readLen);
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
      auto newEnd = std::unique(
          sortStartIt, sortEndIt,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto rcHitsEnd = hits.size();

    // If we had both forward and RC hits, then merge them
    if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
      // Merge the forward and reverse hits
      std::inplace_merge(
          hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd,
          hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid < b.tid;
          });
      // And get rid of duplicate transcript IDs
      auto newEnd = std::unique(
          hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    // Return true if we had any valid hits and false otherwise.
    return foundHit;
    }

  /** Find the mappings for a read **/
  bool oldOperator(std::string& read,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  SASearcher<RapMapIndexT>& saSearcher,
                  rapmap::utils::MateStatus mateStatus,
                  bool strictCheck = false, bool consistentHits = false) {

    using QuasiAlignment = rapmap::utils::QuasiAlignment;
    using MateStatus = rapmap::utils::MateStatus;

    auto& rankDict = rmi_->rankDict;
    auto& txpStarts = rmi_->txpOffsets;
    auto& SA = rmi_->SA;
    auto& khash = rmi_->khash;
    auto& text = rmi_->seq;
    uint32_t sampFactor{1};
    auto salen = SA.size();
    auto hashEnd = khash.end();

    auto readLen = read.length();
    auto maxDist = 1.5 * readLen;
    auto k = rapmap::utils::my_mer::k();
    auto readStartIt = read.begin();
    auto readEndIt = read.end();

    auto readRevStartIt = read.rbegin();
    auto readRevEndIt = read.rend();

    auto rb = read.begin();
    auto re = rb + k;
    // Were the last "rightmost" hit ended on the query
    auto prevHitTail = read.begin();
    
    OffsetT lbLeftFwd = 0, ubLeftFwd = 0;
    OffsetT lbLeftRC = 0, ubLeftRC = 0;
    OffsetT lbRightFwd = 0, ubRightFwd = 0;
    OffsetT lbRightRC = 0, ubRightRC = 0;
    OffsetT matchedLen;

    uint32_t fwdHit{0};
    uint32_t rcHit{0};

    size_t fwdCov{0};
    size_t rcCov{0};
    bool foundHit = false;
    bool isRev = false;
    rapmap::utils::my_mer mer;
    rapmap::utils::my_mer rcMer;

    enum HitStatus { ABSENT = -1, UNTESTED = 0, PRESENT = 1 };
    // Record if k-mers are hits in the
    // fwd direction, rc direction or both
    struct KmerDirScore {
      KmerDirScore(rapmap::utils::my_mer kmerIn, int32_t kposIn,
                   HitStatus fwdScoreIn, HitStatus rcScoreIn)
          : kmer(kmerIn), kpos(kposIn), fwdScore(fwdScoreIn),
            rcScore(rcScoreIn) {}
      KmerDirScore() : kpos(0), fwdScore(UNTESTED), rcScore(UNTESTED) {}
      bool operator==(const KmerDirScore& other) const {
        return kpos == other.kpos;
      }
      bool operator<(const KmerDirScore& other) const {
        return kpos < other.kpos;
      }
      void print() {
        std::cerr << "{ " << kmer.to_str() << ", " << kpos << ", "
                  << ((fwdScore) ? "PRESENT" : "ABSENT") << ", "
                  << ((rcScore) ? "PRESENT" : "ABSENT") << "}\t";
      }
      rapmap::utils::my_mer kmer;
      int32_t kpos;
      HitStatus fwdScore;
      HitStatus rcScore;
    };

    // This allows implementing our heurisic for comparing
    // forward and reverse-complement strand matches
    std::vector<KmerDirScore> kmerScores;

    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;

    std::vector<SAIntervalHit> fwdSAInts;
    std::vector<SAIntervalHit> rcSAInts;

    std::vector<uint32_t> leftTxps, leftTxpsRC;
    std::vector<uint32_t> rightTxps, rightTxpsRC;

    // The number of bases that a new query position (to which
    // we skipped) should overlap the previous extension. A
    // value of 0 means no overlap (the new search begins at the next
    // base) while a value of (k - 1) means that k-1 bases (one less than
    // the k-mer size) must overlap.

    OffsetT skipOverlap = k-1;
    //OffsetT skipOverlap = 0;

    // Number of nucleotides to skip when encountering a homopolymer k-mer.
    OffsetT homoPolymerSkip = k / 2;

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
      if (mer.is_homopolymer()) {
        rb += homoPolymerSkip;
        re += homoPolymerSkip;
        continue;
      }
      rcMer = mer.get_reverse_complement();

      // See if we can find this k-mer in the hash
      auto merIt = khash.find(mer.get_bits(0, 2 * k));
      auto rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));

      // If we can find the k-mer in the hash, get its SA interval
      if (merIt != hashEnd) {
        OffsetT lb = merIt->second.begin();
        OffsetT ub = merIt->second.end();

        // lb must be 1 *less* then the current lb
        auto lbRestart = std::max(static_cast<OffsetT>(0), lb - 1);

        // Extend the SA interval using the read sequence as far as
        // possible
        std::tie(lbLeftFwd, ubLeftFwd, matchedLen) =
            saSearcher.extendSearchNaive(lbRestart, ub, k, rb, readEndIt);

        // If the SA interval is valid, and not too wide, then record
        // the hit.
        OffsetT diff = ubLeftFwd - lbLeftFwd;
        if (ubLeftFwd > lbLeftFwd and diff < maxInterval_) {
          auto queryStart = std::distance(read.begin(), rb);
          fwdSAInts.emplace_back(lbLeftFwd, ubLeftFwd, matchedLen, queryStart,
                                 false);
          // iomutex_.lock();
          // std::cerr << "matchedLen = " << matchedLen << " (qp = " <<
          // queryStart << ")\n";
          // iomutex_.unlock();
          fwdCov += matchedLen;
          if (strictCheck) {
            ++fwdHit;
            // If we also match this k-mer in the rc direction
            if (rcMerIt != hashEnd) {
              ++rcHit;
              kmerScores.emplace_back(mer, pos, PRESENT, PRESENT);
            } else { // Otherwise it doesn't match in the rc direction
              kmerScores.emplace_back(mer, pos, PRESENT, ABSENT);
            }

            // If we didn't end the match b/c we exhausted the query
            // test the mismatching k-mer in the other strand
            // TODO: check for 'N'?
            if (rb + matchedLen < readEndIt) {
              auto kmerPos =
                  std::distance(readStartIt, rb + matchedLen - skipOverlap);
              mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
              kmerScores.emplace_back(mer, kmerPos, ABSENT, UNTESTED);
            }
          } else { // no strict check
            ++fwdHit;
            if (rcMerIt != hashEnd) {
              ++rcHit;
            }
          }
        }
      }

      // See if the reverse complement k-mer is in the hash
      if (rcMerIt != hashEnd) {
        lbLeftRC = rcMerIt->second.begin();
        ubLeftRC = rcMerIt->second.end();
        OffsetT diff = ubLeftRC - lbLeftRC;
        if (ubLeftRC > lbLeftRC) {
          // The original k-mer didn't match in the foward direction
          if (!fwdHit) {
            ++rcHit;
            if (strictCheck) {
              kmerScores.emplace_back(mer, pos, ABSENT, PRESENT);
            }
          }
        }
      }

      // If we had a hit with either k-mer then we can
      // break out of this loop to look for the next informative position
      if (fwdHit + rcHit > 0) {
        foundHit = true;
        break;
      }
      ++rb;
      ++re;
    }

    // If we went the entire length of the read without finding a hit
    // then we can bail.
    if (!foundHit) {
      return false;
    }

    bool lastSearch{false};
    // If we had a hit on the forward strand
    if (fwdHit) {

      // The length of this match
      auto matchLen = fwdSAInts.front().len;
      // The iterator to where this match began
      rb = read.begin() + fwdSAInts.front().queryPos;

      // How much of the sequence is left after the
      // end of our match.
      auto remainingLength = std::distance(rb + matchLen, readEndIt);

      // We haven't done LCE yet, so any match here is "verified".
      // If we're fewer than k bases from the end, don't bother to
      // do any more work in this branch.
      if (remainingLength < k) { goto doneFwd; }
      
      // [lb, ub) is the suffix array interval for the MMP (maximum mappable
      // prefix) of the k-mer we found.  The NIP (next informative position)
      // in the sequence is the position after
      // the LCE (longest common extension) of T[SA[lb]:] and T[SA[ub-1]:]

      // If NIP skipping is disabled, then just treat the NIP as the MMP.
      auto lce = disableNIP_ ? matchLen
                             : saSearcher.lce(lbLeftFwd, ubLeftFwd - 1,
                                              matchLen, remainingLength);

      auto fwdSkip = std::max(static_cast<OffsetT>(matchLen) - skipOverlap,
                              static_cast<OffsetT>(lce) - skipOverlap);

      // The maximum position at which we could start our next search.
      size_t maxPos = static_cast<OffsetT>(std::distance(readStartIt, rb) + fwdSkip);
      size_t nextInformativePosition = maxPos;

      
      // If NIP skipping is *enabled*, and we got to the current position
      // by doing an LCE query, then we allow ourselves to *double check*
      // by querying the last k-mer in the read.
      // Otherwise, we just take the skip we're given.
      OffsetT nipSkip;
      if (!disableNIP_ and (lce > matchLen)) {
	size_t lastKmerPos = std::max(static_cast<OffsetT>(0),
				    static_cast<OffsetT>(readLen) - static_cast<OffsetT>(k));
	nextInformativePosition = std::min(lastKmerPos, nextInformativePosition);
      } 

      /*
      auto nipSkip = disableNIP_ ? static_cast<OffsetT>(readLen)
                                 : std::max(static_cast<OffsetT>(0),
                                            static_cast<OffsetT>(readLen) -
                                                static_cast<OffsetT>(k));
      size_t nextInformativePosition = std::min(
          nipSkip,
          static_cast<OffsetT>(std::distance(readStartIt, rb) + fwdSkip));
      */
      
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
          if (mer.is_homopolymer()) {
            rb += homoPolymerSkip;
            re = rb + k;
            continue;
          }
          auto merIt = khash.find(mer.get_bits(0, 2 * k));

          if (merIt != hashEnd) {
            if (strictCheck) {
              ++fwdHit;
              kmerScores.emplace_back(mer, pos, PRESENT, UNTESTED);
              auto rcMer = mer.get_reverse_complement();
              auto rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));
              if (rcMerIt != hashEnd) {
                ++rcHit;
                kmerScores.back().rcScore = PRESENT;
              }
            }

            lbRightFwd = merIt->second.begin();
            ubRightFwd = merIt->second.end();

            // lb must be 1 *less* then the current lb
            lbRightFwd = std::max(static_cast<OffsetT>(0), lbRightFwd - 1);
            std::tie(lbRightFwd, ubRightFwd, matchedLen) =
                saSearcher.extendSearchNaive(lbRightFwd, ubRightFwd, k, rb,
                                             readEndIt);

            OffsetT diff = ubRightFwd - lbRightFwd;
            if (ubRightFwd > lbRightFwd and diff < maxInterval_) {
              auto queryStart = std::distance(read.begin(), rb);
              fwdSAInts.emplace_back(lbRightFwd, ubRightFwd, matchedLen,
                                     queryStart, false);
              // iomutex_.lock();
              // std::cerr << "matchedLen = " << matchedLen << " (qp = " <<
              // queryStart << ")\n";
              // iomutex_.unlock();
              fwdCov += matchedLen;
              // If we didn't end the match b/c we exhausted the query
              // test the mismatching k-mer in the other strand
              // TODO: check for 'N'?
              if (strictCheck and rb + matchedLen < readEndIt) {
                auto kmerPos =
                    std::distance(readStartIt, rb + matchedLen - skipOverlap);
                mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
                // TODO: 04/11/16
                kmerScores.emplace_back(mer, kmerPos, UNTESTED, UNTESTED);
              }
            }

            if (lastSearch) {
              break;
            }
            auto mismatchIt = rb + matchedLen;
            if (mismatchIt < readEndIt) {
              auto remainingDistance = std::distance(mismatchIt, readEndIt);
              auto lce = disableNIP_
                             ? matchedLen
                             : saSearcher.lce(lbRightFwd, ubRightFwd - 1,
                                              matchedLen, remainingDistance);

	      // Where we would jump if we just used the MMP
	      auto skipMatch = mismatchIt - skipOverlap;
	      // Where we would jump if we used the LCE
	      auto skipLCE = rb + lce - skipOverlap;

	      auto maxSkip = std::max(skipMatch, skipLCE);
	      rb = maxSkip;
	      
	      // If NIP skipping is *enabled*, and we got to the current position
	      // by doing an LCE query, then we allow ourselves to *double check*
	      // by querying the last k-mer in the read.
	      // Otherwise, we just take the skip we're given.
	      if (!disableNIP_ and (lce > matchLen)) {
		if (readLen > k) {
		  rb = std::min(readEndIt - k, rb);
		}
	      } 
	      /* 
              if (rb > (readEndIt - k)) {
                rb = readEndIt - k;
                lastSearch = true;
                if (disableNIP_) {
                  valid = false;
                }
              }
	      */
              re = rb + k;
	      if (re == readEndIt) { lastSearch = true; }

            } else {
	      /*
              lastSearch = true;
              rb = readEndIt - k;
              re = rb + k;
              if (disableNIP_) {
                valid = false;
              }
	      */
	      goto doneFwd;
            }

          } else {
            rb += sampFactor;
            re = rb + k;
          }
        }
      }
    }
  doneFwd:

    lastSearch = false;
    if (rcHit >= fwdHit) {
      size_t pos{read.length() - k};

      auto revReadEndIt = read.rend();

      auto revRB = read.rbegin();
      auto revRE = revRB + k;

      auto invalidPosIt = revRB;
      while (revRE <= revReadEndIt) {

        revRE = revRB + k;
        if (revRE > revReadEndIt) {
          break;
        }

        // See if this k-mer would contain an N
        // only check if we don't yet know that there are no remaining
        // Ns
        if (invalidPosIt != revReadEndIt) {
          invalidPosIt = std::find_if(revRB, revRE, [](const char c) -> bool {
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
        if (mer.is_homopolymer()) {
          revRB += homoPolymerSkip;
          revRE += homoPolymerSkip;
          continue;
        }
        rcMer = mer.get_reverse_complement();
        auto rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));

        // If we found the k-mer
        if (rcMerIt != hashEnd) {
          if (strictCheck) {
            ++rcHit;
            kmerScores.emplace_back(mer, pos, UNTESTED, PRESENT);
            auto merIt = khash.find(mer.get_bits(0, 2 * k));
            if (merIt != hashEnd) {
              ++fwdHit;
              kmerScores.back().fwdScore = PRESENT;
            }
          }

          lbRightRC = rcMerIt->second.begin();
          ubRightRC = rcMerIt->second.end();

          // lb must be 1 *less* then the current lb
          // We can't move any further in the reverse complement direction
          lbRightRC = std::max(static_cast<OffsetT>(0), lbRightRC - 1);
          std::tie(lbRightRC, ubRightRC, matchedLen) =
              saSearcher.extendSearchNaive(lbRightRC, ubRightRC, k, revRB,
                                           revReadEndIt, true);

          OffsetT diff = ubRightRC - lbRightRC;
          if (ubRightRC > lbRightRC and diff < maxInterval_) {
            auto queryStart = std::distance(read.rbegin(), revRB);
            rcSAInts.emplace_back(lbRightRC, ubRightRC, matchedLen, queryStart,
                                  true);
            rcCov += matchedLen;
            // If we didn't end the match b/c we exhausted the query
            // test the mismatching k-mer in the other strand
            // TODO: check for 'N'?
            if (strictCheck and revRB + matchedLen < revReadEndIt) {
              auto kmerPos = std::distance(revRB + matchedLen, revReadEndIt);
              mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
              // TODO: 04/11/16
              kmerScores.emplace_back(mer, kmerPos, UNTESTED, UNTESTED);
            }
          }

          if (lastSearch) {
            break;
          }
          auto mismatchIt = revRB + matchedLen;
          if (mismatchIt < revReadEndIt) {
            auto remainingDistance = std::distance(mismatchIt, revReadEndIt);
            auto lce = disableNIP_
                           ? matchedLen
                           : saSearcher.lce(lbRightRC, ubRightRC - 1,
                                            matchedLen, remainingDistance);


	      // Where we would jump if we just used the MMP
	      auto skipMatch = mismatchIt - skipOverlap;
	      // Where we would jump if we used the LCE
	      auto skipLCE = revRB + lce - skipOverlap;

	      auto maxSkip = std::max(skipMatch, skipLCE);
	      revRB = maxSkip;
	      
	      // If NIP skipping is *enabled*, and we got to the current position
	      // by doing an LCE query, then we allow ourselves to *double check*
	      // by querying the last k-mer in the read.
	      // Otherwise, we just take the skip we're given.
	      if (!disableNIP_ and (lce > matchedLen)) {
		if (readLen > k) {
		  revRB = std::min(revReadEndIt - k, revRB);
		}
	      } 


	      revRE = revRB + k;
	      if (revRE == revReadEndIt) { lastSearch = true; }
	      /*
 
            // Where we would jump if we just used the MMP
            auto skipMatch = mismatchIt - skipOverlap;
            // Where we would jump if we used the lce
            auto skipLCE = revRB + lce - skipOverlap;
            // Choose the larger of the two
            revRB = std::max(skipLCE, skipMatch);
            if (revRB > (revReadEndIt - k)) {
              revRB = revReadEndIt - k;
              lastSearch = true;
              if (disableNIP_) {
                valid = false;
              }
            }
            revRE = revRB + k;
	      */
          } else {
	    /*
            lastSearch = true;
            revRB = revReadEndIt - k;
            revRE = revRB + k;
            if (disableNIP_) {
              valid = false;
            }
	    */
	    goto doneRC;
          }

        } else {
          revRB += sampFactor;
          revRE = revRB + k;
        }
      }
    }
  doneRC:

    if (strictCheck) {
      // The first two conditions shouldn't happen
      // but I'm just being paranoid here
      if (fwdHit > 0 and rcHit == 0) {
        rcSAInts.clear();
      } else if (rcHit > 0 and fwdHit == 0) {
        fwdSAInts.clear();
      } else {
        std::sort(kmerScores.begin(), kmerScores.end());
        auto e = std::unique(kmerScores.begin(), kmerScores.end());
        // Compute the score for the k-mers we need to
        // test in both the forward and rc directions.
        int32_t fwdScore{0};
        int32_t rcScore{0};
        // For every kmer score structure
        // std::cerr << "[\n";
        for (auto kmsIt = kmerScores.begin(); kmsIt != e;
             ++kmsIt) { //: kmerScores) {
          auto& kms = *kmsIt;
          // If the forward k-mer is untested, then test it
          if (kms.fwdScore == UNTESTED) {
            auto merIt = khash.find(kms.kmer.get_bits(0, 2 * k));
            kms.fwdScore = (merIt != hashEnd) ? PRESENT : ABSENT;
          }
          // accumulate the score
          fwdScore += kms.fwdScore;

          // If the rc k-mer is untested, then test it
          if (kms.rcScore == UNTESTED) {
            rcMer = kms.kmer.get_reverse_complement();
            auto rcMerIt = khash.find(rcMer.get_bits(0, 2 * k));
            kms.rcScore = (rcMerIt != hashEnd) ? PRESENT : ABSENT;
          }
          // accumulate the score
          rcScore += kms.rcScore;
          // kms.print();
          // std::cerr << "\n";
        }
        // std::cerr << "]\n";
        // If the forward score is strictly greater
        // then get rid of the rc hits.
        if (fwdScore > rcScore) {
          rcSAInts.clear();
        } else if (rcScore > fwdScore) {
          // If the rc score is strictly greater
          // get rid of the forward hits
          fwdSAInts.clear();
        }
      }
    }

    // Coverage requirements only make sense if
    // we have disabled NIP skipping.
    if (covReq_ > 0.0 and disableNIP_) {
      double fwdFrac{0.0};
      double rcFrac{0.0};
      if (fwdSAInts.size() > 0) {
        fwdFrac = fwdCov / static_cast<double>(readLen);
        if (fwdFrac < covReq_) {
          fwdSAInts.clear();
        }
      }
      if (rcSAInts.size() > 0) {
        rcFrac = rcCov / static_cast<double>(readLen);
        if (rcFrac < covReq_) {
          rcSAInts.clear();
        }
      }
      /*
      if (fwdFrac + rcFrac > 0.0) {
          iomutex_.lock();
          std::cerr << fwdFrac << ", " << rcFrac << "\n";
          iomutex_.unlock();
      }
      */
    }

    auto fwdHitsStart = hits.size();
    // If we had > 1 forward hit
    if (fwdSAInts.size() > 1) {
      auto processedHits = rapmap::hit_manager::intersectSAHits(
          fwdSAInts, *rmi_, readLen, consistentHits);
      rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                               hits, mateStatus);
    } else if (fwdSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = fwdSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        hits.emplace_back(txpID, hitPos, true, readLen);
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
      auto newEnd = std::unique(
          hits.begin() + initialSize, hits.end(),
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto fwdHitsEnd = hits.size();

    auto rcHitsStart = fwdHitsEnd;
    // If we had > 1 rc hit
    if (rcSAInts.size() > 1) {
      auto processedHits = rapmap::hit_manager::intersectSAHits(
          rcSAInts, *rmi_, readLen, consistentHits);
      rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                               hits, mateStatus);
    } else if (rcSAInts.size() == 1) { // only 1 hit!
      auto& saIntervalHit = rcSAInts.front();
      auto initialSize = hits.size();
      for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
        auto globalPos = SA[i];
        auto txpID = rmi_->transcriptAtPosition(globalPos);
        // the offset into this transcript
        auto pos = globalPos - txpStarts[txpID];
        int32_t hitPos = pos - saIntervalHit.queryPos;
        hits.emplace_back(txpID, hitPos, false, readLen);
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
      auto newEnd = std::unique(
          sortStartIt, sortEndIt,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    auto rcHitsEnd = hits.size();

    // If we had both forward and RC hits, then merge them
    if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
      // Merge the forward and reverse hits
      std::inplace_merge(
          hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd,
          hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid < b.tid;
          });
      // And get rid of duplicate transcript IDs
      auto newEnd = std::unique(
          hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid == b.tid;
          });
      hits.resize(std::distance(hits.begin(), newEnd));
    }
    // Return true if we had any valid hits and false otherwise.
    return foundHit;
  }




private:

inline void getSAHits_(
               SASearcher<RapMapIndexT>& saSearcher,
               std::string& read,
               std::string::iterator startIt,
               rapmap::utils::SAInterval<OffsetT>* startInterval,
               size_t& cov,
               uint32_t& strandHits,
               uint32_t& otherStrandHits,
               std::vector<rapmap::utils::SAIntervalHit<OffsetT>>& saInts,
               std::vector<KmerDirScore>& kmerScores,
               bool isRC // true if read is the reverse complement, false otherwise
               ) {
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;
    auto& khash = rmi_->khash;
    auto hashEnd = khash.end();
    auto strictCheck = strictCheck_;
    auto readLen = read.length();
    auto readStartIt = read.begin();
    auto readEndIt = read.end();
    OffsetT matchedLen{0};

    auto k = rapmap::utils::my_mer::k();
    auto skipOverlap = 0;//k - 1;
    OffsetT homoPolymerSkip = k / 2;

    auto rb = readStartIt;
    auto re = rb + k; 
    OffsetT lb, ub;
    size_t invalidPos{0};

    rapmap::utils::my_mer mer, complementMer;
    auto merIt = hashEnd;
    auto complementMerIt = hashEnd;
    size_t pos{0};
    size_t sampFactor{1};
    bool lastSearch{false};

    // If we have some place to start that we have already computed
    // then use it.
    bool canSkipSetup{startIt > readStartIt};

    if (canSkipSetup) {
        rb = startIt;
        invalidPos = std::distance(readStartIt, rb);
        lb = startInterval->begin();
        ub = startInterval->end();
        goto skipSetup;
    }

    while (re <= readEndIt) {
        // The distance from the beginning of the read to the
        // start of the k-mer
        pos = std::distance(readStartIt, rb);

        //re = rb + k;
        if (re > readEndIt) {
          break;
        }

        // See if this k-mer would contain an N
        // only check if we don't yet know that there are no remaining
        // Ns
        if (invalidPos != std::string::npos) {
          invalidPos = read.find_first_of("nN", pos);
        }

        // If the first N is within k bases, then this k-mer is invalid
        if (invalidPos < pos + k) {
          // Skip to the k-mer starting at the next position
          // (i.e. right past the N)
          rb = read.begin() + invalidPos + 1;
          re = rb + k;
          // Go to the next iteration of the while loop
          continue;
        }

        if (re <= readEndIt) {
        // Get the k-mer and query it in the hash
        mer = rapmap::utils::my_mer(read.c_str() + pos);
        if (mer.is_homopolymer()) {
          rb += homoPolymerSkip;
          re += homoPolymerSkip;
          continue;
        }
        complementMer = mer.get_reverse_complement();
        merIt = khash.find(mer.get_bits(0, 2 * k));

        // If we found the k-mer
        if (merIt != hashEnd) {
          if (strictCheck) {
            ++strandHits;
            kmerScores.emplace_back(mer, pos, ABSENT, ABSENT);
            if (isRC) {
                kmerScores.back().rcScore = PRESENT;
            } else {
                kmerScores.back().fwdScore = PRESENT;
            }

            complementMerIt = khash.find(complementMer.get_bits(0, 2 * k));
            if (complementMerIt != hashEnd) {
                ++otherStrandHits;
                if (isRC) { // if we are RC, the other strand is fwd
                    kmerScores.back().fwdScore = PRESENT;
                } else {
                    kmerScores.back().rcScore = PRESENT;
                }
            }
          }
        lb = merIt->second.begin();
        ub = merIt->second.end();
    skipSetup:
          // lb must be 1 *less* then the current lb
          // We can't move any further in the reverse complement direction
          lb = std::max(static_cast<OffsetT>(0), lb - 1);
          std::tie(lb, ub, matchedLen) =
              saSearcher.extendSearchNaive(lb, ub, k, rb, re);

          OffsetT diff = ub - lb;
          if (ub > lb and diff < maxInterval_) {
            uint32_t queryStart = static_cast<uint32_t>(std::distance(readStartIt, rb));
            saInts.emplace_back(lb, ub, matchedLen, queryStart, isRC);

            cov += matchedLen;

            // If we didn't end the match b/c we exhausted the query
            // test the mismatching k-mer in the other strand
            // TODO: check for 'N'?
            if (strictCheck_ and rb + matchedLen < readEndIt) {
              int32_t kmerPos =
                  static_cast<int32_t>(
                                       std::distance(readStartIt, 
                                       rb + matchedLen - skipOverlap));
              mer = rapmap::utils::my_mer(read.c_str() + kmerPos);
              // If we're on the reverse complement strand, then
              // we have to adjust kmerPos to be with respect to the 
              // forward strand.
              if (isRC) {
                  auto kp = kmerPos;
                  kmerPos = readLen - kp - k;
                  mer = mer.get_reverse_complement();
              }
              // TODO: 04/11/16
              // Can we change this yet?
              kmerScores.emplace_back(mer, kmerPos, UNTESTED, UNTESTED);
              if (isRC) {
                  kmerScores.back().rcScore = ABSENT;
              } else {
                  kmerScores.back().fwdScore = ABSENT;
              }
            }
          }

          if (lastSearch) { return; }

          auto mismatchIt = rb + matchedLen;
          if (mismatchIt < readEndIt) {
            auto remainingDistance = std::distance(mismatchIt, readEndIt);
            auto lce = disableNIP_
                           ? matchedLen
                           : saSearcher.lce(lb, ub - 1,
                                            matchedLen, remainingDistance);


            // Where we would jump if we just used the MMP
            auto skipMatch = mismatchIt - skipOverlap;
            // Where we would jump if we used the LCE
            auto skipLCE = rb + lce - skipOverlap;

            auto maxSkip = std::max(skipMatch, skipLCE);
            rb = maxSkip;
     
            // If NIP skipping is *enabled*, and we got to the current position
            // by doing an LCE query, then we allow ourselves to *double check*
            // by querying the last k-mer in the read.
            // Otherwise, we just take the skip we're given.
            if (!disableNIP_ and (lce > matchedLen)) {
                if (readLen > k) {
                    rb = std::min(readEndIt - k, rb);
                }
            } 


            re = rb + k;
            if (re == readEndIt) { lastSearch = true; }
          } else {
              return;
          }
        } else {
          rb += sampFactor;
          re = rb + k;
        }
    }
      }
}   

  RapMapIndexT* rmi_;
  bool disableNIP_{false};
  double covReq_{0.0};
  OffsetT maxInterval_;
  bool strictCheck_;
    std::string rcBuffer_;
  std::mutex iomutex_;
};

#endif // SA_COLLECTOR_HPP
