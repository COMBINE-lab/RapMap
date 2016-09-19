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

  /** Get the current coverage requirement for mappings (0 means no requirement)
   * **/
  double getCoverageRequirement() const { return covReq_; }

  /** If any hit has a suffix array interval of this length or larger, just skip
   * it **/
  void setMaxInterval(OffsetT maxInterval) { maxInterval_ = maxInterval; }

  /** Get the maximum allowable suffix array interval **/
  OffsetT getMaxInterval(OffsetT maxInterval) const { return maxInterval_; }

  /** Get/Set usage of strict-checking **/
  bool getStrictCheck() const { return strictCheck_; };
  void setStrictCheck(bool sc) { strictCheck_ = sc; }

  /** Construct an SACollector given an index **/
  SACollector(RapMapIndexT* rmi) : rmi_(rmi), 
                                   disableNIP_(false),
                                   covReq_(0.0),
                                   maxInterval_(1000), 
                                   strictCheck_(false) {}

  enum HitStatus { ABSENT = -1, UNTESTED = 0, PRESENT = 1 };
  // Record if k-mers are hits in the
  // fwd direction, rc direction or both
  struct KmerDirScore {
    KmerDirScore(rapmap::utils::my_mer kmerIn, int32_t kposIn,
                 HitStatus fwdScoreIn, HitStatus rcScoreIn)
        : kmer(kmerIn), kpos(kposIn), fwdScore(fwdScoreIn), rcScore(rcScoreIn) {
    }
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
                  bool consistentHits = false) {

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

    bool useCoverageCheck{disableNIP_ and strictCheck_};

    // This allows implementing our heurisic for comparing
    // forward and reverse-complement strand matches
    std::vector<KmerDirScore> kmerScores;

    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;

    std::vector<SAIntervalHit> fwdSAInts;
    std::vector<SAIntervalHit> rcSAInts;

    std::vector<uint32_t> leftTxps, leftTxpsRC;
    std::vector<uint32_t> rightTxps, rightTxpsRC;

    // Number of nucleotides to skip when encountering a homopolymer k-mer.
    OffsetT homoPolymerSkip = k / 2;
    auto merIt = hashEnd;
    auto rcMerIt = hashEnd;
    size_t pos{0};
    size_t invalidPos{0};
    // Find a hit within the read
    // While we haven't fallen off the end
    while (re <= readEndIt) {
      // Get the k-mer at the current start position.
      // And make sure that it's valid (contains no Ns).
      pos = std::distance(readStartIt, rb);

      // See if this k-mer would contain an N
      // only check if we don't yet know that there are no remaining
      // Ns
      if (invalidPos != std::string::npos) {
          invalidPos = read.find_first_of("nN", pos);
          if (invalidPos <= pos + k) {
              rb = read.begin() + invalidPos + 1;
              re = rb + k;
              continue;
          }
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

      // If we can find the k-mer in the hash
      if (merIt != hashEnd) {
        
        if (strictCheck_) {
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
          if (strictCheck_) {
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
    
    bool didCheckFwd{false};
    // If we had a hit on the forward strand
    if (fwdHit) {
        didCheckFwd = true;
      getSAHits_(saSearcher,
                 read,             // the read
                 rb,               // where to start the search
                 &(merIt->second), // pointer to the search interval
                 fwdCov, fwdHit, rcHit, fwdSAInts, kmerScores, false);
    }

    // Good enough?
    //bool checkRC = (rcHit >= fwdHit);
    // otherwise use this one
    bool checkRC = useCoverageCheck ? (rcHit > 0) : (rcHit >= fwdHit);

    // If we had a hit on the reverse complement strand
    if (checkRC) {
      rapmap::utils::reverseRead(read, rcBuffer_);
      getSAHits_(saSearcher,
                 rcBuffer_,         // the read
                 rcBuffer_.begin(), // where to start the search
                 nullptr,           // pointer to the search interval
                 rcCov, rcHit, fwdHit, rcSAInts, kmerScores, true);
    }
    
    // Now, if we *didn't* check the forward strand at first, but we encountered fwd hits 
    // while looking at the RC strand, then check the fwd strand now
    
    // Good enough?
    //bool checkFwd = (fwdHit >= rcHit);
    // otherwise use this one
    bool checkFwd = useCoverageCheck ? (fwdHit > 0) : (fwdHit >= rcHit);
    if (!didCheckFwd and checkFwd) {
        didCheckFwd = true;
      getSAHits_(saSearcher,
                 read,             // the read
                 read.begin(),               // where to start the search
                 nullptr, // pointer to the search interval
                 fwdCov, fwdHit, rcHit, fwdSAInts, kmerScores, false);
    }

    if (strictCheck_) {
      // If we're computing coverage, then we can make use of that info here
      if (useCoverageCheck) {
        if (fwdCov > rcCov) {
          rcSAInts.clear();
        } else if (rcCov > fwdCov) {
          fwdSAInts.clear();
        }
      } else { // use the k-mer "spot check"
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
      SASearcher<RapMapIndexT>& saSearcher, std::string& read,
      std::string::iterator startIt,
      rapmap::utils::SAInterval<OffsetT>* startInterval, size_t& cov,
      uint32_t& strandHits, uint32_t& otherStrandHits,
      std::vector<rapmap::utils::SAIntervalHit<OffsetT>>& saInts,
      std::vector<KmerDirScore>& kmerScores,
      bool isRC // true if read is the reverse complement, false otherwise
      ) {
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT>;
    auto& khash = rmi_->khash;

    auto hashEnd = khash.end();
    auto readLen = read.length();
    auto readStartIt = read.begin();
    auto readEndIt = read.end();
    OffsetT matchedLen{0};

    auto k = rapmap::utils::my_mer::k();
    auto skipOverlapMMP = k - 1; 
    auto skipOverlapNIP = k - 1;
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
    size_t prevMMPEnd{0};

    // If we have some place to start that we have already computed
    // then use it.
    bool canSkipSetup{startInterval != nullptr};

    if (canSkipSetup) {
      rb = startIt;
      re = rb + k;
      pos = std::distance(readStartIt, rb);
      invalidPos = pos;
      lb = startInterval->begin();
      ub = startInterval->end();
      goto skipSetup;
    }

    while (re <= readEndIt) {
      // The distance from the beginning of the read to the
      // start of the k-mer
      pos = std::distance(readStartIt, rb);

      // See if this k-mer would contain an N
      // only check if we don't yet know that there are no remaining
      // Ns
      if (invalidPos != std::string::npos) {
        invalidPos = read.find_first_of("nN", pos);
        // If the first N is within k bases, then this k-mer is invalid
        if (invalidPos < pos + k) {
            // Skip to the k-mer starting at the next position
            // (i.e. right past the N)
            rb = read.begin() + invalidPos + 1;
            re = rb + k;
            // Go to the next iteration of the while loop
            continue;
        }
      }

      // Get the k-mer
      mer = rapmap::utils::my_mer(read.c_str() + pos);
      // If this is a homopolymer, then skip it
      if (mer.is_homopolymer()) {
        char nuc = std::toupper(*rb);
        while (std::toupper(*re) == nuc) {
          ++rb;
          ++re;
          // if we walk off the end, then we're done
          if (re > readEndIt) {
            return;
          }
        }
        // we found a different character
        pos = std::distance(readStartIt, rb);
        /*
      rb += homoPolymerSkip;
      re += homoPolymerSkip;
      // If the default skip jumps us off the end of the read
      // then try to check the last k-mer
      if (re >= readEndIt and !lastSearch) {
        rb = readEndIt - k;
        re = rb + k;
        // but give up if that's still a homopolymer
        lastSearch = true;
      }
      continue;
        */
      }
      // If it's not a homopolymer, then get the complement
      // k-mer and query both in the hash.
      complementMer = mer.get_reverse_complement();
      merIt = khash.find(mer.get_bits(0, 2 * k));

      // If we found the k-mer
      if (merIt != hashEnd) {
        if (strictCheck_) {
          ++strandHits;
          // If we're on the reverse complement strand, then
          // we have to adjust kmerPos to be with respect to the
          // forward strand.
          if (isRC) {
            auto kp = pos;
            pos = readLen - kp - k;
            mer = mer.get_reverse_complement();
          }
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
            saSearcher.extendSearchNaive(lb, ub, k, rb, readEndIt);

        OffsetT diff = ub - lb;
        if (ub > lb and diff < maxInterval_) {
          uint32_t queryStart =
              static_cast<uint32_t>(std::distance(readStartIt, rb));
          saInts.emplace_back(lb, ub, matchedLen, queryStart, isRC);

          size_t matchOffset = std::distance(readStartIt, rb);
          size_t correction = 0;
          
          // NOTE: prevMMPEnd points 1 position past the last *match* of the
          // previous MMP (i.e. it points to the *first mismatch*).  This is 
          // why we ignore the case where prevMMPEnd == matchOffset, and why 
          // we don't have to add 1 to correction.
          if (prevMMPEnd > matchOffset) {
            correction = prevMMPEnd - matchOffset;
          }
          // Update the coverage and position of the last MMP match
          cov += (matchedLen - correction);
          prevMMPEnd = matchOffset + matchedLen;

          // If we didn't end the match b/c we exhausted the query
          // test the mismatching k-mer in the other strand
          if (strictCheck_ and rb + matchedLen < readEndIt) {
            int32_t kmerPos = static_cast<int32_t>(
                std::distance(readStartIt, rb + matchedLen - skipOverlapMMP));

            // validNucs is true if mer contained no 'Ns'
            bool validNucs = mer.from_chars(read.c_str() + kmerPos);
            if (validNucs) {
                // no hit on the current strand

                // If we're on the reverse complement strand, then
                // we have to adjust kmerPos to be with respect to the
                // forward strand.
                if (isRC) {
                  auto kp = kmerPos;
                  kmerPos = readLen - kp - k;
                  mer = mer.get_reverse_complement();
                }
                kmerScores.emplace_back(mer, kmerPos, ABSENT, ABSENT);
                
                // Test the complementary strand
                complementMer = mer.get_reverse_complement();
                complementMerIt = khash.find(complementMer.get_bits(0, 2 * k));
                if (complementMerIt != hashEnd) {
                    ++otherStrandHits;
                    if (isRC) { // if we are RC, the other strand is fwd
                        kmerScores.back().fwdScore = PRESENT;
                    } else {
                        kmerScores.back().rcScore = PRESENT;
                    }
                }
            } // The k-mer was valid (had no Ns)
          } // we're performing a strict check
        } // This hit was worth recording --- occurred fewer then maxInterval_ times

        // If we've previously declared that the search that just occurred was our last, then we're done!
        if (lastSearch) {
          return;
        }

        // Otherwise, figure out how we should continue the search.
        auto mismatchIt = rb + matchedLen;
        // If we reached the end of the read, then we're done.
        if (mismatchIt >= readEndIt) {
            return;
        }

        auto remainingDistance = std::distance(mismatchIt, readEndIt);
        auto lce = disableNIP_ ? matchedLen
                               : saSearcher.lce(lb, ub - 1, matchedLen,
                                                remainingDistance);

        // Where we would jump if we just used the MMP
        auto skipMatch = mismatchIt - skipOverlapMMP;
        //if (skipMatch + k )
        // Where we would jump if we used the LCE
        auto skipLCE = rb + lce - skipOverlapNIP;
        // Pick the maximum of the two
        auto maxSkip = std::max(skipMatch, skipLCE);
        // And that's where our new search will start
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

        // If the search ends at the end of the read, then
        // set the flag that we need not try after this.
        if (re == readEndIt) {
          lastSearch = true;
        }

      } else { // If we couldn't match this k-mer, move on to the next.
        rb += sampFactor;
        re = rb + k;
      }
    }
  }

  RapMapIndexT* rmi_;
  bool disableNIP_;
  double covReq_;
  OffsetT maxInterval_;
  bool strictCheck_;
  std::string rcBuffer_;
};

#endif // SA_COLLECTOR_HPP
