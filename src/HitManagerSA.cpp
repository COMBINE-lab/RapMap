#include "HitManagerSA.hpp"
#include "BooMap.hpp"
#include "FastxParser.hpp"
#include "FrugalBooMap.hpp"
#include "HitManager.hpp"
#include "SelectiveAlignment.hpp"

#include "edlib.h"
//#include "ssw_cpp.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>

#include <type_traits>

#define FRAG_LEN 1000 //assume fragment length of 1000

namespace rapmap {
namespace hit_manager_sa {

/*
template <typename RapMapIndexT, typename OffsetT>
int32_t alignRead(std::string& read,
        RapMapIndexT& rmi,
        OffsetT& pos){

    AlignerEngine ae_ ;
    auto& txpStarts = rmi_->txpOffsets ;
    auto& SA = rmi_->SA;
    auto& concatText = rmi_->seq;
    auto& txpLens = rmi_->txpLens ;
    auto& khash = rmi_->khash ;


}*/

std::string printOrient(rapmap::hit_manager_sa::MateSAStatus str){
  switch(str){

  case rapmap::hit_manager_sa::MateSAStatus::LEFT_END_FWD:
    return "LEFT_END_FWD" ;
  case rapmap::hit_manager_sa::MateSAStatus::LEFT_END_RC:
    return "LEFT_END_RC" ;
  case rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_FWD:
    return "RIGHT_END_FWD" ;
  case rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_RC:
    return "RIGHT_END_RC" ;
  
  }
  return "don't know" ;
}

template <typename RapMapIndexT>
bool rescueOrphan(fastx_parser::ReadPair& rpair,
                  SAHitMap& alignedMap,
                  rapmap::hit_manager_sa::MateSAStatus str,
                  std::vector<QuasiAlignment>& jointHits,
                  RapMapIndexT& rmi){

  bool flag = (rpair.first.name == "SRR2912443.678362") ;
  //std::cerr << rpair.first.name << "\n" ;

  using OffsetT = typename RapMapIndexT::IndexType;
  SAHitMap mateMap ;

  if(flag) std::cerr << "\nWe are in Orphan Rescue map size: "<<alignedMap.size() <<"\n" ;

  constexpr const int32_t signedZero{0};
  //constexpr const int32_t COV_THRESHOLD = rpair.first.seq.length() + rpair.second.seq.length() - 31 ;

  typedef struct covInfo {
  public:
    covInfo() {
      cov = 0;
      endPos = 0;
      rc = false;
      mmp = false;
    };
    // covInfo(SATxpQueryPos s,int c, int e,bool r) : saTxpQp(s), cov(c),
    // endPos(e) {rc = r;};
    covInfo(int c, int e, bool r, bool m) : cov(c), endPos(e), rc(r), mmp(m){};
    // SATxpQueryPos saTxpQp ;
    void operator=(const covInfo& CV) {
      cov = CV.cov;
      endPos = CV.endPos;
      rc = CV.rc;
      mmp = CV.mmp;
    }

    int cov;
    int endPos;
    bool rc;
    bool mmp;
    bool orphanActive{false} ;
    int32_t matePos{0} ;
    bool mateDir;

  };


  bool orphanRes{false} ;
  //AlignerEngine ae ;
  //find out transcripts from alignedMap end



  for(auto& ph : alignedMap){
    auto tid = ph.first ;
    if(alignedMap[tid].active){

      if(flag)std::cerr << "\nThis transcript is active\n" ;

      auto& tqvec = alignedMap[tid].tqvec ;
      //take the leftmost position within transcript
      std::sort(tqvec.begin(), tqvec.end(),
                [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                  return a.pos < b.pos;
                });
      //find the covered one
      std::map<int32_t, covInfo> hitCov ;
      for(auto& tq : tqvec){
        int numCov = 0 ;
        int32_t hitKey = tq.pos - tq.queryPos ;
        //not in dict
        if(hitCov.count(hitKey) == 0){
          covInfo cInfo(tq.matchedLen, tq.pos + tq.matchedLen,
                        tq.queryRC,
                        tq.matchedLen == rpair.first.seq.length());
          hitCov[hitKey] = cInfo ;
        }else{
          auto& curr = hitCov[hitKey] ;
          if(curr.endPos > tq.pos){
            curr.cov += tq.pos + tq.matchedLen - curr.endPos ;
          }else{
            curr.cov += tq.matchedLen ;
            curr.endPos = tq.pos + tq.matchedLen ;
          }
        }
      }//hash values for this transcript, now align
      using pair_type = typename decltype(hitCov)::value_type  ;
      //left most hit with most coverage
      /*
      auto pr = std::max_element
        (
         std::begin(hitCov), std::end(hitCov),
         [] (const pair_type & p1, const pair_type & p2) {
           return ((p1.second.cov < p2.second.cov) || ((p1.second.cov == p2.second.cov) && (p1.first < p2.first)));
         }
         );*/


      //test without string_view
      for(auto& pr : hitCov){

        auto aligneStartPos = pr.first ;

        auto mateStartPos = pr.first + 1 ;
        auto mateEndPos = mateStartPos + FRAG_LEN ;

        auto& txpStarts = rmi.txpOffsets ;

        auto& concatText = rmi.seq;
        //stx::string_view concatTextView(concatText);
        auto& txpLens = rmi.txpLens ;

        OffsetT globalPos = txpStarts[tid] ;
        OffsetT thisTxpLen = txpLens[tid] ;


        if(flag) std::cout << "Before arithmatic "<< rpair.first.name << "\t" << rpair.first.seq.length() << "\t" << rpair.second.seq.length() <<"\t mateEndPos: "
                  << mateEndPos << "\t mateStartPos: "
                  << mateStartPos  << "\t" << thisTxpLen <<"\t" << printOrient(str) <<"\n" ;


        std::string readToAlign;
        bool mateDir{false} ;
        // If LEFT_END_RC is already hanging
        // skip it
        // If RIGHT_END_FWD is hanging
        bool skip{false};

        switch(str){
        case rapmap::hit_manager_sa::MateSAStatus::LEFT_END_FWD :
          readToAlign = rapmap::utils::reverseComplement(rpair.second.seq) ;
          mateStartPos = std::max(0, mateStartPos) ;
          mateEndPos = std::min(static_cast<int32_t>(thisTxpLen) , mateEndPos) ;
          hitCov[aligneStartPos].mateDir = false ;
          break ;
        case rapmap::hit_manager_sa::MateSAStatus::LEFT_END_RC :
          mateStartPos = std::max(aligneStartPos - FRAG_LEN + static_cast<int32_t>(rpair.second.seq.length()), 0) ;
          mateEndPos = std::min(aligneStartPos + static_cast<int32_t>(rpair.second.seq.length()) - 1, static_cast<int32_t>(thisTxpLen)) ;
          readToAlign = rpair.second.seq ;
          hitCov[aligneStartPos].mateDir = true;
          break ;
        case rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_FWD :
          mateStartPos = std::max(0,mateStartPos) ;
          mateEndPos = std::min(static_cast<int32_t>(thisTxpLen), mateEndPos) ;
          readToAlign = rapmap::utils::reverseComplement(rpair.first.seq) ;
          hitCov[aligneStartPos].mateDir = true ;
          break ;
        case rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_RC :
          mateStartPos = std::max(aligneStartPos - FRAG_LEN + static_cast<int32_t>(rpair.first.seq.length()), 0) ;
          mateEndPos = std::min(aligneStartPos + static_cast<int32_t>(rpair.first.seq.length()) - 1, static_cast<int32_t>(thisTxpLen)) ;
          readToAlign = rpair.first.seq ;
          hitCov[aligneStartPos].mateDir = false ;
          break ;

        }

        if(mateEndPos - mateStartPos + 1 < static_cast<int32_t>(readToAlign.length()/2)){
          continue ;
        }

        const char* thisTxpSeq = concatText.data() + globalPos + mateStartPos;
        bool localOrphanActive{false} ;

        if(flag) std::cout << "After arithmatic \t"<< readToAlign.length() << "\t mateEndPos: " << mateEndPos << "\t mateStartPos: "
                           << mateStartPos  <<"\t" <<(mateEndPos - mateStartPos + 1) << "\t" << thisTxpLen << "\t skip " << skip << "\n" ;
        auto thisEdlibResult = edlibAlign(readToAlign.c_str(), readToAlign.length(), thisTxpSeq,
        (mateEndPos-mateStartPos+1) , edlibNewAlignConfig(10, EDLIB_MODE_HW, EDLIB_TASK_LOC)); // hardcoded for now

        auto thisEditDistance = thisEdlibResult.editDistance ;
        int* startLocations = thisEdlibResult.startLocations ;

        if(startLocations != NULL and thisEditDistance != -1){
          hitCov[aligneStartPos].cov += (readToAlign.length() - thisEditDistance) ;
          hitCov[aligneStartPos].matePos = startLocations[0] + mateStartPos ;
          hitCov[aligneStartPos].orphanActive = true ;
          localOrphanActive = true ;
        }

        /*
        if(thisEditDistance == -1 and startLocations == NULL){
          std::cerr << "\n What the hell\n" ;
        }else{
          std::cerr << "\n working " << thisEditDistance <<"\n" ;
          }*/
        


        /*
        //Using https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/
        //instead of edlib

        int32_t maskLen = rpair.first.seq.length()/2 ;
        maskLen = maskLen < 15 ? 15 : maskLen ;
        // Declares a default Aligner
        StripedSmithWaterman::Aligner aligner;
        // Declares a default filter
        StripedSmithWaterman::Filter filter;
        // Declares an alignment that stores the result
        StripedSmithWaterman::Alignment alignment;
        // Aligns the query to the ref
        bool success = aligner.Align(readToAlign.c_str(),
                                    thisTxpSeq,
                                    (mateEndPos - mateStartPos + 1),
                                    filter,
                                    &alignment,
                                    maskLen);
        */


        if(flag) std::cerr << "\nAbout to rescue orphan\t"<< thisEditDistance << "\n" ;
       
        if(flag and localOrphanActive) std::cerr << "\naligned end maps to\ttranscript: "<<tid<<"\tAligned Pos: "
                          <<pr.first<<"\tmate maps from: "
                          <<mateStartPos<<"\t"
                          <<"ends at: "<<mateEndPos
                          <<"\ttxpLen "<<thisTxpLen<<"\n"
                          <<rpair.first.seq<<"\n"
                          <<rpair.second.seq<<"\n"
                          <<readToAlign<<"\n"
                          <<std::string(thisTxpSeq,1000)<<"\n";



      }

      //take out the active and most covered candidate, if there is a tie take the left most
      //first erase the parts that are non active
      auto maxIt = hitCov.begin() ;
      for(auto it = hitCov.begin() ; it != hitCov.end() ; ++it){
        if(!maxIt->second.orphanActive){
          maxIt = it ;
        }else if(it->second.orphanActive and it->second.cov > maxIt->second.cov){
          maxIt = it ;
        }
      }

      //fill in jointHits
      if(maxIt->second.orphanActive){
        if(flag) std::cerr << "\nOrphan rescued\n" ;
        //orphan rescued
        orphanRes = true ;

        auto hitPos1 = maxIt->first  ;
        auto hitPos2 = maxIt->second.matePos ;
        int32_t startRead1 = std::max(hitPos1, signedZero);
        int32_t startRead2 = std::max(hitPos2, signedZero);

        bool read1First{(startRead1 < startRead2)};
        int32_t fragStartPos = read1First ? startRead1 : startRead2;
        int32_t fragEndPos = read1First
                                  ? (startRead2 + rpair.second.seq.length())
                                  : (startRead1 + rpair.first.seq.length());
        uint32_t fragLen = fragEndPos - fragStartPos;
        if(str == rapmap::hit_manager_sa::MateSAStatus::LEFT_END_FWD || str == rapmap::hit_manager_sa::MateSAStatus::LEFT_END_RC){
          jointHits.emplace_back(tid, hitPos1, !hitCov[hitPos1].rc,
                                rpair.first.seq.length(),
                                0, // leftHitCov[vp.first].saTxpQp.lcpLength,
                                fragLen, true);
          auto& qaln = jointHits.back();

          qaln.mateLen = rpair.second.seq.length();
          qaln.matePos = hitPos2;

          qaln.mateIsFwd = !(hitCov[hitPos1].mateDir); //! vp.second.queryRC ;
          qaln.mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
          qaln.editD = 0;
          qaln.mateEditD = 0;
          qaln.resOrphan = true ;

        }else{
          jointHits.emplace_back(tid, hitPos2, !hitCov[hitPos1].mateDir,
                                  rpair.first.seq.length(),
                                  0, // leftHitCov[vp.first].saTxpQp.lcpLength,
                                  fragLen, true);

          auto& qaln = jointHits.back();

          qaln.mateLen = rpair.second.seq.length();
          qaln.matePos = hitPos2;

          qaln.mateIsFwd = !(hitCov[hitPos1].rc); //! vp.second.queryRC ;
          qaln.mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
          qaln.editD = 0;
          qaln.mateEditD = 0;
          qaln.resOrphan = true ;
        }

      }//edit distance checked


    } //end active check
    else{
      if(flag) std::cerr << "\nThis transcript is not active\n" ;
    }

  } // for loop

  if(jointHits.size() > 0){
    std::sort(jointHits.begin(), jointHits.end(),
              [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                return a.tid < b.tid;
              } );
  }

  return orphanRes ;

}

template <typename RapMapIndexT>
bool mergeLeftRightMap(fastx_parser::ReadPair& rpair, SAHitMap& leftMap,
                       SAHitMap& rightMap,
                       std::vector<QuasiAlignment>& jointHits,
                       uint32_t editDistance, RapMapIndexT& rmi, uint32_t maxInsertSize_) {
  // using OffsetT = typename RapMapIndexT::IndexType;
  typedef struct covInfo {
  public:
    covInfo() {
      cov = 0;
      endPos = 0;
      rc = false;
      mmp = false;
    };
    // covInfo(SATxpQueryPos s,int c, int e,bool r) : saTxpQp(s), cov(c),
    // endPos(e) {rc = r;};
    covInfo(int c, int e, bool r, bool m) : cov(c), endPos(e), rc(r), mmp(m){};
    // SATxpQueryPos saTxpQp ;
    void operator=(const covInfo& CV) {
      cov = CV.cov;
      endPos = CV.endPos;
      rc = CV.rc;
      mmp = CV.mmp;
    }

    int cov;
    int endPos;
    bool rc;
    bool mmp;
    std::vector<int32_t> gapsBegin;
    std::vector<int32_t> gapsEnd;

 };


  //debug
  //std::string testRead = "GCGAGAGTCTGTCTCAAAAACAACAACCACCACCACCACCACCACAACAACAACGAAACTATCTGCAAACTTTGC";

  constexpr const int32_t signedZero{0};
  //uint32_t maxInsertSize_{1000};
  bool foundHit{false};

//  SECollector<RapMapIndexT> hitSECollector(&rmi);
  std::vector<QuasiAlignment> hits;

  std::set<uint32_t> commonTids;
  std::set<uint32_t> leftTids;
  for (auto& phLeft : leftMap) {
    auto leftTid = phLeft.first;
    leftTids.insert(leftTid);
  }

  for (auto& phRight : rightMap) {
    auto rightTid = phRight.first;
    if (leftTids.find(rightTid) != leftTids.end()) {
      commonTids.insert(rightTid);
    }
  }

  // This seems redundant --- if there are
  // no commonTids, then foundHit will remain the
  // default value (false) listed above.
  //if (commonTids.size() > 0) {


    // for each common tid check
    // if both the end tell us
    // something
    for (auto tid : commonTids) {

      if (leftMap[tid].active and rightMap[tid].active) {

        auto& lefttqvec = leftMap[tid].tqvec;
        auto& righttqvec = rightMap[tid].tqvec;

        std::sort(lefttqvec.begin(), lefttqvec.end(),
                  [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                    return a.pos < b.pos;
                  });

        std::sort(righttqvec.begin(), righttqvec.end(),
                  [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                    return a.pos < b.pos;
                  });

        // look at the things in pair
        // Make possible pairs
        // Discard pairs which are
        // impossible in the sense
        // that they are far(1000 bp) apart

        std::map<int32_t, covInfo> leftHitCov;
        {
          auto it = lefttqvec.begin();
          while (it != lefttqvec.end()) {
            int numOfCov = 0;
            int32_t hitKey = it->pos - it->queryPos;
            if (leftHitCov.count(hitKey) == 0) {
              covInfo cInfo(it->matchedLen, it->pos + it->matchedLen,
                            it->queryRC,
                            it->matchedLen == rpair.first.seq.length());
              leftHitCov[hitKey] = cInfo;
	      if(it->queryPos>0){
		leftHitCov[hitKey].gapsBegin.push_back(0);
	        leftHitCov[hitKey].gapsEnd.push_back(it->queryPos);
	      }
            } else {
              auto& curr = leftHitCov[hitKey];
              if (curr.endPos > it->pos) {
                curr.cov += it->pos + it->matchedLen - curr.endPos;
                curr.endPos = it->pos + it->matchedLen;
              } else {
 	        if(it->pos > curr.endPos){
		  leftHitCov[hitKey].gapsBegin.push_back(curr.endPos-hitKey);
	          leftHitCov[hitKey].gapsEnd.push_back(it->pos-hitKey);
	        }
                curr.cov += it->matchedLen;
                curr.endPos = it->pos + it->matchedLen;
              }
            }

            ++it;
          }
	  for(auto it = leftHitCov.begin(); it!=leftHitCov.end(); ++it){
	    if(it->second.endPos - it->first < rpair.first.seq.length()){
	      it->second.gapsBegin.push_back(it->second.endPos - it->first);
	      it->second.gapsEnd.push_back(rpair.first.seq.length());
	    }
	  }
        }


        std::map<int32_t, covInfo> rightHitCov;
        {
          auto it = righttqvec.begin();
          while (it != righttqvec.end()) {
            int numOfCov = 0;

            int32_t hitKey = it->pos - it->queryPos;
            if (rightHitCov.count(hitKey) == 0) {
              covInfo cInfo(it->matchedLen, it->pos + it->matchedLen,
                            it->queryRC,
                            it->matchedLen == rpair.second.seq.length());
              rightHitCov[hitKey] = cInfo;
	      if(it->queryPos>0){
		rightHitCov[hitKey].gapsBegin.push_back(0);
	        rightHitCov[hitKey].gapsEnd.push_back(it->queryPos);
	      }
            } else {
              auto& curr = rightHitCov[hitKey];
              if (curr.endPos > it->pos) {
                curr.cov += (it->pos + it->matchedLen - curr.endPos);
                curr.endPos = it->pos + it->matchedLen;
              } else {
                if(it->pos > curr.endPos){
		  rightHitCov[hitKey].gapsBegin.push_back(curr.endPos-hitKey);
	          rightHitCov[hitKey].gapsEnd.push_back(it->pos-hitKey);
	        }
                curr.cov += it->matchedLen;
                curr.endPos = it->pos + it->matchedLen;

              }
            }
            ++it;
          }
	  for(auto it = rightHitCov.begin(); it!=rightHitCov.end(); ++it){

	    if(it->second.endPos - it->first < rpair.second.seq.length()){
	      it->second.gapsBegin.push_back(it->second.endPos - it->first);
	      it->second.gapsEnd.push_back(rpair.second.seq.length());
	    }
	  }
        }

        int maxCovScore = 0;
        std::vector<std::pair<int32_t, int32_t>> validPairs;
        bool leftFirst =
            leftHitCov.cbegin()->first - rightHitCov.cbegin()->first < 0;
        if (leftFirst) {
          for (auto& leftTQPos : leftHitCov) {
            int32_t hitPos1 = leftTQPos.first;  // pos - leftTQPos.queryPos ;
            if (std::abs(hitPos1 - rightHitCov.cbegin()->first) >
                maxInsertSize_) {
              continue;
            }
            for (auto& rightTQPos : rightHitCov) {

              // Check for the distance
              // go with it otherwise
              // go with something else
              int32_t hitPos2 = rightTQPos.first; // pos - rightTQPos.queryPos ;

              if (std::abs(hitPos1 - hitPos2) < maxInsertSize_) {
                validPairs.push_back(std::make_pair(hitPos1, hitPos2));
                auto score = leftTQPos.second.cov + rightTQPos.second.cov;
                if (maxCovScore < score) { maxCovScore = score; }
              } else {
                break;
              }
            }
          } // end for*/

        } else {

          for (auto& rightTQPos : rightHitCov) {
            int32_t hitPos2 = rightTQPos.first; // pos - rightTQPos.queryPos ;
            if (std::abs(hitPos2 - leftHitCov.cbegin()->first) >
                maxInsertSize_) {
              continue;
            }
            for (auto& leftTQPos : leftHitCov) {
              // Check for the distance
              // go with it otherwise
              // go with something else
              int32_t hitPos1 = leftTQPos.first;  // pos - leftTQPos.queryPos ;

              if (std::abs(hitPos1 - hitPos2) < maxInsertSize_) {
                validPairs.push_back(std::make_pair(hitPos1, hitPos2));
                auto score = leftTQPos.second.cov + rightTQPos.second.cov;
                if (maxCovScore < score) { maxCovScore = score; }
              } else {
                break;
              }
            }
          } // end for*/
        }

        // Now I have valid pairs so I can directly make QuasiAlignment vectors
        //
        if (validPairs.size() > 0) {
          foundHit = true;
        }

        hits.clear();
        for (auto& vp : validPairs) {
          // I will try to make QuasiAlignment objects
          // real alignments and edit distance are not present
          // yet


          int32_t hitPos1 = vp.first;  // vp.first.pos - vp.first.queryPos ;
          int32_t hitPos2 = vp.second; // vp.second.pos - vp.second.queryPos ;



          if (leftHitCov[hitPos1].cov + rightHitCov[hitPos2].cov ==
              maxCovScore) {
            int32_t startRead1 = std::max(hitPos1, signedZero);
            int32_t startRead2 = std::max(hitPos2, signedZero);

            bool read1First{(startRead1 < startRead2)};
            int32_t fragStartPos = read1First ? startRead1 : startRead2;
            int32_t fragEndPos = read1First
                                     ? (startRead2 + rpair.second.seq.length())
                                     : (startRead1 + rpair.first.seq.length());
            uint32_t fragLen = fragEndPos - fragStartPos;
            jointHits.emplace_back(tid, hitPos1, !leftHitCov[vp.first].rc,
                                   rpair.first.seq.length(),
                                   0, // leftHitCov[vp.first].saTxpQp.lcpLength,
                                   fragLen, true);
            auto& qaln = jointHits.back();
            // qaln.fwd = !leftHitCov[vp.first].rc;
            qaln.mateLen = rpair.second.seq.length();
            qaln.matePos = hitPos2;
            qaln.mateIsFwd =
                !(rightHitCov[vp.second].rc); //! vp.second.queryRC ;
            qaln.mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
            if (leftHitCov[hitPos1].mmp) {
              qaln.toAlign = true;
              qaln.editD = 0;
            }
            if (rightHitCov[hitPos2].mmp) {
              qaln.mateToAlign = true;
              qaln.mateEditD = 0;
            }

	    qaln.gapsBegin = leftHitCov[vp.first].gapsBegin;
	    qaln.mateGapsBegin = rightHitCov[vp.second].gapsBegin;
	    qaln.gapsEnd = leftHitCov[vp.first].gapsEnd;
	    qaln.mateGapsEnd = rightHitCov[vp.second].gapsEnd;
            break;
          }
        }

      } // end if
    }   // end for

    // This seems redundant
    //} else {
    //foundHit = false;
    //}

  // since we use an ordered set, these should be sorted.
  /*
  if(jointHits.size()>0){
    std::sort(jointHits.begin(), jointHits.end(),
              [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                return a.tid < b.tid;
              } );
  }
  */
  return foundHit;
}


//template <typename RapMapIndexT>
bool mergeMap(fastx_parser::ReadSeq& rp, SAHitMap& leftMap,
                       std::vector<QuasiAlignment>& jointHits) {

	for (auto& ph : leftMap) {
		auto tid = ph.first;
		//sort transcripts with respect to
		//positions
		std::sort(ph.second.tqvec.begin(),
			ph.second.tqvec.end(),
			[](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
				return a.pos < b.pos;
			} );

		std::map<int32_t , int> hitCov ;
		std::map<int32_t , int> hitEndPos ;
		auto it = ph.second.tqvec.begin();
		while(it != ph.second.tqvec.end()){
		    int numOfCov = 0;
		    int32_t hitKey = it->pos - it->queryPos ;
		    if(hitCov.count(hitKey) == 0){
			hitCov[hitKey] = it->matchedLen;
			hitEndPos[hitKey] = it->pos + it->matchedLen;
		    }
		    else{
			if(hitEndPos[hitKey] > it->pos){
			    hitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );

			}else{
			    hitCov[hitKey] += it->matchedLen;
			    hitEndPos[hitKey] = it->pos + it->matchedLen ;
			}
		    }
		    ++it ;
		}
		//std::cout << "\n outside while \n";

		auto minPosIt = std::min_element(ph.second.tqvec.begin(),
				ph.second.tqvec.end(),
				[](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
				    return a.pos < b.pos;
				});

		using pair_type = decltype(hitCov)::value_type;
		auto maxCov = std::max_element(hitCov.begin(),hitCov.end(),
			[](const pair_type& p1, const pair_type& p2) -> bool {
				return p1.second < p2.second ;
			});

		bool hitRC = ph.second.tqvec[0].queryRC;
		int32_t hitPos = maxCov->first;
		bool isFwd = !hitRC;
		jointHits.emplace_back(tid, hitPos, isFwd, rp.seq.length(), 0,0,true);

                auto& qaln = jointHits.back();
                // qaln.fwd = !leftHitCov[vp.first].rc;
                qaln.mateStatus = rapmap::utils::MateStatus::SINGLE_END;
                if (maxCov->second== rp.seq.length()) {
                  qaln.toAlign = true;
                  qaln.editD = 0;
                }
	}
	if(jointHits.size()>0)
		return true;
	else
		return false;
}



template <typename RapMapIndexT>
bool mergeLeftRightSAInts(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>&
        rightFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, RapMapIndexT& rmi, bool maxNumHits,
    bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_) {

  using OffsetT = typename RapMapIndexT::IndexType;
  AlignerEngine ae_;

  auto& SA = rmi.SA;
  auto& txpStarts = rmi.txpOffsets;
  auto& txpIDs = rmi.positionIDs;

  //uint32_t maxInsertSize_{1000};

  bool foundHit{false};

  SAHitMap leftFwdMap;
  SAHitMap leftRcMap;
  SAHitMap rightFwdMap;
  SAHitMap rightRcMap;

  if (leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0) {
    if (leftFwdSAInts.size() > 1) {
      leftFwdMap = rapmap::hit_manager::unionSAHits(
          leftFwdSAInts, rmi, rpair.first.seq.length(), consistentHits);
    } else if (leftFwdSAInts.size() > 0) {
      auto inHit = leftFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        leftFwdMap[tid].active = true;
        leftFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }

    if (leftRcSAInts.size() > 1) {
      leftRcMap = rapmap::hit_manager::unionSAHits(
          leftRcSAInts, rmi, rpair.first.seq.length(), consistentHits);
    } else if (leftRcSAInts.size() > 0) {
      auto inHit = leftRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                          inHit.lcpLength, inHit.queryRC);
        leftRcMap[tid].active = true;
        leftRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  if (rightFwdSAInts.size() > 0 or rightRcSAInts.size() > 0) {
    if (rightFwdSAInts.size() > 1) {
      rightFwdMap = rapmap::hit_manager::unionSAHits(
          rightFwdSAInts, rmi, rpair.second.seq.length(), consistentHits);
    } else if (rightFwdSAInts.size() > 0) {
      auto inHit = rightFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        rightFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                            inHit.lcpLength, inHit.queryRC);
        rightFwdMap[tid].active = true;
        rightFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
    if (rightRcSAInts.size() > 1) {
      rightRcMap = rapmap::hit_manager::unionSAHits(
          rightRcSAInts, rmi, rpair.second.seq.length(), consistentHits);
    } else if (rightRcSAInts.size() > 0) {
      auto inHit = rightRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        rightRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        rightRcMap[tid].active = true;
        rightRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  //jointHits in two different orientations
  //std::vector<QuasiAlignment> fwdRcHits;
  //std::vector<QuasiAlignment> rcFwdHits;

  // there are two possibilities
  bool fwdRc{false};
  bool rcFwd{false};
  // 1. The left is from forward and right is from reverse
  if (leftFwdMap.size() > 0 and rightRcMap.size() > 0) {
    // only consider transcripts that are common between both
    fwdRc = mergeLeftRightMap(rpair, leftFwdMap, rightRcMap, jointHits,//fwdRcHits,
                              editDistance, rmi,maxInsertSize_);
    foundHit = true;
  }
  size_t fwdRcOffset = jointHits.size();

  if (leftRcMap.size() > 0 and rightFwdMap.size() > 0) {
    rcFwd = mergeLeftRightMap(rpair, leftRcMap, rightFwdMap, jointHits,//rcFwdHits,
                              editDistance, rmi, maxInsertSize_);
    foundHit = true;
  }

  // merge two sets if jointHits in different orientations instead of sort
  // use std::merge for this
  std::inplace_merge(jointHits.begin(), jointHits.begin() + fwdRcOffset,
                     jointHits.end(),
                     [](const QuasiAlignment& l, const QuasiAlignment& r)->bool {
                       return l.tid < r.tid;
                     });


  //manage orphans here
  bool orphanRes{false} ;

  bool flag = (rpair.first.name == "10000007_1_4452_4695_155/1") ;

  if(flag)std::cerr << "\n Status of foundHit "<<leftFwdMap.size() << "\t" << leftRcMap.size() << "\t"
            << rightFwdMap.size() << "\t" << rightRcMap.size() << "\n" ;

  if(!foundHit){


    if(leftFwdMap.size() != 0 or rightRcMap.size() != 0){
      //rescue right end
      if(leftFwdMap.size() > 0){
        //code here
        orphanRes = rescueOrphan(rpair, leftFwdMap,
                                rapmap::hit_manager_sa::MateSAStatus::LEFT_END_FWD,
                                 jointHits,
                                rmi) ;
      }else{
        //code here
        orphanRes = rescueOrphan(rpair,
                                rightRcMap,
                                rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_RC,
                                jointHits,
                                rmi) ;
      }
    }

    if(leftRcMap.size() != 0 or rightFwdMap.size() != 0){
      //rescue right end
      if(leftRcMap.size() > 0){
        //code here
        orphanRes = rescueOrphan(rpair, leftRcMap,
                                rapmap::hit_manager_sa::MateSAStatus::LEFT_END_RC,
                                 jointHits,
                                rmi) ;
      }else{
        //code here
        orphanRes = rescueOrphan(rpair, rightFwdMap,
                                rapmap::hit_manager_sa::MateSAStatus::RIGHT_END_FWD,
                                 jointHits,
                                rmi) ;
      }
    }
  }

  if(orphanRes){
    foundHit = true ;
  }


  /*
  uint32_t i = 0;
  uint32_t j = 0;
  while(i<fwdRcHits.size() and j<rcFwdHits.size()){
    if(fwdRcHits[i].tid < rcFwdHits[j].tid){
      jointHits.push_back(fwdRcHits[i]);
      i++;
    } else {
      jointHits.push_back(rcFwdHits[j]);
      j++;
    }
  }
  while (i<fwdRcHits.size()) {
    jointHits.push_back(fwdRcHits[i]);
    i++;
  }
  while (j<rcFwdHits.size()) {
    jointHits.push_back(rcFwdHits[j]);
    j++;
  }
  */

  /*if(fwdRc and rcFwd and jointHits.size()>1){
    std::sort(jointHits.begin(), jointHits.end(),
              [](const QuasiAlignment& a, const QuasiAlignment& b)-> bool {
                return a.tid < b.tid;
              } );
  }*/
  return foundHit;
}


template <typename RapMapIndexT>
bool mergeSAInts(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
    std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, RapMapIndexT& rmi, bool maxNumHits,
    bool consistentHits, rapmap::utils::HitCounters& hctr) {

  //if(leftFwdSAInts.size()==0 and leftRcSAInts.size()==0)
  //	std::cout<<"err\n";

  using OffsetT = typename RapMapIndexT::IndexType;
  AlignerEngine ae_;

  auto& SA = rmi.SA;
  auto& txpStarts = rmi.txpOffsets;
  auto& txpIDs = rmi.positionIDs;

  //uint32_t maxInsertSize_{1000};

  bool foundHit{false};

  SAHitMap leftFwdMap;
  SAHitMap leftRcMap;

  if (leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0) {
    if (leftFwdSAInts.size() > 1) {
      leftFwdMap = rapmap::hit_manager::unionSAHits(
          leftFwdSAInts, rmi, rp.seq.length(), consistentHits);
    } else if (leftFwdSAInts.size() > 0) {
      auto inHit = leftFwdSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                           inHit.lcpLength, inHit.queryRC);
        leftFwdMap[tid].active = true;
        leftFwdMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }

    if (leftRcSAInts.size() > 1) {
      leftRcMap = rapmap::hit_manager::unionSAHits(
          leftRcSAInts, rmi, rp.seq.length(), consistentHits);
    } else if (leftRcSAInts.size() > 0) {
      auto inHit = leftRcSAInts.front();
      for (auto i = inHit.begin; i < inHit.end; ++i) {
        auto globalPos = SA[i];
        auto tid = rmi.transcriptAtPosition(globalPos);
        auto txpPos = globalPos - txpStarts[tid];
        leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos,
                                          inHit.lcpLength, inHit.queryRC);
        leftRcMap[tid].active = true;
        leftRcMap[tid].tqvec.back().matchedLen = inHit.len;
      }
    }
  }

  bool fwdRc{false};
  bool rcFwd{false};

  if (leftFwdMap.size() > 0) {
    // only consider transcripts that are common between both
    fwdRc = mergeMap(rp, leftFwdMap, jointHits);
    foundHit = true;
  }

  size_t fwdRcOffset = jointHits.size();

  if (leftRcMap.size() > 0) {
    rcFwd = mergeMap(rp, leftRcMap, jointHits);
    foundHit = true;
  }

  // merge two sets if jointHits in different orientations instead of sort
  // use std::merge for this
  std::inplace_merge(jointHits.begin(), jointHits.begin() + fwdRcOffset,
                     jointHits.end(),
                     [](const QuasiAlignment& l, const QuasiAlignment& r)->bool {
                       return l.tid < r.tid;
                     });

  return foundHit;

}



using SAIndex32BitDense =
    RapMapSAIndex<int32_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int32_t>,
                                    rapmap::utils::KmerKeyHasher>>;
using SAIndex64BitDense =
    RapMapSAIndex<int64_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int64_t>,
                                    rapmap::utils::KmerKeyHasher>>;
using SAIndex32BitPerfect =
    RapMapSAIndex<int32_t,
                  PerfectHashT<uint64_t, rapmap::utils::kmerVal<int32_t>>>;
using SAIndex64BitPerfect =
    RapMapSAIndex<int64_t,
                  PerfectHashT<uint64_t, rapmap::utils::kmerVal<int64_t>>>;

template bool mergeLeftRightSAInts<SAIndex32BitDense>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex64BitDense>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex32BitPerfect>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);

template bool mergeLeftRightSAInts<SAIndex64BitPerfect>(
    fastx_parser::ReadPair& rpair, bool lhp, bool rhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr,
    uint32_t editDistance, uint32_t maxInsertSize_);



template bool mergeSAInts<SAIndex32BitDense>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex64BitDense>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitDense& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex32BitPerfect>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex32BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);

template bool mergeSAInts<SAIndex64BitPerfect>(
    fastx_parser::ReadSeq& rp, bool lhp,
    std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
    std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
    std::vector<QuasiAlignment>& jointHits, SAIndex64BitPerfect& rmi,
    bool maxNumHits, bool consistentHits, rapmap::utils::HitCounters& hctr);







template bool mergeLeftRightMap<SAIndex32BitDense>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex32BitDense& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex64BitDense>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex64BitDense& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex32BitPerfect>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex32BitPerfect& rmi, uint32_t maxInsertSize_);

template bool mergeLeftRightMap<SAIndex64BitPerfect>(
    fastx_parser::ReadPair& rpair, SAHitMap& leftMap, SAHitMap& rightMap,
    std::vector<QuasiAlignment>& jointHits, uint32_t editDistance,
    SAIndex64BitPerfect& rmi, uint32_t maxInsertSize_);




  template bool rescueOrphan<SAIndex32BitDense>(
                                                fastx_parser::ReadPair& rpair,
                                                SAHitMap& alignedMap,
                                                rapmap::hit_manager_sa::MateSAStatus str,
                                                std::vector<QuasiAlignment>& jointHits,
                                                SAIndex32BitDense& rmi) ;

  template bool rescueOrphan<SAIndex64BitDense>(
                                                     fastx_parser::ReadPair& rpair,
                                                     SAHitMap& alignedMap,
                                                     rapmap::hit_manager_sa::MateSAStatus str,
                                                     std::vector<QuasiAlignment>& jointHits,
                                                     SAIndex64BitDense& rmi) ;

  template bool rescueOrphan<SAIndex32BitPerfect>(
                                                       fastx_parser::ReadPair& rpair,
                                                       SAHitMap& alignedMap,
                                                       rapmap::hit_manager_sa::MateSAStatus str,
                                                       std::vector<QuasiAlignment>& jointHits,
                                                       SAIndex32BitPerfect& rmi) ;

  template bool rescueOrphan<SAIndex64BitPerfect>(
                                                       fastx_parser::ReadPair& rpair,
                                                       SAHitMap& alignedMap,
                                                       rapmap::hit_manager_sa::MateSAStatus str,
                                                       std::vector<QuasiAlignment>& jointHits,
                                                       SAIndex64BitPerfect& rmi) ;






}
}
