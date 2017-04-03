//jello

#ifndef __SELECTIVE_ALIGNMENT_HPP__
#define __SELECTIVE_ALIGNMENT_HPP__

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"
#include "EditDistance.hpp"
#include "edlib.h"

#include <algorithm>
#include <iostream>
#include <iterator>

//#include "dtl/dtl.hpp"

static constexpr int8_t rc_table[128] = {
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 15
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 31
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 787
      78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 63
      78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 79
      78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 95
      78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 101
      78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // 127
  };

/*
template <typename OffsetT>
int32_t SECollector::hammingDist(QuasiAlignment& qa, std::string& read, std::string& txpSeq,  OffsetT maxOffset, int maxDist){
       int32_t hamming = 0;
       auto minOffset = (qa.pos < 0) ? -(qa.pos) : 0;
       hamming += minOffset;
       //int maxOffset = std::min(static_cast<int>(read.size()), static_cast<int>(transcript.RefLength - qa.pos));
       hamming += (read.size() - maxOffset);
       //const char* txpSeq = transcript.Sequence();
       //const char* readSeq{nullptr};
       auto readLen = read.length();
       if (qa.fwd) {
       for (int j=minOffset; j<maxOffset; ++j) {
            hamming += (read[j]!=txpSeq[qa.pos+j]);
                if (hamming > maxDist) { return hamming;  }
            }
        } else {
            for (int j=minOffset; j<maxOffset; ++j) {
             hamming += ((char)rc_table[(uint8_t)read[(readLen-j)-1]] != txpSeq[qa.pos+j]);
                if (hamming > maxDist) { return hamming;  }
            }
        }
        return hamming;

}
*/



template <typename RapMapIndexT> class SECollector{
public:
  using OffsetT = typename RapMapIndexT::IndexType;

  /** Construct an SECollector given an index **/
    SECollector(RapMapIndexT* rmi)
        : rmi_(rmi) {}

    //int32_t hammingDist(QuasiAlignment& q, std::string& read, std::string& seq,  Offset trancriptLen, int maxDist);


   template <typename ReadStructT>
       bool operator()(ReadStructT& readT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  int32_t editThreshold
                  ) {

	  //go through the lcp length of each SA Interval for each hit
	  //figure out if it needs to be changed
          //
           if(hits.size() == 0){
               return false ;
           }

          auto& read = readT.seq;
	  auto& txpStarts = rmi_->txpOffsets ;
	  auto& SA = rmi_->SA;
	  auto& concatText = rmi_->seq;
	  auto& txpLens = rmi_->txpLens ;




	  auto& startHit = hits.front();
	  auto lcpLength = startHit.lcpLength ;
	  auto readLen = read.length();

          //uint8_t editThreshold = readLen/2 ;
	  int startEditDistance = 0;

            auto& readName = readT.name;
            // If the read name contains multiple space-separated parts,
            // print only the first
            size_t splitPos = readName.find(' ');
            if (splitPos < readName.length()) {
                readName[splitPos] = '\0';
            } else {
                splitPos = readName.length();
            }

            // trim /1 from the pe read
            if (splitPos > 2 and readName[splitPos - 2] == '/') {
                readName[splitPos - 2] = '\0';
            }

	  // There are two possible ways to go
	  // from here
	  // We have a large lcp so we will calculate alignment only once
	  // and gonna copy for rest of the things
	  // Nevertheless we have to align the first read to the transcript

	  // Align the startHit any way
	  // get transcript id
	  // sequence and read sequence
	  // and position
	  {
                  uint32_t txpID = startHit.tid ;
                  int32_t pos = startHit.pos;
                  int32_t startOffset, endOffset ;
                  OffsetT globalPos = txpStarts[txpID];
                  OffsetT thisTxpLen = txpLens[txpID];

                  // start and end would
                  // ideally be -5 and +5
                  // of the start position

                  // If the start is before start
                  // of transcript
                  /*
                  if(pos >= 5){
                      globalPos = globalPos - 5;
                  }else{
                      globalPos = txpStarts[txpID];
                  }

                  if(pos + readLen + 5 <= thisTxpLen){
                      endOffset = readLen + 5;
                  }else{
                      endOffset = thisTxpLen - pos ;
                  }
                  */

                  auto overHangLeft = (pos < 0)?-(pos):0;
                  auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                  globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;
                  auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                  extend = (overHangRight > 0)?(extend-overHangRight):extend;

                  auto thisTxpSeq = concatText.substr(globalPos, extend);



                  EdlibAlignResult startEdlibResult;
                  if(startHit.fwd){
                    startEdlibResult = edlibAlign(read.c_str(), read.length(), thisTxpSeq.c_str(), thisTxpSeq.length(), edlibNewAlignConfig(editThreshold, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                  }else{
                    auto revRead = rapmap::utils::reverseComplement(read);
                    startEdlibResult = edlibAlign(revRead.c_str(), read.length(), thisTxpSeq.c_str(), thisTxpSeq.length(), edlibNewAlignConfig(editThreshold, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                  }


                /*
                  if(readName == "SRR1265495.29995"){
                    std::cout << "\nRead Header: "<<readName ;
                    std::cout << "\nLCP Length: "<< (int)lcpLength ;
                    std::cout << "\nI am at the beginning of hit list, the edit distance is: " << startEditDistance << "\n" ;
                  }
                */

                  startHit.editD = startEdlibResult.editDistance;
                  startEditDistance = startEdlibResult.editDistance;

                  if(startEditDistance != -1){
                    startHit.toAlign = true;
                    char* cigar_ = edlibAlignmentToCigar(startEdlibResult.alignment, startEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                    std::string cigar(cigar_);
                    startHit.cigar = cigar ;
                 }else{
                    startHit.toAlign = false ;

                 }
          }


          if(hits.size() > 1){
            if ( false &&  lcpLength >= readLen){
                      if(startEditDistance != -1){
                          for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                                      //selectedHits.emplace_back(hitsIt->tid,hitsIt->pos,hitsIt->fwd,hitsItreadLen->readLen,startEditDistance,"II");
                                      hitsIt->editD = startEditDistance;
                                      hitsIt->toAlign = true;
                                      hitsIt->cigar = startHit.cigar;
                              }
                      }else{
                          for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                                    hitsIt->toAlign = false ;
                          }
                      }
              }else{

                  for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                              uint32_t txpID = hitsIt->tid ;
                              int32_t pos = hitsIt->pos;
                              int32_t startOffset, endOffset ;
                              OffsetT globalPos = txpStarts[txpID];
                              OffsetT thisTxpLen = txpLens[txpID];

                              /*
                              if(pos >= 5){
                                  globalPos = globalPos - 5;
                              }else{
                                  globalPos = txpStarts[txpID];
                              }

                              if(pos + readLen + 5 <= thisTxpLen){
                                  endOffset = readLen + 5;
                              }else{
                                  endOffset = thisTxpLen - pos ;
                              }
                              */

                              auto overHangLeft = (pos < 0)?-(pos):0;
                              auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                              globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;
                              auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                              extend = (overHangRight > 0)?(extend-overHangRight):extend;

                              auto thisTxpSeq = concatText.substr(globalPos, extend);

                              //auto thisEditDistance = edit_distance(read, thisTxpSeq, 50) ;
                              EdlibAlignResult thisEdlibResult;

                              if(hitsIt->fwd){

                                thisEdlibResult = edlibAlign(read.c_str(), read.length(), thisTxpSeq.c_str(), thisTxpSeq.length(), edlibNewAlignConfig( editThreshold, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                              }else{
                                auto revRead = rapmap::utils::reverseComplement(read);
                                thisEdlibResult = edlibAlign(revRead.c_str(), read.length(), thisTxpSeq.c_str(), thisTxpSeq.length(), edlibNewAlignConfig(editThreshold, EDLIB_MODE_NW, EDLIB_TASK_PATH));
                              }
                              auto thisEditDistance = thisEdlibResult.editDistance ;


                              if(thisEditDistance != -1){
                                      //selectedHits.emplace_back(txpID,pos,startHit.fwd,hitsIt->readLen, thisEditDistance,"II");
                                      hitsIt->editD = thisEditDistance;
                                      char* cigar_ = edlibAlignmentToCigar(thisEdlibResult.alignment, thisEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                                      std::string cigar(cigar_);
                                      hitsIt->toAlign = true;
                                      hitsIt->cigar = cigar;
                              }
                            }

            }
        }

          return true ;
   }
private:
  RapMapIndexT* rmi_ ;

};

#endif
