#ifndef __SELECTIVE_ALIGNMENT_HPP__
#define __SELECTIVE_ALIGNMENT_HPP__

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"
#include "EditDistance.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>






template <typename RapMapIndexT> class SECollector{
public:
  using OffsetT = typename RapMapIndexT::IndexType;

  /** Construct an SECollector given an index **/
    SECollector(RapMapIndexT* rmi)
        : rmi_(rmi) {}




   template <typename ReadStructT>
       bool operator()(ReadStructT& readT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits
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

          uint8_t editThreshold = readLen/2 ;
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
                  OffsetT globalPos = pos + txpStarts[txpID];
                  OffsetT thisTxpLen = txpLens[txpID];

                  // start and end would
                  // ideally be -5 and +5
                  // of the start position

                  // If the start is before start
                  // of transcript
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

                  //error checking
                  if(globalPos > concatText.length()){
                      std::cout << "What is wrong ! global position " << globalPos << " bigger than " << concatText.length() << " Length of transcript " << thisTxpLen  << "\n" ;
                  }

                  auto thisTxpSeq = concatText.substr(globalPos, endOffset);
		  //compute edit distance
		  startEditDistance = edit_distance(read, thisTxpSeq, 50) ;


                  if(readName == "SRR1265495.29995"){
                    std::cout << "\nRead Header: "<<readName ;
                    std::cout << "\nLCP Length: "<< (int)lcpLength ;
                    std::cout << "\nI am at the beginning of hit list, the edit distance is: " << startEditDistance << "\n" ;
                  }

                  startHit.editD = startEditDistance;
                  startHit.toAlign = true;
                  startHit.cigar = "II";

	  }


          if(hits.size() > 1){
            if (lcpLength > readLen/2){
                      if(startEditDistance < editThreshold){
                          for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                                      //selectedHits.emplace_back(hitsIt->tid,hitsIt->pos,hitsIt->fwd,hitsItreadLen->readLen,startEditDistance,"II");
                                      hitsIt->editD = startEditDistance;
                                      hitsIt->toAlign = true;
                                      hitsIt->cigar = "II";
                              }
                      }
              }else{

                  for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                              uint32_t txpID = hitsIt->tid ;
                              int32_t pos = hitsIt->pos;
                              int32_t startOffset, endOffset ;
                              OffsetT globalPos = pos + txpStarts[txpID];
                              OffsetT thisTxpLen = txpLens[txpID];
                      //error checkingo
                      //
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

                              if(globalPos > concatText.length()){
                                  std::cout << " What is wrong ! global position " << globalPos << " bigger than " << concatText.length() << " for tid: " << txpID  << "\n" ;
                                  std::cout << " transcript start position  "<< txpStarts[txpID] << " transcript Length: " << thisTxpLen << "\n" ;
                                  std::cout << " Number of transcripts " << txpLens.size() << "\n" ;
                              }

                              auto thisTxpSeq = concatText.substr(globalPos, endOffset);
                              //compute edit distance


                              auto thisEditDistance = edit_distance(read, thisTxpSeq, 50) ;

                             if(readName == "SRR1265495.29995"){
                                std::cout << "\nRead Header: "<<readName ;
                                std::cout << "\nLCP Length: "<< (int)lcpLength ;
                                std::cout << "\nI am at the middle of hit list, the edit distance is: " << thisEditDistance << "\n" ;
                              }

                              if(thisEditDistance < editThreshold){
                                      //selectedHits.emplace_back(txpID,pos,startHit.fwd,hitsIt->readLen, thisEditDistance,"II");
                                      hitsIt->editD = thisEditDistance;
                                      hitsIt->toAlign = true;
                                      hitsIt->cigar = "II";
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
