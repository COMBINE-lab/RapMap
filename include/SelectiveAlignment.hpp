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





   void operator()(std::string& read,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  SASearcher<RapMapIndexT>& saSearcher,
                  rapmap::utils::MateStatus mateStatus,
                  bool consistentHits = false) {

	  //go through the lcp length of each SA Interval for each hit
	  //figure out if it needs to be changed
	  auto& txpStarts = rmi_->txpOffsets ;
	  auto& SA = rmi_->SA;
	  auto& concatText = rmi_->seq;
	  auto& txpLens = rmi_->txpLens ;




	  auto& startHit = hits.front();
	  auto lcpLength = startHit.lcpLength ;
	  auto readLen = read.length();

          uint8_t editThreshold = readLen/2 ;
	  int startEditDistance = 0;

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
		  auto& txpID = startHit.tid ;
		  auto pos = startHit.pos;
		  auto globalPos = pos + txpStarts[txpID];
		  auto thisTxpLen = txpLens[txpID];

                  //error checking
                  if(thisTxpLen < pos){
                      std::cout << "What is wrong ! " << thisTxpLen << " " << pos << "\n" ;
                  }

		  auto thisTxpSeq = concatText.substr(globalPos,thisTxpLen-pos);
		  //compute edit distance
		  startEditDistance = edit_distance(read,thisTxpSeq,1000) ;
                  startHit.editD = startEditDistance;
                  startHit.toAlign = true;
                  startHit.cigar = "II";

	  }

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
			  auto& txpID = hitsIt->tid ;
			  auto pos = hitsIt->pos;
			  auto globalPos = pos + txpStarts[txpID];
			  auto thisTxpLen = txpLens[txpID];
                  //error checking
                          if(thisTxpLen < pos){
                              std::cout << "What is wrong ! " << thisTxpLen << " " << pos << "\n" ;
                          }
			  auto thisTxpSeq = concatText.substr(globalPos,thisTxpLen-pos);
			  //compute edit distance


			  auto thisEditDistance = edit_distance(read,thisTxpSeq,1000) ;

			  if(thisEditDistance < editThreshold){
				  //selectedHits.emplace_back(txpID,pos,startHit.fwd,hitsIt->readLen, thisEditDistance,"II");
	  			  hitsIt->editD = thisEditDistance;
	  			  hitsIt->toAlign = true;
	  			  hitsIt->cigar = "II";
			  }
			}

	  }



  }
private:
  RapMapIndexT* rmi_ ;

};

#endif
