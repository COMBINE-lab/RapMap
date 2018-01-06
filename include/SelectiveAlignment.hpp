//jello

#ifndef __SELECTIVE_ALIGNMENT_HPP__
#define __SELECTIVE_ALIGNMENT_HPP__

#include "RapMapSAIndex.hpp"
#include "RapMapUtils.hpp"
#include "SASearcher.hpp"
#include "EditDistance.hpp"
#include "edlib.h"
#include "sparsepp/spp.h"
#include "string_view.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>

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
  std::atomic<uint64_t> numCacheHits{0};

  /** Construct an SECollector given an index **/
    SECollector(RapMapIndexT* rmi)
        : rmi_(rmi) {}

   /* 
    ~SECollector() {
        if (numCacheHits > 0) { std::cerr << "\n\n\n numCacheHits : " << numCacheHits << "\n\n\n"; }
        else {
            std::cerr  << "\n\n in SECollector destructor --- no cached hits!\n\n";
        }
    }
    */

    //int32_t hammingDist(QuasiAlignment& q, std::string& read, std::string& seq,  Offset trancriptLen, int maxDist);
    class SubAlignmentKey {
        public:
        //SubAlignmentKey(const SubAlignmentKey&& o) = default;
        //SubAlignmentKey& operator=(const SubAlignmentKey&& o) = default;
        uint32_t qstart;
        uint32_t qlen;
        bool fwd;
        uint64_t refhash;
        std::size_t getHash() const { return XXH64(reinterpret_cast<const void*>(this), sizeof(this), 314); }
        bool operator==(const SubAlignmentKey& o) const {
            return o.qstart == qstart and o.qlen == qlen and o.fwd == fwd and o.refhash == refhash;
        }
    };

    class SubKeyHasher {
        public:
        size_t operator()(const SubAlignmentKey& k) const { return k.getHash(); }
    };

   template <typename ReadStructT>
       bool compute(ReadStructT& readT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  int32_t editThreshold
                  ) {
           using QuasiAlignment = rapmap::utils::QuasiAlignment;
           using MateStatus = rapmap::utils::MateStatus;
           using kmerVal = rapmap::utils::kmerVal<OffsetT>;
	  //go through the lcp length of each SA Interval for each hit
	  //figure out if it needs to be changed


           if(hits.size() == 0){
               return false ;
           }

          //Mapping information
          //Loaed
          //std::cout << "\n HERE \n" ;
          spp::sparse_hash_set<int32_t> tidset ;

	  auto& txpStarts = rmi_->txpOffsets ;
	  auto& SA = rmi_->SA;
	  auto& concatText = rmi_->seq;
      stx::string_view concatTextView(concatText);
	  auto& txpLens = rmi_->txpLens ;
          auto& khash = rmi_->khash ;

          auto k = rapmap::utils::my_mer::k();


          auto& read = readT.seq;
          stx::string_view readView(read);

          rapmap::utils::my_mer mer;

          //std::vector<SAIntervalHit> fwdSAInts ;
          //std::vector<SAIntervalHit> rcSAInts ;
          kmerVal tSAInt ;


          //start hit information
	  auto& startHit = hits.front();
	  auto lcpLength = startHit.lcpLength ;
	  auto readLen = readView.length();

          //uint8_t editThreshold = readLen/2 ;
	  int32_t startEditDistance = 0;
	  uint32_t eDistScale = 2;

          bool skipLCPOpt{false};
          bool currentValid{false};
          stx::string_view currentKmer;

        //TODO check if lcp is really returning same sequences
          //stx::string_view firsttidString ;


	  //debug
	  //std::string testRead = "xGCGAGAGTCTGTCTCAAAAACAACAACCACCACCACCACCACCACAACAACAACGAAACTATCTGCAAACTTTGC";

    /*
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
    */

	  // There are two possible ways to go
	  // from here
	  // We have a large lcp so we will calculate alignment only once
	  // and gonna copy for rest of the things
	  // Nevertheless we have to align the first read to the transcript

	  // Align the startHit any way
	  // get transcript id
	  // sequence and read sequence
	  // and position
      bool multiHit = hits.size() > 1;
      
      spp::sparse_hash_map<SubAlignmentKey, int, SubKeyHasher> edmap;
            auto cacheResult = [&edmap](uint32_t gapStart, uint32_t refGapLen, bool readFW, const char* refSeq, int edist) -> void  {
                auto seqhash = XXH64(reinterpret_cast<const void*>(refSeq), refGapLen, 314);
                SubAlignmentKey k{gapStart, refGapLen, readFW, seqhash};
                edmap[k] = edist;
            };

          if(!startHit.toAlign){
                  uint32_t txpID = startHit.tid ;
                  int32_t pos = startHit.pos;
                  int32_t startOffset, endOffset ;
                  OffsetT globalPos = txpStarts[txpID];
                  OffsetT thisTxpLen = txpLens[txpID];

                  /*auto overHangLeft = (pos < 0)?-(pos):0;
                  auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                  globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;
                  auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                  extend = (overHangRight > 0)?(extend-overHangRight):extend;*/

                  /* ROB: get rid of the copy */
                  //auto thisTxpSeq = concatText.substr(globalPos, extend);
                  //const char* thisTxpSeq = concatText.data() + globalPos;
                  //int thisTargetLen = extend;


                  /* ROB : slight interface change */
                  //EdlibAlignResult startEdlibResult;

                  /*if(startHit.fwd){
                    ae_(read.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                  }else{
                    auto revRead = rapmap::utils::reverseComplement(read);
                    ae_(revRead.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                  }
                  auto& startEdlibResult = ae_.result();

                  startHit.editD = startEdlibResult.editDistance;

                  int32_t startEditDistance = startEdlibResult.editDistance;*/

                  //@debug purpose
                  //auto readName = rapmap::utils::getReadName(readT) ;
		  startEditDistance  = 0;
                  startHit.editD = 0; ;
		  for(int i=0;i<startHit.gapsBegin.size();++i){
                    auto overHangLeft = (pos+startHit.gapsBegin[i] < 0)?-(pos+startHit.gapsBegin[i]):0;
		    uint32_t gapLen = startHit.gapsEnd[i] - startHit.gapsBegin[i];
                    auto overHangRight = (pos+startHit.gapsEnd[i] > thisTxpLen)?(pos+startHit.gapsEnd[i]-thisTxpLen):0;

                    globalPos = (overHangLeft == 0)?(pos+startHit.gapsBegin[i]+txpStarts[txpID]):txpStarts[txpID];
                    auto extend = (overHangLeft > 0)?(gapLen - overHangLeft):gapLen ;
                    extend = (overHangRight > 0)?(extend-overHangRight):extend;

		    const char* thisTxpSeq = concatText.data() + globalPos;
		    uint32_t thisTargetLen = extend;
		    if(startHit.fwd){
                      ae_(readView.substr(startHit.gapsBegin[i],gapLen).data(), gapLen, thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
		    } else  {
		      auto revRead = rapmap::utils::reverseComplement(read);
              stx::string_view revReadView(revRead);
                      ae_(revReadView.substr(startHit.gapsBegin[i],gapLen).data(), gapLen, thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
		    }
		    
                    auto& startEdlibResult = ae_.result();


		    if(startEdlibResult.editDistance != -1 and startEditDistance!=-1) {
                      startEditDistance += startEdlibResult.editDistance;
                      startHit.editD += startEdlibResult.editDistance;
                      if (multiHit) { cacheResult(startHit.gapsBegin[i], thisTargetLen, startHit.fwd, thisTxpSeq, startEdlibResult.editDistance); }
	 	    } else {
		      startEditDistance = -1;
		      startHit.editD = -1;
		    }
	
		    

		  }


		  //debug
		  /*if (read == testRead and startEditDistance1!=-1 and startEditDistance != startEditDistance1){
		    for(int i=0;i<startHit.gapsBegin.size();++i){
		      std::cout<<startEditDistance1 <<"\t" << startEditDistance << "\t" <<startHit.gapsBegin[i]<<"\t"<<startHit.gapsEnd[i]<<"\n"<<startHit.fwd<<"\t"<<read<<"\t"<<globalPos<<"\n";
                      auto txpID = rmi_->transcriptAtPosition(globalPos);auto& txpNames = rmi_->txpNames;
                      std::cout << rmi_->txpNames[txpID] << "\n";
		    }
		  }*/
                  /*
                  if(readName == "15006_1_6022_314_181"){
                      std::cout << "\nIn selective Alignment start  edit distance: "<< startEditDistance << " length: "<< thisTxpLen << " pos " << pos << " fwd "<<  startHit.fwd << "\n";
                  }*/


                  if(startEditDistance <= editThreshold and startEditDistance>=0){
                    startHit.toAlign = true;
                    /* ROB: No CIGAR right now */
                    /*
                    char* cigar_ = edlibAlignmentToCigar(startEdlibResult.alignment, startEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                    std::string cigar(cigar_);
                    startHit.cigar = cigar ;
                    */
                 }else{
                    startHit.toAlign = false ;
                 }

                /* Take the kmer from the transcript */
                currentKmer = concatTextView.substr(globalPos, k);
                //firsttidString = concatTextView.substr(globalPos, 75);


                if(currentKmer.length() == k and currentKmer.find_first_of('$') == std::string::npos){
                    mer = currentKmer;
                    currentValid= true ;
                }else{
                    currentValid = false ;
                }
          }

	  spp::sparse_hash_map<int,OffsetT> tidPos; 
       
          if(currentValid){
              auto bits = mer.word(0);
              auto hashIt = khash.find(bits);
              if(hashIt != khash.end()){
                  //std::cout << "\n Found \n";
                  tSAInt = hashIt->second;
                  lcpLength = hashIt->second.lcpLength ;
                  for(OffsetT i = tSAInt.interval.begin() ; i != tSAInt.interval.end(); ++i ){
                      auto iGlobalPos = SA[i];
                      auto txpID = rmi_->transcriptAtPosition(iGlobalPos);
                      tidset.insert(txpID);
                      tidPos[txpID] = iGlobalPos;
                }
              }else{
                  skipLCPOpt = true ;
              }

          }else{
              skipLCPOpt = true ;
          }




          if(multiHit){
            for(auto hitsIt= hits.begin()+1 ; hitsIt != hits.end() ; ++hitsIt){
                uint32_t txpID = hitsIt->tid ;
                auto search = tidset.find(txpID);
                int32_t pos = hitsIt->pos;
                int32_t startOffset, endOffset ;
                OffsetT globalPos = txpStarts[txpID];
                OffsetT thisTxpLen = txpLens[txpID];
                //auto overHangLeft = (pos < 0)?-(pos):0;
                //auto overHangRight = (pos+readLen > thisTxpLen)?(pos+readLen-thisTxpLen):0;

                //globalPos = (overHangLeft == 0)?(pos+globalPos):globalPos;

                if(hitsIt->toAlign){
                    continue;
                }
                if ( (!skipLCPOpt) && (search != tidset.end()) && (lcpLength >= readLen) && (tidPos[txpID] == globalPos)){
                    if(startEditDistance <= editThreshold and startEditDistance>=0){
                                hitsIt->editD = startEditDistance;
                                hitsIt->toAlign = true;
                                //std::cout << " LCP Length " << lcpLength << "\n";

                                /* ROB: No CIGAR right now */
                                //hitsIt->cigar = startHit.cigar;
                          }else{
                               hitsIt->editD = startEditDistance;
                                hitsIt->toAlign = false ;
                          }
                }else{
                      //auto extend = (overHangLeft > 0)?(readLen - overHangLeft):readLen ;
                      //extend = (overHangRight > 0)?(extend-overHangRight):extend;

                      /* ROB: get rid of the copy */
                      //auto thisTxpSeq = concatText.substr(globalPos, extend);
                      //const char* thisTxpSeq = concatText.data() + globalPos;
                      //int thisTargetLen = extend;

                      //auto thisEditDistance = edit_distance(read, thisTxpSeq, 50) ;
                      /* ROB : slight interface change */
                      //EdlibAlignResult thisEdlibResult;
                      /*if(hitsIt->fwd){
                        ae_(read.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold*+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                      }else{
                        auto revRead = rapmap::utils::reverseComplement(read);
                        ae_(revRead.c_str(), read.length(), thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                      }
                      auto& thisEdlibResult = ae_.result();
                      auto thisEditDistance = thisEdlibResult.editDistance;*/

		      int32_t thisEditDistance = 0;
		      for(int i=0;i<hitsIt->gapsBegin.size();++i){
                   
   		        int overHangLeft = (pos+hitsIt->gapsBegin[i] < 0)?-(pos+hitsIt->gapsBegin[i]):0;
		        uint32_t gapLen = hitsIt->gapsEnd[i] - hitsIt->gapsBegin[i];
                        int32_t overHangRight = (pos+hitsIt->gapsEnd[i] > thisTxpLen)?(pos+hitsIt->gapsEnd[i]-thisTxpLen):0;

                        globalPos = (overHangLeft == 0)?(pos+hitsIt->gapsBegin[i]+txpStarts[txpID]):txpStarts[txpID];
                        int32_t extend = (overHangLeft > 0)?(gapLen - overHangLeft):gapLen ;
                        extend = (overHangRight > 0)?(extend-overHangRight):extend;

                        const char* thisTxpSeq = concatText.data() + globalPos;
                        uint32_t thisTargetLen = extend;

                        int edist = -1;
                        bool useCached = false;
		        if(hitsIt->fwd){
                          auto seqhash = XXH64(reinterpret_cast<const void*>(thisTxpSeq), thisTargetLen, 314);
                          SubAlignmentKey k{hitsIt->gapsBegin[i], thisTargetLen, true, seqhash};
                          auto edistIt = edmap.find(k);
                          if (edistIt == edmap.end())  {
                            ae_(readView.substr(hitsIt->gapsBegin[i],gapLen).data(), gapLen, thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                          } else {
                            useCached = true;
                            edist = edistIt->second;
                            ++numCacheHits;
                          }
		        } else  {
		          auto revRead = rapmap::utils::reverseComplement(read);
                  stx::string_view revReadView(revRead);
                          auto seqhash = XXH64(reinterpret_cast<const void*>(thisTxpSeq), thisTargetLen, 314);
                          SubAlignmentKey k{hitsIt->gapsBegin[i], thisTargetLen, false, seqhash};
                          auto edistIt = edmap.find(k);
                         if (edistIt == edmap.end())  {
                           ae_(revReadView.substr(hitsIt->gapsBegin[i],gapLen).data(), gapLen, thisTxpSeq, thisTargetLen, edlibNewAlignConfig((editThreshold+1)*eDistScale, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE));
                         } else {
                           useCached = true;
                           edist = edistIt->second;
                           ++numCacheHits;
                          }
		        }

                        if (!useCached) {
                          auto& thisEdlibResult = ae_.result();
                          edist = thisEdlibResult.editDistance;
                        }


		        if(edist!=-1 and thisEditDistance != -1){
                          thisEditDistance += edist;
                          if (!useCached) {
                              cacheResult(hitsIt->gapsBegin[i], thisTargetLen, hitsIt->fwd, thisTxpSeq, edist); 
                          }
			} else {
			   hitsIt->editD = -1;
			   thisEditDistance = -1;
		        }

		      }

		      //debug
		      /*if(read==testRead and thisEditDistance1!=-1 and thisEditDistance1 != thisEditDistance ){
		        for(int i=0;i<hitsIt->gapsBegin.size();++i){
		          std::cout<<thisEditDistance1 <<"\t" << thisEditDistance << "\t" <<hitsIt->gapsBegin[i]<<"\t"<<hitsIt->gapsEnd[i]<<"\n"<<hitsIt->fwd<<"\t"<<read<<"\t"<<globalPos<<"\n";
                          auto txpID = rmi_->transcriptAtPosition(globalPos);auto& txpNames = rmi_->txpNames;
                          std::cout << rmi_->txpNames[txpID] << "\n";
		        }
		      }*/     


                      if(thisEditDistance <= editThreshold and thisEditDistance>=0){
                        //selectedHits.emplace_back(txpID,pos,startHit.fwd,hitsIt->readLen, thisEditDistance,"II");
                              hitsIt->editD = thisEditDistance;
                              hitsIt->toAlign = true;
                              /* ROB: No CIGAR right now */
                              /*
                              char* cigar_ = edlibAlignmentToCigar(thisEdlibResult.alignment, thisEdlibResult.alignmentLength, EDLIB_CIGAR_STANDARD);
                              std::string cigar(cigar_);
                              hitsIt->cigar = cigar;
                              */
                      } else {
                        hitsIt->editD = thisEditDistance;
                        hitsIt->toAlign = false;
                      }
                }

            }
        }

          //std::cout << "\n DONE \n" ;
          return true ;
   }

   template <typename ReadStructT>
       void operator()(ReadStructT& leftReadT,ReadStructT& rightReadT,
                  std::vector<rapmap::utils::QuasiAlignment>& hits,
                  int32_t editThreshold
               ) {
           auto leftHits = hits;
           auto rightHits = hits;
           //compute(leftReadT,hits,editThreshold);
           //for(auto hit: hits)
             //  std::cout<< hit.editD << "\n";
           for(auto& hit: rightHits){
                hit.pos = hit.matePos;
                hit.fwd = hit.mateIsFwd;
                hit.toAlign = hit.mateToAlign;
                hit.editD = hit.mateEditD;
		hit.gapsBegin = hit.mateGapsBegin;
		hit.gapsEnd = hit.mateGapsEnd;
           }

           compute(rightReadT,rightHits,editThreshold);
           compute(leftReadT,leftHits,editThreshold);

           bool edit = false; bool edit_r = true;
           for(int i=0; i<leftHits.size(); i++) {
               if(leftHits[i].editD==-1)
                    leftHits[i].editD = (editThreshold+1)*3;
                if(rightHits[i].editD==-1)
                    rightHits[i].editD = (editThreshold+1)*3;

                hits[i].editD = leftHits[i].editD + rightHits[i].editD;
                hits[i].toAlign =  hits[i].editD <= 2*editThreshold;

           }
       }

    template <typename ReadStructT>
	void operator()(ReadStructT& leftReadT,
			std::vector<rapmap::utils::QuasiAlignment>& hits,
			int32_t editThreshold) { 
			compute(leftReadT,hits,editThreshold);
	}


private:
  RapMapIndexT* rmi_ ;
  AlignerEngine ae_;

};

#endif
