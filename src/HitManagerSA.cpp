#include "HitManager.hpp"
#include "HitManagerSA.hpp"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "FastxParser.hpp"

#include "edlib.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>

#include <type_traits>

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

    //template <typename readT>
    bool mergeLeftRightMap(fastx_parser::ReadPair& rpair,
    		              SAHitMap& leftMap,
						  SAHitMap& rightMap,
						  std::vector<QuasiAlignment>& jointHits
    		               ){
    	//using OffsetT = typename RapMapIndexT::IndexType;

    	constexpr const int32_t signedZero{0};
    	uint32_t maxInsertSize_{1000} ;
    	bool foundHit{false} ;


    		std::set<int> commonTids ;
    		for(auto& phLeft : leftMap){
    			auto leftTid = phLeft.first ;
    			for(auto& phRight: rightMap){
    				auto rightTid = phRight.first ;
    				if(leftTid == rightTid){
    					commonTids.insert(leftTid);
    				}
    			}
    		}
    		if(commonTids.size() > 0){
    			//for each common tid check
    			//if both the end tell us
    			//something
    			for(auto tid:commonTids){
    				if(leftMap[tid].active and rightMap[tid].active){
						auto& lefttqvec = leftMap[tid].tqvec ;
						auto& righttqvec = rightMap[tid].tqvec ;
						std::sort(lefttqvec.begin(),
								  lefttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						std::sort(righttqvec.begin(),
								  righttqvec.end(),
								  [](const SATxpQueryPos& a, const SATxpQueryPos& b)-> bool {
															return a.pos < b.pos;
													} );

						//look at the things in pair
						//Make possible pairs
						//Discard pairs which are
						//impossible in the sense
						//that they are far(1000 bp) apart
						std::map<int32_t , int> leftHitCov ;
						{
							std::map<int32_t , int> hitEndPos ;
							auto it = lefttqvec.begin();
							while(it != lefttqvec.end()){
								int numOfCov = 0;
								int32_t hitKey = it->pos - it->queryPos ;
								if(leftHitCov.count(hitKey) == 0){
									leftHitCov[hitKey] = it->matchedLen;
									hitEndPos[hitKey] = it->pos + it->matchedLen;
								}
								else{
									if(hitEndPos[hitKey] > it->pos){
										leftHitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}else{
										leftHitCov[hitKey] += it->matchedLen;
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}
								}
								++it ;
							}
						}

						std::map<int32_t , int> rightHitCov ;
						{
							std::map<int32_t , int> hitEndPos ;
							auto it = righttqvec.begin();
							while(it != righttqvec.end()){
								int numOfCov = 0;
								int32_t hitKey = it->pos - it->queryPos ;
								if(rightHitCov.count(hitKey) == 0){
									rightHitCov[hitKey] = it->matchedLen;
									hitEndPos[hitKey] = it->pos + it->matchedLen;
								}
								else{
									if(hitEndPos[hitKey] > it->pos){
										rightHitCov[hitKey] += ( it->pos + it->matchedLen - hitEndPos[hitKey] );
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}else{
										rightHitCov[hitKey] += it->matchedLen;
										hitEndPos[hitKey] = it->pos + it->matchedLen ;
									}
								}
								++it ;
							}
						}


						std::vector<std::pair<SATxpQueryPos, SATxpQueryPos>> validPairs ;

						for(auto& leftTQPos : lefttqvec){
							for(auto& rightTQPos: righttqvec){
								// Check for the distance
								// go with it otherwise
								// go with something else
								int32_t hitPos1 = leftTQPos.pos - leftTQPos.queryPos ;
								int32_t hitPos2 = rightTQPos.pos - rightTQPos.queryPos ;

								if(std::abs(hitPos1 - hitPos2) < maxInsertSize_){
									validPairs.push_back(std::make_pair(leftTQPos, rightTQPos)) ;
								}
							}
						}//end for

						//Now I have valid pairs so I can directly make QuasiAlignment vectors
						//
						if(validPairs.size() > 0){
							foundHit = true;
						}

						for(auto& vp: validPairs){
							//I will try to make QuasiAlignment objects
							//real alignments and edit distance are not present
							//yet

							int32_t hitPos1 = vp.first.pos - vp.first.queryPos ;
							int32_t hitPos2 = vp.second.pos - vp.second.queryPos ;

							int32_t startRead1 = std::max(hitPos1, signedZero);
							int32_t startRead2 = std::max(hitPos2, signedZero);

							bool read1First{(startRead1 < startRead2)};
							int32_t fragStartPos = read1First ? startRead1 : startRead2;
							int32_t fragEndPos = read1First ?
							       (startRead2 + rpair.second.seq.length()) : (startRead1 + rpair.first.seq.length());
							uint32_t fragLen = fragEndPos - fragStartPos;

							jointHits.emplace_back(
									tid,
									hitPos1,
									!vp.first.queryRC,
									rpair.first.seq.length(),
									vp.first.lcpLength,
									fragLen,
									true
									);
							auto& qaln = jointHits.back();

							qaln.mateLen = rpair.second.seq.length();
							qaln.matePos = hitPos2 ;
							qaln.mateIsFwd = !vp.second.queryRC ;
						}
    				}//end if
    			}//end for
    		}else{
    			foundHit = false;
    		}
    	return foundHit ;
    }


    template <typename RapMapIndexT>
	bool mergeLeftRightSAInts(
						fastx_parser::ReadPair& rpair,
						bool lhp,
						bool rhp,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightFwdSAInts,
						std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightRcSAInts,
						std::vector<QuasiAlignment>& jointHits,
						RapMapIndexT& rmi,
						bool maxNumHits,
						bool consistentHits,
						rapmap::utils::HitCounters& hctr
						){

    	using OffsetT = typename RapMapIndexT::IndexType;

    	AlignerEngine ae_ ;

    	auto& SA = rmi.SA;
    	auto& txpStarts = rmi.txpOffsets;
    	auto& txpIDs = rmi.positionIDs;

    	OffsetT maxInsertSize_{1000} ;

    	bool foundHit{false} ;




    	SAHitMap leftFwdMap ;
    	SAHitMap leftRcMap ;
    	SAHitMap rightFwdMap ;
    	SAHitMap rightRcMap ;

    	if(leftFwdSAInts.size() > 0 or leftRcSAInts.size() > 0){
    		if(leftFwdSAInts.size() > 1){
    			leftFwdMap = rapmap::hit_manager::unionSAHits(leftFwdSAInts, rmi, rpair.first.seq.length(),consistentHits);
    		}else if(leftFwdSAInts.size() > 0){
    			auto inHit = leftFwdSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				leftFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				leftFwdMap[tid].active = true;
    				leftFwdMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}

    		if(leftRcSAInts.size() > 1){
    			leftRcMap = rapmap::hit_manager::unionSAHits(leftRcSAInts, rmi, rpair.first.seq.length(), consistentHits);
    		}else if(leftRcSAInts.size() > 0){
    			auto inHit = leftRcSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				leftRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				leftRcMap[tid].active = true;
    				leftRcMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    	}

    	if(rightFwdSAInts.size() > 0 or rightRcSAInts.size() > 0){
    		if(rightFwdSAInts.size() > 1){
    			rightFwdMap = rapmap::hit_manager::unionSAHits(rightFwdSAInts, rmi, rpair.second.seq.length(),consistentHits);
    		}else if(rightFwdSAInts.size() > 0){
    			auto inHit = rightFwdSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				rightFwdMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				rightFwdMap[tid].active = true;
    				rightFwdMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    		if(rightRcSAInts.size() > 1){
    			rightRcMap = rapmap::hit_manager::unionSAHits(rightRcSAInts, rmi, rpair.second.seq.length(), consistentHits);
    		}else if(rightRcSAInts.size() > 0){
    			auto inHit = rightRcSAInts.front() ;
    			for(auto i = inHit.begin ; i < inHit.end; ++i){
    				auto globalPos = SA[i] ;
    				auto tid = rmi.transcriptAtPosition(globalPos) ;
    				auto txpPos = globalPos - txpStarts[tid] ;
    				rightRcMap[tid].tqvec.emplace_back(txpPos, inHit.queryPos, inHit.lcpLength, inHit.queryRC);
    				rightRcMap[tid].active = true;
    				rightRcMap[tid].tqvec.back().matchedLen = inHit.len ;
    			}
    		}
    	}

    	//there are two possibilities
    	bool fwdRc{false} ;
    	bool rcFwd{false} ;
    	//1. The left is from forward and right is from reverse
    	if(leftFwdMap.size() > 0 and rightRcMap.size() > 0){
    		// only consider transcripts that are common between both
    		fwdRc = mergeLeftRightMap(rpair, leftFwdMap, rightRcMap, jointHits) ;
    		foundHit = true ;
    	}

    	if(leftRcMap.size() > 0 and rightFwdMap.size() > 0){
    		rcFwd = mergeLeftRightMap(rpair, leftRcMap, rightFwdMap, jointHits) ;
    		foundHit = true ;
    	}


    	return foundHit ;
    }

    using SAIndex32BitDense = RapMapSAIndex<int32_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int32_t>,
                   					     rapmap::utils::KmerKeyHasher>>;
    using SAIndex64BitDense = RapMapSAIndex<int64_t, RegHashT<uint64_t, rapmap::utils::kmerVal<int64_t>,
   									     rapmap::utils::KmerKeyHasher>>;
    using SAIndex32BitPerfect = RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int32_t>>>;
    using SAIndex64BitPerfect = RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int64_t>>>;

    template
    bool mergeLeftRightSAInts<SAIndex32BitDense>(
         						fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex32BitDense& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex64BitDense>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex64BitDense& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex32BitPerfect>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int32_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int32_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex32BitPerfect& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);
    template
    bool mergeLeftRightSAInts<SAIndex64BitPerfect>(
    							fastx_parser::ReadPair& rpair,
         						bool lhp,
         						bool rhp,
         						std::vector<SAIntervalHit<int64_t>>& leftFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& leftRcSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightFwdSAInts,
         						std::vector<SAIntervalHit<int64_t>>& rightRcSAInts,
         						std::vector<QuasiAlignment>& jointHits,
         						SAIndex64BitPerfect& rmi,
         						bool maxNumHits,
         						bool consistentHits,
         						rapmap::utils::HitCounters& hctr
         						);



    }
}
