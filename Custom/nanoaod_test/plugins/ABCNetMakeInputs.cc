#include <iostream>
#include <memory>
#include <cmath>
#include <numeric>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //for pat
#include "DataFormats/Candidate/interface/Candidate.h" // for reco
#include "ABCNetMakeInputs.h"

using namespace abcnet;

template <typename A, typename B> void zip(const std::vector<A> &a, const std::vector<B> &b, std::vector<std::pair<A,B>> &zipped) {
  for(size_t i=0; i<a.size(); ++i) {
    zipped.push_back(std::make_pair(a[i], b[i]));
  }
}

template <typename A, typename B> void unzip( const std::vector<std::pair<A, B>> &zipped, std::vector<A> &a, std::vector<B> &b) {
  for(size_t i=0; i<a.size(); i++) {
    a[i] = zipped[i].first;
    b[i] = zipped[i].second;
  }
}

std::unordered_map<std::string, std::vector<float>> ABCNetMakeInputs::makeFeatureMap (const reco::CandidateView * pfCol, bool debug) {
  
  //store PF candidates and their pts into vectors
  //sort the PF candidates based on their pt
  std::vector<pat::PackedCandidate> PFCands;
  std::vector<float> Pts;
  
  for (auto const & aPF : *pfCol ) {
    const pat::PackedCandidate *lPack = dynamic_cast<const pat::PackedCandidate*>(&aPF);    
    PFCands.push_back(*lPack);
    Pts.push_back(lPack->pt());
  } // end loop over PF candidates

  //zip the vectors
  std::vector<std::pair<pat::PackedCandidate,float>> zipped_vec;
  zip(PFCands, Pts, zipped_vec);
  // Sort the vector of pairs
  std::sort(std::begin(zipped_vec), std::end(zipped_vec), [&](const auto& a, const auto& b) { return a.second > b.second; });
  // Write the sorted pairs back to the original vectors
  unzip(zipped_vec, PFCands, Pts);
  
  //fill feature map
  std::unordered_map<std::string, std::vector<float>> fts;
  for (auto const & aPF : PFCands) {
      
    fts["PFCandEta"].push_back(aPF.eta()); //f0
    fts["PFCandPhi"].push_back(aPF.phi()); //f1
    fts["PFCandLogPt"].push_back(log(aPF.pt())); //f2
    fts["PFCandLogE"].push_back(log(aPF.energy())); //f3c
    fts["PFCandCharge"].push_back(aPF.charge()); //f4
    if (abs(aPF.pdgId()) == 11) fts["PFCandIsEle"].push_back(1.0); else fts["PFCandIsEle"].push_back(0.0); //f5
    if (abs(aPF.pdgId()) == 13) fts["PFCandIsMu"].push_back(1.0); else fts["PFCandIsMu"].push_back(0.0); //f6
    if (abs(aPF.pdgId()) == 211 || abs(aPF.pdgId()) == 2 || abs(aPF.pdgId()) == 1) fts["PFCandIsHad"].push_back(1.0); else fts["PFCandIsHad"].push_back(0.0); //f7
    if (abs(aPF.pdgId()) == 130) fts["PFCandIsNeutralHad"].push_back(1.0); else fts["PFCandIsNeutralHad"].push_back(0.0); //f8
    if (abs(aPF.pdgId()) == 22) fts["PFCandIsPhoton"].push_back(1.0); else fts["PFCandIsPhoton"].push_back(0.0); //f9
    fts["PFCandNumHits"].push_back(aPF.numberOfHits()); //f14
    fts["PFCandNumLayersHit"].push_back(aPF.trackerLayersWithMeasurement()); //f15
    fts["PFCandFromPV"].push_back(aPF.fromPV()); //f16
    fts["PFCandTrackHighPurity"].push_back(aPF.trackHighPurity()); //f17
    if (aPF.bestTrack()) {
      if ( isinf(aPF.dz()/aPF.dzError()) ) fts["PFCandDZSig"].push_back(999.0); else fts["PFCandDZSig"].push_back(aPF.dz()/aPF.dzError()); //f10
      if ( isinf(aPF.dxy()/aPF.dxyError()) ) fts["PFCandDXYSig"].push_back(999.0); else fts["PFCandDXYSig"].push_back(aPF.dxy()/aPF.dxyError()); //f11
      const auto *trk = aPF.bestTrack();
      fts["PFCandNormChi2"].push_back(trk->normalizedChi2()); //f13
      fts["PFCandQuality"].push_back(trk->qualityMask()); //f18
    }
    else {
      fts["PFCandDZSig"].push_back(0.0);
      fts["PFCandDXYSig"].push_back(0.0);
      fts["PFCandNormChi2"].push_back(999.0);
      fts["PFCandQuality"].push_back(0.0);
    }
    if (aPF.pdgId() == 130 || aPF.pdgId() == 1) fts["PFCandHCalFrac"].push_back(aPF.hcalFraction()); else if (aPF.isIsolatedChargedHadron()) fts["PFCandHCalFrac"].push_back(aPF.rawHcalFraction()); else fts["PFCandHCalFrac"].push_back(0.0); //f12
    
  }

  return fts;

};