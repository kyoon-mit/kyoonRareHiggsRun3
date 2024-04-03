#ifndef JPsiCC_GenAnalyzer_h
#define JPsiCC_GenAnalyzer_h

#include "ROOT/RVec.hxx"
#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"

using namespace ROOT::VecOps;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecF = ROOT::VecOps::RVec<float>;

RVecF SelectByIdx(const RVecF& FloatVector,
                  const RVecI& Idx);
int HiggsIdx(const RVecI& GenPart_pdgId,
             const RVecI& GenPart_genPartIdxMother);
RVecI HiggsDaughtersIdx(const RVecI& GenPart_pdgId,
                        const RVecI& GenPart_genPartIdxMother);
RVecI HiggsDaughtersPDG(const RVecI& GenPart_pdgId,
                        const RVecI& GenPart_genPartIdxMother);
RVecI GenericDaughtersIdx(const RVecI& GenPart_pdgId,
                          const RVecI& GenPart_genPartIdxMother,
                          const RVecI& mother_idx);
RVecI GenericDaughtersPDG(const RVecI& GenPart_pdgId,
                          const RVecI& GenPart_genPartIdxMother,
                          const RVecI& mother_idx);
float DeltaR(
    const float& eta1, const float& phi1,
    const float& eta2, const float& phi2);
int IndexFindPDG(
    const RVecI& vec_pdgId, const int& PDGID);
std::array<int, 2> SortByPt(
    const int& index1, const int& index2, RVecF& vec_pt); 
ROOT::Math::PxPyPzEVector SumPxPyPzE(
    const float& px1, const float& py1, const float& pz1, const float& E1,
    const float& px2, const float& py2, const float& pz2, const float& E2);
ROOT::Math::PxPyPzMVector SumPxPyPzM(
    const float& px1, const float& py1, const float& pz1, const float& m1,
    const float& px2, const float& py2, const float& pz2, const float& m2);
ROOT::Math::PtEtaPhiEVector SumPtEtaPhiM(
    const float& pt1, const float& eta1, const float& phi1, const float& E1,
    const float& pt2, const float& eta2, const float& phi2, const float& E2);

#endif