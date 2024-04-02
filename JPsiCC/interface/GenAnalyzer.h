#ifndef JPsiCC_GenAnalyzer_h
#define JPsiCC_GenAnalyzer_h

#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;
using RVecI = ROOT::VecOps::RVec<int>;

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

#endif