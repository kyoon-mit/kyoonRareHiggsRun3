#include "../interface/GenAnalyzer.h"

using namespace ROOT;

/*
* Finds gen-level indices of the Higgs daughter particles
*/
RVecI HiggsDaughtersIdx(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother
) {
    RVecI daughters_idx;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (int i=idx_high; i>=idx_low; i--) {
        if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]]==25 && // if mother is Higgs
            GenPart_pdgId[i]!=25 // but itself is not the Higgs
            ) {
            daughters_idx.push_back(i);
        }
    }
    return daughters_idx;
}

/*
* Finds PDG ID of the Higgs daughter particles
*/
RVecI HiggsDaughtersPDG(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother
) {
    RVecI daughters_pdg;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (int i=idx_high; i>=idx_low; i--) {
        if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]]==25 && // if mother is Higgs
            GenPart_pdgId[i]!=25 // but itself is not the Higgs
            ) {
            daughters_pdg.push_back(GenPart_pdgId[i]);
        }
    }
    return daughters_pdg;
}

/*
* Finds gen-level indices of all the particles given in the input.
*
* The returned vector will list the decay particle after each mother particle.
*
* @param mother_idx Gen-level indices of the particle whose daughters are to be found.
*/
RVecI GenericDaughtersIdx(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother,
    const RVecI& mother_idx // indices of the mothers whom you want to check the daughters of
) {
    RVecI daughters_idx;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (const int& idx: mother_idx) {
        daughters_idx.push_back(idx); // insert mother idx for reference
        for (int i=idx_high; i>=idx_low; i--) {
            if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
            if (GenPart_genPartIdxMother[i]==idx) {
                daughters_idx.push_back(i);
            }
        }
    }
    return daughters_idx;
}

/*
* Finds PDG ID of all the particles given in the input.
*
* The returned vector will list the decay particle after each mother particle.
*
* @param mother_idx Gen-level indices of the particle whose daughters are to be found.
*/
RVecI GenericDaughtersPDG(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother,
    const RVecI& mother_idx // indices of the mothers whom you want to check the daughters of
) {
    RVecI daughters_pdg;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (const int& idx: mother_idx) {
        daughters_pdg.push_back(GenPart_pdgId[idx]); // insert mother PDG for reference
        for (int i=idx_high; i>=idx_low; i--) {
            if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
            if (GenPart_genPartIdxMother[i]==idx) {
                daughters_pdg.push_back(GenPart_pdgId[i]);
            }
        }
    }
    return daughters_pdg;
}