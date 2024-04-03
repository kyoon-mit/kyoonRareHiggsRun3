#include "../interface/GenAnalyzer.h"
#include <math.h> 

using namespace ROOT;

/*
* Select items in a vector of float given a vector of indices
*/
RVecF SelectByIdx(
    const RVecF& FloatVector,
    const RVecI& Idx
) {
    RVecF return_vals;
    for (const int& idx: Idx) {
        return_vals.push_back(FloatVector[idx]);
    }
    return return_vals;
}

/*
* Finds the index of the Higgs
*/
int HiggsIdx(
    const RVecI& GenPart_pdgId,
    const RVecI& GenPart_genPartIdxMother
) {
    int higgs_idx = -1;
    int idx_low = 0;
    int idx_high = GenPart_pdgId.size();
    for (int i=idx_high; i>=idx_low; i--) {
        if (GenPart_genPartIdxMother[i] < 0 || GenPart_genPartIdxMother[i] > idx_high) continue; // TODO: why these nonsensical indices?
        if (GenPart_pdgId[GenPart_genPartIdxMother[i]]==25 && // if mother is Higgs
            GenPart_pdgId[i]!=25 // but itself is not the Higgs
            ) {
            higgs_idx = GenPart_genPartIdxMother[i]; // mother as the Higgs
            break;
        }
    }
    return higgs_idx;
}

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

/*
* Computes the dR between two particles.
*/
float DeltaR(
    const float& eta1, const float& phi1,
    const float& eta2, const float& phi2
) {
    return sqrt(pow(eta1-eta2, 2.) + pow(phi1-phi2, 2.));
}

/*
* Find index based on PDG.
*/
int IndexFindPDG(
    const RVecI& vec_pdgId, const int& PDGID
) {
    int idx = 999999999;
    for (long unsigned int i=0; i<vec_pdgId.size(); i++) {
        if (vec_pdgId[i] == PDGID) {
            idx = i;
            break;
        }
    }
    return idx;
}

/*
* Sort the leading and sub-leading particle indices according to their transverse momenta.
*/
std::array<int, 2> SortByPt(
    const int& index1, const int& index2, RVecF& vec_pt
) {
    std::array<int, 2> sorted_indices;
    if (vec_pt[index1] >= vec_pt[index2]) {
        sorted_indices[0] = index1;
        sorted_indices[1] = index2;
    } else {
        sorted_indices[0] = index2;
        sorted_indices[1] = index1;
    }
    return sorted_indices;
}

/*
* Computes the PxPyPzE four-vector sum of two particles.
*/
ROOT::Math::PxPyPzEVector SumPxPyPzE(
    const float& px1, const float& py1, const float& pz1, const float& E1,
    const float& px2, const float& py2, const float& pz2, const float& E2
) {
    ROOT::Math::PxPyPzEVector particle1(px1, py1, pz1, E1);
    ROOT::Math::PxPyPzEVector particle2(px2, py2, pz2, E2);

    ROOT::Math::PxPyPzEVector sum_vector = particle1 + particle2;
    return sum_vector;
}

/*
* Computes the PxPyPzM four-vector sum of two particles.
*/
ROOT::Math::PxPyPzMVector SumPxPyPzM(
    const float& px1, const float& py1, const float& pz1, const float& m1,
    const float& px2, const float& py2, const float& pz2, const float& m2
) {
    ROOT::Math::PxPyPzMVector particle1(px1, py1, pz1, m1);
    ROOT::Math::PxPyPzMVector particle2(px2, py2, pz2, m2);

    ROOT::Math::PxPyPzMVector sum_vector = particle1 + particle2;
    return sum_vector;
}

/*
* Computes the PtEtaPhiE four-vector sum of two particles.
*/
ROOT::Math::PtEtaPhiEVector SumPtEtaPhiE(
    const float& pt1, const float& eta1, const float& phi1, const float& E1,
    const float& pt2, const float& eta2, const float& phi2, const float& E2
) {
    ROOT::Math::PtEtaPhiEVector particle1(pt1, eta1, phi1, E1);
    ROOT::Math::PtEtaPhiEVector particle2(pt2, eta2, phi2, E2);

    ROOT::Math::PtEtaPhiEVector sum_vector = particle1 + particle2;
    return sum_vector;
}