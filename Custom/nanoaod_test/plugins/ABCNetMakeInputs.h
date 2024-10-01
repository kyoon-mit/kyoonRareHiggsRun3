#ifndef ABCNetMakeInputs_h
#define ABCNetMakeInputs_h

#include <string>
#include <vector>
#include <unordered_map>

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"

namespace abcnet {

  class ABCNetMakeInputs {
  public:
    ABCNetMakeInputs() {}
    ~ABCNetMakeInputs() {}

    static std::unordered_map<std::string, std::vector<float>> makeFeatureMap( const reco::CandidateView * PFCol, bool debug = false);
    constexpr static unsigned max_num_PFCandidates = 4000;

  };
}  // namespace abcnet

#endif  // ABCNetMakeInputs_h
