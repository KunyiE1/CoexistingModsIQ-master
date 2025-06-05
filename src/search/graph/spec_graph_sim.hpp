#ifndef TOPPIC_SEARCH_GRAPH_SPEC_GRAPH_SIM_HPP_
#define TOPPIC_SEARCH_GRAPH_SPEC_GRAPH_SIM_HPP_

#include <vector>

#include "ms/spec/prm_peak.hpp"
#include "ms/spec/spectrum_set.hpp"
#include "search/graph/dist.hpp"
#include "search/graph/graph.hpp"

namespace toppic {

class SpecGraph_sim {
 public:

  SpecGraph_sim(std::vector<std::pair<PeakPtr, std::string>> peak_vec,
                MassGraphPtr mass_graph_ptr, double convert_ratio);

  MassGraphPtr getMassGraphPtr() {return graph_ptr_;}

  DistVec getDistVec() {return dist_;}

  std::pair<PeakPtr, std::string> getPeakPtr(int i) {return peak_vec_[i];}

  const std::vector<std::pair<PeakPtr, std::string>>& getPeakPtrVec() {return peak_vec_;}


  

  int getPeakDist(int v1, int v2);

 private:

  int node_num_;

  int pair_num_;

  std::vector<int> peak_dists_;

  MassGraphPtr graph_ptr_;

  std::vector<std::pair<PeakPtr, std::string>> peak_vec_;

  DistVec dist_;

  int getVecIndex(int v1, int v2);

  void compSpecDistances(double convert_ratio);
};

typedef std::shared_ptr<SpecGraph_sim> SpecGraphPtr_sim;
typedef std::vector<SpecGraphPtr_sim> SpecGraphPtrVec_sim;

}  // namespace toppic

#endif /* SPEC_GRAPH_HPP_ */