//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_SEARCH_GRAPH_ALIGN_GRAPH_ALIGN_HPP_
#define TOPPIC_SEARCH_GRAPH_ALIGN_GRAPH_ALIGN_HPP_

#include "prsm/prsm.hpp"
#include "search/oneptmsearch/diagonal_header.hpp"
#include "search/graph/graph.hpp"
#include "search/graph/proteo_graph.hpp"
#include "search/graph/spec_graph.hpp"
#include "search/graphalign/graph_dp_node.hpp"
#include "search/graphalign/graph_result_node.hpp"
#include "search/graphalign/graph_align_mng.hpp"

namespace toppic {


typedef std::vector<std::vector<std::vector<std::vector<std::pair<int, int>>>>> ConsistentPairs;
typedef std::vector<std::vector<std::vector<std::pair<unsigned int, std::vector<std::pair<std::pair<unsigned short, unsigned short>,std::vector<std::pair<unsigned short,unsigned short>>>>>>>> NewConsPairs;
typedef std::pair<std::pair<int,Vertex_AGraph>,std::pair<int,Vertex_AGraph>> PrePathInfo;

struct hashKey_tuple{
  template <class T1, class T2, class T3>
  size_t operator()(const std::tuple<T1, T2, T3>& p) const
  {
    auto hash1 = std::hash<T1>{}(std::get<0>(p));
    auto hash2 = std::hash<T2>{}(std::get<1>(p));
    auto hash3 = std::hash<T3>{}(std::get<2>(p));
    return hash1 ^ hash2 ^ hash3;
  } 
};


typedef std::pair<std::pair<unsigned short,unsigned short>,std::tuple<short int, unsigned int, std::vector<std::pair<unsigned short,unsigned short>>>> prePosition;




class GraphAlign {
 public:
  GraphAlign(GraphAlignMngPtr mng_ptr, ProteoGraphPtr proteo_graph_ptr,
             SpecGraphPtr spec_graph_ptr);

  //void process(short prot_start, short spec_start);

  void TopMGFast();

  void process_v2();

  //void oriProcess();

  //void diagonal_v2(short prot_start, short spec_start);

  //void diagonal_v3(short prot_start, short spec_start);

  PrsmPtr geneResult(int s, int m);

  PrsmPtr geneResult(int s);

 private:
  GraphAlignMngPtr mng_ptr_;

  ProteoGraphPtr proteo_graph_ptr_;

  MassGraphPtr pg_;

  int proteo_ver_num_;

  SpecGraphPtr spec_graph_ptr_;

  MassGraphPtr sg_;

  int spec_ver_num_;

  int n_unknown_shift_;

  DistVec spec_dist_;

  DistVec2D dist_vec_;

  ConsistentPairs cons_pairs_;

  GraphDpNodePtrVec2D table_;

  GraphResultNodePtrVec3D result_nodes_;

  GraphResultNodePtrVec2D nodes_2d_;

  DiagonalHeaderPtrVec diag_headers_; 

  DiagonalHeaderPtrVec2D diag_headers_2d_; 
  

  NewConsPairs new_cons_pairs_;

  std::vector<int> delta;

  std::vector<int> deltaL;

  std::vector<int> deltaR;

  std::vector<int> spectrumMass;

  int maxDelta; //For spectrum

  std::vector<int> originalMass;
  std::vector<std::pair<int,int>> MaxRed;
  std::vector<std::pair<short int, short int>> specRange;


  //void getConsistentPairs();

  //void addToConsistentPairs(int m, const std::vector<std::pair<int, int>> & sp_pair_ij,
  //                          const std::vector<std::pair<int, int>> & pg_pair_ij);

  void initTable();


  GraphDpNodePtr compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num);

  GraphDpNodePtr compBestShiftNode(int i, int j, int s, int m);

  void updateBestShiftNode(int i, int j, int s, int m);

  void dp();

  GraphResultNodePtrVec backtrace(int s, int m);

  void backtrace();

  void getNodeDiagonals(int s, int m);

  void geneHeaders();


  void getNewConsPair();

  //void getNewConsPairOri();


  void addToNewConsistentPairs(int m, const std::vector<std::pair<std::pair<int, int>, std::vector<std::pair<unsigned short,unsigned short>>>> & sp_pair_ij,
                            const std::vector<std::pair<std::pair<int, int>, std::vector<std::pair<unsigned short,unsigned short>>>> & pg_pair_ij);


  void getDelta();

  void computeT_s();

  void computeT_v2(bool case1, double ptm_mass);

  short int binarySearch(std::vector<int> spectrumMass, int key, bool lower);

  void findSpecRange(short prot_start, short spec_start);

  void findModMass(short prot_start);

  //void smallTij(short prot_start, short spec_start);

  void cleanMemory();


  void deleteOverlap();

  void deleteOverlap_v2();

  void printDelta();

  void quantification_case1(std::vector<std::tuple<int, int, int>> endNodes, std::vector<std::vector<std::vector<short int>>> T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> E);

  void quantification_case2(std::vector<std::tuple<int, int, int>> endNodes, std::vector<std::tuple<int, int, int>> endNodes2, std::vector<std::vector<std::vector<short int>>> T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> E);
//  void add_all_info(int i, int j, int k,
//                  std::vector<std::vector<std::vector<short int>>> T,
//                  std::vector<std::vector<std::vector<std::vector<prePosition>>>> E,
//                  AlignmentGraphPtr alignGraph_ptr,
//                  std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple>> vertexMapPtr,
//                  int ori_index, int exactMass, int blackMass, std::vector<std::pair<unsigned short, unsigned short>> modInfo);

  std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> buildAdditionPath(AlignmentGraphPtr alignGraph_ptr, int T1, int T2, Vertex_AGraph beginNode1, std::vector<std::vector<Vertex_AGraph>> layers1, double a_inten);

  std::pair<std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>>, std::vector<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> storePath(AlignmentGraphPtr alignGraph_ptr, PrePathInfo endPath, PrePathInfo prePath, std::vector<std::vector<std::vector<std::pair<double, PrePathInfo>>>> D, std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> additionPath, int T1Larger);


  std::vector<std::vector<std::pair<double, PrePathInfo>>> countError(AlignmentGraphPtr alignGraph_ptr, std::vector<std::vector<Vertex_AGraph>> layers, int layer, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, double abund_a, double abund_b);
  std::vector<std::vector<std::pair<double, PrePathInfo>>> countError_case2(AlignmentGraphPtr alignGraph_ptr, std::vector<std::vector<Vertex_AGraph>> layers1, std::vector<std::vector<Vertex_AGraph>> layers2, int layer1, int layer2, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, double a_inten, double b_inten);


  std::pair<double, PrePathInfo> getMinPreError(AlignmentGraphPtr alignGraph_ptr, std::vector<Vertex_AGraph> pre_layer, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, Vertex_AGraph node_a, Vertex_AGraph node_b);
  std::pair<double, PrePathInfo> getMinPreError_case2(AlignmentGraphPtr alignGraph_ptr, std::vector<Vertex_AGraph> pre_layer1, std::vector<Vertex_AGraph> pre_layer2, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, Vertex_AGraph node_a, Vertex_AGraph node_b);


  void add_all_info(int i, int j, int k, std::vector<std::vector<std::vector<short>>> T,
                    std::vector<std::vector<std::vector<std::vector<prePosition>>>> E,
                    const AlignmentGraphPtr &alignGraph_ptr,
                    const std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple>> &vertexMapPtr,
                    int ori_index, int exactMass, int blackMass,
                    std::vector<std::pair<unsigned short, unsigned short>> modInfo);
  void BFS(std::queue<std::tuple<int,int,int>> node_queue,
           std::vector<std::vector<std::vector<short int>>> T,
           std::vector<std::vector<std::vector<std::vector<prePosition>>>> E,
           AlignmentGraphPtr &alignGraph_ptr,
           std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple>> &vertexMapPtr);

  std::pair<double, std::vector<std::vector<Vertex_AGraph>>> findLargestIntensity(AlignmentGraphPtr &alignGraph_ptr, Vertex_AGraph beginNode);

};


typedef std::shared_ptr<GraphAlign> GraphAlignPtr;
typedef boost::graph_traits<MassGraph>::out_edge_iterator out_edge_iter;
typedef boost::graph_traits<AlignmentGraph>::out_edge_iterator out_alignEdge_iter;
typedef boost::graph_traits<AlignmentGraph>::in_edge_iterator in_alignEdge_iter;




}  // namespace toppic

#endif
