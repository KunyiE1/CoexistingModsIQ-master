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

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <limits>
#include <unordered_map>
#include <cmath>
#include <list>
#include <ctime>
#include <queue>
#include <iomanip>

#include "seq/proteoform_factory.hpp"
#include "ms/spec/extend_ms_factory.hpp"
#include "search/oneptmsearch/diagonal.hpp"
#include "search/oneptmsearch/diagonal_header.hpp"
#include "search/graph/graph.hpp"
#include "search/graphalign/graph_align_a.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "prsm/prsm_algo.hpp"

#define USER_DEFINE_MAX_MOD 5
#define USER_DEFINE_PARAMETER 5


namespace toppic {

std::vector<int> getMinMaxProtDist2(DistVec2D dist_vec) {
  std::vector<int> min_dist, max_dist;
  for (size_t i = 0; i < dist_vec.size(); i++) {
    if (dist_vec[i].size() == 0) continue;
    min_dist.push_back(dist_vec[i][0].dist_);
    max_dist.push_back(dist_vec[i][dist_vec[i].size() - 1].dist_);
  }
  std::sort(min_dist.begin(), min_dist.end());
  std::sort(max_dist.begin(), max_dist.end());
  std::vector<int> res;
  res.push_back(min_dist[0]);
  res.push_back(max_dist[max_dist.size() - 1]);
  return res;
}

GraphAlignSim::GraphAlignSim(GraphAlignMngPtr mng_ptr,
                       ProteoGraphPtr proteo_graph_ptr,
                       SpecGraphPtr_sim spec_graph_ptr) {
  LOG_DEBUG("Graph constructor start");
  mng_ptr_ = mng_ptr;
  proteo_graph_ptr_ = proteo_graph_ptr;
  spec_graph_ptr_ = spec_graph_ptr;

  dist_vec_ = proteo_graph_ptr_->getDistVec2D();
  spec_dist_ = spec_graph_ptr_->getDistVec();

  std::vector<int> cutoff = getMinMaxProtDist2(dist_vec_);
  cutoff[0] -= mng_ptr_->getIntTolerance();
  cutoff[1] += mng_ptr_->getIntTolerance();
  //std::cout << "cutoff[0]: " << cutoff[0] << ", cutoff[1]: " << cutoff[0] << std::endl;

  spec_dist_.erase(
      std::remove_if(spec_dist_.begin(), spec_dist_.end(),
                     [cutoff](Dist d){return d.dist_ < cutoff[0] || d.dist_ > cutoff[1];}),
      spec_dist_.end());

  std::sort(spec_dist_.begin(), spec_dist_.end(), distVecUp);

  for (int i = 0; i < mng_ptr->max_known_mods_ + 1; i++) {
    std::sort(dist_vec_[i].begin(), dist_vec_[i].end(), distVecUp);
  }

  pg_ = proteo_graph_ptr_->getMassGraphPtr();
  sg_ = spec_graph_ptr_->getMassGraphPtr();
  proteo_ver_num_ = num_vertices(*pg_.get());
  spec_ver_num_ = num_vertices(*sg_.get());
  LOG_DEBUG("Graph constructor end");
}


/*
void GraphAlignSim::getConsistentPairs() {
  LOG_DEBUG("consistent pair start");
  int tole = mng_ptr_->getIntTolerance();
  //std::cout << "tole in getConsistentPairs is: " << tole << std::endl;
  LOG_DEBUG("Integer error tolerance " << tole);
  std::vector<std::pair<int, int>> empty_list;
  std::vector<std::vector<std::pair<int, int>>> empty_vec(mng_ptr_->max_known_mods_ + 1, empty_list);
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::vector<std::pair<int, int>>>> empty_vec_2d;
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    cons_pairs_.push_back(empty_vec_2d);
  }

  int min_dist = mng_ptr_->getIntMinConsistentDist();
  //std::cout << "min_dist: " << min_dist << std::endl;
  for (size_t m = 0; m < dist_vec_.size(); m++) {
    if (dist_vec_[m].size() == 0) continue;

    size_t prot_idx_min = 0, prot_idx_max = 0;
    for (size_t spec_idx = 0; spec_idx < spec_dist_.size(); spec_idx++) {
      Dist distance = spec_dist_[spec_idx];
      int sp_dist = distance.dist_;
      if (sp_dist < min_dist) continue;

      bool flag = true;
      while (prot_idx_min < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_min].dist_ >= sp_dist - tole) {
          flag = false;
        } else {
          prot_idx_min++;
        }
      }

      if (prot_idx_min >= dist_vec_[m].size()) continue;

      prot_idx_max = std::max(prot_idx_min, prot_idx_max);

      flag = true;

      while (prot_idx_max < dist_vec_[m].size() && flag) {
        if (dist_vec_[m][prot_idx_max].dist_ > sp_dist + tole) {
          flag = false;
        } else {
          prot_idx_max++;
        }
      }

      for (size_t t = prot_idx_min; t < prot_idx_max; t++) {
        Dist new_distance = dist_vec_[m][t];
        if (std::abs(sp_dist - new_distance.dist_) <= tole){
          //std::cout << "m: " << m << "spec pair: " << distance.pair_ij_.size()<< ";proteo pair: " << new_distance.pair_ij_.size() << std::endl;
          addToConsistentPairs(m, distance.pair_ij_, new_distance.pair_ij_);
        }
      }
    }
    //std::cout << "dist_vec_[" << m << "].size() = " << a.size() << "; Pair list size is: " << totalSize << std::endl;
    dist_vec_[m].clear();
    //for (int a = 0; a < cons_pairs_[4][29][m].size(); a++){
    //  std::cout << "(4, 29): (" << cons_pairs_[4][29][m][a].first << ", " << cons_pairs_[4][29][m][a].second << ")" << std::endl;
    //}
  }
  dist_vec_.clear();
  LOG_DEBUG("consistent pair end");
}



void GraphAlignSim::addToConsistentPairs(int m, const std::vector<std::pair<int, int>> & sp_pair_ij,
                                      const std::vector<std::pair<int, int>> & pg_pair_ij) {
  for (size_t k = 0; k < pg_pair_ij.size(); k++) {
    for (size_t sp = 0; sp < sp_pair_ij.size(); sp++) {
      int pr_v1 = pg_pair_ij[k].first;
      int pr_v2 = pg_pair_ij[k].second;
      int sp_v1 = sp_pair_ij[sp].first;
      int sp_v2 = sp_pair_ij[sp].second;
      std::pair<int, int> pre_pair(pr_v1, sp_v1);
      cons_pairs_[pr_v2][sp_v2][m].push_back(pre_pair);
    }
  }
}
*/


void GraphAlignSim::initTable() {
  LOG_DEBUG("init table start");
  double node_score = 1.0;
  for (int i = 0; i < proteo_ver_num_; i++) {
    GraphDpNodePtrVec node_vec;
    for (int j = 0; j < spec_ver_num_; j++) {
      GraphDpNodePtr node_ptr
          = std::make_shared<GraphDpNode>(i, j, node_score, mng_ptr_->n_unknown_shift_,
                                          mng_ptr_->max_known_mods_);
      node_vec.push_back(node_ptr);
    }
    table_.push_back(node_vec);
  }
  LOG_DEBUG("init table step 1");
  if (mng_ptr_->whole_protein_only_) {
    int init_aa_num = 2;
    if (init_aa_num > proteo_ver_num_) {
      init_aa_num = proteo_ver_num_;
    }
    for (int i = 0; i < init_aa_num; i++) {
      table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, node_score);
      table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
      // LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
    }
  }
  else{
    for (int i = 0; i < proteo_ver_num_; i++) {
      table_[i][0]->updateTable(0, 0, GRAPH_ALIGN_TYPE_NULL, 0,  nullptr, node_score);
      table_[i][0]->updateBestShiftNode(0, 0, 0, table_[i][0]);
      // LOG_DEBUG("type " << table_[i][0]->getType(0) << " first index " << table_[i][0]->getFirstIdx() << " second index " << table_[i][0]->getSecondIdx());
    }
  } 
  LOG_DEBUG("init table end");
}


GraphDpNodePtr GraphAlignSim::compBestVariableNode(int i, int j, int s, int m, int &best_edge_mod_num) {
  int best_prev_score = -1;
  best_edge_mod_num = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  for (int p = 0; p <= m; p++) {
    for (size_t q = 0; q < cons_pairs_[i][j][p].size(); q++) {
      std::pair<int, int> pair = cons_pairs_[i][j][p][q];
      int pi = pair.first;
      int pj = pair.second;
      // LOG_DEBUG("pi " << pi << " pj " << pj);
      int score = table_[pi][pj]->getBestScore(s, m-p);
      if (score > best_prev_score) {
        best_prev_score = score;
        best_prev_node = table_[pi][pj];
        best_edge_mod_num = p;
      }
    }
  }
  //if(i == 4 && j == 29 && m == 0){
  //  std::cout << "score: " << best_prev_score << "; (" << best_prev_node->getFirstIdx() << ", " << best_prev_node->getSecondIdx() << "); mod num: " << best_edge_mod_num << std::endl;
  //}
  return best_prev_node;
}


GraphDpNodePtr GraphAlignSim::compBestShiftNode(int i, int j, int s, int m) {
  if (s == 0) {
    return nullptr;
  }
  int best_prev_score = -1;
  GraphDpNodePtr best_prev_node = nullptr;
  GraphDpNodePtr up_node = table_[i-1][j];
  int score = up_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = up_node->getBestShiftNodePtr(s-1, m);
  }
  GraphDpNodePtr left_node = table_[i][j-1];
  score = left_node->getBestShiftScore(s-1, m);
  if (score > best_prev_score) {
    best_prev_score = score;
    best_prev_node = left_node->getBestShiftNodePtr(s-1, m);
  }
  return best_prev_node;
}


void GraphAlignSim::updateBestShiftNode(int i, int j, int s, int m) {
  // update best node
  if (table_[i-1][j]->getBestShiftScore(s, m) > table_[i][j-1]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i-1][j]->getBestShiftScore(s, m), table_[i-1][j]->getBestShiftNodePtr(s, m));
  } else {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j-1]->getBestShiftScore(s, m), table_[i][j-1]->getBestShiftNodePtr(s, m));
  }

  if (table_[i][j]->getBestScore(s, m) > table_[i][j]->getBestShiftScore(s, m)) {
    table_[i][j]->updateBestShiftNode(s, m, table_[i][j]->getBestScore(s, m), table_[i][j]);
  }
}

void GraphAlignSim::dp() {
  LOG_DEBUG("dp start");
  //std::cout << "n_unknown_shift_: " << mng_ptr_->n_unknown_shift_ << std::endl;
  //std::cout << "max_known_mods_: " << mng_ptr_->max_known_mods_ << std::endl;
  //std::cout << "- std::numeric_limits<double>::max(): " << - std::numeric_limits<double>::max() << std::endl;
  for (int i = 1; i < proteo_ver_num_; i++) {
    for (int j = 1; j < spec_ver_num_; j++) {
      // compute for zero shift
      for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
        for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
          int edge_mod_num;
          GraphDpNodePtr best_var_node = compBestVariableNode(i, j, s, m, edge_mod_num);
          double var_score;
          if (best_var_node == nullptr) {
            var_score = - std::numeric_limits<double>::max();
          } else {
            var_score = best_var_node->getBestScore(s, m-edge_mod_num);
          }
          GraphDpNodePtr best_shift_node = compBestShiftNode(i, j, s, m);
          double shift_score;
          if (best_shift_node != nullptr) {
            shift_score = best_shift_node->getBestScore(s-1, m);
          } else {
            shift_score = - std::numeric_limits<double>::max();
          }
          double new_score = table_[i][j]->getNodeScore();
          // LOG_DEBUG("new score " << new_score);
          
          if (var_score >= shift_score) {
            if (var_score ==  - std::numeric_limits<double>::max()) {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_NULL,
                                        0, nullptr, -std::numeric_limits<int>::max());
            } else {
              table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_VARIABLE,
                                        edge_mod_num, best_var_node, var_score + new_score);
            }
          } else {
            table_[i][j]->updateTable(s, m, GRAPH_ALIGN_TYPE_UNEXPECTED,
                                      0, best_shift_node, shift_score + new_score);
          }
          updateBestShiftNode(i, j, s, m);
        }
      }
    }
  }
  LOG_DEBUG("dp end");
}

ReturnStruct GraphAlignSim::backtrace(int s, int m) {
  // find the best score;
  int best_score = -1;
  GraphDpNodePtr best_node_ptr = nullptr;

  /*
  for (int i = 0; i < proteo_ver_num_; i++) {
    int score = table_[i][spec_ver_num_-1]->getBestScore(s, m);
    if (score > best_score) {
      best_score = score;
      best_node_ptr = table_[i][spec_ver_num_-1];
    }
  }
  */

  for (int i = 0; i < proteo_ver_num_; i++) {
    for(int j = 0; j < spec_ver_num_; j++){
      int score = table_[i][j]->getBestScore(s, m);
      if (score > best_score) {
        best_score = score;
        best_node_ptr = table_[i][j];
      }
    }
  }

  
  int shift = s;
  int mod = m;
  LOG_DEBUG("best score " << best_score);
  GraphResultNodePtrVec results;
  short prot_start, spec_start;
  if (best_score > 0) {
    GraphDpNodePtr cur_node_ptr = best_node_ptr;
    //std::cout << "---------------------current node---------------------" << std::endl;
    //std::cout << "current best score: " << best_score << std::endl;
    while (cur_node_ptr != nullptr) {
      LOG_DEBUG("cur node " << cur_node_ptr);
      results.push_back(std::make_shared<GraphResultNode>(cur_node_ptr, shift, mod));
      int type = cur_node_ptr->getPrevEdgeType(shift, mod);
      prot_start = cur_node_ptr->getFirstIdx();
      spec_start = cur_node_ptr->getSecondIdx();
      //std::cout << "(" << cur_node_ptr->getFirstIdx() << ", " << cur_node_ptr->getSecondIdx() << ") previous edge type: " << type << "; PrevEdgeModNum: " << cur_node_ptr->getPrevEdgeModNum(shift,mod) << "; get node score: " << cur_node_ptr->getNodeScore() << "; get best score: " << cur_node_ptr->getBestScore(shift,mod) << std::endl;
      LOG_DEBUG("type " << type << " shift " << shift << " first index " << cur_node_ptr->getFirstIdx() << " second index " << cur_node_ptr->getSecondIdx());
      int prev_edge_mod_num = cur_node_ptr->getPrevEdgeModNum(shift, mod);
      cur_node_ptr = cur_node_ptr->getPrevNodePtr(shift, mod);
      LOG_DEBUG("get prev node ");
      if (type == GRAPH_ALIGN_TYPE_UNEXPECTED) {
        shift--;
      }
      mod = mod - prev_edge_mod_num;
    }
  }
  LOG_DEBUG("obtained result node ptr vec");
  std::reverse(results.begin(), results.end());
  ReturnStruct para = std::make_tuple(best_score,prot_start,spec_start,results);
  return para;
}

ReturnStruct GraphAlignSim::backtrace() {
  LOG_DEBUG("start back trace");
  //std::cout << "*******************start backtrace*********************" << std::endl;
  result_nodes_.clear();
  //std::cout << "n_unknow_shift: " << mng_ptr_->n_unknown_shift_ << std::endl;
  //std::cout << "max_knwon_mods: " << mng_ptr_->max_known_mods_ << std::endl;
  int best_score = 0;
  ReturnStruct best_para;
  for (int s = 0; s <= mng_ptr_->n_unknown_shift_; s++) {
    LOG_DEBUG("shift num " << s);
    GraphResultNodePtrVec2D nodes_2d;
    for (int m = 0; m <= mng_ptr_->max_known_mods_; m++) {
      //std::cout << "--------------------------------" << std::endl;
      //std::cout << "parameter m = " << m << std::endl;
      ReturnStruct para = backtrace(s,m);
      int current_score = std::get<0>(para);
      if(current_score > best_score){
        best_score = current_score;
        best_para = para;
      }
      nodes_2d.push_back(std::get<3>(para));
    }
    result_nodes_.push_back(nodes_2d);
  }
  // print the size of result_nodes_;
  return best_para;
  
  LOG_DEBUG("end back trace");
}

/*
ReturnStruct GraphAlignSim::topMG(){
  clock_t startTime, endTime;
  startTime = clock();
  getConsistentPairs();
  initTable();
  dp();
  ReturnStruct best_para = backtrace();
  endTime = clock();
  double runningTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
  std::cout << "Totle running time is: " << runningTime * 1000 << " ms" << std::endl;
  return best_para;
  
}
*/


void GraphAlignSim::TopMGFast(){
  //std::cout << "here start" << std::endl;
  clock_t startTime, midTime, endTime;
  
  startTime = clock();
  //std::cout << "start1" << std::endl;
  getDelta_ori();
  //std::cout << "start12" << std::endl;
  deleteOverlap_v2();
  //std::cout << "start34" << std::endl;
  getNewConsPair();
  //std::cout << "end" << std::endl;
  
  midTime = clock();
  //computeT_s();
  bool case1 = false;
  double ptm_mass = 79.966331;
  std::cout<<"bgn computTv2"<<std::endl;
  computeT_v2(case1, ptm_mass);

  endTime = clock();
  double cijTime = (double)(midTime - startTime) / CLOCKS_PER_SEC;
  double runningTime = (double)(endTime - midTime) / CLOCKS_PER_SEC;
  std::cout << "Running time of creating cij is: " << cijTime * 1000 << " ms" << std::endl;
  std::cout << "Running time of dp is: " << runningTime * 1000 << " ms" << std::endl;
  std::cout << "Totle running time is: " << (runningTime + cijTime) * 1000 << " ms" << std::endl; 
}



void GraphAlignSim::getDelta_complex(double totalMass){

  spectrumMass.clear();
  auto aa = spec_graph_ptr_->getPeakPtrVec();
  for(auto peakIter = aa.begin(); peakIter != aa.end(); peakIter++){
    int mass = static_cast<int>(std::round((*peakIter).first->getPosition() * mng_ptr_->convert_ratio_));
    spectrumMass.push_back(mass);
  }

  delta.clear();


  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  maxDelta = 0;
  delta.push_back(0);
  for(int i = 1; i < spec_ver_num_; i++){
    auto a = aa[i].second;
    //auto a = spec_graph_ptr_->getPrmPeakPtr(i)->getBaseTypePtr()->getName();
    if(a == "Original"){
      int temp = static_cast<int>(std::round(spectrumMass[i] * ppo + 0.1 * mng_ptr_->convert_ratio_));
      delta.push_back(temp);
      if(temp > maxDelta) maxDelta = temp;
    }
    else if(a == "Reversed") {
      int temp = static_cast<int>(std::round(spectrumMass[i] * ppo + totalMass * mng_ptr_->convert_ratio_ * ppo + 0.1 * mng_ptr_->convert_ratio_));
      delta.push_back(temp);
      if(temp > maxDelta) maxDelta = temp;
    }    
  }
}

void GraphAlignSim::getDelta_ori(){

  spectrumMass.clear();
  auto aa = spec_graph_ptr_->getPeakPtrVec();
  for(auto peakIter = aa.begin(); peakIter != aa.end(); peakIter++){
    int mass = static_cast<int>(std::round((*peakIter).first->getPosition() * mng_ptr_->convert_ratio_));
    spectrumMass.push_back(mass);
  }

  delta.clear();
  delta.push_back(0);
  for(int i = 1; i < spec_ver_num_; i++){
    delta.push_back(27);
  }
}


void GraphAlignSim::deleteOverlap(){
  deltaL.clear();
  deltaR.clear();
  deltaL.push_back(0);
  deltaR.push_back(0);
  for(int i = 1; i < spec_ver_num_; i++){
    if((spectrumMass[i-1] + deltaR[i-1]) < (spectrumMass[i] - delta[i])) {
      deltaL.push_back(delta[i]);
      deltaR.push_back(delta[i]);
      continue;
    }
    else{  //overlap
      if(spectrumMass[i] > (spectrumMass[i-1] + deltaR[i-1])){   //peak i out of the delta range of peak i-1;
        if(spectrumMass[i] - delta[i] < spectrumMass[i-1] - deltaL[i-1]){
          int iter = i-1;
          while(spectrumMass[iter] + deltaR[iter] > spectrumMass[i] - delta[i]){
            iter = iter - 1;
          }
          if(iter < i-1){
            deltaL[iter+1] = spectrumMass[iter+1] - std::max(spectrumMass[i]-delta[i], spectrumMass[iter]+deltaR[iter]+1);
          }
          else{
            deltaL[iter] = spectrumMass[iter] - (spectrumMass[i] - delta[i]);
          }
        }
        deltaL.push_back(std::max(0, spectrumMass[i]-(spectrumMass[i-1]+deltaR[i-1]+1)));
        deltaR.push_back(delta[i]);
        continue;

      }
      else{  // peak i in the delta range of peak i-1;
        if(spectrumMass[i]-delta[i] < spectrumMass[i-1] - deltaL[i-1]){
          int iter = i-1;
          int current  = spectrumMass[i] - delta[i];
          while(spectrumMass[iter] + deltaR[iter] > current){
            iter = iter - 1;
          }
          if(iter < i-1){
            if(spectrumMass[iter] + deltaR[iter] + 1 < spectrumMass[iter+1] - deltaL[iter+1]){
              deltaL[iter+1] = spectrumMass[iter+1] - std::min(spectrumMass[iter+1]-deltaL[iter+1], std::max(current, spectrumMass[iter]+deltaR[iter]+1));
            }
          }
          else{
            deltaL[iter] = spectrumMass[iter] - (spectrumMass[i] - delta[i]);
          }
        }
        int temp = deltaR[i-1];
        deltaR[i-1] = std::max(0, spectrumMass[i] - 1 - spectrumMass[i-1]);
        deltaL.push_back(0);
        deltaR.push_back(std::max(delta[i], spectrumMass[i-1] + temp - spectrumMass[i]));
        continue;             
      }
    }
  }
}


void GraphAlignSim::deleteOverlap_v2(){
  deltaL.clear();
  deltaR.clear();
  deltaL.push_back(0);
  deltaR.push_back(0);
  for(int i = 1; i < spec_ver_num_; i++){
    if((spectrumMass[i-1] + deltaR[i-1]) < (spectrumMass[i] - delta[i])) {
      deltaL.push_back(delta[i]);
      deltaR.push_back(delta[i]);
      continue;
    }
    else{  //overlap
      if(spectrumMass[i] > (spectrumMass[i-1] + deltaR[i-1])){   //peak i out of the delta range of peak i-1;
        if(spectrumMass[i] - delta[i] < spectrumMass[i-1] - deltaL[i-1]){
          int iter = i-1;
          while(spectrumMass[iter] + deltaR[iter] > spectrumMass[i] - delta[i]){
            iter = iter - 1;
          }
          if(iter < i-1){
            deltaL[iter+1] = spectrumMass[iter+1] - std::max(spectrumMass[i]-delta[i], spectrumMass[iter]+deltaR[iter]+1);
          }
          else{
            deltaL[iter] = spectrumMass[iter] - (spectrumMass[i] - delta[i]);
          }
        }
        deltaL.push_back(std::max(0, spectrumMass[i]-(spectrumMass[i-1]+deltaR[i-1]+1)));
        deltaR.push_back(delta[i]);
        continue;

      }
      else{  // peak i in the delta range of peak i-1;
        if(spectrumMass[i]-delta[i] < spectrumMass[i-1] - deltaL[i-1]){
          int iter = i-1;
          int current  = spectrumMass[i] - delta[i];
          while(spectrumMass[iter] + deltaR[iter] > current){
            iter = iter - 1;
          }
          if(iter < i-1){
            if(spectrumMass[iter] + deltaR[iter] + 1 < spectrumMass[iter+1] - deltaL[iter+1]){
              deltaL[iter+1] = spectrumMass[iter+1] - std::min(spectrumMass[iter+1]-deltaL[iter+1], std::max(current, spectrumMass[iter]+deltaR[iter]+1));
            }
          }
          else{
            deltaL[iter] = spectrumMass[iter] - (spectrumMass[i] - delta[i]);
          }
        }
        int temp = deltaR[i-1];
        deltaR[i-1] = std::max(0, spectrumMass[i] - 1 - spectrumMass[i-1]);
        deltaL.push_back(0);
        deltaR.push_back(std::max(delta[i], spectrumMass[i-1] + temp - spectrumMass[i]));
        continue;             
      }
    }
  }
}



void GraphAlignSim::getNewConsPair() {

  //int tole = maxDelta * 2;
  //std::cout << "maxDelta: " << maxDelta << std::endl;
  //int totalTole = 0;

  std::vector<std::pair<unsigned int, std::vector<std::pair<std::pair<unsigned short, unsigned short>,std::vector<std::pair<unsigned short,unsigned short>>>>>> empty_vec;
  for (int i = 0; i < proteo_ver_num_; i++) {
    std::vector<std::vector<std::pair<unsigned int, std::vector<std::pair<std::pair<unsigned short, unsigned short>,std::vector<std::pair<unsigned short,unsigned short>>>>>>> empty_vec_2d;
    for (int j = 0; j < spec_ver_num_; j++) {
      empty_vec_2d.push_back(empty_vec);
    }
    new_cons_pairs_.push_back(empty_vec_2d);
  }

  int min_dist = mng_ptr_->getIntMinConsistentDist();
  for (size_t m = 0; m < dist_vec_.size(); m++) {
    if (dist_vec_[m].size() == 0) continue;
    size_t spec_idx_min = 0, specidx = 0;
    size_t maxIdx = spec_dist_.size();
    for (size_t prot_idx = 0; prot_idx < dist_vec_[m].size(); prot_idx++){
      int pr_dist = dist_vec_[m][prot_idx].dist_;
      if (pr_dist < min_dist) continue;
      //(pr_dist == 321397 || pr_dist == 107033 || pr_dist == 267325 || pr_dist == 270616 || pr_dist == 515198)
      bool flag1 = true, flag2 = true;
      maxIdx = spec_dist_.size();
      specidx = spec_idx_min;
      //std::cout << "prot_idx: " << prot_idx << ", " << spec_idx_min << ", " << specidx << ", " << maxIdx << std::endl;
      while (specidx < spec_dist_.size() && specidx < maxIdx){
        auto specPair = spec_dist_[specidx].pair_ij_;
        int specDist = spec_dist_[specidx].dist_;
        //std::cout << "specidx: " << specidx << "; ( "<< specDist << ", " << pr_dist;
        int maxD = 0;
        int minD = 0;
        for(auto iter = specPair.begin(); iter != specPair.end(); iter++){
          int tempMax = deltaL[(*iter).first.first] + deltaR[(*iter).first.second];
          int tempMin = deltaR[(*iter).first.first] + deltaL[(*iter).first.second];
          if(tempMax > maxD) maxD = tempMax;
          if(tempMin > minD) minD = tempMin;
        }
        // if(specDist == 35949 && pr_dist == 35949){
        //   std::cout << "specPair: " << std::endl;
        //   for(int i11 = 0; i11 < specPair.size(); i11++){
        //     std::cout << "(" << specPair[i11].first.first << ", " << specPair[i11].first.second << "),";
        //   }
        //   std::cout << std::endl;
        //   std::cout << "protPair: " << std::endl;
        //   for(int i22 = 0; i22 < dist_vec_[m][prot_idx].pair_ij_.size(); i22++){
        //     std::cout << "(" << dist_vec_[m][prot_idx].pair_ij_[i22].first.first << ", " << dist_vec_[m][prot_idx].pair_ij_[i22].first.second << "),";
        //   }
        //   std::cout << std::endl;
        // }

        if(specDist <= pr_dist){
          if(specDist + maxD < pr_dist) {
            //std::cout << "no match1;" << std::endl;
            specidx++;
          }
          else{
            if(flag1 == true){
              spec_idx_min = specidx;
              flag1 = false;
            }
            //std::cout << "match1 with spec_idx_min: " << spec_idx_min << std::endl;
            addToNewConsistentPairs(pr_dist, specPair, dist_vec_[m][prot_idx].pair_ij_);
            specidx++;
          }
        }
        else{
          if(specDist - minD > pr_dist){
            if(flag2 == true){
              maxIdx = specidx + USER_DEFINE_PARAMETER;
              flag2 = false;
            }
            //std::cout << "no match2 with maxIdx: " << maxIdx << std::endl;
            specidx++;
          }
          else{
            //std::cout << "match2" << std::endl;
            addToNewConsistentPairs(pr_dist, specPair, dist_vec_[m][prot_idx].pair_ij_);
            flag2 = true;
            maxIdx++;
            specidx++;
          }
        }
        
      }
    }
    dist_vec_[m].clear();
  }
  dist_vec_.clear();
  LOG_DEBUG("consistent pair end");
  //std::cout << "End getting new Cons_Pairs." << std::endl;
}




void GraphAlignSim::addToNewConsistentPairs(int mass, const std::vector<std::pair<std::pair<int, int>, std::vector<std::pair<unsigned short,unsigned short>>>> & sp_pair_ij,
                                      const std::vector<std::pair<std::pair<int, int>, std::vector<std::pair<unsigned short,unsigned short>>>> & pg_pair_ij) {
  
  std::vector<std::pair<std::pair<unsigned short, unsigned short>,std::vector<std::pair<unsigned short,unsigned short>>>> prePairList;
  for (size_t k = 0; k < pg_pair_ij.size(); k++) {
    for (size_t sp = 0; sp < sp_pair_ij.size(); sp++) {
      prePairList.clear();
      short int pr_v1 = pg_pair_ij[k].first.first;
      short int pr_v2 = pg_pair_ij[k].first.second;
      short int sp_v1 = sp_pair_ij[sp].first.first;
      short int sp_v2 = sp_pair_ij[sp].first.second;
      int maxD = deltaL[sp_v1] + deltaR[sp_v2];
      int minD = deltaR[sp_v1] + deltaL[sp_v2];
      int specMass = spectrumMass[sp_v2] - spectrumMass[sp_v1];
      if(specMass - mass <= maxD || mass - specMass <= minD){
        short int cSize = new_cons_pairs_[pr_v2][sp_v2].size();
        if(cSize > 0){
          if(new_cons_pairs_[pr_v2][sp_v2][cSize - 1].first == mass){
            prePairList = new_cons_pairs_[pr_v2][sp_v2][cSize - 1].second;
            bool find = false;
            for(int iter = 0; iter < prePairList.size(); iter++){
              if(prePairList[iter].first.first == pr_v1 && prePairList[iter].first.second == sp_v1){
                find = true;
              }
            }
            if(find == false){
              prePairList.insert(prePairList.begin(), std::make_pair(std::make_pair(pr_v1, sp_v1), pg_pair_ij[k].second));
              new_cons_pairs_[pr_v2][sp_v2][cSize - 1].second = prePairList;
            }
          } else{
            prePairList.emplace_back(std::make_pair(std::make_pair(pr_v1, sp_v1), pg_pair_ij[k].second));
            new_cons_pairs_[pr_v2][sp_v2].emplace_back(std::make_pair(mass, prePairList));
          }
        }else{
          prePairList.emplace_back(std::make_pair(std::make_pair(pr_v1, sp_v1), pg_pair_ij[k].second));
          new_cons_pairs_[pr_v2][sp_v2].emplace_back(std::make_pair(mass, prePairList));
        }
        // if(pr_v2 == 141 && sp_v2 == 18){
        //   std::cout << "new_cons_pairs_[pr_v2][sp_v2]:" << new_cons_pairs_[pr_v2][sp_v2].size()<<std::endl;
        // }
      }
      /*
      std::cout << "Update c[" << pr_v2 << "][" << sp_v2 << "] as {";
      auto li = new_cons_pairs_[pr_v2][sp_v2]; 
      for(int iter = 0; iter < li.size(); iter++){
        std::cout << "(" << li[iter].first << ",";
        auto l2 = li[iter].second;
        for(int i2 = 0; i2 < l2.size(); i2++){
          std::cout << "<" << l2[i2].first.first << "," << l2[i2].first.second << ",";
          auto l3 = l2[i2].second;
          for(int i3 = 0; i3 < l3.size(); i3++){
            std::cout << "$" << l3[i3].first << "," << l3[i3].second << "$,";
          }
          std::cout << ">,";
        }
        std::cout << "),";
      }
      std::cout << "}" << std::endl;
      */
    }  
  }
}




void GraphAlignSim:: computeT_v2(bool case1, double ptm_mass){ //whole spectrum and surfix of protein

  //**********Build and initialize T[i,j,k] and E[i,j,k]******************
  
  std::vector<std::vector<std::vector<short int>>> T(proteo_ver_num_);
  std::vector<std::vector<std::vector<std::vector<prePosition>>>> E(proteo_ver_num_);
  for(int i = 0; i < proteo_ver_num_; i++){
    T[i].resize(spec_ver_num_);
    E[i].resize(spec_ver_num_);
    for(int j = 0; j < spec_ver_num_; j++){
      int size_k = deltaL[j] + deltaR[j] + 1;
      T[i][j].resize(size_k);
      E[i][j].resize(size_k);
    }
  }

  //initialize T[i,j,k] to be 1;
  for(size_t i = 0; i < T.size(); i++){
    for(size_t j = 0; j < T[i].size(); j++){
      for(size_t k = 0; k < T[i][j].size(); k++){
//        if(j == 0){
        if(j == 0 && i == 0 && k == 0){
          T[i][j][k] = 1;
        }
        else T[i][j][k] = -1;

      }
    }
  }


  //std::cout << "ini" << std::endl;

  for(int i = 0; i < proteo_ver_num_; i++){
    for(int j = 0; j < spec_ver_num_; j++){
      auto cij = new_cons_pairs_[i][j];
      for(auto m_iter = cij.begin(); m_iter != cij.end(); m_iter++){
        int current_pointer = 0;
        // Computing T[i,j,0];
        auto list = m_iter->second;
        unsigned int exactM = m_iter->first;
        int listSize = list.size();
        bool preFound = false;
        for(int q = 0; q < listSize; q++){
          unsigned short i_pre = list[q].first.first;
          unsigned short j_pre = list[q].first.second;
          auto mod_info = list[q].second;
          if(spectrumMass[j_pre] >= (spectrumMass[j] - deltaL[j] - exactM - deltaR[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] - deltaL[j] - exactM + deltaL[j_pre])){
            short int k_pre = spectrumMass[j] - deltaL[j] - exactM - spectrumMass[j_pre] + deltaL[j_pre];
            current_pointer = q;
            preFound = true;
            if(T[i_pre][j_pre][k_pre] + 1 > 0){
              if(T[i][j][0] == T[i_pre][j_pre][k_pre] + 1){
                E[i][j][0].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
              }
              if(T[i][j][0] < T[i_pre][j_pre][k_pre] + 1){
                T[i][j][0] = T[i_pre][j_pre][k_pre] + 1;
//                E[i][j][0].clear();
                E[i][j][0].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
              }
            }
            break;
          }
        }
        //std::cout << "current i: " << i << " and j: " << j << std::endl;

        //Computing T[i,j,k] to T[i,j,deltaL[j] + deltaR[j]];
        for(int k = 1; k <= deltaL[j] + deltaR[j]; k++){
          bool update = false;
          short int offsent = k - deltaL[j];
          short int i_pre = list[current_pointer].first.first;
          short int j_pre = list[current_pointer].first.second;
          auto mod_info = list[current_pointer].second;
          if(preFound == true){
            if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - deltaR[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + deltaL[j_pre])){
              update = true;
              short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + deltaL[j_pre];
              //int i_star = E[i][j][k-1].first.first;
              //int j_star = E[i][j][k-1].first.second;
              if(T[i_pre][j_pre][k_pre] + 1 > 0){
                if(T[i][j][k] == T[i_pre][j_pre][k_pre] + 1){
                  E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                }
                if(T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                  T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
//                  E[i][j][k].clear();
                  E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                }
              }
//              if(i_star == list[current_pointer].first && j_star == list[current_pointer].second){
//                if(T[i][j][k-1] > T[i][j][k] && T[i][j][k-1] >= T[i_star][j_star][k_star + 1] + 1){
//                  T[i][j][k] = T[i][j][k-1];
//                  E[i][j][k] = std::make_pair(E[i][j][k-1].first, E[i][j][k-1].second + 1);
//                }
//                else if(T[i_star][j_star][k_star + 1] + 1 > T[i][j][k] && T[i_star][j_star][k_star + 1] + 1 > T[i][j][k-1]){
//                  T[i][j][k] = T[i_star][j_star][k_star + 1] + 1;
//                  E[i][j][k] = std::make_pair(E[i][j][k-1].first, E[i][j][k-1].second + 1);
//                }
//              }
            }
            else{
              int iter_point = current_pointer;
              while(iter_point <= listSize - 1){
                short int i_pre = list[iter_point].first.first;
                short int j_pre = list[iter_point].first.second;
                auto mod_info = list[iter_point].second;
                if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - deltaR[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + deltaL[j_pre])){
                  short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + deltaL[j_pre];
                  update = true;
                  current_pointer = iter_point;
                  if(T[i_pre][j_pre][k_pre] + 1 > 0){
                    if(T[i][j][k] == T[i_pre][j_pre][k_pre] + 1){
                      E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                    }
                    if(T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                      T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
//                      E[i][j][k].clear();
                      E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                    }
                  }
                  break;
                }
                iter_point++;
              }
            }
            if(current_pointer == listSize -1 && update == false && spectrumMass[list[current_pointer].first.second] < (spectrumMass[j] + offsent - exactM - deltaR[list[current_pointer].first.second])) break;  
            if(update == false) preFound = false;
          }
          else{
            int iter_point = current_pointer;
            while(iter_point <= listSize - 1){
                short int i_pre = list[iter_point].first.first;
                short int j_pre = list[iter_point].first.second;
                auto mod_info = list[iter_point].second;
              if(spectrumMass[j_pre] >= (spectrumMass[j] + offsent - exactM - deltaR[j_pre]) && spectrumMass[j_pre] <= (spectrumMass[j] + offsent - exactM + deltaL[j_pre])){
                short int k_pre = spectrumMass[j] + offsent - exactM - spectrumMass[j_pre] + deltaL[j_pre];
                preFound = true;
                current_pointer = iter_point;
                if(T[i_pre][j_pre][k_pre] + 1 > 0){
                  if(T[i][j][k] == T[i_pre][j_pre][k_pre] + 1){
                    E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                  }
                  if(T[i][j][k] < T[i_pre][j_pre][k_pre] + 1){
                    T[i][j][k] = T[i_pre][j_pre][k_pre] + 1;
//                    E[i][j][k].clear();
                    E[i][j][k].emplace_back(std::make_pair(std::make_pair(i_pre, j_pre), std::make_tuple(k_pre - deltaL[j_pre], exactM, mod_info)));
                  }
                }
                break;
              }
              iter_point++;
            }
          }
          /*
          if(T[i][j][k]>0){
            std::cout << "T[" << i << "," << j << "," << k << "] = " << T[i][j][k] << "; preNode: {";
            for(int piter = 0; piter < E[i][j][k].size(); piter++){
              std::cout << "[" << E[i][j][k][piter].first.first << ", " << E[i][j][k][piter].first.second << ", " << std::get<0>(E[i][j][k][piter].second) << "]," << std::get<1>(E[i][j][k][piter].second) << ";";
            }
            std::cout << "}" << std::endl;
          }
          */
          
        }
      }
    }
    new_cons_pairs_[i].clear();
    std::vector<std::vector<std::pair<unsigned int, std::vector<std::pair<std::pair<unsigned short, unsigned short>,std::vector<std::pair<unsigned short,unsigned short>>>>>>>().swap(new_cons_pairs_[i]);
  }
  new_cons_pairs_.clear();
  NewConsPairs().swap(new_cons_pairs_);

  Quantification(T,E);




//
//  //std::cout << "herereeee" << std::endl;
//
//  bool right = false;
//  bool left = false;
//  int rightPeak;
//  int leftPeak;
//  int modDelta = 27;
//  if(case1 == false){
//    int mass2 = spectrumMass[spectrumMass.size()-1] - static_cast<int>(std::round(ptm_mass * mng_ptr_->convert_ratio_));
//    //std::cout << mass2 << std::endl;
//    for(int backIter = spectrumMass.size()-1; backIter >= 0; backIter--){
//      if(spectrumMass[backIter] <= mass2 + modDelta && right == false && spectrumMass[backIter] >= mass2 - modDelta){
//        //std::cout << backIter << std::endl;
//        right = true;
//        rightPeak = backIter;
//        continue;
//      }
//      if(spectrumMass[backIter] < mass2 - modDelta && left == false){
//        left = true;
//        leftPeak = backIter + 1;
//        break;
//      }
//    }
//  }
//  //std::cout << spectrumMass[leftPeak-1] << ", " << spectrumMass[leftPeak] << "," << spectrumMass[leftPeak+1] << std::endl;
//  //std::cout << "(" << leftPeak << ", " << rightPeak << ")" << std::endl;
//
//  // std::cout << "test:---" << std::endl;
//  // std::cout << "peak mass 1-0 : " << spectrumMass[2] - spectrumMass[0] << ", " << (spectrumMass[2] - spectrumMass[0])/mng_ptr_->convert_ratio_ << std::endl;
//  // pg_;
//  // auto v11 = vertex(144, *pg_.get());
//  // auto v22 = vertex(145, *pg_.get());
//  // auto bedge = edge(v11, v22, *(pg_.get()));
//  // int mm = (*pg_)[bedge.first].int_mass_;
//  // std::cout << mm << "? " << std::endl;
//  // std::cout << new_cons_pairs_[141][18].size() << std::endl;
//  // int aa1 = 949;
//  // int aa2 = 20;
//  // auto mu = T[aa1][aa2];
//  // int big = 0;
//  // int nn = 0;
//  // for(int ff = 0; ff < mu.size(); ff++){
//  //   if(T[aa1][aa2][ff] > big){
//  //     big = T[aa1][aa2][ff];
//  //     nn = ff;
//  //   }
//  // }
//  // std::cout << "T: " << T[aa1][aa2][nn] << std::endl;
//  // std::cout << new_cons_pairs_[142][18].size() << std::endl;
//  // for(int qq = 0; qq < new_cons_pairs_[141][18].size(); qq++){
//  //   if(new_cons_pairs_[141][18][qq].first == 35949){
//  //     auto bbb = new_cons_pairs_[141][18][qq].second;
//  //     for(int pp = 0; pp < bbb.size(); pp++){
//  //       std::cout << bbb[pp].first.first << ", " << bbb[pp].first.second << std::endl;
//  //     }
//  //   }
//  // }
//
//
//
//  //***************protein: anywhere ~ end  &  spectrum: 0 ~ end***********************
//
//  /*
//  int bigT = 0;
//  int i_final = 0;
//  int j_final = 0;
//  int k_final = 0;
//  for(size_t k = 0; k < T[proteo_ver_num_ - 1][spec_ver_num_-1].size(); k++){
//    if(T[proteo_ver_num_ - 1][spec_ver_num_-1][k] > bigT){
//      bigT = T[proteo_ver_num_ - 1][spec_ver_num_-1][k];
//      i_final = proteo_ver_num_ - 1;
//      j_final = spec_ver_num_-1;
//      k_final = k;
//    }
//  }
//
//
//  std::cout << "The biggest T[i][j][k] is: T[" << i_final << "][" << j_final << "][" << k_final << "] = " << bigT << std::endl;
//  */
//
//  //***************protein: anywhere  &  spectrum: 0 ~ end***********************
//  int bigT = 0;
//  std::vector<std::tuple<int, int, int>> endNodes;
//  for(int iMax = 0; iMax < proteo_ver_num_; iMax++){
//    for(size_t k = 0; k < T[iMax][spec_ver_num_-1].size(); k++){
//      if(T[iMax][spec_ver_num_-1][k] >= bigT){
//        if(T[iMax][spec_ver_num_-1][k] == bigT){
//          endNodes.push_back(std::make_tuple(iMax, spec_ver_num_-1, k));
//        }
//        if(T[iMax][spec_ver_num_-1][k] > bigT){
//          bigT = T[iMax][spec_ver_num_-1][k];
//          endNodes.clear();
//          endNodes.push_back(std::make_tuple(iMax, spec_ver_num_-1, k));
//          //std::cout << "??: " << endNodes.size() << std::endl;
//        }
//      }
//    }
//  }
//  int bigT2 = 0;
//  std::vector<std::tuple<int, int, int>> endNodes2;
//  if(case1 == false){
//    if(right == true && left == true){
//      std::cout << "We can find case2 for this spectrum." << std::endl;
//      for(int iMax = 0; iMax < proteo_ver_num_; iMax++){
//        for(int pIter = leftPeak; pIter <= rightPeak; pIter++){
//          for(size_t k = 0; k < T[iMax][pIter].size(); k++){
//            if(T[iMax][pIter][k] >= bigT2){
//              if(T[iMax][pIter][k] == bigT2){
//                endNodes2.push_back(std::make_tuple(iMax, pIter, k));
//              }
//              if(T[iMax][pIter][k] > bigT2){
//                bigT2 = T[iMax][pIter][k];
//                endNodes2.clear();
//                endNodes2.push_back(std::make_tuple(iMax, pIter, k));
//              }
//            }
//          }
//        }
//      }
//    }
//    else if(right == false && left == true){
//      std::cout << "We can not find case2 for this spectrum." << std::endl;
//    }
//  }
//  std::cout << "The biggest T[i][j][k] is: " << std::max(bigT,bigT2) << std::endl;
//
//  //std::cout << "??? " << endNodes2.size() << std::endl;
//
//  //std::cout << endNodes.size() << std::endl;
//
//  //*****************************Backtrace********************************
//  if(endNodes.size() > 0){
//    short modification = 0;
//    int i_final = std::get<0>(endNodes[0]);
//    int j_final = std::get<1>(endNodes[0]);
//    int k_final = std::get<2>(endNodes[0]);
//    std::cout << "-----------------p: " << proteo_ver_num_ << "---------------------s: " << spec_ver_num_ << "-------------------------------- " << std::endl;
//    std::cout << "Protein_index  |  Peak_index  |  Peak_oriPosition  |  Peak_modPosition  |  Shifting" << std::endl;
//    std::cout << "       " << i_final << "      |      " << j_final << "      |      " << spectrumMass[j_final] << "     |     " << spectrumMass[j_final] + k_final - deltaL[j_final] << "      |      " << k_final - deltaL[j_final] << std::endl;
//    short int i = i_final;
//    short int j = j_final;
//    short int k = k_final;
//    while(T[i][j][k] > 1) {
//      short int i_star = E[i][j][k][0].first.first;
//      short int j_star = E[i][j][k][0].first.second;
//      short int k_value = std::get<0>(E[i][j][k][0].second);
//      short int k_star = k_value + deltaL[j_star];
//      //modification = modification + std::get<2>(E[i][j][k][0].second);
//      i = i_star;
//      j = j_star;
//      k = k_star;
//
//      std::cout << "       " << i << "      |      " << j << "      |     " << spectrumMass[j] << "     |     " << spectrumMass[j] + k - deltaL[j] << "      |      " << k - deltaL[j] << std::endl;
//    }
//
//
//  //std::cout << "The total modification of this alignment is: " << modification << std::endl;
//
//    for(int iter0 = endNodes.size()-1; iter0 > 0; iter0--){
//      for(int iterBack = iter0-1; iterBack >= 0; iterBack--){
//        if(std::get<0>(endNodes[iter0]) == std::get<0>(endNodes[iterBack]) && std::get<1>(endNodes[iter0]) == std::get<1>(endNodes[iterBack])){
//          endNodes.erase(endNodes.begin()+iter0);
//          break;
//        }
//      }
//    }
//  }
//  if(case1 == false){
//    for(int iter_case2 = endNodes2.size()-1; iter_case2 > 0; iter_case2--){
//      for(int iterBack2 = iter_case2-1; iterBack2 >= 0; iterBack2--){
//        if(std::get<0>(endNodes2[iter_case2]) == std::get<0>(endNodes2[iterBack2]) && std::get<1>(endNodes2[iter_case2]) == std::get<1>(endNodes2[iterBack2])){
//          endNodes2.erase(endNodes2.begin()+iter_case2);
//          break;
//        }
//      }
//    }
//  }
//
//
//  //*****************************Backtrace********************************
//
//  clock_t startTime, endTime;
//
//  startTime = clock();
//  //std::cout << "??? " << endNodes2.size() << std::endl;
//
//  if(case1 == true){
//    quantification_case1(endNodes, T, E);
//  }
//  else{
//    if(endNodes2.size() == 0){
//      quantification_case1(endNodes, T, E);
//    }
//    else{
//      quantification_case2(endNodes, endNodes2, T, E);
//    }
//  }
//
//
//
//  endTime = clock();
//  double runningTime = (double)(endTime - startTime) / CLOCKS_PER_SEC;
//  std::cout << "Running time of quantification is: " << runningTime * 1000 << " ms" << std::endl;
//
//
//  T.clear();
//  E.clear();
//  //spectrumMass.clear();
//  //delta.clear();
//  std::vector<std::vector<std::vector<short int>>>().swap(T);
//  std::vector<std::vector<std::vector<std::vector<prePosition>>>>().swap(E);
//  //std::vector<int>().swap(spectrumMass);
//  //std::vector<int>().swap(delta);
}

    static void outputProgressBar(float progress) {
        // 进度条的长度
        int barWidth = 50;

        // 计算已完成的进度比例
        int pos = progress * barWidth;

        // 输出进度条
        std::cout << "Progress [" << std::setw(barWidth) << std::setfill('=') << "" << "] " << std::fixed
                  << std::setprecision(1) << (progress * 100) << "%\r";
        std::cout.flush(); // 刷新输出缓冲区，确保立即显示

        // 如果进度已经达到100%，输出换行
        if (progress >= 1)
            std::cout << std::endl;
    }


void GraphAlignSim::Quantification(std::vector<std::vector<std::vector<short int>>> & T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> & E){
  auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple2>>();
  AlignmentGraphPtr alignGraph_ptr = std::make_shared<AlignmentGraph>();
  double max_inten = ConstructBacktrackingGraph(T, E, alignGraph_ptr);

  if(max_inten < 0){
      std::ofstream outFile;
      std::string result_file = mng_ptr_->align_result_file_;
      outFile.open(result_file, std::ios::app);
      outFile << "alignment not found" << std::endl;
      outFile.close();

      result_file = mng_ptr_->protoform_result_file_;
      outFile.open(result_file,std::ios::app);
      outFile << "alignment not found" << std::endl;
      outFile << std::endl;
      outFile.close();

      result_file = mng_ptr_->result_file_;
      outFile.open(result_file, std::ios::app);
      outFile  <<  std::endl;
      outFile.close();

  }else {

      std::cout << "backtrackinggraph built" << std::endl;

      std::unordered_map<int, std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>> pair_map;
      std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> Q2 = ConstructQ2(alignGraph_ptr, pair_map);

      double min_error = 999999999999;
      std::vector<std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>> best_B;
      double best_perc_a;
      double best_perc_b;
      double best_q1;
      double best_q2;
      double best_a;
      double best_b;
      for (int abund_a = 1; abund_a < 100; abund_a++) {
//    std::cout<<"abund_a = "<<abund_a<<std::endl;
          outputProgressBar(float(abund_a) / 100.0);
          for (double abund_b = 1; abund_b < 100; abund_b++) {
              double q1 = 0.01 * abund_a * max_inten;
              double q2 = 0.01 * abund_b * max_inten;
//      std::cout<<"abund_a = "<<abund_a<<std::endl;
              int ver_num = num_vertices(*alignGraph_ptr);
              std::vector<std::vector<double>> D(ver_num, std::vector<double>(ver_num, -1));
              std::vector<std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>> B(ver_num,
                                                                                  std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>(
                                                                                          ver_num));
//      if(abund_a==11 && abund_b==2) {
              FindTwoPath(D, B, Q2, alignGraph_ptr, pair_map, q1, q2);
//      }
//      std::cout<<"error: "<<D[0][0]<<" with abund a,b = "<<abund_a<<","<<abund_b<<std::endl;
              if ((abund_a == 1 && abund_b == 1) || D[0][0] < min_error) {
                  min_error = D[0][0];
                  best_B = B;
                  best_perc_a = abund_a / (abund_a + abund_b);
                  best_perc_b = abund_b / (abund_a + abund_b);
                  best_q1 = q1;
                  best_q2 = q2;
                  best_a = abund_a;
                  best_b = abund_b;
              }
              D.clear();
              B.clear();
          }
      }
      std::cout << std::endl << "min error: " << min_error << ", ab: (" << best_a << "," << best_b << "), abundance: ("
                << std::fixed << std::setprecision(4) << best_perc_a << "," << best_perc_b << "), intensity: ("
                << best_q1 << "," << best_q2 << ")" << std::endl;
      std::ofstream outFile;
       std::string result_file = mng_ptr_->align_result_file_;
      outFile.open(result_file, std::ios::app);
      outFile << "min error: " << min_error << ", abundance: (" << std::fixed << std::setprecision(4) << best_perc_a
              << "," << best_perc_b << "), intensity: (" << best_q1 << "," << best_q2 << ")" << std::endl;
      outFile.close();

      result_file = mng_ptr_->result_file_;
      outFile.open(result_file, std::ios::app);
      outFile <<  min_error << "\t";
      outFile.close();

      backtracking(best_B, E, alignGraph_ptr);

      result_file = mng_ptr_->result_file_;
      outFile.open(result_file, std::ios::app);
      outFile  << std::fixed << std::setprecision(4) << best_perc_a
              << "\t" << best_perc_b << "\t" << best_q1 << "\t" << best_q2 << std::endl;
      outFile.close();

  }
  std::cout<<"Quantification Finish"<<std::endl;
}



void GraphAlignSim::FindTwoPath(std::vector<std::vector<double>> & D,
                                std::vector<std::vector<std::pair<Vertex_AGraph,Vertex_AGraph>>> & B,
                                std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> & Q2,
                                AlignmentGraphPtr & alignGraph_ptr,
                                std::unordered_map<int, std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>> & pair_map,
                                double q1, double q2){

  int ver_num = num_vertices(*alignGraph_ptr);
  std::vector<std::vector<bool>> flag_table(ver_num, std::vector<bool>(ver_num, false));
  std::vector<std::vector<bool>> cut_table(ver_num, std::vector<bool>(ver_num, false));
  for(int n = Q2.size() - 1; n >= 0; n--){
//    std::cout<<"n="<<n<<std::endl;
    std::pair<Vertex_AGraph, Vertex_AGraph> node_pair = Q2[n];
    int j1 = (*alignGraph_ptr)[Q2[n].first].j_;
    int j2 = (*alignGraph_ptr)[Q2[n].second].j_;
    int i1 = (*alignGraph_ptr)[Q2[n].first].i_;
    int i2 = (*alignGraph_ptr)[Q2[n].second].i_;
    int idx1 = (*alignGraph_ptr)[Q2[n].first].idx_;
    int idx2 = (*alignGraph_ptr)[Q2[n].second].idx_;
    int pair_idx = idx1 * ver_num + idx2;

    if(flag_table[idx1][idx2]){
            continue;
    }

    if(i1 == proteo_ver_num_ - 1 && i2 == proteo_ver_num_ - 1){
      D[idx1][idx2] = 0;
//      if(j1 != j2){
//        D[idx1][idx2] = abs(q1 - spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity()) + abs(q2 - spec_graph_ptr_->getPeakPtrVec()[j2].first->getIntensity());
//      }else{
//        D[idx1][idx2] = abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1 - q2);
//      }
    }else{
      std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> pair_vec = pair_map[pair_idx];
      if(pair_vec.empty()){
        cut_table[idx1][idx2] = true;
        continue;
      }

      double min_error = 999999999999;
      double local_error;
      std::pair<Vertex_AGraph , Vertex_AGraph> source_pair;
      bool cur_cut = true;

      if(n == 0) {
          for (int ii = 0; ii < pair_vec.size(); ii++) {
              Vertex_AGraph cur_node1 = pair_vec[ii].first;
              Vertex_AGraph cur_node2 = pair_vec[ii].second;
              int cur_idx1 = (*alignGraph_ptr)[cur_node1].idx_;
              int cur_idx2 = (*alignGraph_ptr)[cur_node2].idx_;
              if(!cut_table[cur_idx1][cur_idx2]) {
                  cur_cut = false;
                  local_error = D[cur_idx1][cur_idx2];
              }
              if (ii == 0 || local_error < min_error) {
                  min_error = local_error;
                  source_pair = pair_vec[ii];
              }
          }
      }else {
          if (i1 == proteo_ver_num_ - 1 && i2 < proteo_ver_num_ - 1) {
              for (int ii = 0; ii < pair_vec.size(); ii++) {
                  Vertex_AGraph cur_node = pair_vec[ii].second;
                  int cur_idx = (*alignGraph_ptr)[cur_node].idx_;
                  int cur_i = (*alignGraph_ptr)[cur_node].i_;
                  if (!cut_table[idx1][cur_idx]) {
                      cur_cut = false;
                      if (j1 != j2) {
                          local_error = D[idx1][cur_idx] +
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j2].first->getIntensity() - q2) +
                                        q2 * (cur_i - i2 - 1);
//                              0 * (cur_i - i2 - 1);
                      } else {
                          local_error = D[idx1][cur_idx] -
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1) +
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1 - q2) +
                                        q2 * (cur_i - i2 - 1);
//                              0 * (cur_i - i2 - 1);
                      }
                      if (ii == 0 || local_error < min_error) {
                          min_error = local_error;
                          source_pair = pair_vec[ii];
//              source_pair = std::make_pair((*alignGraph_ptr)[pair_vec[ii].first].idx_,
//                                           (*alignGraph_ptr)[pair_vec[ii].second].idx_);
                      }
                  }

              }
          }
          if (i1 < proteo_ver_num_ - 1 && i2 == proteo_ver_num_ - 1) {
              for (int ii = 0; ii < pair_vec.size(); ii++) {
                  Vertex_AGraph cur_node = pair_vec[ii].first;
                  int cur_idx = (*alignGraph_ptr)[cur_node].idx_;
                  int cur_i = (*alignGraph_ptr)[cur_node].i_;
                  if (!cut_table[cur_idx][idx2]) {
                      cur_cut = false;
                      if (j1 != j2) {
                          local_error = D[cur_idx][idx2] +
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1) +
                                        q1 * (cur_i - i1 - 1);
//                            0 * (cur_i - i1 - 1);
                      } else {
                          local_error = D[cur_idx][idx2] -
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j2].first->getIntensity() - q2) +
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1 - q2) +
                                        q1 * (cur_i - i1 - 1);
//                                0 * (cur_i - i1 - 1);
                      }
                      if (ii == 0 || local_error < min_error) {
                          min_error = local_error;
                          source_pair = pair_vec[ii];
//              source_pair = std::make_pair((*alignGraph_ptr)[pair_vec[ii].first].idx_,
//                                           (*alignGraph_ptr)[pair_vec[ii].second].idx_);
                      }
                  }
              }
          }
          if (i1 < proteo_ver_num_ - 1 && i2 < proteo_ver_num_ - 1) {
              for (int ii = 0; ii < pair_vec.size(); ii++) {
                  Vertex_AGraph cur_node1 = pair_vec[ii].first;
                  Vertex_AGraph cur_node2 = pair_vec[ii].second;
                  int cur_idx1 = (*alignGraph_ptr)[cur_node1].idx_;
                  int cur_idx2 = (*alignGraph_ptr)[cur_node2].idx_;
                  int cur_i1 = (*alignGraph_ptr)[cur_node1].i_;
                  int cur_i2 = (*alignGraph_ptr)[cur_node2].i_;
                  if (!cut_table[cur_idx1][cur_idx2]) {
//              std::cout<<ii<<std::endl;
                      cur_cut = false;


                      if (j1 == j2) {
                          local_error = D[cur_idx1][cur_idx2] +
                                        abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1 - q2) +
                                        q1 * (cur_i1 - i1 - 1) + q2 * (cur_i2 - i2 - 1);
//                            0 * (cur_i1 - i1 - 1) + 0 * (cur_i2 - i2 - 1);

//            std::cout<<local_error<<std::endl;
                      } else {
                          if (j1 < j2) {
                              local_error = D[cur_idx1][idx2] +
                                            abs(spec_graph_ptr_->getPeakPtrVec()[j1].first->getIntensity() - q1) +
                                            q1 * (cur_i1 - i1 - 1);
//                              0 * (cur_i1 - i1 - 1);

                          } else {
                              local_error = D[idx1][cur_idx2] +
                                            abs(spec_graph_ptr_->getPeakPtrVec()[j2].first->getIntensity() - q2) +
                                            q2 * (cur_i2 - i2 - 1);
//                                0 * (cur_i2 - i2 - 1);
                          }
                      }

                      if (ii == 0 || local_error < min_error) {
                          min_error = local_error;
                          source_pair = pair_vec[ii];
                      }
                  }
              }
          }
      }
      cut_table[idx1][idx2] = cur_cut;
      D[idx1][idx2] = min_error;
      B[idx1][idx2] = source_pair;
    }
    flag_table[idx1][idx2] = true;
  }

}

void GraphAlignSim::backtracking(std::vector<std::vector<std::pair<Vertex_AGraph,Vertex_AGraph>>> & B, std::vector<std::vector<std::vector<std::vector<prePosition>>>> & E, AlignmentGraphPtr & alignGraph_ptr){
    std::vector<Vertex_AGraph> path1;
    std::vector<Vertex_AGraph> path2;

    std::pair<int,int> cur_pair = B[0][0];
    while(!((*alignGraph_ptr)[cur_pair.first].i_ == proteo_ver_num_ - 1 && (*alignGraph_ptr)[cur_pair.second].i_ == proteo_ver_num_ - 1)){
        path1.push_back(cur_pair.first);
        path2.push_back(cur_pair.second);
        int idx1 = (*alignGraph_ptr)[cur_pair.first].idx_;
        int idx2 = (*alignGraph_ptr)[cur_pair.second].idx_;
        std::pair<Vertex_AGraph, Vertex_AGraph> source_pair = B[idx1][idx2];
//        std::pair<Vertex_AGraph,Vertex_AGraph> source_pair = B[cur_pair.first][cur_pair.second];
        cur_pair = source_pair;
    }
    path1.push_back(cur_pair.first);
    path2.push_back(cur_pair.second);

//    std::cout<<"======1251624124="<<std::endl;

    std::ofstream outFile;
    std::string result_file = mng_ptr_->align_result_file_;
    outFile.open(result_file,std::ios::app);
    outFile<<"===Path1==="<<std::endl;
    outFile<<std::setw(5)<< "Peak" << std::setw(5) << "Node" << "\t" << "Mods" << std::endl;
    outFile.close();
    rebuildPath(path1, E, alignGraph_ptr);

    outFile.open(result_file,std::ios::app);
    outFile<<"===Path2==="<<std::endl;
    outFile<<std::setw(5)<< "Peak" << std::setw(5) << "Node" << "\t" << "Mods" << std::endl;
    outFile.close();
    rebuildPath(path2, E, alignGraph_ptr);

}

void GraphAlignSim::rebuildPath(std::vector<Vertex_AGraph> & path, std::vector<std::vector<std::vector<std::vector<prePosition>>>> & E, AlignmentGraphPtr & alignGraph_ptr){
  std::ofstream outFile;
  std::string result_file = mng_ptr_->align_result_file_;
  outFile.open(result_file,std::ios::app);
  std::vector<std::pair<unsigned short, unsigned short>> total_mods;
  for(int n = 0; n < path.size(); n ++) {
    int pre_i, pre_j, pre_k;
    if (n == 0) {
      pre_i = 0;
      pre_j = 0;
      pre_k = 0;
    } else {
      pre_i = (*alignGraph_ptr)[path[n - 1]].i_;
      pre_j = (*alignGraph_ptr)[path[n - 1]].j_;
      pre_k = (*alignGraph_ptr)[path[n - 1]].k_;
    }
    int cur_i = (*alignGraph_ptr)[path[n]].i_;
    int cur_j = (*alignGraph_ptr)[path[n]].j_;
    int cur_k = (*alignGraph_ptr)[path[n]].k_;
    int cur_k_value = cur_k - deltaL[cur_j];
    std::vector<prePosition> preNodes = E[cur_i][cur_j][cur_k];
    prePosition prepos;
    for (int m = 0; m < preNodes.size(); m++) {
      int i_pre = preNodes[m].first.first;
      int j_pre = preNodes[m].first.second;
      int k_value = std::get<0>(preNodes[m].second);
      int k_pre = k_value + deltaL[j_pre];
      if (i_pre == pre_i && j_pre == pre_j && k_pre == pre_k) {
        prepos = preNodes[m];
        break;
      }
    }
    std::vector<std::pair<unsigned short, unsigned short>> modinfo = std::get<2>(prepos.second);
    outFile<<std::setw(5)<< cur_j << std::setw(5) << cur_i << "\t";
    for(int l = 0; l < modinfo.size(); l++){
      outFile << "<" << mng_ptr_->ptm_map_[modinfo[l].first]->getName() << "," << modinfo[l].second << ">";
      total_mods.push_back(modinfo[l]);
    }
    outFile << std::endl;
    if(cur_i == proteo_ver_num_ - 1){
      break;
    }
  }
  outFile.close();

  result_file = mng_ptr_->protoform_result_file_;
  outFile.open(result_file,std::ios::app);
  for(int l = 0; l < total_mods.size(); l++){
    outFile<< "<" << mng_ptr_->ptm_map_[total_mods[l].first]->getName() << "," << total_mods[l].second << ">";
  }
  outFile << std::endl;
  outFile.close();


    result_file = mng_ptr_->result_file_;
    outFile.open(result_file,std::ios::app);
    for(int l = 0; l < total_mods.size(); l++){
        outFile<< "<" << mng_ptr_->ptm_map_[total_mods[l].first]->getName()<<","<< total_mods[l].second << ">";
    }
    outFile << "\t";
    outFile.close();

}


std::vector<Vertex_AGraph> GraphAlignSim::getSourceNodes(Vertex_AGraph & v, AlignmentGraphPtr & alignGraph_ptr){
  std::vector<Vertex_AGraph> node_list;
  in_alignEdge_iter in_i, in_end;
  for(tie(in_i, in_end) = boost::in_edges(v, *alignGraph_ptr); in_i != in_end; ++in_i){
    Edge_AGraph e = *in_i;
    Vertex_AGraph source = boost::source(e, *alignGraph_ptr);
    node_list.emplace_back(source);
  }
  return node_list;
}




std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> GraphAlignSim::ConstructQ2(AlignmentGraphPtr & alignGraph_ptr, std::unordered_map<int, std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>>> & pair_map){
//  std::pair<vertex_iter, vertex_iter> vp;
//  Vertex_AGraph v_start;
//  for(vp = boost::vertices(*alignGraph_ptr); vp.first != vp.second; ++vp.first){
//      std::cout<<(*alignGraph_ptr)[*vp.first].i_<<std::endl;
//    (*alignGraph_ptr)[*vp.first].idx_ = (*alignGraph_ptr)[*vp.first].i_ * spec_ver_num_ + (*alignGraph_ptr)[*vp.first].j_;
//    if((*alignGraph_ptr)[*vp.first].i_ == 0 && (*alignGraph_ptr)[*vp.first].j_ == 0 && (*alignGraph_ptr)[*vp.first].k_ == 0){
//      v_start = *vp.first;
//      break;
//    }
//  }

  Vertex_AGraph v_start = vertex(0, *alignGraph_ptr.get());

//  std::cout<<(*alignGraph_ptr)[v_start].i_<<std::endl;

  std::queue<std::pair<Vertex_AGraph, Vertex_AGraph>> Q1;
  std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> Q2;
  Q1.push(std::make_pair(v_start,v_start));
  while(!Q1.empty()) {
    std::pair<Vertex_AGraph, Vertex_AGraph> cur_pair = Q1.front();
    std::vector<std::pair<Vertex_AGraph, Vertex_AGraph>> pair_vec;

    Q1.pop();
    Q2.push_back(cur_pair);
    std::vector<Vertex_AGraph> node_list_1 = getSourceNodes(cur_pair.first, alignGraph_ptr);
    std::vector<Vertex_AGraph> node_list_2 = getSourceNodes(cur_pair.second, alignGraph_ptr);

    if((!node_list_1.empty()) && (!node_list_2.empty())) {
      if ((*alignGraph_ptr)[cur_pair.first].j_ == (*alignGraph_ptr)[cur_pair.second].j_) {
        for (int x = 0; x < node_list_1.size(); x++) {
          for (int y = 0; y < node_list_2.size(); y++) {
            if (abs((*alignGraph_ptr)[node_list_1[x]].i_ - (*alignGraph_ptr)[node_list_2[y]].i_) <=
                mng_ptr_->max_head_diff) {
              Q1.push(std::make_pair(node_list_1[x], node_list_2[y]));
              pair_vec.push_back(std::make_pair(node_list_1[x], node_list_2[y]));
            }
          }
        }
      } else {
        if ((*alignGraph_ptr)[cur_pair.first].j_ > (*alignGraph_ptr)[cur_pair.second].j_) {
          for (int y = 0; y < node_list_2.size(); y++) {
            if (abs((*alignGraph_ptr)[cur_pair.first].i_ - (*alignGraph_ptr)[node_list_2[y]].i_) <=
                mng_ptr_->max_head_diff) {
              Q1.push(std::make_pair(cur_pair.first, node_list_2[y]));
              pair_vec.push_back(std::make_pair(cur_pair.first, node_list_2[y]));
            }
          }
        } else {
          for (int x = 0; x < node_list_1.size(); x++) {
            if (abs((*alignGraph_ptr)[node_list_1[x]].i_ - (*alignGraph_ptr)[cur_pair.second].i_) <=
                mng_ptr_->max_head_diff) {
              Q1.push(std::make_pair(node_list_1[x], cur_pair.second));
              pair_vec.push_back(std::make_pair(node_list_1[x], cur_pair.second));
            }
          }
        }
      }
    }else{
      if(node_list_1.empty() && (!node_list_2.empty())){
        for (int y = 0; y < node_list_2.size(); y++){
          Q1.push(std::make_pair(cur_pair.first, node_list_2[y]));
          pair_vec.push_back(std::make_pair(cur_pair.first, node_list_2[y]));
        }
      }
      if((!node_list_1.empty()) && node_list_2.empty()){
        for (int x = 0; x < node_list_1.size(); x++){
          Q1.push(std::make_pair(node_list_1[x], cur_pair.second));
          pair_vec.push_back(std::make_pair(node_list_1[x], cur_pair.second));
        }
      }
    }
    int pair_idx = (*alignGraph_ptr)[cur_pair.first].idx_ * num_vertices(*alignGraph_ptr) + (*alignGraph_ptr)[cur_pair.second].idx_;
    pair_map.insert(std::make_pair(pair_idx, pair_vec));
  }
  std::cout<<"Q2 built"<<std::endl;
  return Q2;
}


double GraphAlignSim::ConstructBacktrackingGraph(std::vector<std::vector<std::vector<short int>>> & T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> & E, AlignmentGraphPtr & alignGraph_ptr) {
  int i_n = proteo_ver_num_ - 1;
  int y_m = spectrumMass.size() - 1;
  std::vector<std::tuple<int, int, int>> endNodes;
//  for (int y = 0; y < spectrumMass.size(); y++)
    for (size_t k = 0; k < T[i_n][y_m].size(); k++) {
      if (T[i_n][y_m][k] >= mng_ptr_->alignment_thresh) {
        endNodes.push_back(std::make_tuple(i_n, y_m, k));
      }
    }

    if(endNodes.empty()){
        return -1;
    }

//  }
  auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int, int, int>, int, hashKey_tuple2>>();
  VertexInfo_AGraph v_0(T[0][0][0], 0, 0, 0, 0);
    int v_0_index = add_vertex(v_0, *alignGraph_ptr.get());
//    std::cout<<"v_index "<<v_0_index<<", mapindex"<<vertexMapPtr->size()<<std::endl;
  (*vertexMapPtr).insert(std::make_pair(std::make_tuple(0, 0, 0), vertexMapPtr->size()));



  std::queue<std::tuple<int, int, int>> node_queue;
//  int max_T = 0;
//  int best_node_idx = 0;
//  for (int i = 0; i < endNodes.size(); i++){
//    int i1 = std::get<0>(endNodes[i]);
//    int j1 = std::get<1>(endNodes[i]);
//    int k1 = std::get<2>(endNodes[i]);
//    if(T[i1][j1][k1] > max_T){
//      max_T = T[i1][j1][k1];
//      best_node_idx = i;
//    }
//  }
//  node_queue.push(endNodes[best_node_idx]);
//  int i1 = std::get<0>(endNodes[best_node_idx]);
//  int j1 = std::get<1>(endNodes[best_node_idx]);
//  int k1 = std::get<2>(endNodes[best_node_idx]);
//  VertexInfo_AGraph v1(T[i1][j1][k1], i1, j1, k1);
//  int v_index = add_vertex(v1, *alignGraph_ptr.get());
//  Vertex_AGraph end_v = vertex(v_index, *alignGraph_ptr.get());
//  (*alignGraph_ptr)[end_v].idx_ = v_index;
//  std::cout<<"v_index "<<v_index<<", mapindex"<<vertexMapPtr->size()<<std::endl;
//  (*vertexMapPtr).insert(std::make_pair(std::make_tuple(i1, j1, k1), vertexMapPtr->size()));

  for (int i = 0; i < endNodes.size(); i++) {
    node_queue.push(endNodes[i]);
    int i1 = std::get<0>(endNodes[i]);
    int j1 = std::get<1>(endNodes[i]);
    int k1 = std::get<2>(endNodes[i]);
    VertexInfo_AGraph v1(T[i1][j1][k1], i1, j1, k1);
      int v_index = add_vertex(v1, *alignGraph_ptr.get());
      Vertex_AGraph end_v = vertex(v_index, *alignGraph_ptr.get());
      (*alignGraph_ptr)[end_v].idx_ = v_index;
//      std::cout<<"v_index "<<v_index<<", mapindex"<<vertexMapPtr->size()<<std::endl;
    (*vertexMapPtr).insert(std::make_pair(std::make_tuple(i1, j1, k1), vertexMapPtr->size()));
  }

  double max_inten = BFS_max_inten(node_queue, T, E, alignGraph_ptr, vertexMapPtr);

//    //----------------print out the backtracking graph-------------------------------------
//
//    std::cout << "---Backtracking Graph---" << std::endl;
//    int node_num_ = num_vertices(*alignGraph_ptr.get());
//    std::cout << "The number of nodes is: " << node_num_ << std::endl;
//    for(int v = node_num_ - 1; v > 0; v--){
//
//      Vertex_AGraph vtx = vertex(v, *alignGraph_ptr.get());
//      std::cout<<"Node idx: "<< (*alignGraph_ptr)[vtx].idx_<< " T(" << (*alignGraph_ptr)[vtx].i_ << ", " << (*alignGraph_ptr)[vtx].j_ << ", " << (*alignGraph_ptr)[vtx].k_ << ") = " << (*alignGraph_ptr)[vtx].T_<<": ";
//      Edge_AGraph edge_next;
//      out_alignEdge_iter out_i, out_end;
//      for(tie(out_i, out_end) = boost::out_edges(vtx, *alignGraph_ptr); out_i != out_end; ++out_i){
//        Edge_AGraph e = *out_i;
//        auto modInfo = (*alignGraph_ptr)[e].ModInfo_;
//
//
//        Vertex_AGraph target = boost::target(e, *alignGraph_ptr);
//        std::cout<<"Node: "<< (*alignGraph_ptr)[target].idx_<<"with mod: ";
//          for(int index3 = 0; index3 < modInfo.size(); index3++){
//              std::cout << "[" << mng_ptr_->ptm_map_[modInfo[index3].first]->getName() << ", " << modInfo[index3].second << "], ";
//          }
//      }
//      std::cout<<std::endl;
//
//    }
//
//    //----------------print out the backtracking graph-------------------------------------

  return max_inten;
}



double GraphAlignSim::BFS_max_inten(std::queue<std::tuple<int,int,int>> node_queue,
                        std::vector<std::vector<std::vector<short int>>> T,
                        std::vector<std::vector<std::vector<std::vector<prePosition>>>> E,
                        AlignmentGraphPtr &alignGraph_ptr,
                        std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple2>> &vertexMapPtr){
  double max_inten = 0;
  while(!node_queue.empty()){
    auto node = node_queue.front();
    node_queue.pop();

    int i = std::get<0>(node);
    int j = std::get<1>(node);
    int k = std::get<2>(node);
    double cur_inten = (spec_graph_ptr_->getPeakPtrVec())[j].first->getIntensity();
    if(cur_inten > max_inten){
      max_inten = cur_inten;
    }
    if(T[i][j][k] > 1){
      auto preNodes = E[i][j][k];
      //std::cout << "1: " << preNodes.size() << std::endl;
      int ori_index = vertexMapPtr->find(node)->second;
      for(auto & preNode : preNodes){
        int i_pre = preNode.first.first;
        int j_pre = preNode.first.second;
        int k_value = std::get<0>(preNode.second);
        int k_pre = k_value + deltaL[j_pre];

        int v_index;
        if (vertexMapPtr->find(std::make_tuple(i_pre, j_pre, k_pre)) == vertexMapPtr->end())
        {
          VertexInfo_AGraph v_info(T[i_pre][j_pre][k_pre], i_pre, j_pre, k_pre);
          v_index = add_vertex(v_info, *alignGraph_ptr.get());
          Vertex_AGraph cur_v = vertex(v_index, *alignGraph_ptr.get());
          (*alignGraph_ptr)[cur_v].idx_ = v_index;

          vertexMapPtr->insert(std::make_pair(std::make_tuple(i_pre, j_pre, k_pre), v_index));
          node_queue.emplace(i_pre, j_pre, k_pre);
          //verIndex = verIndex + 1;
          //should return 'verIndex';
        }
        else
        {
          v_index = vertexMapPtr->find(std::make_tuple(i_pre, j_pre, k_pre))->second;
        }


        int exactMass = std::get<1>(preNode.second);
        int blackMass = proteo_graph_ptr_->getSeqMass(i_pre,i); //need to be confirmed
        auto modInfo = std::get<2>(preNode.second);
        //std::cout << i_pre << ", " << j_pre << ", " << k_pre << std::endl;

        Vertex_AGraph v_pre = vertex(v_index, *alignGraph_ptr.get());
        Vertex_AGraph v2 = vertex(ori_index, *alignGraph_ptr.get());
        EdgeInfo_AGraph edge_info(blackMass, exactMass, std::move(modInfo));
        add_edge(v2, v_pre, edge_info, *alignGraph_ptr.get());

      }
    }

  }
  return max_inten;
}



void GraphAlignSim::quantification_case1(std::vector<std::tuple<int, int, int>> endNodes, std::vector<std::vector<std::vector<short int>>> T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> E){
  std::cout << "Start case 1." << std::endl;
  int num = endNodes.size();
  double minimumError;
  bool firstCase = true;
  double abundancy_a, abundancy_b;
  std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> resultPath;
  //std::cout << "endNodes.size(): " << num << std::endl;
  //AlignmentGraphPtrVec aGraphVec;
  //std::vector<std::tuple<double, double, double, std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>>>> ResultVec;
  //AlignmentGraphPtr alignGraph_ptr = std::make_shared<AlignmentGraph>();
  //auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple>>();
  for(int num_iter = 0; num_iter < num; num_iter++){
    //std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple> vertexMap;
    auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple2>>();
    //int verIndex = 0;
    //std::cout << "ccc " << std::endl;
    AlignmentGraphPtr alignGraph_ptr = std::make_shared<AlignmentGraph>();
    std::queue<std::tuple<int,int,int>> node_queue;
    std::tuple<int, int, int> endNode = endNodes[num_iter];
    int i = std::get<0>(endNode);
    int j = std::get<1>(endNode);
    int k = std::get<2>(endNode);
    node_queue.push(endNode);
    //add_all_info(i, j, k, T, E, alignGraph_ptr, vertexMapPtr, 0, 0);
    VertexInfo_AGraph v(T[i][j][k], i, j, k);
    add_vertex(v, *alignGraph_ptr.get());
    (*vertexMapPtr).insert(std::make_pair(std::make_tuple(i,j,k), vertexMapPtr->size()));
    BFS(node_queue, T, E, alignGraph_ptr, vertexMapPtr);
    //int ori_index = verIndex;
    //verIndex = verIndex + 1;
//    if(T[i][j][k] > 1){
//      auto preNodes = E[i][j][k];
//      //std::cout << "1: " << preNodes.size() << std::endl;
//      for(int preNode_iter = 0; preNode_iter < preNodes.size(); preNode_iter++){
//        int i_pre = preNodes[preNode_iter].first.first;
//        int j_pre = preNodes[preNode_iter].first.second;
//        int k_value = std::get<0>(preNodes[preNode_iter].second);
//        int k_pre = k_value + deltaL[j_pre];
//        int exactMass = std::get<1>(preNodes[preNode_iter].second);
//        int blackMass = proteo_graph_ptr_->getSeqMass(i_pre,i); //need to be confirmed
//        auto modInfo = std::get<2>(preNodes[preNode_iter].second);
//        //std::cout << i_pre << ", " << j_pre << ", " << k_pre << std::endl;
//        add_all_info(i_pre, j_pre, k_pre, T, E, alignGraph_ptr, vertexMapPtr, 0, exactMass, blackMass, modInfo);
//
//      }

    //aGraphVec.push_back(alignGraph_ptrt);

    //----------------print out the backtracking graph-------------------------------------
    /*
    std::cout << "---Backtracking Graph---" << std::endl;
    int node_num_ = num_vertices(*alignGraph_ptr.get());
    std::cout << "The number of nodes is: " << node_num_ << std::endl;
    Vertex_AGraph v_start = *boost::vertices(*alignGraph_ptr).first;
    Vertex_AGraph v_next = v_start;
    for(int index = 0; index < node_num_; index++){
      std::cout << "T(" << (*alignGraph_ptr)[v_next].i_ << ", " << (*alignGraph_ptr)[v_next].j_ << ", " << (*alignGraph_ptr)[v_next].k_ << ") = " << (*alignGraph_ptr)[v_next].T_ << std::endl;
      std::pair<out_alignEdge_iter, out_alignEdge_iter> outEdge = boost::out_edges(v_next, *alignGraph_ptr);
      out_alignEdge_iter edgeIter_start = outEdge.first;
      out_alignEdge_iter edgeIter_next = edgeIter_start;
      int outEdgeNum = boost::out_degree(v_next, *alignGraph_ptr);
      Edge_AGraph edge_next;
      for(int index2 = 0; index2 < outEdgeNum; index2++){
        edge_next = *edgeIter_next;
        Vertex_AGraph target = boost::target(edge_next, *alignGraph_ptr);
        std::cout << "   Edge " << index2+1 << ": to T(" << (*alignGraph_ptr)[target].i_ << ", " << (*alignGraph_ptr)[target].j_ << ", " << (*alignGraph_ptr)[target].k_ << ") = " << (*alignGraph_ptr)[target].T_ << "; | ";
        std::cout << "Exact mass is " << (*alignGraph_ptr)[edge_next].exactMass_ << "; | Black mass is " << (*alignGraph_ptr)[edge_next].blackMass_ << "; | Mod infor are ";
        auto modInfo = (*alignGraph_ptr)[edge_next].ModInfo_;
        int sum = 0;
        for(int index3 = 0; index3 < modInfo.size(); index3++){
          std::cout << "[" << modInfo[index3].first << ", " << modInfo[index3].second << "], ";
          sum = sum + std::round(modInfo[index3].first * mng_ptr_->convert_ratio_);
        }
        std::cout << "; | Exact mass - Black mass = " << (*alignGraph_ptr)[edge_next].exactMass_ - (*alignGraph_ptr)[edge_next].blackMass_ << "; | Mod Sum = " << sum << std::endl;

        ++edgeIter_next;
      }
      ++v_next;
    }
    */
    //----------------print out the backtracking graph-------------------------------------
    
    //AlignmentGraphPtr alignGraph_ptr = aGraphVec[u];
    int node_num_ = num_vertices(*alignGraph_ptr.get());
    Vertex_AGraph v_start = *boost::vertices(*alignGraph_ptr).first;
    if(node_num_ == (*alignGraph_ptr)[v_start].T_){
      std::cout << "For the alignment graph end at T[" << (*alignGraph_ptr)[v_start].i_ << "][" << (*alignGraph_ptr)[v_start].j_ << "][" << (*alignGraph_ptr)[v_start].k_ << "] = " << (*alignGraph_ptr)[v_start].T_ << ", there is only one path." << std::endl;
      continue;
    }
    // If there exists more than one path
    else if(node_num_ > (*alignGraph_ptr)[v_start].T_){
      //double minimumError = -1.0;
      //bool firstCase = true;
      //double abundancy_a, abundancy_b;
      Vertex_AGraph beginNode = v_start;
      //std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> resultPath;
      auto res = findLargestIntensity(alignGraph_ptr,beginNode);
      double biggestInten = res.first;
      auto layers = res.second;
      for(double abund_a = 0.01; abund_a < 1.0; abund_a = abund_a + 0.01){
        for(double abund_b = 0.01; abund_b < 1.0; abund_b = abund_b + 0.01){
          //std::cout << abund_a << ", " << abund_b << std::endl;
          std::vector<std::vector<std::vector<std::pair<double, PrePathInfo>>>> D;
          std::vector<std::vector<std::pair<double, PrePathInfo>>> path_a;
          std::vector<std::pair<double, PrePathInfo>> path_b;
          //std::vector<Vertex_AGraph> pre_layer;
          
          //double inten_a, inten_b;
          //double inten_totle = (spec_graph_ptr_->getPrmPeakPtrVec())[(*(aGraphVec[graph_iter]))[v_start].j_]->getBasePeakPtr()->getIntensity();
          
          //std::cout << "The end peak intensity: " << inten_totle << std::endl;
          //inten_a = inten_totle * abund_a;
          //inten_b = inten_totle * abund_b;
          //std::cout << "inten_a: " << inten_a << ", b: " << inten_b << std::endl; 
          auto endInten = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[beginNode].j_].first->getIntensity();
          double endError = abs(endInten - abund_a * biggestInten - abund_b * biggestInten);
          path_b.emplace_back(std::make_pair(endError, std::make_pair(std::make_pair(0,beginNode),std::make_pair(0,beginNode))));
          path_a.emplace_back(path_b);
          D.emplace_back(path_a);
          //pre_layer.emplace_back(beginNode);
          for(int layer = 1; layer < layers.size(); layer++){
            //std::cout << "D[layer-1]=D[" << layer << "-1].size = " << D[layer-1].size() << std::endl;
            auto result = countError(alignGraph_ptr, layers, layer, D[layer-1], abund_a * biggestInten, abund_b * biggestInten);
            //pre_layer.clear();
            //pre_layer = result.first;
            D.emplace_back(result);
          }
          

          //Find two paths with the minimum peak intensity errors
          //std::cout << "1111" << std::endl;
          //double minError = -1;
          //bool first = true;
          PrePathInfo endPath, prePath;
          int size = D.size();
          //std::cout << "size: " << D[size-1].size() << std::endl;
          for(int iter_a = 0; iter_a < D[size-1].size(); iter_a++){
            for(int iter_b = 0; iter_b < D[size-1][iter_a].size(); iter_b++){
              if(firstCase == true){
                abundancy_a = abund_a * biggestInten;
                abundancy_b = abund_b * biggestInten;
                minimumError = D[size-1][iter_a][iter_b].first;
                prePath = D[size-1][iter_a][iter_b].second;
                endPath = std::make_pair(std::make_pair(iter_a, layers[layers.size()-1][iter_a]), std::make_pair(iter_b, layers[layers.size()-1][iter_b]));
                firstCase = false;
                resultPath.clear();
                //std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> returnValue;
                std::vector<std::pair<unsigned short, unsigned short>> info1, info2;
                resultPath.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[endPath.first.second].i_, (*alignGraph_ptr)[endPath.first.second].j_, (*alignGraph_ptr)[endPath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[endPath.second.second].i_,(*alignGraph_ptr)[endPath.second.second].j_, (*alignGraph_ptr)[endPath.second.second].k_, info2)));
                for(int p = D.size()-2; p >= 0; p--){
                  //std::cout << "    " << (*(aGraphVec[graph_iter]))[prePath.first.second].j_ << "    |    " << (*(aGraphVec[graph_iter]))[prePath.second.second].j_ << std::endl;
                  info1.clear();
                  info2.clear();
                  auto a = edge(prePath.first.second, endPath.first.second, *(alignGraph_ptr.get()));
                  auto b = edge(prePath.second.second, endPath.second.second, *(alignGraph_ptr.get()));
                  info1 = (*alignGraph_ptr)[a.first].ModInfo_;
                  info2 = (*alignGraph_ptr)[b.first].ModInfo_;
                  //std::cout << "666" << std::endl;
                  resultPath.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[prePath.first.second].i_, (*alignGraph_ptr)[prePath.first.second].j_, (*alignGraph_ptr)[prePath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[prePath.second.second].i_, (*alignGraph_ptr)[prePath.second.second].j_, (*alignGraph_ptr)[prePath.second.second].k_, info2)));
                  //std::cout << "777" << std::endl;
                  int m = prePath.first.first;
                  int n = prePath.second.first;
                  endPath = prePath;                  
                  prePath = D[p][m][n].second;
                       
                }
              }
              else{
                if(minimumError > D[size-1][iter_a][iter_b].first){
                  abundancy_a = abund_a * biggestInten;
                  abundancy_b = abund_b * biggestInten;
                  minimumError = D[size-1][iter_a][iter_b].first;
                  prePath = D[size-1][iter_a][iter_b].second;
                  endPath = std::make_pair(std::make_pair(iter_a, layers[layers.size()-1][iter_a]), std::make_pair(iter_b, layers[layers.size()-1][iter_b]));
                  resultPath.clear();
                  //std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> returnValue;
                  std::vector<std::pair<unsigned short, unsigned short>> info1, info2;
                  resultPath.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[endPath.first.second].i_, (*alignGraph_ptr)[endPath.first.second].j_, (*alignGraph_ptr)[endPath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[endPath.second.second].i_,(*alignGraph_ptr)[endPath.second.second].j_, (*alignGraph_ptr)[endPath.second.second].k_, info2)));
                  
                  for(int p = D.size()-2; p >= 0; p--){
                    //std::cout << "    " << (*(aGraphVec[graph_iter]))[prePath.first.second].j_ << "    |    " << (*(aGraphVec[graph_iter]))[prePath.second.second].j_ << std::endl;
                    info1.clear();
                    info2.clear();
                    auto a = edge(prePath.first.second, endPath.first.second, *(alignGraph_ptr.get()));
                    auto b = edge(prePath.second.second, endPath.second.second, *(alignGraph_ptr.get()));
                    //std::cout << a.second << ", " << b.second << std::endl;
                    info1 = (*alignGraph_ptr)[a.first].ModInfo_;
                    info2 = (*alignGraph_ptr)[b.first].ModInfo_;
                    //std::cout << "666" << std::endl;
                    resultPath.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[prePath.first.second].i_, (*alignGraph_ptr)[prePath.first.second].j_, (*alignGraph_ptr)[prePath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[prePath.second.second].i_, (*alignGraph_ptr)[prePath.second.second].j_, (*alignGraph_ptr)[prePath.second.second].k_, info2)));
                    //std::cout << "777" << std::endl;
                    int m = prePath.first.first;
                    int n = prePath.second.first;
                    endPath = prePath;
                    //std::cout << "555" << std::endl;
                    //std::cout << "--------" << m << ", " << n << ", " << p << std::endl;
                  
                    prePath = D[p][m][n].second;
                         
                  }
                  
                }
              }
            }
          }
        }
      }
      //std::cout << "prePath.first.first: " << prePath.first.first << ", " << prePath.second.first << std::endl;     
    }
  }

  //----------------print out the unordered_map-------------------------------------
  // for(auto map_iter = vertexMapPtr->begin(); map_iter != vertexMapPtr->end(); map_iter++){
  //   std::cout << "[" << std::get<0>(map_iter->first) << ", " << std::get<1>(map_iter->first) << ", " << std::get<2>(map_iter->first) << "] = " << map_iter->second << std::endl;
  // }
  //----------------print out the unordered_map-------------------------------------
  //std::cout << "14543" << std::endl;
  if(resultPath.size()>0){
    int count = 0;
    double p_totle = abundancy_a + abundancy_b;
    double p1 = abundancy_a / p_totle;
    double p2 = abundancy_b / p_totle;
    std::cout << "The minError is: " << minimumError << " with the abundancy (" << p1 << ", " << p2 << ")." << std::endl;
    std::cout << "Two paths are: " << std::endl;
    std::cout << "        Path A           |          Path B        " << std::endl; 
    for(int iter = 0; iter < resultPath.size(); iter++){
      std::cout << "    peak:" << std::get<1>(resultPath[iter].first) << " (Res: " << std::get<0>(resultPath[iter].first) << ", Err: " << std::get<2>(resultPath[iter].first) << ")";
      auto mod1 = std::get<3>(resultPath[iter].first);
      std::cout << "    Mod: " << mod1.size() << "; ";
      for(int iter2 = 0; iter2 < mod1.size(); iter2++){
        std::cout << "<" << mod1[iter2].first << "," << mod1[iter2].second << ">";
      }
      std::cout << "     |peak:" << std::get<1>(resultPath[iter].second) << " (Res: " << std::get<0>(resultPath[iter].second) << ", Err: " << std::get<2>(resultPath[iter].second) << ")";
      auto mod2 = std::get<3>(resultPath[iter].second);
      std::cout << "    Mod: " << mod2.size() << "; ";
      for(int iter3 = 0; iter3 < mod2.size(); iter3++){
        std::cout << "<" << mod2[iter3].first << "," << mod2[iter3].second << ">";
      }
      std::cout << std::endl;
      if(std::get<1>(resultPath[iter].first) != std::get<1>(resultPath[iter].second) || std::get<0>(resultPath[iter].first) != std::get<0>(resultPath[iter].second)){
        count = count + 1;
      }
    }
    if(count > 0){
      std::cout << "This case reports two different paths. The number of different peaks is: " << count << std::endl;
    }
    else{
      std::cout << "This case reports two same paths." << std::endl;
    }
  }
  

  //std::cout << "For the alignment graph end at T[" << (*alignGraph_ptr)[v_start].i_ << "][" << (*alignGraph_ptr)[v_start].j_ << "][" << (*alignGraph_ptr)[v_start].k_ << "] = " << (*alignGraph_ptr)[v_start].T_ << ", we find two paths: " << std::endl;

  
}


void GraphAlignSim::quantification_case2(std::vector<std::tuple<int, int, int>> endNodes1, std::vector<std::tuple<int, int, int>> endNodes2, std::vector<std::vector<std::vector<short int>>> T, std::vector<std::vector<std::vector<std::vector<prePosition>>>> E){
  std::cout << "Start case 2." << std::endl;
  double minimumError;
  bool firstCase = true;
  double abundancy_a, abundancy_b;
  int T1_1_T2_2_T12_0;
  std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> resultPath;
  std::vector<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>> singlePath;
  
  int num1 = endNodes1.size();
  int num2 = endNodes2.size();
  //std::cout << "endNodes1.size(): " << num1 << std::endl;
  //std::cout << "endNodes2.size(): " << num2 << std::endl;
  //AlignmentGraphPtrVec aGraphVec;
  //std::vector<std::tuple<double, double, double, std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>>>> ResultVec;
  //AlignmentGraphPtr alignGraph_ptr = std::make_shared<AlignmentGraph>();
  //auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple>>();
  for(int num_iter1 = 0; num_iter1 < num1; num_iter1++){
    for(int num_iter2 = 0; num_iter2 < num2; num_iter2++){
      //std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple> vertexMap;
      auto vertexMapPtr = std::make_shared<std::unordered_map<std::tuple<int,int,int>, int, hashKey_tuple2>>();
      //int verIndex = 0;
      //std::cout << "ccc " << std::endl;
      AlignmentGraphPtr alignGraph_ptr = std::make_shared<AlignmentGraph>();
      std::queue<std::tuple<int,int,int>> node_queue;
      std::tuple<int, int, int> endNode1 = endNodes1[num_iter1];
      std::tuple<int, int, int> endNode2 = endNodes2[num_iter2];
      int i1 = std::get<0>(endNode1);
      int j1 = std::get<1>(endNode1);
      int k1 = std::get<2>(endNode1);
      node_queue.push(endNode1);
      //add_all_info(i, j, k, T, E, alignGraph_ptr, vertexMapPtr, 0, 0);
      VertexInfo_AGraph v1(T[i1][j1][k1], i1, j1, k1);
      add_vertex(v1, *alignGraph_ptr.get());
      (*vertexMapPtr).insert(std::make_pair(std::make_tuple(i1,j1,k1), vertexMapPtr->size()));
      int i2 = std::get<0>(endNode2);
      int j2 = std::get<1>(endNode2);
      int k2 = std::get<2>(endNode2);
      node_queue.push(endNode2);
      //add_all_info(i, j, k, T, E, alignGraph_ptr, vertexMapPtr, 0, 0);
      VertexInfo_AGraph v2(T[i2][j2][k2], i2, j2, k2);
      add_vertex(v2, *alignGraph_ptr.get());
      (*vertexMapPtr).insert(std::make_pair(std::make_tuple(i2,j2,k2), vertexMapPtr->size()));
      BFS(node_queue, T, E, alignGraph_ptr, vertexMapPtr);
      //int ori_index = verIndex;
      //verIndex = verIndex + 1;
      //std::cout << "ddd " << std::endl;
  //    if(T[i][j][k] > 1){
  //      auto preNodes = E[i][j][k];
  //      //std::cout << "1: " << preNodes.size() << std::endl;
  //      for(int preNode_iter = 0; preNode_iter < preNodes.size(); preNode_iter++){
  //        int i_pre = preNodes[preNode_iter].first.first;
  //        int j_pre = preNodes[preNode_iter].first.second;
  //        int k_value = std::get<0>(preNodes[preNode_iter].second);
  //        int k_pre = k_value + deltaL[j_pre];
  //        int exactMass = std::get<1>(preNodes[preNode_iter].second);
  //        int blackMass = proteo_graph_ptr_->getSeqMass(i_pre,i); //need to be confirmed
  //        auto modInfo = std::get<2>(preNodes[preNode_iter].second);
  //        //std::cout << i_pre << ", " << j_pre << ", " << k_pre << std::endl;
  //        add_all_info(i_pre, j_pre, k_pre, T, E, alignGraph_ptr, vertexMapPtr, 0, exactMass, blackMass, modInfo);
  //
  //      }

      //aGraphVec.push_back(alignGraph_ptrt);

      //----------------print out the backtracking graph-------------------------------------
      /*
      std::cout << "---Backtracking Graph---" << std::endl;
      int node_num_ = num_vertices(*alignGraph_ptr.get());
      std::cout << "The number of nodes is: " << node_num_ << std::endl;
      Vertex_AGraph v_start = *boost::vertices(*alignGraph_ptr).first;
      Vertex_AGraph v_next = v_start;
      for(int index = 0; index < node_num_; index++){
        std::cout << "T(" << (*alignGraph_ptr)[v_next].i_ << ", " << (*alignGraph_ptr)[v_next].j_ << ", " << (*alignGraph_ptr)[v_next].k_ << ") = " << (*alignGraph_ptr)[v_next].T_ << std::endl;
        std::pair<out_alignEdge_iter, out_alignEdge_iter> outEdge = boost::out_edges(v_next, *alignGraph_ptr);
        out_alignEdge_iter edgeIter_start = outEdge.first;
        out_alignEdge_iter edgeIter_next = edgeIter_start;
        int outEdgeNum = boost::out_degree(v_next, *alignGraph_ptr);
        Edge_AGraph edge_next;
        for(int index2 = 0; index2 < outEdgeNum; index2++){
          edge_next = *edgeIter_next;
          Vertex_AGraph target = boost::target(edge_next, *alignGraph_ptr);
          std::cout << "   Edge " << index2+1 << ": to T(" << (*alignGraph_ptr)[target].i_ << ", " << (*alignGraph_ptr)[target].j_ << ", " << (*alignGraph_ptr)[target].k_ << ") = " << (*alignGraph_ptr)[target].T_ << "; | ";
          std::cout << "Exact mass is " << (*alignGraph_ptr)[edge_next].exactMass_ << "; | Black mass is " << (*alignGraph_ptr)[edge_next].blackMass_ << "; | Mod infor are ";
          auto modInfo = (*alignGraph_ptr)[edge_next].ModInfo_;
          int sum = 0;
          for(int index3 = 0; index3 < modInfo.size(); index3++){
            std::cout << "[" << modInfo[index3].first << ", " << modInfo[index3].second << "], ";
            sum = sum + std::round(modInfo[index3].first * mng_ptr_->convert_ratio_);
          }
          std::cout << "; | Exact mass - Black mass = " << (*alignGraph_ptr)[edge_next].exactMass_ - (*alignGraph_ptr)[edge_next].blackMass_ << "; | Mod Sum = " << sum << std::endl;

          ++edgeIter_next;
        }
        ++v_next;
      }
      */
      //----------------print out the backtracking graph-------------------------------------
      
      //AlignmentGraphPtr alignGraph_ptr = aGraphVec[u];
      //int node_num_ = num_vertices(*alignGraph_ptr.get());
      //std::cout << "1" << std::endl;
      int v_index1 = vertexMapPtr->find(std::make_tuple(i1, j1, k1))->second;
      int v_index2 = vertexMapPtr->find(std::make_tuple(i2, j2, k2))->second;
      Vertex_AGraph beginNode1 = vertex(v_index1, *alignGraph_ptr.get());
      Vertex_AGraph beginNode2 = vertex(v_index2, *alignGraph_ptr.get());
      // if(node_num_ == (*alignGraph_ptr)[v_start].T_){
      //   std::cout << "For the alignment graph end at T[" << (*alignGraph_ptr)[v_start].i_ << "][" << (*alignGraph_ptr)[v_start].j_ << "][" << (*alignGraph_ptr)[v_start].k_ << "] = " << (*alignGraph_ptr)[v_start].T_ << ", there is only one path." << std::endl;
      //   continue;
      // }
      // If there exists more than one path


      //std::cout << "2" << std::endl;
      auto res1 = findLargestIntensity(alignGraph_ptr,beginNode1);
      auto res2 = findLargestIntensity(alignGraph_ptr,beginNode2);
      double biggestInten = std::max(res1.first, res2.first);
      auto layers1 = res1.second;
      auto layers2 = res2.second;
      int T1 = T[i1][j1][k1];
      int T2 = T[i2][j2][k2];
      //std::cout << "3" << std::endl;
      for(double abund_a = 0.01; abund_a < 1.0; abund_a = abund_a + 0.01){
        for(double abund_b = 0.01; abund_b < 1.0; abund_b = abund_b + 0.01){
          std::vector<std::vector<std::vector<std::pair<double, PrePathInfo>>>> D;
          std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> additionPath;
          //std::cout << T1 << ", " << T2 << std::endl;
          //double preError;
          if(T1 > T2){
            additionPath = buildAdditionPath(alignGraph_ptr, T1, T2, beginNode1, layers1, abund_a * biggestInten);
          
            //----------------------normal---------------------------------
            std::vector<std::vector<std::pair<double, PrePathInfo>>> index1;
            auto specialLayer = additionPath[additionPath.size()-1];
            for(int layer1Iter = 0; layer1Iter < layers1[T1-T2].size(); layer1Iter++){
              std::vector<std::pair<double, PrePathInfo>> index2;
              //auto NodeInf = std::make_pair(specialLayer[layer1Iter].first, std::make_pair(specialLayer[layer1Iter].second, specialLayer[layer1Iter].second));
              for(int layer2Iter = 0; layer2Iter < layers2[0].size(); layer2Iter++){
                double combinError;
                if(layers1[T1-T2][layer1Iter] == layers2[0][layer2Iter]){
                  combinError = specialLayer[layer1Iter].first + abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers2[0][layer2Iter]].j_].first->getIntensity()- abund_b * biggestInten - abund_a * biggestInten);
                }
                else{
                  combinError = specialLayer[layer1Iter].first + abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers2[0][layer2Iter]].j_].first->getIntensity()-abund_b * biggestInten);
                }
                //double combinError = specialLayer[layer1Iter].first + abs((spec_graph_ptr_->getPrmPeakPtrVec())[(*alignGraph_ptr)[layers2[0][layer2Iter]].j_]->getBasePeakPtr()->getIntensity()-abund_b * biggestInten);
                auto NodeInf = std::make_pair(combinError, std::make_pair(std::make_pair(specialLayer[layer1Iter].second.first, specialLayer[layer1Iter].second.second),std::make_pair(layer2Iter, layers2[0][layer2Iter])));
                index2.emplace_back(NodeInf);
              }
              index1.emplace_back(index2);
            }
            D.emplace_back(index1);

            for(int commonIter = 1; commonIter < T2; commonIter++){
              auto result = countError_case2(alignGraph_ptr, layers1, layers2, commonIter+(T1-T2), commonIter, D[commonIter-1], abund_a * biggestInten, abund_b * biggestInten);
              D.emplace_back(result);
            }
          }
          else if(T2 > T1){
            additionPath = buildAdditionPath(alignGraph_ptr, T2, T1, beginNode2, layers2, abund_b * biggestInten);
            //std::cout << "endF" << std::endl;
            std::vector<std::vector<std::pair<double, PrePathInfo>>> index1;
            auto specialLayer = additionPath[additionPath.size()-1];
            for(int layer1Iter = 0; layer1Iter < layers1[0].size(); layer1Iter++){
              std::vector<std::pair<double, PrePathInfo>> index2;
              for(int layer2Iter = 0; layer2Iter < layers2[T2-T1].size(); layer2Iter++){
                //std::cout << "ww: " << layers2[T2-T1].size() << std::endl;
                double combinError;
                if(layers2[T2-T1][layer2Iter] == layers1[0][layer1Iter]){
                  combinError = specialLayer[layer2Iter].first + abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers1[0][layer1Iter]].j_].first->getIntensity()- abund_b * biggestInten - abund_a * biggestInten);
                }
                else{
                  combinError = specialLayer[layer2Iter].first + abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers1[0][layer1Iter]].j_].first->getIntensity()-abund_a * biggestInten);
                }
                //double combinError = specialLayer[layer1Iter].first + abs((spec_graph_ptr_->getPrmPeakPtrVec())[(*alignGraph_ptr)[layers1[0][layer2Iter]].j_]->getBasePeakPtr()->getIntensity()-abund_b * biggestInten);
                auto NodeInf = std::make_pair(combinError, std::make_pair(std::make_pair(layer1Iter, layers1[0][layer1Iter]),std::make_pair(specialLayer[layer2Iter].second.first, specialLayer[layer2Iter].second.second)));
                index2.emplace_back(NodeInf);
              }
              index1.emplace_back(index2);
            }
            D.emplace_back(index1);
            //std::cout << "endF2" << std::endl;

            for(int commonIter = 1; commonIter < T1; commonIter++){
              auto result = countError_case2(alignGraph_ptr, layers1, layers2, commonIter, commonIter+(T2-T1), D[commonIter-1], abund_a * biggestInten, abund_b * biggestInten);
              D.emplace_back(result);
            }
          }
          else if(T2 == T1){
            std::vector<std::vector<std::pair<double, PrePathInfo>>> index1;
            //auto specialLayer = additionPath[additionPath.size()-1];
            for(int layer1Iter = 0; layer1Iter < layers1[0].size(); layer1Iter++){
              std::vector<std::pair<double, PrePathInfo>> index2;
              //auto NodeInf = std::make_pair(specialLayer[layer1Iter].first, std::make_pair(specialLayer[layer1Iter].second, specialLayer[layer1Iter].second));
              for(int layer2Iter = 0; layer2Iter < layers2[0].size(); layer2Iter++){
                double combinError;
                if(layers1[0][layer1Iter] == layers2[0][layer2Iter]){
                  combinError = abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers1[0][layer1Iter]].j_].first->getIntensity()- abund_b * biggestInten - abund_a * biggestInten);
                }
                else{
                  combinError = abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers1[0][layer1Iter]].j_].first->getIntensity()- abund_a * biggestInten) + abs((spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[layers2[0][layer2Iter]].j_].first->getIntensity()- abund_b * biggestInten);
                }
                auto NodeInf = std::make_pair(combinError, std::make_pair(std::make_pair(layer1Iter, layers1[0][layer1Iter]),std::make_pair(layer2Iter, layers2[0][layer2Iter])));
                index2.emplace_back(NodeInf);
              }
              index1.emplace_back(index2);
            }
            D.emplace_back(index1);

            for(int commonIter = 1; commonIter < T1; commonIter++){
              auto result = countError_case2(alignGraph_ptr, layers1, layers2, commonIter, commonIter, D[commonIter-1], abund_a * biggestInten, abund_b * biggestInten);
              D.emplace_back(result);
            }
          }

          
          //Find two paths with the minimum peak intensity errors
          //std::cout << "1111" << std::endl;
          //double minError = -1;
          //bool first = true;
          PrePathInfo endPath, prePath;
          int size = D.size();
          //std::cout << "size: " << D[size-1].size() << std::endl;
          for(int iter_a = 0; iter_a < D[size-1].size(); iter_a++){
            for(int iter_b = 0; iter_b < D[size-1][iter_a].size(); iter_b++){
              if(firstCase == true){
                
                abundancy_a = abund_a * biggestInten;
                abundancy_b = abund_b * biggestInten;
                minimumError = D[size-1][iter_a][iter_b].first;
                prePath = D[size-1][iter_a][iter_b].second;
                endPath = std::make_pair(std::make_pair(iter_a, layers1[layers1.size()-1][iter_a]), std::make_pair(iter_b, layers2[layers2.size()-1][iter_b]));
                firstCase = false;
                //std::cout << "here?" << std::endl;                
                if(T1>T2){
                  T1_1_T2_2_T12_0 = 1;
                  auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                  resultPath = twoResultPath.first;
                  singlePath = twoResultPath.second;
                }
                else if(T2>T1){
                  T1_1_T2_2_T12_0 = 2;
                  auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                  resultPath = twoResultPath.first;
                  singlePath = twoResultPath.second;
                  //std::cout << "here2?" << std::endl;
                }
                else if(T2 == T1){
                  T1_1_T2_2_T12_0 = 0;
                  auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                  resultPath = twoResultPath.first;
                  singlePath = twoResultPath.second;
                }
              }
              else{
                if(minimumError > D[size-1][iter_a][iter_b].first){
                  //std::cout << "here" << std::endl;
                  abundancy_a = abund_a * biggestInten;
                  abundancy_b = abund_b * biggestInten;
                  minimumError = D[size-1][iter_a][iter_b].first;
                  prePath = D[size-1][iter_a][iter_b].second;
                  endPath = std::make_pair(std::make_pair(iter_a, layers1[layers1.size()-1][iter_a]), std::make_pair(iter_b, layers2[layers2.size()-1][iter_b]));
                  if(T1>T2){
                    T1_1_T2_2_T12_0 = 1;
                    auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                    resultPath = twoResultPath.first;
                    singlePath = twoResultPath.second;
                  }
                  else if(T2>T1){
                    T1_1_T2_2_T12_0 = 2;
                    auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                    resultPath = twoResultPath.first;
                    singlePath = twoResultPath.second;
                  }
                  else if(T2 == T1){
                    T1_1_T2_2_T12_0 = 0;
                    auto twoResultPath = storePath(alignGraph_ptr, endPath, prePath, D, additionPath, T1_1_T2_2_T12_0);
                    resultPath = twoResultPath.first;
                    singlePath = twoResultPath.second;
                  }
                }
              }
            }
          }
        }
      }      
    }    
  }

  //----------------print out the unordered_map-------------------------------------
  // for(auto map_iter = vertexMapPtr->begin(); map_iter != vertexMapPtr->end(); map_iter++){
  //   std::cout << "[" << std::get<0>(map_iter->first) << ", " << std::get<1>(map_iter->first) << ", " << std::get<2>(map_iter->first) << "] = " << map_iter->second << std::endl;
  // }
  //----------------print out the unordered_map-------------------------------------

  int count = 0;
  double p_totle = abundancy_a + abundancy_b;
  double p1 = abundancy_a / p_totle;
  double p2 = abundancy_b / p_totle;
  std::cout << "T1_1_T2_2_T12_0 = " << T1_1_T2_2_T12_0 << std::endl;
  std::cout << "The minError is: " << minimumError << " with the abundancy (" << p1 << ", " << p2 << ")." << std::endl;
  std::cout << "Two paths are: " << std::endl;
  std::cout << "        Path A           |          Path B        " << std::endl; 
  for(int iter = 0; iter < resultPath.size(); iter++){
    std::cout << "    peak:" << std::get<1>(resultPath[iter].first) << " (Res: " << std::get<0>(resultPath[iter].first) << ", Err: " << std::get<2>(resultPath[iter].first) << ")";
    auto mod1 = std::get<3>(resultPath[iter].first);
    std::cout << "    Mod: " << mod1.size() << "; ";
    for(int iter2 = 0; iter2 < mod1.size(); iter2++){
      std::cout << "<" << mod1[iter2].first << "," << mod1[iter2].second << ">";
    }
    std::cout << "     |peak:" << std::get<1>(resultPath[iter].second) << " (Res: " << std::get<0>(resultPath[iter].second) << ", Err: " << std::get<2>(resultPath[iter].second) << ")";
    auto mod2 = std::get<3>(resultPath[iter].second);
    std::cout << "    Mod: " << mod2.size() << "; ";
    for(int iter3 = 0; iter3 < mod2.size(); iter3++){
      std::cout << "<" << mod2[iter3].first << "," << mod2[iter3].second << ">";
    }
    std::cout << std::endl;
    if(std::get<1>(resultPath[iter].first) != std::get<1>(resultPath[iter].second) || std::get<0>(resultPath[iter].first) != std::get<0>(resultPath[iter].second)){
      count = count + 1;
    }
  }
  for(int iter_s = 0; iter_s < singlePath.size(); iter_s++){
    if(T1_1_T2_2_T12_0 == 2){
      std::cout << "           /           |    peak: " << std::get<1>(singlePath[iter_s]) << " (Res: " << std::get<0>(singlePath[iter_s]) << ", Err: " << std::get<2>(singlePath[iter_s]) << ")";
      auto mod1 = std::get<3>(singlePath[iter_s]);
      std::cout << "    Mod: " << mod1.size() << "; ";
      for(int iter2 = 0; iter2 < mod1.size(); iter2++){
        std::cout << "<" << mod1[iter2].first << "," << mod1[iter2].second << ">";
      }
      std::cout << std::endl;
    }
    if(T1_1_T2_2_T12_0 == 1){
      std::cout << "    peak: " << std::get<1>(singlePath[iter_s]) << " (Res: " << std::get<0>(singlePath[iter_s]) << ", Err: " << std::get<2>(singlePath[iter_s]) << ")";
      auto mod1 = std::get<3>(singlePath[iter_s]);
      std::cout << "    Mod: " << mod1.size() << "; ";
      for(int iter2 = 0; iter2 < mod1.size(); iter2++){
        std::cout << "<" << mod1[iter2].first << "," << mod1[iter2].second << ">         |           /           ";
      }
      std::cout << std::endl;
    }
  }
  if(count > 0){
    std::cout << "This case reports two different paths. The number of different peaks is: " << count << std::endl;
  }
  else{
    std::cout << "This case reports two same paths." << std::endl;
  }
  //std::cout << "For the alignment graph end at T[" << (*alignGraph_ptr)[v_start].i_ << "][" << (*alignGraph_ptr)[v_start].j_ << "][" << (*alignGraph_ptr)[v_start].k_ << "] = " << (*alignGraph_ptr)[v_start].T_ << ", we find two paths: " << std::endl;

  
}

std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> GraphAlignSim::buildAdditionPath(AlignmentGraphPtr alignGraph_ptr, int T1, int T2, Vertex_AGraph beginNode1, std::vector<std::vector<Vertex_AGraph>> layers1, double a_inten){
  //std::cout << "aa" << std::endl;
  auto nodePointer1 = beginNode1;
  double capacity1 = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[beginNode1].j_].first->getIntensity();
  double preError = abs(capacity1 - a_inten);
  auto endInfor1 = std::make_pair(preError, std::make_pair(0, beginNode1));
  std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>> layer0;
  std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> returnResult;
  layer0.emplace_back(endInfor1);
  returnResult.emplace_back(layer0);
  //std::cout << "b" << std::endl;
  for(int pointer1 = T1-1; pointer1 >= T2; pointer1--){
    //std::cout << "l1" << std::endl;
    auto currentLayer1 = layers1[layers1.size()-pointer1];
    //std::cout << "l2" << std::endl;
    auto preLayer1 = layers1[layers1.size()-pointer1-1];
    //std::cout << "l3" << std::endl;
    
    //std::cout << "c" << std::endl;
    std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>> layerInfo;
    for(int ppointer1 = 0; ppointer1 < currentLayer1.size(); ppointer1++){
      //std::cout << "size: " << currentLayer1.size() << std::endl;
      bool FirstNodeinLayer = true;
      double capacity_temp = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[currentLayer1[ppointer1]].j_].first->getIntensity();
      double curError1 = abs(capacity_temp - a_inten);

      std::pair<in_alignEdge_iter, in_alignEdge_iter> inEdge1 = boost::in_edges(currentLayer1[ppointer1], *alignGraph_ptr);
      in_alignEdge_iter in_edgeIter1 = inEdge1.first;
      int inEdgeNum1 = boost::in_degree(currentLayer1[ppointer1], *alignGraph_ptr);
      Edge_AGraph in_edge1;
      //std::vector<int> path1;
      double minPreError;
      std::pair<double, std::pair<int, Vertex_AGraph>> NodeInfo;
      //std::cout << "d" << std::endl;
      for(int j = 0; j < inEdgeNum1; j++){
        
        //std::cout << "j " << inEdgeNum1 << std::endl;
        in_edge1 = *in_edgeIter1;
        Vertex_AGraph source = boost::source(in_edge1, *alignGraph_ptr);
        
        //std::cout << "outdegree: " << boost::out_degree(source, *alignGraph_ptr)<<std::endl;
        int q = find(preLayer1.begin(), preLayer1.end(), source) - preLayer1.begin();
        //std::cout << "q:" << preLayer1.size() << std::endl;
        if(q == preLayer1.size()){
          //std::cout << "e1" << std::endl;
          ++in_edgeIter1;
        }
        else{
          //std::cout << "e2" << std::endl;
          double preError = returnResult[layers1.size()-pointer1-1][q].first;
          if(FirstNodeinLayer == true){
            //std::cout << "e3" << std::endl;
            minPreError = preError;
            NodeInfo = std::make_pair(minPreError + curError1, std::make_pair(q, source));
            FirstNodeinLayer = false;
            //std::cout << "e3end" << std::endl;
          }
          else{
            //std::cout << "e4" << std::endl;
            if(minPreError > preError){
              //std::cout << "e5" << std::endl;
              minPreError = preError;
              NodeInfo = std::make_pair(minPreError + curError1, std::make_pair(q, source));
            }
          }
          //std::cout << "2222" << std::endl;
          ++in_edgeIter1; 
          //std::cout << "222224" << std::endl;
        }                        
      }
      layerInfo.emplace_back(NodeInfo);  
      //std::cout << "here1" << std::endl;              
    }
    returnResult.emplace_back(layerInfo);
    //std::cout << "here2" << std::endl; 
  }
  return returnResult;
}




std::pair<std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>>, std::vector<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> GraphAlignSim::storePath(AlignmentGraphPtr alignGraph_ptr, PrePathInfo endPath, PrePathInfo prePath, std::vector<std::vector<std::vector<std::pair<double, PrePathInfo>>>> D, std::vector<std::vector<std::pair<double, std::pair<int, Vertex_AGraph>>>> additionPath, int T1Larger){
  //abundancy_a = abund_a * biggestInten;
  //abundancy_b = abund_b * biggestInten;
  //resultPath.clear();
  std::vector<std::pair<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>, std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>>> returnValue;
  std::vector<std::pair<unsigned short, unsigned short>> info1, info2;
  returnValue.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[endPath.first.second].i_, (*alignGraph_ptr)[endPath.first.second].j_, (*alignGraph_ptr)[endPath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[endPath.second.second].i_,(*alignGraph_ptr)[endPath.second.second].j_, (*alignGraph_ptr)[endPath.second.second].k_, info2)));
  //std::cout << "222" << std::endl;
  for(int p = D.size()-2; p >= 0; p--){
    //std::cout << "    " << (*(aGraphVec[graph_iter]))[prePath.first.second].j_ << "    |    " << (*(aGraphVec[graph_iter]))[prePath.second.second].j_ << std::endl;
    info1.clear();
    info2.clear();
    auto a = edge(prePath.first.second, endPath.first.second, *(alignGraph_ptr.get()));
    auto b = edge(prePath.second.second, endPath.second.second, *(alignGraph_ptr.get()));
    //std::cout << "p: " << p << std::endl;
    //std::cout << a.second << ", " << b.second << std::endl;
    info1 = (*alignGraph_ptr)[a.first].ModInfo_;
    info2 = (*alignGraph_ptr)[b.first].ModInfo_;
    //std::cout << "666" << std::endl;
    returnValue.emplace_back(std::make_pair(std::make_tuple((*alignGraph_ptr)[prePath.first.second].i_, (*alignGraph_ptr)[prePath.first.second].j_, (*alignGraph_ptr)[prePath.first.second].k_, info1), std::make_tuple((*alignGraph_ptr)[prePath.second.second].i_, (*alignGraph_ptr)[prePath.second.second].j_, (*alignGraph_ptr)[prePath.second.second].k_, info2)));
    //std::cout << "777" << std::endl;
    int m = prePath.first.first;
    int n = prePath.second.first;
    endPath = prePath;
    //std::cout << "555" << std::endl;
    //std::cout << "--------" << m << ", " << n << ", " << p << std::endl;
    
    prePath = D[p][m][n].second;
    
        
  }

  std::vector<std::tuple<int, int, int, std::vector<std::pair<unsigned short, unsigned short>>>> singlePath;
  int position;
  Vertex_AGraph node;
  if(T1Larger == 1){
    position = endPath.first.first;
    node = endPath.first.second;
  }
  else if(T1Larger == 2){
    position = endPath.second.first;
    node = endPath.second.second;
  }
  //std::cout << "555" << std::endl;
  //int position = std::max(endPath.first.first, endPath.second.first);
  //auto preNode;
  int temp = position;
  auto tempNode = node;
  std::pair<int ,Vertex_AGraph> preNode;
  int inEdgeNum_a = boost::in_degree(tempNode, *alignGraph_ptr);
  //std::cout << "qq: " << inEdgeNum_a << std::endl;
  for(int iter = additionPath.size()-1; iter >= 1; iter--){
    preNode = additionPath[iter][temp].second;
    int outEdgeNum_a = boost::out_degree(preNode.second, *alignGraph_ptr);
    //std::cout << "pp: " << outEdgeNum_a << std::endl;
    //std::cout << "temp: " << temp << std::endl;
    info1.clear();
    //info2.clear();
    //std::cout << "iter:" << iter << std::endl;
    auto a = edge(preNode.second, tempNode, *(alignGraph_ptr.get()));
    //auto b = edge(prePath.second.second, endPath.second.second, *(alignGraph_ptr.get()));
    //std::cout << a.second << std::endl;
    info1 = (*alignGraph_ptr)[a.first].ModInfo_;
    //info2 = (*alignGraph_ptr)[b.first].ModInfo_;
    singlePath.emplace_back((*alignGraph_ptr)[preNode.second].i_, (*alignGraph_ptr)[preNode.second].j_, (*alignGraph_ptr)[preNode.second].k_, info1);
    tempNode = preNode.second;
    temp = preNode.first;
    //std::cout << "??" << std::endl;
  }
  return std::make_pair(returnValue, singlePath);
}





std::vector<std::vector<std::pair<double, PrePathInfo>>> GraphAlignSim::countError(AlignmentGraphPtr alignGraph_ptr, std::vector<std::vector<Vertex_AGraph>> layers, int layer, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, double a_inten, double b_inten){
  auto pre_layer = layers[layer-1];
  auto current_layer = layers[layer];
  //std::cout << "current size: " << current_layer.size() << ", " << pre_layer.size();
  //std::cout << "peak: " << (*alignGraph_ptr)[current_layer[0]].j_ << std::endl;
  //if(current_layer.size() == 2){
  //  std::cout << "peak: " << (*alignGraph_ptr)[current_layer[1]].j_ << std::endl;
  //}
  std::vector<std::vector<std::pair<double, PrePathInfo>>> errorVecVec;
  std::vector<double> current_inten;
    // double total_inten = 0.0;
    // for(int p = 0; p < current_layer.size(); p++) {
    //   double temp_inten = (spec_graph_ptr_->getPrmPeakPtrVec())[(*alignGraph_ptr)[current_layer[p]].j_]->getBasePeakPtr()->getIntensity();
    //   total_inten = total_inten + temp_inten;
    // }
  for(int a = 0; a < current_layer.size(); a++){
    std::vector<std::pair<double, PrePathInfo>> errorVec;
    for(int b = 0; b < current_layer.size(); b++){
      double error;
      //std::cout << a << ", " << b << std::endl;
      if(a == b){
        double capacity = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer[a]].j_].first->getIntensity();
        error = abs(capacity - a_inten - b_inten);
        //std::cout << "capacity: " << capacity << ", error: " << error << std::endl;
      }
      else if(a != b){
        double cap_a = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer[a]].j_].first->getIntensity();
        double cap_b = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer[b]].j_].first->getIntensity();
        error = abs(cap_a - a_inten) + abs(cap_b - b_inten);
        //std::cout << "error_a: " << cap_a << ", error - inten_a: " << abs(cap_a - abund_a) << std::endl;
        //std::cout << "error_b: " << cap_b << ", error - inten_b: " << abs(cap_b - (1-abund_a)) << std::endl;
      }
      //std::cout << "before" << std::endl;
      auto pre_position = getMinPreError(alignGraph_ptr, pre_layer, D_layer, current_layer[a], current_layer[b]);
      //std::cout << "after" << std::endl;
      double new_error = error + pre_position.first;
      //std::cout << "new_error: " << new_error << ", pre_error: " << pre_position.first << std::endl;
      errorVec.emplace_back(std::make_pair(new_error, pre_position.second));
    }
    errorVecVec.emplace_back(errorVec);
  }
  return errorVecVec;
}

std::vector<std::vector<std::pair<double, PrePathInfo>>> GraphAlignSim::countError_case2(AlignmentGraphPtr alignGraph_ptr, std::vector<std::vector<Vertex_AGraph>> layers1, std::vector<std::vector<Vertex_AGraph>> layers2, int layer1, int layer2, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, double a_inten, double b_inten){
  auto pre_layer1 = layers1[layer1-1];
  auto pre_layer2 = layers2[layer2-1];
  auto current_layer1 = layers1[layer1];
  auto current_layer2 = layers2[layer2];
  //std::cout << "current size: " << current_layer.size() << ", " << pre_layer.size();
  //std::cout << "peak: " << (*alignGraph_ptr)[current_layer[0]].j_ << std::endl;
  //if(current_layer.size() == 2){
  //  std::cout << "peak: " << (*alignGraph_ptr)[current_layer[1]].j_ << std::endl;
  //}
  std::vector<std::vector<std::pair<double, PrePathInfo>>> errorVecVec;
  std::vector<double> current_inten;
    // double total_inten = 0.0;
    // for(int p = 0; p < current_layer.size(); p++) {
    //   double temp_inten = (spec_graph_ptr_->getPrmPeakPtrVec())[(*alignGraph_ptr)[current_layer[p]].j_]->getBasePeakPtr()->getIntensity();
    //   total_inten = total_inten + temp_inten;
    // }

  for(int a = 0; a < current_layer1.size(); a++){
    std::vector<std::pair<double, PrePathInfo>> errorVec;
    for(int b = 0; b < current_layer2.size(); b++){
      double error;
      //std::cout << a << ", " << b << std::endl;
      if(current_layer1[a] == current_layer2[b]){
        double capacity = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer1[a]].j_].first->getIntensity();
        error = abs(capacity - a_inten - b_inten);
        //std::cout << "capacity: " << capacity << ", error: " << error << std::endl;
      }
      else if(current_layer1[a] != current_layer2[b]){
        double cap_a = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer1[a]].j_].first->getIntensity();
        double cap_b = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[current_layer2[b]].j_].first->getIntensity();
        error = abs(cap_a - a_inten) + abs(cap_b - b_inten);
        //std::cout << "error_a: " << cap_a << ", error - inten_a: " << abs(cap_a - abund_a) << std::endl;
        //std::cout << "error_b: " << cap_b << ", error - inten_b: " << abs(cap_b - (1-abund_a)) << std::endl;
      }
      //std::cout << "before" << std::endl;
      auto pre_position = getMinPreError_case2(alignGraph_ptr, pre_layer1, pre_layer2, D_layer, current_layer1[a], current_layer2[b]);
      //std::cout << "after" << std::endl;
      double new_error = error + pre_position.first;
      //std::cout << "new_error: " << new_error << ", pre_error: " << pre_position.first << std::endl;
      errorVec.emplace_back(std::make_pair(new_error, pre_position.second));
    }
    errorVecVec.emplace_back(errorVec);
  }
  return errorVecVec;
}



std::pair<double, PrePathInfo> GraphAlignSim::getMinPreError(AlignmentGraphPtr alignGraph_ptr, std::vector<Vertex_AGraph> pre_layer, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, Vertex_AGraph node_a, Vertex_AGraph node_b){
  std::pair<in_alignEdge_iter, in_alignEdge_iter> inEdge_a = boost::in_edges(node_a, *alignGraph_ptr);
  in_alignEdge_iter in_edgeIter_a = inEdge_a.first;
  int inEdgeNum_a = boost::in_degree(node_a, *alignGraph_ptr);
  Edge_AGraph in_edge_a;
  std::vector<int> path_a;
  //std::cout << pre_layer.size() << " = size" << std::endl;
  for(int i = 0; i < inEdgeNum_a; i++){
    in_edge_a = *in_edgeIter_a;
    Vertex_AGraph source = boost::source(in_edge_a, *alignGraph_ptr);
    int p = find(pre_layer.begin(), pre_layer.end(), source) - pre_layer.begin();
    if(p == pre_layer.size()){
      ++in_edgeIter_a;
    }
    else{
      path_a.emplace_back(p);   
      ++in_edgeIter_a;
    }
     
  }

  std::pair<in_alignEdge_iter, in_alignEdge_iter> inEdge_b = boost::in_edges(node_b, *alignGraph_ptr);
  in_alignEdge_iter in_edgeIter_b = inEdge_b.first;
  int inEdgeNum_b = boost::in_degree(node_b, *alignGraph_ptr);
  Edge_AGraph in_edge_b;
  std::vector<int> path_b;
  for(int j = 0; j < inEdgeNum_b; j++){
    in_edge_b = *in_edgeIter_b;
    Vertex_AGraph source = boost::source(in_edge_b, *alignGraph_ptr);
    int q = find(pre_layer.begin(), pre_layer.end(), source) - pre_layer.begin();
    if(q == pre_layer.size()){
      ++in_edgeIter_b;
    }
    else{
      path_b.emplace_back(q);   
      ++in_edgeIter_b; 
    }
    
  }

  //std::cout << "path_a size: " << path_a.size() << ", " << path_b.size() << std::endl;
  //std::cout << path_a[0] << ", " << path_b[0] << std::endl;
  //if(path_a.size() == 2) {
  //  std::cout << "a: " << path_a[1] << std::endl;
  //}
  //if(path_b.size() == 2) {
  //  std::cout << "b: " << path_b[1] << std::endl;
  //}

  double min_error;
  bool first = true;
  PrePathInfo pre_pair;
  //Vertex_AGraph node1, node2;
  for(int m = 0; m < path_a.size(); m++){
    for(int n = 0; n < path_b.size(); n++){
      //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
      if(first == true){
        min_error = D_layer[path_a[m]][path_b[n]].first;
        //find(pre_layer.begin(), pre_layer.end(), path_a[m]);
        pre_pair = std::make_pair(std::make_pair(path_a[m], pre_layer[path_a[m]]), std::make_pair(path_b[n], pre_layer[path_b[n]]));
        //node1 = pre_layer[path_a[m]];
        //node2 = pre_layer[path_b[n]];
        //std::cout << "1PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
        first = false;
      }
      else{
        //std::cout << "---" << std::endl;
        //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
        if(min_error > D_layer[path_a[m]][path_b[n]].first){
          //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
          min_error = D_layer[path_a[m]][path_b[n]].first;
          pre_pair = std::make_pair(std::make_pair(path_a[m], pre_layer[path_a[m]]), std::make_pair(path_b[n], pre_layer[path_b[n]]));
          //node1 = pre_layer[path_a[m]];
          //node2 = pre_layer[path_b[n]];
          //std::cout << "2PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
        }
      }
    }
  }
  
  //std::cout << "PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
  return std::make_pair(min_error, pre_pair);

}

std::pair<double, PrePathInfo> GraphAlignSim::getMinPreError_case2(AlignmentGraphPtr alignGraph_ptr, std::vector<Vertex_AGraph> pre_layer1, std::vector<Vertex_AGraph> pre_layer2, std::vector<std::vector<std::pair<double, PrePathInfo>>> D_layer, Vertex_AGraph node_a, Vertex_AGraph node_b){
  std::pair<in_alignEdge_iter, in_alignEdge_iter> inEdge_a = boost::in_edges(node_a, *alignGraph_ptr);
  in_alignEdge_iter in_edgeIter_a = inEdge_a.first;
  int inEdgeNum_a = boost::in_degree(node_a, *alignGraph_ptr);
  Edge_AGraph in_edge_a;
  std::vector<int> path_a;
  //std::cout << pre_layer.size() << " = size" << std::endl;
  for(int i = 0; i < inEdgeNum_a; i++){
    in_edge_a = *in_edgeIter_a;
    Vertex_AGraph source = boost::source(in_edge_a, *alignGraph_ptr);
    int p = find(pre_layer1.begin(), pre_layer1.end(), source) - pre_layer1.begin();
    if(p == pre_layer1.size()){
      ++in_edgeIter_a;
    }
    else{
      path_a.emplace_back(p);   
      ++in_edgeIter_a;
    }
     
  }

  std::pair<in_alignEdge_iter, in_alignEdge_iter> inEdge_b = boost::in_edges(node_b, *alignGraph_ptr);
  in_alignEdge_iter in_edgeIter_b = inEdge_b.first;
  int inEdgeNum_b = boost::in_degree(node_b, *alignGraph_ptr);
  Edge_AGraph in_edge_b;
  std::vector<int> path_b;
  for(int j = 0; j < inEdgeNum_b; j++){
    in_edge_b = *in_edgeIter_b;
    Vertex_AGraph source = boost::source(in_edge_b, *alignGraph_ptr);
    int q = find(pre_layer2.begin(), pre_layer2.end(), source) - pre_layer2.begin();
    if(q == pre_layer2.size()){
      ++in_edgeIter_b;
    }
    else{
      path_b.emplace_back(q);   
      ++in_edgeIter_b; 
    }
    
  }

  //std::cout << "path_a size: " << path_a.size() << ", " << path_b.size() << std::endl;
  //std::cout << path_a[0] << ", " << path_b[0] << std::endl;
  //if(path_a.size() == 2) {
  //  std::cout << "a: " << path_a[1] << std::endl;
  //}
  //if(path_b.size() == 2) {
  //  std::cout << "b: " << path_b[1] << std::endl;
  //}

  double min_error;
  bool first = true;
  PrePathInfo pre_pair;
  //Vertex_AGraph node1, node2;
  for(int m = 0; m < path_a.size(); m++){
    for(int n = 0; n < path_b.size(); n++){
      //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
      if(first == true){
        min_error = D_layer[path_a[m]][path_b[n]].first;
        //find(pre_layer.begin(), pre_layer.end(), path_a[m]);
        pre_pair = std::make_pair(std::make_pair(path_a[m], pre_layer1[path_a[m]]), std::make_pair(path_b[n], pre_layer2[path_b[n]]));
        //node1 = pre_layer[path_a[m]];
        //node2 = pre_layer[path_b[n]];
        //std::cout << "1PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
        first = false;
      }
      else{
        //std::cout << "---" << std::endl;
        //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
        if(min_error > D_layer[path_a[m]][path_b[n]].first){
          //std::cout << path_a[m] << ", " << path_b[n] << std::endl;
          min_error = D_layer[path_a[m]][path_b[n]].first;
          pre_pair = std::make_pair(std::make_pair(path_a[m], pre_layer1[path_a[m]]), std::make_pair(path_b[n], pre_layer2[path_b[n]]));
          //node1 = pre_layer[path_a[m]];
          //node2 = pre_layer[path_b[n]];
          //std::cout << "2PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
        }
      }
    }
  }
  
  //std::cout << "PrePair: " << (*alignGraph_ptr)[node1].j_ << ", " << (*alignGraph_ptr)[node2].j_ << std::endl;
  return std::make_pair(min_error, pre_pair);

}




void GraphAlignSim::add_all_info(int i, int j, int k,
                              std::vector<std::vector<std::vector<short int>>> T, 
                              std::vector<std::vector<std::vector<std::vector<prePosition>>>> E, 
                              const AlignmentGraphPtr &alignGraph_ptr,
                              const std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple2>> &vertexMapPtr,
                              int ori_index, int exactMass, int blackMass, std::vector<std::pair<unsigned short, unsigned short>> modInfo){
  int v_index;
  if (vertexMapPtr->find(std::make_tuple(i, j, k)) == vertexMapPtr->end())
  {
      VertexInfo_AGraph v_info(T[i][j][k], i, j, k);
      v_index = add_vertex(v_info, *alignGraph_ptr.get());
      vertexMapPtr->insert(std::make_pair(std::make_tuple(i, j, k), v_index));
      //verIndex = verIndex + 1;
      //should return 'verIndex';
  }
  else
  {
      v_index = vertexMapPtr->find(std::make_tuple(i, j, k))->second;
  }
  Vertex_AGraph v_pre = vertex(v_index, *alignGraph_ptr.get());
  Vertex_AGraph v2 = vertex(ori_index, *alignGraph_ptr.get());
  EdgeInfo_AGraph edge_info(blackMass, exactMass, std::move(modInfo));
  add_edge(v2, v_pre, edge_info, *alignGraph_ptr.get());

  if(T[i][j][k]>1){
    //std::cout << "T[" << i << "][" << j << "][" << k << "]=" << T[i][j][k]<< std::endl;
    auto preNodes = E[i][j][k];
    //std::cout << "preNodes.size(): " << preNodes.size() << std::endl;
    for(auto & preNode : preNodes){
      int i_pre = preNode.first.first;
      int j_pre = preNode.first.second;
      int k_value = std::get<0>(preNode.second);
      int k_pre = k_value + deltaL[j_pre];
      int exactMass2 = std::get<1>(preNode.second);
      int blackMass2 = proteo_graph_ptr_->getSeqMass(i_pre, i); // need to be confirmed
      auto modInfo2 = std::get<2>(preNode.second);
      //std::cout << "node_num: " << num_vertices(*alignGraph_ptr.get()) << std::endl;
      add_all_info(i_pre,j_pre,k_pre, T, E, alignGraph_ptr, vertexMapPtr, v_index, exactMass2, blackMass2, modInfo2);
    }
  }
}





void GraphAlignSim::BFS(std::queue<std::tuple<int,int,int>> node_queue,
                     std::vector<std::vector<std::vector<short int>>> T,
                     std::vector<std::vector<std::vector<std::vector<prePosition>>>> E,
                     AlignmentGraphPtr &alignGraph_ptr,
                     std::shared_ptr<std::unordered_map<std::tuple<int, int, int>, int, toppic::hashKey_tuple2>> &vertexMapPtr){
    double max_inten = 0;
    while(!node_queue.empty()){
        auto node = node_queue.front();
        node_queue.pop();

        int i = std::get<0>(node);
        int j = std::get<1>(node);
        int k = std::get<2>(node);
        double cur_inten = (spec_graph_ptr_->getPeakPtrVec())[j].first->getIntensity();
        if(cur_inten > max_inten){
          max_inten = cur_inten;
        }
        if(T[i][j][k] > 1){
            auto preNodes = E[i][j][k];
            //std::cout << "1: " << preNodes.size() << std::endl;
            int ori_index = vertexMapPtr->find(node)->second;
            for(auto & preNode : preNodes){
                int i_pre = preNode.first.first;
                int j_pre = preNode.first.second;
                int k_value = std::get<0>(preNode.second);
                int k_pre = k_value + deltaL[j_pre];

                int v_index;
                if (vertexMapPtr->find(std::make_tuple(i_pre, j_pre, k_pre)) == vertexMapPtr->end())
                {
                    VertexInfo_AGraph v_info(T[i_pre][j_pre][k_pre], i_pre, j_pre, k_pre);
                    v_index = add_vertex(v_info, *alignGraph_ptr.get());
                    vertexMapPtr->insert(std::make_pair(std::make_tuple(i_pre, j_pre, k_pre), v_index));
                    node_queue.emplace(i_pre, j_pre, k_pre);
                    //verIndex = verIndex + 1;
                    //should return 'verIndex';
                }
                else
                {
                    v_index = vertexMapPtr->find(std::make_tuple(i_pre, j_pre, k_pre))->second;
                }


                int exactMass = std::get<1>(preNode.second);
                int blackMass = proteo_graph_ptr_->getSeqMass(i_pre,i); //need to be confirmed
                auto modInfo = std::get<2>(preNode.second);
                //std::cout << i_pre << ", " << j_pre << ", " << k_pre << std::endl;

                Vertex_AGraph v_pre = vertex(v_index, *alignGraph_ptr.get());
                Vertex_AGraph v2 = vertex(ori_index, *alignGraph_ptr.get());
                EdgeInfo_AGraph edge_info(blackMass, exactMass, std::move(modInfo));
                add_edge(v2, v_pre, edge_info, *alignGraph_ptr.get());

            }
        }

    }

}


std::pair<double, std::vector<std::vector<Vertex_AGraph>>> GraphAlignSim::findLargestIntensity(AlignmentGraphPtr &alignGraph_ptr, Vertex_AGraph beginNode){
  double biggestInten = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[beginNode].j_].first->getIntensity();
  std::vector<std::vector<Vertex_AGraph>> nodeLayers;
  std::vector<Vertex_AGraph> current_layer;
  std::vector<Vertex_AGraph> pre_layer;
  pre_layer.emplace_back(beginNode);
  nodeLayers.emplace_back(pre_layer);
  for(int layer_index = (*alignGraph_ptr)[beginNode].T_; layer_index > 1; layer_index--){
    for(int i = 0; i < pre_layer.size(); i++){
      std::pair<out_alignEdge_iter, out_alignEdge_iter> outEdge = boost::out_edges(pre_layer[i], *alignGraph_ptr);
      out_alignEdge_iter edgeIter_start = outEdge.first;
      out_alignEdge_iter edgeIter_next = edgeIter_start;
      int outEdgeNum = boost::out_degree(pre_layer[i], *alignGraph_ptr);
      Edge_AGraph edge_next;
      for(int j = 0; j < outEdgeNum; j++){
        edge_next = *edgeIter_next;
        Vertex_AGraph target = boost::target(edge_next, *alignGraph_ptr);
        if(find(current_layer.begin(), current_layer.end(), target) == current_layer.end()){
          current_layer.emplace_back(target);
          double inten = (spec_graph_ptr_->getPeakPtrVec())[(*alignGraph_ptr)[target].j_].first->getIntensity();
          if(inten > biggestInten){
            biggestInten = inten;
          }
        }      
        ++edgeIter_next;
      }
    }
    nodeLayers.emplace_back(current_layer);
    pre_layer = current_layer;
    current_layer.clear();
  }
  return std::make_pair(biggestInten, nodeLayers);  
}



}  // namespace toppic

