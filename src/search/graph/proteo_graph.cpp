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

#include <set>
#include <vector>

#include "common/util/logger.hpp"
#include "common/base/mod_base.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/residue_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/residue_seq.hpp"
#include "seq/proteoform_util.hpp"

#include "search/graph/proteo_graph.hpp"

namespace toppic {

ProteoGraph::ProteoGraph(FastaSubSeqPtr fasta_seq_ptr, ModPtrVec fix_mod_ptr_vec,
                         MassGraphPtr graph_ptr, bool is_nme,
                         double convert_ratio, int max_mod_num,
                         int max_ptm_sum_mass, int proteo_graph_gap,
                         int var_ptm_in_gap):
    is_nme_(is_nme),
    proteo_graph_gap_(proteo_graph_gap),
    var_ptm_in_gap_(var_ptm_in_gap) {
      //std::cout << "geneDbProteoformPtr" << std::endl;
      db_proteo_ptr_ = proteoform_util::geneDbProteoformPtr(fasta_seq_ptr, fix_mod_ptr_vec,
                                                            fasta_seq_ptr->getSubSeqStart());
      graph_ptr_ = graph_ptr;

      node_num_ = num_vertices(*graph_ptr.get());
      LOG_DEBUG("node num " << node_num_);
      pair_num_ = node_num_ * (proteo_graph_gap_ + 1);
      //std::cout << "compSeqMasses" << std::endl;
      compSeqMasses(convert_ratio);
      //sstd::cout << "compDistances" << std::endl;
//    compAddedPrecs(max_mod_num, max_ptm_sum_mass);
      compDistances(max_mod_num, max_ptm_sum_mass);
    }

int ProteoGraph::getVecIndex(int v1, int v2) {
  int index =  (proteo_graph_gap_ + 1) * v1 + (v2 - v1);
  return index;
}



int ProteoGraph::getSeqMass(int v1, int v2) {
  int index = getVecIndex(v1, v2);
  return seq_masses_[index];
}

void ProteoGraph::compSeqMasses(double convert_ratio) {
  ResSeqPtr res_seq_ptr = db_proteo_ptr_->getResSeqPtr();
  seq_masses_ = std::vector<int>(pair_num_, 0);
  for (int i = 0; i < node_num_; i ++) {
    int mass = 0;
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      int cur_mass = std::round(res_seq_ptr->getResiduePtr(j-1)->getMass() * convert_ratio);
      mass += cur_mass;
      int index = getVecIndex(i, j);
      seq_masses_[index] = mass;
    }
  }
}

void ProteoGraph::compAddedPrecs(int max_mod_num, int max_ptm_sum_mass){
    MassGraph *g_p = graph_ptr_.get();
    // get mass without ptms

    std::vector<std::vector<std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>>>> dist_vecs;
    for (int i = 0; i < node_num_; i++) {
        //std::vector<std::pair<double,int>> empty_mod;
        std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>> empty_set; //needs to be conformed
        std::vector<std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>> > one_pair_vec;
        for (int j = 0; j < max_mod_num + 1; j ++) {
            one_pair_vec.push_back(empty_set);
        }
        dist_vecs.push_back(one_pair_vec);
    }
    //initialize 0,0
    std::vector<std::pair<unsigned short,unsigned short>> empty_mod;
    dist_vecs[0][0].insert(std::make_pair(0, empty_mod));

    for (int j =  1; j < node_num_; j++) {
        Vertex v2 = vertex(j, *g_p);
        Vertex pre_v2 = vertex(j-1, *g_p);
        int index = j;
        int pre_index = j - 1;
        boost::graph_traits<MassGraph>::out_edge_iterator ei, ei_end;
        boost::tie(ei, ei_end) = out_edges(pre_v2, *g_p);
        for ( ; ei != ei_end; ++ei) {
            if (target(*ei, *g_p) == v2) {
                MassGraph::edge_descriptor e = *ei;
                int d =(*g_p)[e].int_mass_;
                int change = (*g_p)[e].alter_type_;
                if(index == 6){
                    std::cout<<"===="<<d<<std::endl;
                }
                for (int k = 0; k < var_ptm_in_gap_ + 1; k++) {
                    if (k == max_mod_num &&
                        (change == AlterType::PROTEIN_VARIABLE->getId()
                         || change == AlterType::VARIABLE->getId())) {
                        continue;
                    }
                    for (auto it=dist_vecs[pre_index][k].begin();
                         it != dist_vecs[pre_index][k].end(); it++) {
                        int new_d = d + (*it).first;
                        if(index == 6){
                            std::cout<<"newd "<<new_d<<std::endl;
                        }
                        if (std::abs(new_d - seq_masses_[index]) <= max_ptm_sum_mass) {
                            if (change == AlterType::PROTEIN_VARIABLE->getId()
                                || change == AlterType::VARIABLE->getId()) {
                                //dist_vecs[index][k+1].insert(new_d);
                                bool add1 = true;
                                if(dist_vecs[index][k+1].size() > 0){
                                    for(auto iter_a = dist_vecs[index][k+1].begin(); iter_a != dist_vecs[index][k+1].end(); iter_a++){
                                        if((*iter_a).first == new_d){
                                            add1 = false;
                                            break;
                                        }
                                    }
                                }
                                if(add1 == true){
                                    unsigned short unID = (*g_p)[e].res_ptr_->getPtmPtr()->getUnimodId();
                                    auto oldMod = (*it).second;
                                    oldMod.push_back(std::make_pair(unID, j));
                                    dist_vecs[index][k+1].insert(std::make_pair(new_d, oldMod));
                                    //std::cout << "add: dist_vecs[" << index << "][" <<  k+1 << "]: (" << new_d << "," << oldMod.size() << std::endl;
                                }
                            } else {
                                bool add2 = true;
                                if(dist_vecs[index][k].size() > 0){
                                    for(auto iter_b = dist_vecs[index][k].begin(); iter_b != dist_vecs[index][k].end(); iter_b++){
                                        if((*iter_b).first == new_d){
                                            add2 = false;
                                            break;
                                        }
                                    }
                                }
                                if(add2 == true){
                                    if(index==6){
                                        std::cout<<"+_+_+_"<<new_d<<std::endl;
                                    }
                                    dist_vecs[index][k].insert(std::make_pair(new_d, (*it).second));
                                    //std::cout << "add: dist_vecs[" << index << "][" <<  k << "]: (" << new_d << "," << (*it).second.size() << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
        //For one index, here is the ending for building it.
        /*
        auto a = dist_vecs[index];
        for(int i1 = 0; i1 < a.size(); i1++){
          std::cout << "dist_vecs[" << index << "][" << i1 << "]: {";
          auto b = a[i1];
          for(auto i2 = b.begin(); i2 != b.end(); i2++){
            std::cout << "(" << (*i2).first << ",";
            auto c = (*i2).second;
            for(int i3 = 0; i3 < c.size(); i3++){
              std::cout << "<" << c[i3].first << "," << c[i3].second << ">";
            }
            std::cout << ")";
          }
          std::cout << "}" << std::endl;
        }
        */
    }
    int prec_index = node_num_ - 1;
    std::cout<<"prec_index"<<prec_index<<std::endl;
    std::set<int> prec_mass_set;
    for(int k = 0; k < max_mod_num + 1; k ++) {
        for (auto it = dist_vecs[prec_index][k].begin();
             it != dist_vecs[prec_index][k].end(); it++) {
            if ((*it).first != 0) {
                prec_mass_set.insert((*it).first);
            }
        }
    }
    prec_mass_set_ = prec_mass_set;
}


void ProteoGraph::compDistances(int max_mod_num, int max_ptm_sum_mass) {
  MassGraph *g_p = graph_ptr_.get();
  // get mass without ptms

  std::vector<std::vector<std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>>>> dist_vecs;
  for (int i = 0; i < pair_num_; i++) {
    //std::vector<std::pair<double,int>> empty_mod;
    std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>> empty_set; //needs to be conformed
    std::vector<std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>> > one_pair_vec;
    for (int j = 0; j < max_mod_num + 1; j ++) {
      one_pair_vec.push_back(empty_set);
    }
    dist_vecs.push_back(one_pair_vec);
  }


  // initialize pair (i, i)
  for (int i = 0; i < node_num_; i++) {
    int index = getVecIndex(i, i);
        std::vector<std::pair<unsigned short,unsigned short>> empty_mod;
        dist_vecs[index][0].insert(std::make_pair(0, empty_mod));
  }
  //std::cout << "if here" << std::endl;
  for (int i = 0; i < node_num_ - 1; i++) {
    for (int j = i + 1; j < node_num_ && j <= i + proteo_graph_gap_; j++) {
      Vertex v2 = vertex(j, *g_p);
      Vertex pre_v2 = vertex(j-1, *g_p);
      int index = getVecIndex(i, j);
      int pre_index = getVecIndex(i, j-1);
      boost::graph_traits<MassGraph>::out_edge_iterator ei, ei_end;
      boost::tie(ei, ei_end) = out_edges(pre_v2, *g_p);
      for ( ; ei != ei_end; ++ei) {
        if (target(*ei, *g_p) == v2) {
          MassGraph::edge_descriptor e = *ei;
          int d =(*g_p)[e].int_mass_;
          int change = (*g_p)[e].alter_type_;
//          if(index == 6){
//              std::cout<<"===="<<d<<std::endl;
//          }
          for (int k = 0; k < var_ptm_in_gap_ + 1; k++) {
            if (k == max_mod_num &&
                (change == AlterType::PROTEIN_VARIABLE->getId()
                 || change == AlterType::VARIABLE->getId())) {
              continue;
            }
            for (auto it=dist_vecs[pre_index][k].begin();
                 it != dist_vecs[pre_index][k].end(); it++) {
              int new_d = d + (*it).first;
//              if(index == 6){
//                  std::cout<<"newd "<<new_d<<std::endl;
//              }
              if (std::abs(new_d - seq_masses_[index]) <= max_ptm_sum_mass) {
                if (change == AlterType::PROTEIN_VARIABLE->getId()
                    || change == AlterType::VARIABLE->getId()) {
                  //dist_vecs[index][k+1].insert(new_d);
                  bool add1 = true;
                  if(dist_vecs[index][k+1].size() > 0){
                    for(auto iter_a = dist_vecs[index][k+1].begin(); iter_a != dist_vecs[index][k+1].end(); iter_a++){
                      if((*iter_a).first == new_d){
                        add1 = false;
                        break;
                      }
                    }
                  }
                  if(add1 == true){
                    unsigned short unID = (*g_p)[e].res_ptr_->getPtmPtr()->getUnimodId();
                    auto oldMod = (*it).second;
                    oldMod.push_back(std::make_pair(unID, j));
                    dist_vecs[index][k+1].insert(std::make_pair(new_d, oldMod));
                    //std::cout << "add: dist_vecs[" << index << "][" <<  k+1 << "]: (" << new_d << "," << oldMod.size() << std::endl;
                  }
                } else {
                  bool add2 = true;
                  if(dist_vecs[index][k].size() > 0){
                    for(auto iter_b = dist_vecs[index][k].begin(); iter_b != dist_vecs[index][k].end(); iter_b++){
                      if((*iter_b).first == new_d){
                        add2 = false;
                        break;
                      }
                    }
                  }
                  if(add2 == true){
//                    if(index==6){
//                        std::cout<<"+_+_+_"<<new_d<<std::endl;
//                    }
                    dist_vecs[index][k].insert(std::make_pair(new_d, (*it).second));
                    //std::cout << "add: dist_vecs[" << index << "][" <<  k << "]: (" << new_d << "," << (*it).second.size() << std::endl;
                  }
                }
              }
            }
          }
        }
      }
      //For one index, here is the ending for building it.
      /*
      auto a = dist_vecs[index];
      for(int i1 = 0; i1 < a.size(); i1++){
        std::cout << "dist_vecs[" << index << "][" << i1 << "]: {";
        auto b = a[i1];
        for(auto i2 = b.begin(); i2 != b.end(); i2++){
          std::cout << "(" << (*i2).first << ",";
          auto c = (*i2).second;
          for(int i3 = 0; i3 < c.size(); i3++){
            std::cout << "<" << c[i3].first << "," << c[i3].second << ">";
          }
          std::cout << ")";
        }
        std::cout << "}" << std::endl;
      }
      */
    }
  }
  //std::cout << "end here" << std::endl;




  dist_vec_.reserve(max_mod_num);
  DistVec tmp;
  for (int k = 0; k < max_mod_num + 1; k++) {
    dist_vec_.push_back(tmp);
  }

  for (int k = 0; k < max_mod_num + 1; k++) {
    addToDistVec(graph_ptr_, dist_vecs, node_num_, k, dist_vec_[k], proteo_graph_gap_);
  }
  //std::cout << "end" << std::endl;
}

}  // namespace toppic

