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

#include <algorithm>

#include "search/graph/dist.hpp"

namespace toppic {

int getVecIndex(int v1, int v2, int gap) {
  int index = (gap + 1) * v1 + (v2 - v1);
  return index;
}

void addToDistVec(MassGraphPtr graph_ptr, const std::vector<std::vector<std::set<std::pair<int, std::vector<std::pair<unsigned short,unsigned short>>>>>> & dist_vecs,
                  int node_num, int mod_num, DistVec & dist_vec, int gap) {
  std::set<Dist> dist_set;
  for (int i = 0; i < node_num - 1; i++) {
    for (int j = i + 1; j < node_num && j <= i + gap; j++) {
      int index = getVecIndex(i, j, gap);
      for (auto it=dist_vecs[index][mod_num].begin();
           it != dist_vecs[index][mod_num].end(); it++) {
        if ((*it).first == 0) continue;
        Dist tmp = Dist(graph_ptr, (*it).first, i, j, (*it).second);
        auto search = dist_set.find(tmp);
        if (search != dist_set.end()) {
          search->pair_ij_.push_back(std::make_pair(std::pair<int, int>(i, j),(*it).second));
        } else {
          dist_set.insert(tmp);
        }
      }
    }
  }

  std::copy(dist_set.begin(), dist_set.end(), std::back_inserter(dist_vec));
}

}  // namespace toppic
