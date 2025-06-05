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

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

//#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
//#include "common/thread/simple_thread_pool.hpp"
//#else
//#include <sys/wait.h>
//#endif

#include "common/thread/simple_thread_pool.hpp"

#include "common/util/file_util.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_sub_util.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_xml_writer_util.hpp"
#include "search/graph/proteo_graph_reader.hpp"
#include "search/graph/spec_graph_reader.hpp"
#include "search/graphalign/graph_align.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "search/graph/spec_graph_sim.hpp"
#include "search/graphalign/graph_align_a.hpp"


namespace toppic {


std::function<void()> geneTask(GraphAlignMngPtr mng_ptr,
                               ModPtrVec var_mod_ptr_vec,
                               ModPtrVec N_mod_ptr_vec,
                               int spectrum_num, int idx) {
  return [mng_ptr, var_mod_ptr_vec, N_mod_ptr_vec, spectrum_num, idx]() {




    std::unordered_map<int,PtmPtr> ptm_map;
    for(int mm = 0; mm < var_mod_ptr_vec.size(); mm++){
        PtmPtr ptm = var_mod_ptr_vec[mm]->getModResiduePtr()->getPtmPtr();
        int mod_id = ptm->getUnimodId();
        ptm_map.insert(std::make_pair(mod_id, ptm));
    }
    mng_ptr->ptm_map_ = ptm_map;

    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

    SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
      int group_spec_num = prsm_para_ptr->getGroupSpecNum();
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

      MsAlignReader sp_reader(sp_file_name, group_spec_num,
                              sp_para_ptr->getActivationPtr(),
                              sp_para_ptr->getSkipList());

      SpecGraphReader spec_reader(sp_file_name,
                                  prsm_para_ptr->getGroupSpecNum(),
                                  mng_ptr->convert_ratio_,
                                  sp_para_ptr);

      double nterm_fix_shift;
      if (prsm_para_ptr->getNtermFixModPtrVec().empty()){
          nterm_fix_shift = 0;
      }else{
          nterm_fix_shift = prsm_para_ptr->getNtermFixModPtrVec()[0]->getShift();
      }


      SpectrumSetPtr spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
      ProteoAnnoPtr proteo_anno_ptr
              = std::make_shared<ProteoAnno>(prsm_para_ptr->getFixModPtrVec(),
//                                             prsm_para_ptr->getNtermFixModPtrVec(),
                                             prsm_para_ptr->getProtModPtrVec(),
                                             var_mod_ptr_vec,
                                             N_mod_ptr_vec);

      std::string sp_directory = sp_file_name;
      size_t pos = sp_directory.find_last_of('/');
      if (pos != std::string::npos) {
          sp_directory.erase(pos + 1);
      }
      mng_ptr->resultpath = sp_directory;

      std::ofstream bgnoutFile;
      std::string bgn_file = mng_ptr->resultpath + "results.txt";
      mng_ptr->result_file_ = bgn_file;
      bgnoutFile.open(bgn_file, std::ios::app);
//                      outFile << "===spectrum id: " << spec_id << "===" << std::endl;
      bgnoutFile <<  "ID\toffset\tError\tMod1\tMod2\tAbund1\tAbund2\tq1\tq2\n";
      bgnoutFile.close();


      std::fstream f(sp_directory + "ref_peptide.txt");
      std::string line;
      std::vector<std::pair<int, std::string>> ref_peps;
      while(getline(f,line)){
          std::vector<std::string> l = str_util::split(line, "\t");
          int sp_id = std::stoi(l[0]);
          std::string pep = l[1];
          pep = pep.substr(0, pep.length() - 1);
          ref_peps.push_back(std::make_pair(sp_id, pep));
      }
      f.close();

      int cnt = 0;
      clock_t start,end;
      start = clock();
      while (spec_set_ptr != nullptr) {
          cnt += group_spec_num;
//          if (spec_set_ptr->isValid() && spec_set_ptr->getSpectrumId() == 71) {
              if (spec_set_ptr->isValid()) {
              //***** read spectrum
              int spec_id = spec_set_ptr->getSpectrumId();
              DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
              double prec_mono_mass = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
//              SpecGraphPtrVec spec_ptr_vec
//                      = spec_reader.getNextSpecGraphPtrVec(spec_set_ptr, mng_ptr->prec_error_);
              std::vector<double> prec_error_vec;
              double IM = mass_constant::getIsotopeMass();
              if(mng_ptr->use_adjusted_precmass){
                  prec_error_vec = {0, -IM, IM};
              }else{
                  prec_error_vec = {0};;
              }
              for(int ee = 0; ee < prec_error_vec.size(); ee ++) {

                  double adjusted_prec_mono_mass = prec_mono_mass + prec_error_vec[ee];

                  SpectrumSetPtr adjusted_spec_set_ptr
                          = std::make_shared<SpectrumSet>(deconv_ms_ptr_vec, sp_para_ptr, adjusted_prec_mono_mass);
                  PrmMsPtrVec ms_two_vec = adjusted_spec_set_ptr->getMsTwoPtrVec();


                  std::cout << "===spectrum id: " << spec_id << "===offset: " << prec_error_vec[ee] <<"==="<<std::endl;

                  //***** find ref peptide and build PMG
                  std::string ref_raw_seq;
                  bool ref_found = false;
                  for (int p = 0; p < ref_peps.size(); p++) {
                      if (ref_peps[p].first == spec_id) {
                          ref_raw_seq = ref_peps[p].second;
                          ref_found = true;
                          break;
                      }
                  }
                  if (ref_found) {
                      FastaSeqPtr seq_ptr = std::make_shared<FastaSeq>("sp" + std::to_string(spec_id) + "-RefPeptide",
                                                                       ref_raw_seq);
                      FastaSubSeqPtr sub_seq_ptr = std::make_shared<FastaSubSeq>(seq_ptr, 0,
                                                                                 seq_ptr->getAcidPtmPairLen());
                      proteo_anno_ptr->anno(ref_raw_seq, true);
                      MassGraphPtr graph_ptr = getMassGraphPtr(proteo_anno_ptr, mng_ptr->convert_ratio_);
                      ProteoGraphPtr proteo_ptr = std::make_shared<ProteoGraph>(sub_seq_ptr,
                                                                                prsm_para_ptr->getFixModPtrVec(),
                                                                                graph_ptr,
                                                                                proteo_anno_ptr->isNme(),
                                                                                mng_ptr->convert_ratio_,
                                                                                mng_ptr->max_known_mods_,
                                                                                mng_ptr->getIntMaxPtmSumMass(),
                                                                                mng_ptr->proteo_graph_gap_,
                                                                                mng_ptr->var_ptm_in_gap_);
                      std::cout << "PMG built" << std::endl;

                      // add mod prec peak
//                  std::set<int> prec_mass_set = proteo_ptr->getPrecMassSet();
                      PrmPeakPtrVec peak_vec = prm_ms::getPrmPeakPtrs(ms_two_vec, sp_para_ptr->getPeakTolerancePtr());
//                  PrmPeakPtrVec peak_vec = prm_ms::getPrmPeakPtrs_addPrecs(ms_two_vec, sp_para_ptr->getPeakTolerancePtr(),prec_mass_set,mng_ptr->convert_ratio_);


                      PrmPeakPtrVec no_nterm_peak_vec;
                      no_nterm_peak_vec.push_back(peak_vec[0]);
                      for(int v = 1; v < peak_vec.size();  v++){
                          double ori_mass = peak_vec[v]->getMonoMass();
                          double new_mass = ori_mass - nterm_fix_shift;
                          if(new_mass > mng_ptr->peak_min_mass_) {
                              peak_vec[v]->setMonoMass(ori_mass - nterm_fix_shift);
                              no_nterm_peak_vec.push_back(peak_vec[v]);
                          }
                      }

                      peak_vec = no_nterm_peak_vec;


//                      std::ofstream outFile_peaks;
//                      std::string peaks_file = mng_ptr->resultpath + "processed_spectrum.txt";
//                      mng_ptr->protoform_result_file_ = peaks_file;
//                      outFile_peaks.open(peaks_file, std::ios::app);
////                      outFile << "===spectrum id: " << spec_id << "===" << std::endl;
//                      outFile_peaks << "===spectrum id: " << spec_id << "===offset: " << prec_error_vec[ee] <<"==="<<std::endl;
//                      for(int pe = 0; pe < peak_vec.size();pe++){
//                          outFile_peaks << pe << '\t' << peak_vec[pe]->getMonoMass() << "\t" << peak_vec[pe]->getIntensity() <<std::endl;
//                      }
//                      outFile_peaks.close();



                      //***** build SMG
                      std::vector<std::pair<PeakPtr, std::string>> peakVec;
                      for (int p = 0; p < peak_vec.size(); p++) {
                          double peakMass = peak_vec[p]->getMonoMass();
                          double peakInten = peak_vec[p]->getIntensity();
                          std::string type = peak_vec[p]->getBaseTypePtr()->getName();
                          Peak cur_peak(peakMass, peakInten);
                          peakVec.push_back(std::make_pair(std::make_shared<Peak>(cur_peak), type));
                      }
                      MassGraphPtr sp_graph_ptr = std::make_shared<MassGraph>();
                      // add mass 0/start nod
                      VertexInfo v(0);
                      add_vertex(v, *sp_graph_ptr.get());
                      for (size_t i = 1; i < peakVec.size(); i++) {
                          // add a new node for the prm
                          VertexInfo cur_v(i);
                          add_vertex(cur_v, *sp_graph_ptr.get());

                          Vertex v1, v2;
                          v1 = vertex(i - 1, *sp_graph_ptr.get());
                          v2 = vertex(i, *sp_graph_ptr.get());

                          double dist = peakVec[i].first->getPosition() - peakVec[i - 1].first->getPosition();
                          EdgeInfo edge_info(dist, mng_ptr->convert_ratio_);
                          add_edge(v1, v2, edge_info, *sp_graph_ptr.get());
                      }
                      SpecGraphPtr_sim spec_graph_ptr = std::make_shared<SpecGraph_sim>(peakVec, sp_graph_ptr,
                                                                                        mng_ptr->convert_ratio_);


                      std::cout << "SMG built" << std::endl;

                      GraphAlignPtr_sim graph_align
                              = std::make_shared<GraphAlignSim>(mng_ptr, proteo_ptr, spec_graph_ptr);



                      std::ofstream outFile;
                      std::string result_file = mng_ptr->resultpath + "align_results_gap" + std::to_string(mng_ptr->proteo_graph_gap_) + "_head" + std::to_string(mng_ptr->max_head_diff) +  ".txt";
                      mng_ptr->align_result_file_ = result_file;
                      outFile.open(result_file, std::ios::app);
//                      outFile << "===spectrum id: " << spec_id << "===" << std::endl;
                      outFile << "===spectrum id: " << spec_id << "===offset: " << prec_error_vec[ee] <<"==="<<std::endl;
                      outFile.close();

                      result_file = mng_ptr->resultpath + "proteoforms_results_gap" + std::to_string(mng_ptr->proteo_graph_gap_) + "_head" + std::to_string(mng_ptr->max_head_diff) +  ".txt";
                      mng_ptr->protoform_result_file_ = result_file;
                      outFile.open(result_file, std::ios::app);
//                      outFile << "===spectrum id: " << spec_id << "===" << std::endl;
                      outFile << "===spectrum id: " << spec_id << "===offset: " << prec_error_vec[ee] <<"==="<<std::endl;
                      outFile.close();


                      result_file = mng_ptr->resultpath + "results.txt";
                      mng_ptr->result_file_ = result_file;
                      outFile.open(result_file, std::ios::app);
//                      outFile << "===spectrum id: " << spec_id << "===" << std::endl;
                      outFile <<  spec_id << "\t" << prec_error_vec[ee] <<"\t";
                      outFile.close();



                      graph_align->TopMGFast();


                      //for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
                      //  PrsmPtr prsm_ptr = graph_align->geneResult(shift);
                      //  if (prsm_ptr != nullptr) {
                      //    writer_ptr->write(prsm_ptr);
                      //  }
                      //}
                      graph_align = nullptr;
                  }
              }
          }
          spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
      }//end while
      end = clock();
      std::cout<<"total time: "<<(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;

      std::ofstream outFile;
      std::string result_file = mng_ptr->resultpath+"align_results.txt";
      outFile.open(result_file,std::ios::app);
      outFile<<"total time: "<<(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;
      outFile.close();
  };
}

std::function<void()> geneTask_old(GraphAlignMngPtr mng_ptr,
                               ModPtrVec var_mod_ptr_vec,
                               ModPtrVec N_mod_ptr_vec,
                               int spectrum_num, int idx) {
    return [mng_ptr, var_mod_ptr_vec, N_mod_ptr_vec, spectrum_num, idx]() {


        PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

        SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
    //*********read spectrum information into this program************
    std::string file_name = "./realData/yeast_400_500-Sample003.ms2";
    //std::cout << file_name << std::endl;
    //std::string file_name = "ourDelta27.txt";
    std::ifstream infile(file_name.c_str());
    if (!infile.is_open()) {
        std::cout << "error!" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<std::vector<std::pair<PeakPtr, std::string>>> specInfoVec;
    std::string line;
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    getline(infile, line);
    //std::cout << line << std::endl;
    std::vector<std::pair<PeakPtr, std::string>> peak_vec;
    peak_vec.push_back(std::make_pair(std::make_shared<Peak>(0.0,0.0), "Original"));
    while (!infile.eof()){
      if(isdigit(line[0])){
        //std::cout << "here" << std::endl;
        int j = line.find(' ');
        std::string peakMassStr = line.substr(0,j);
        std::string peakIntenStr = line.substr(j+1, line.size());
        //std::cout << peakMassStr << ", " << peakIntenStr << std::endl;
        double peakMass = stod(peakMassStr);
        int peakInten = stod(peakIntenStr);
        if(peakInten == 0){
          getline(infile, line);
          continue;
        }
        Peak p(peakMass, peakInten);
        peak_vec.push_back(std::make_pair(std::make_shared<Peak>(p), "Original"));
        getline(infile, line);
      }
      else{
        specInfoVec.emplace_back(peak_vec);
        peak_vec.clear();
        peak_vec.push_back(std::make_pair(std::make_shared<Peak>(0.0,0.0), "Original"));
        getline(infile, line);
        getline(infile, line);
        getline(infile, line);
        getline(infile, line);
        getline(infile, line);
        getline(infile, line);
      }

    }
    specInfoVec.emplace_back(peak_vec);
    peak_vec.clear();



    infile.close();
    //std::cout << "ee";
    //SpecGraphPtr spec_graph_ptr
  //        = std::make_shared<SpecGraph>(adjusted_spec_set_ptr, peak_vec, graph_ptr, convert_ratio_);

  //*********read spectrum information into this program************

  //If wanna save running time
  //1.read protein information into this program;
    std::string database = "./realData/yeast_3.fasta";
    std::ifstream infile2(database.c_str());
    if (!infile2.is_open()) {
        std::cout << "error!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<FastaSeqPtr> proteinInfo;
    std::string line2;
    std::getline(infile2, line2);
    while (!infile2.eof()){
      if (line2.size() == 0) std::getline(infile2, line2);
      if (line2.substr(0,1) == ">"){
        int iter = 1;
        std::string seq_name;
        while(line2[iter] != ' '){
          seq_name.append(line2,iter,1);
          iter++;
        }
        std::string seq_desc = line2.substr(iter+1, line2.size());
        std::string ori_seq;
        std::getline(infile2, line2);
        while(line2.substr(0,1) != ">" && !infile2.eof()){
          ori_seq += line2;
          std::getline(infile2, line2);
        }
        FastaSeqPtr seq_ptr = std::make_shared<FastaSeq>(seq_name, seq_desc, ori_seq);
        proteinInfo.push_back(seq_ptr);
        //proteo_anno_ptr->anno(seq_ptr->getRawSeq(), 1);
        //MassGraphPtr graph_ptr = getMassGraphPtr(proteo_anno_ptr, mng_ptr->convert_ratio_);
      }
    }
    infile2.close();
    ProteoAnnoPtr proteo_anno_ptr
        = std::make_shared<ProteoAnno>(prsm_para_ptr->getFixModPtrVec(),
                                       prsm_para_ptr->getProtModPtrVec(),
                                       var_mod_ptr_vec,
                                       N_mod_ptr_vec);
    SpecGraphPtrVec_sim specGraphSimPtrVec;
    ///////////
    //std::cout << "spec:" << specInfoVec.size()<< std::endl;
    for(int j = 0; j < specInfoVec.size(); j++){
      MassGraphPtr sp_graph_ptr = std::make_shared<MassGraph>();
      // add mass 0/start nod
      VertexInfo v(0);
      add_vertex(v, *sp_graph_ptr.get());
      auto peakVec = specInfoVec[j];
      for (size_t i = 1; i < peakVec.size(); i++) {
        // add a new node for the prm
        VertexInfo cur_v(i);
        add_vertex(cur_v, *sp_graph_ptr.get());

        Vertex v1, v2;
        v1 = vertex(i-1, *sp_graph_ptr.get());
        v2 = vertex(i, *sp_graph_ptr.get());

        double dist = peakVec[i].first->getPosition() - peakVec[i-1].first->getPosition();
        //std::cout << "dist[" << i-1 << "," << i <<"]: " << peakVec[i].first->getPosition() << " - " << peakVec[i-1].first->getPosition() << " = " << dist << std::endl;

        EdgeInfo edge_info(dist, mng_ptr->convert_ratio_);
        add_edge(v1, v2, edge_info , *sp_graph_ptr.get());
      }
      SpecGraphPtr_sim spec_graph_ptr = std::make_shared<SpecGraph_sim>(peakVec, sp_graph_ptr, mng_ptr->convert_ratio_);
      //specGraphSimPtrVec.push_back(spec_graph_ptr);
      std::cout << "Processing No. " << j+1 << " spectrum." << std::endl;
        //ReturnStruct bestAlign;
      //std::pair<short, TopMGFastResultPtrVec> bestAlign;
      std::string bestProt;
      int protSize;
      for(int p = 0;p < proteinInfo.size();p++){
        //std::cout << "protein name: " << proteinInfo[p]->getName() << std::endl;
        std::string prot_name = proteinInfo[p]->getName();
        std::cout << "prot name: " << prot_name << std::endl;
        std::vector<FastaSubSeqPtr> seq_ptr_vec = fasta_sub_util::breakSeq(proteinInfo[p]);
        //std::cout << "1" << std::endl;
        bool pName = true;
        if(pName == true){
            for (size_t q = 0; q < seq_ptr_vec.size(); q++){
            proteo_anno_ptr->anno(seq_ptr_vec[q]->getRawSeq(), seq_ptr_vec[q]->isNTerm());
            //std::cout << "2" << std::endl;
            MassGraphPtr graph_ptr = getMassGraphPtr(proteo_anno_ptr, mng_ptr->convert_ratio_);
            //createHash(graph_ptr,seq_name);
            //std::cout << "3" << std::endl;
            ProteoGraphPtr proteo_ptr = std::make_shared<ProteoGraph>(seq_ptr_vec[q],
                                                                      prsm_para_ptr->getFixModPtrVec(),
                                                                      graph_ptr,
                                                                      proteo_anno_ptr->isNme(),
                                                                      mng_ptr->convert_ratio_,
                                                                      mng_ptr->max_known_mods_,
                                                                      mng_ptr->getIntMaxPtmSumMass(),
                                                                      mng_ptr->proteo_graph_gap_,
                                                                      mng_ptr->var_ptm_in_gap_);
            //std::cout << "333" << std::endl;
            GraphAlignPtr_sim graph_align
                      = std::make_shared<GraphAlignSim>(mng_ptr, proteo_ptr, spec_graph_ptr);
            //graph_align->process_v2();
            //graph_align->diagonal_v2(813, 0);
            //graph_align->diagonal_v3(813, 0);
            //std::cout << "start align" << std::endl;
            graph_align->TopMGFast();
            //graph_align->oriProcess();


            //for (int shift = 0; shift <= mng_ptr->n_unknown_shift_; shift++) {
            //  PrsmPtr prsm_ptr = graph_align->geneResult(shift);
            //  if (prsm_ptr != nullptr) {
            //    writer_ptr->write(prsm_ptr);
            //  }
            //}
            graph_align = nullptr;

          }

        }
      }
      //if (idx == 0) {
      //  std::cout << std::flush << "Mass graph alignment - processing " << cnt
      //      << " of " << spectrum_num << " spectra.\r";
      //}
      //spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
    } // end while;



  };
}

void SimplePrsmFilter(SimplePrsmPtrVec & selected_prsm_ptrs) {
  if (selected_prsm_ptrs.size() == 0) return;
  if (selected_prsm_ptrs.size() == 1) {
    if (selected_prsm_ptrs[0]->getScore() < 10) selected_prsm_ptrs.clear();
    return;
  }
  std::sort(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(), SimplePrsm::cmpScoreDec);
  selected_prsm_ptrs.erase(std::remove_if(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                                          [] (const SimplePrsmPtr & p) {return p->getScore() < 10;}),
                           selected_prsm_ptrs.end());
  selected_prsm_ptrs.erase(std::unique(selected_prsm_ptrs.begin(), selected_prsm_ptrs.end(),
                                       [] (const SimplePrsmPtr & a, const SimplePrsmPtr & b) {
                                         return a->getSpectrumScan() == b->getSpectrumScan()
                                           && a->getSeqName() == b->getSeqName();
                                       }),
                           selected_prsm_ptrs.end());
}

void GraphAlignProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  sp_para_ptr->prec_error_ = 0;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;
  LOG_DEBUG("start reading " << var_mod_file_name);
  ModPtrVec var_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[2];
  ModPtrVec N_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[3];
  LOG_DEBUG("end reading " << var_mod_file_name);

  int spectrum_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());


  //--------filtering start--------------
  /*
  std::string input_file_name
      = file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  //std::cout << "Input_file_name: " << sp_file_name << "." << mng_ptr_->input_file_ext_ << std::endl;

  SimplePrsmReaderPtr simple_prsm_reader = std::make_shared<toppic::SimplePrsmReader>(input_file_name);
  std::shared_ptr<SimplePrsmXmlWriter> graph_filter_writer
      = std::make_shared<SimplePrsmXmlWriter>(file_util::basename(sp_file_name) + ".topmg_graph");

  //std::cout << "graph_filter_writer: " << sp_file_name << ".topmg_graph" << std::endl;
  SimplePrsmPtr prsm_ptr = simple_prsm_reader->readOnePrsm();
  int spec_id = -1;
  SimplePrsmPtrVec selected_prsm_ptrs;
  while (prsm_ptr != nullptr) {
    if (prsm_ptr->getSpectrumId() == spec_id) {
      selected_prsm_ptrs.push_back(prsm_ptr); 
    } else {
      SimplePrsmFilter(selected_prsm_ptrs);
      graph_filter_writer->write(selected_prsm_ptrs); 
      selected_prsm_ptrs.clear();
      spec_id = prsm_ptr->getSpectrumId();
      selected_prsm_ptrs.push_back(prsm_ptr);
    }
    prsm_ptr = simple_prsm_reader->readOnePrsm(); 
  }
  simple_prsm_reader->close();
  SimplePrsmFilter(selected_prsm_ptrs);
  graph_filter_writer->write(selected_prsm_ptrs); 
  graph_filter_writer->close();

  SimplePrsmXmlWriterPtrVec simple_prsm_writer_vec = 
      simple_prsm_xml_writer_util::geneWriterPtrVec(input_file_name, mng_ptr_->thread_num_);

  int cnt = 0;
  simple_prsm_reader
      = std::make_shared<toppic::SimplePrsmReader>(file_util::basename(sp_file_name) + ".topmg_graph");
  prsm_ptr = simple_prsm_reader->readOnePrsm();
  while (prsm_ptr != nullptr) {
    cnt = cnt % mng_ptr_->thread_num_;
    simple_prsm_writer_vec[cnt]->write(prsm_ptr);
    cnt++;
    prsm_ptr = simple_prsm_reader->readOnePrsm();
  }
  simple_prsm_reader->close();
  simple_prsm_xml_writer_util::closeWriterPtrVec(simple_prsm_writer_vec);
  */
  //--------filtering end--------------


//#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  std::vector<ThreadPtr> thread_vec;
  //std::cout << "thread_num_ : " << mng_ptr_->thread_num_ << std::endl;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    ThreadPtr thread_ptr = std::make_shared<boost::thread>(geneTask(mng_ptr_, var_mod_ptr_vec, N_mod_ptr_vec, spectrum_num, i));
    thread_vec.push_back(thread_ptr);
  }

  for (size_t i = 0; i < thread_vec.size(); i++) {
    if (thread_vec[i]->joinable()) thread_vec[i]->join();
  }

  time_t now = time(0);
  tm *ltm = localtime(&now);
  std::cout << "The ending time of the program: " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << std::endl;

  exit(0);

/*
#else
  int n = mng_ptr_->thread_num_;

  pid_t pids;

  for (int i = 0; i < n; i++) {
    pids = fork();
    if (pids < 0) {
      std::abort();
    } else if (pids == 0) {
      auto task = geneTask(reader_ptr, mng_ptr_, var_mod_ptr_vec, spectrum_num, i);
      task();
      exit(0);
    }
  }

  int status;
  while (n > 0) {
    wait(&status);
    --n;
  }
#endif
*/
  std::cout << std::flush << "Mass graph alignment - processing " << spectrum_num
      << " of " << spectrum_num << " spectra." << std::endl;

  // combine result files
  std::vector<std::string> input_exts;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) {
    std::string fname = mng_ptr_->output_file_ext_ + "_" + str_util::toString(i);
    input_exts.push_back(fname);
  }

  int top_num = (mng_ptr_->n_unknown_shift_ + 1) * 4;
  PrsmStrMergePtr merge_ptr
      = std::make_shared<PrsmStrMerge>(sp_file_name, input_exts,
                                       mng_ptr_->output_file_ext_, top_num);
  bool normalization = true;
  merge_ptr->process(normalization);
  merge_ptr = nullptr;

  // remove temporary files
  
  for (int t = 0; t < mng_ptr_->thread_num_; t++) {
    file_util::cleanTempFiles(sp_file_name, mng_ptr_->input_file_ext_ + "_" + str_util::toString(t));
    file_util::cleanTempFiles(sp_file_name, mng_ptr_->output_file_ext_ + "_" + str_util::toString(t));
  }
  
}


}  // namespace toppic

