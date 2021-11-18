#ifndef hhAnalysis_DelphesAnalysis_DelphesLeptonReader_h
#define hhAnalysis_DelphesAnalysis_DelphesLeptonReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h" // DelphesLepton

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"     // ReaderBase

// forward declarations
class TTree;

#include <map>    // std::map
#include <string> // std::string

class DelphesLeptonReader : public ReaderBase
{
 public:
  DelphesLeptonReader();
  DelphesLeptonReader(DelphesLepton::Type type, const std::string & branchName_obj);
  ~DelphesLeptonReader();

  /**
   * @brief Call tree->SetBranchAddress for all lepton branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of lepton objects
   * @return Collection of lepton objects
   */
  std::vector<DelphesLepton>
  read() const;

 protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const int max_nLeptons_;
  std::string branchName_num_;
  std::string branchName_obj_;
  DelphesLepton::Type type_;

  std::string branchName_pt_;
  std::string branchName_eta_;
  std::string branchName_phi_;
  std::string branchName_mass_;
  std::string branchName_charge_;
  std::string branchName_dz_;
  std::string branchName_relIso_;
  std::string branchName_relIsoRhoCorr_;
  std::string branchName_sumPt_;
  std::string branchName_sumPtCh_;
  std::string branchName_sumPtNeu_;
  std::string branchName_sumPtCPU_;

  Int_t nLeptons_;
  Float_t * pt_;
  Float_t * eta_;
  Float_t * phi_;
  Float_t * mass_;
  Int_t   * charge_;
  Float_t * dz_;
  Float_t * relIso_;
  Float_t * relIsoRhoCorr_;
  Float_t * sumPt_;
  Float_t * sumPtCh_;
  Float_t * sumPtNeu_;
  Float_t * sumPtCPU_;

  // CV: make sure that only one DelphesLeptonReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, DelphesLeptonReader *> instances_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesLeptonReader_h
