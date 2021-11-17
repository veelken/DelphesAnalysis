#ifndef hhAnalysis_DelphesAnalysis_DelphesJetReader_h
#define hhAnalysis_DelphesAnalysis_DelphesJetReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h" // DelphesJet

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"  // ReaderBase

// forward declarations
class TTree;

#include <map>    // std::map
#include <string> // std::string

class DelphesJetReader : public ReaderBase
{
 public:
  DelphesJetReader();
  DelphesJetReader(const std::string & branchName_obj);
  ~DelphesJetReader() override;

  /**
   * @brief Call tree->SetBranchAddress for all jet branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of jet objects
   * @return Collection of jet objects
   */
  std::vector<DelphesJet>
  read() const;

 protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const int max_nJets_;
  std::string branchName_num_;
  std::string branchName_obj_;

  std::string branchName_pt_;
  std::string branchName_eta_;
  std::string branchName_phi_;
  std::string branchName_mass_;
  std::string branchName_btag_;
  std::string branchName_btagAlgo_;
  std::string branchName_btagPhys_;
  std::string branchName_flavor_;
  std::string branchName_flavorAlgo_;
  std::string branchName_flavorPhys_;

  Int_t nJets_;
  Float_t * pt_;
  Float_t * eta_;
  Float_t * phi_;
  Float_t * mass_;
  Float_t * charge_;
  Int_t   * btag_;
  Int_t   * btagAlgo_;
  Int_t   * btagPhys_;
  Int_t   * flavor_;
  Int_t   * flavorAlgo_;
  Int_t   * flavorPhys_;

  // CV: make sure that only one DelphesJetReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, DelphesJetReader *> instances_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesJetReader_h
