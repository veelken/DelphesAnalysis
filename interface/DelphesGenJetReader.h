#ifndef hhAnalysis_DelphesAnalysis_DelphesGenJetReader_h
#define hhAnalysis_DelphesAnalysis_DelphesGenJetReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h" // DelphesGenJet

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"     // ReaderBase

// forward declarations
class TTree;

#include <map>    // std::map
#include <string> // std::string

class DelphesGenJetReader : public ReaderBase
{
 public:
  DelphesGenJetReader();
  DelphesGenJetReader(const std::string & branchName_obj);
  ~DelphesGenJetReader();

  /**
   * @brief Call tree->SetBranchAddress for all lepton branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of generator-level particles
   * @return Collection of generator-level particles
   */
  std::vector<DelphesGenJet>
  read() const;

 protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const int max_nGenJets_;
  std::string branchName_num_;
  std::string branchName_obj_;

  std::string branchName_pt_;
  std::string branchName_eta_;
  std::string branchName_phi_;
  std::string branchName_mass_;

  Int_t nGenJets_;
  Float_t * pt_;
  Float_t * eta_;
  Float_t * phi_;
  Float_t * mass_;

  // CV: make sure that only one DelphesGenJetReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, DelphesGenJetReader *> instances_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesGenJetReader_h
