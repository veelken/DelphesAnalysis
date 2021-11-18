#ifndef hhAnalysis_DelphesAnalysis_DelphesGenParticleReader_h
#define hhAnalysis_DelphesAnalysis_DelphesGenParticleReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h" // DelphesGenParticle

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"          // ReaderBase

// forward declarations
class TTree;

#include <map>    // std::map
#include <string> // std::string

class DelphesGenParticleReader : public ReaderBase
{
 public:
  DelphesGenParticleReader();
  DelphesGenParticleReader(const std::string & branchName_obj);
  ~DelphesGenParticleReader();

  /**
   * @brief Call tree->SetBranchAddress for all lepton branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill collection of generator-level particles
   * @return Collection of generator-level particles
   */
  std::vector<DelphesGenParticle>
  read() const;

 protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  const int max_nGenParticles_;
  std::string branchName_num_;
  std::string branchName_obj_;

  std::string branchName_pt_;
  std::string branchName_eta_;
  std::string branchName_phi_;
  std::string branchName_mass_;
  std::string branchName_pdgId_;
  std::string branchName_status_;
  std::string branchName_idxMother_;
  std::string branchName_idxDaughter1_;
  std::string branchName_idxDaughter2_;

  Int_t nGenParticles_;
  Float_t * pt_;
  Float_t * eta_;
  Float_t * phi_;
  Float_t * mass_;
  Int_t   * pdgId_;
  Int_t   * status_;
  Int_t   * idxMother_;
  Int_t   * idxDaughter1_;
  Int_t   * idxDaughter2_;

  // CV: make sure that only one DelphesGenParticleReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, DelphesGenParticleReader *> instances_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesGenParticleReader_h
