#ifndef hhAnalysis_DelphesAnalysis_DelphesMEtReader_h
#define hhAnalysis_DelphesAnalysis_DelphesMEtReader_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEt.h" // DelphesMEt

#include "tthAnalysis/HiggsToTauTau/interface/ReaderBase.h"  // ReaderBase

// forward declarations
class TTree;

#include <map>    // std::map
#include <string> // std::string

class DelphesMEtReader : public ReaderBase
{
public:
  DelphesMEtReader();
  DelphesMEtReader(const std::string & branchName_obj);
  ~DelphesMEtReader();

  /**
   * @brief Call tree->SetBranchAddress for all MEt branches
   */
  std::vector<std::string>
  setBranchAddresses(TTree * tree) override;

  /**
   * @brief Read branches from tree and use information to fill MEt object
   * @return RecoMEt object
   */
  DelphesMEt
  read() const;

 protected:
 /**
   * @brief Initialize names of branches to be read from tree
   */
  void
  setBranchNames();

  std::string branchName_obj_;

  std::string branchName_pt_;
  std::string branchName_phi_;

  DelphesMEt met_;

  // CV: make sure that only one DelphesMEtReader instance exists for a given branchName,
  //     as ROOT cannot handle multiple TTree::SetBranchAddress calls for the same branch.
  static std::map<std::string, int> numInstances_;
  static std::map<std::string, DelphesMEtReader *> instances_;
};

#endif // hhAnalysis_DelphesAnalysis_DelphesMEtReader_h
