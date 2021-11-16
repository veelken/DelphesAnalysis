#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEtReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> DelphesMEtReader::numInstances_;
std::map<std::string, DelphesMEtReader *> DelphesMEtReader::instances_;

DelphesMEtReader::DelphesMEtReader()
  : DelphesMEtReader("MET")
{}

DelphesMEtReader::DelphesMEtReader(const std::string & branchName_obj)
  : branchName_obj_(branchName_obj)
{
  setBranchNames();
}

DelphesMEtReader::~DelphesMEtReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    DelphesMEtReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    instances_[branchName_obj_] = nullptr;
  }
}

void
DelphesMEtReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    branchName_pt_  = Form("%s_%s", branchName_obj_.data(), "pt");
    branchName_phi_ = Form("%s_%s", branchName_obj_.data(), "phi");
    instances_[branchName_obj_] = this;
  }
  ++numInstances_[branchName_obj_];
}

std::vector<std::string>
DelphesMEtReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree);
    bai.setBranchAddress(met_.pt_,  branchName_pt_);
    bai.setBranchAddress(met_.phi_, branchName_phi_);

    return bai.getBoundBranchNames();
  }
  return {};
}

DelphesMEt
DelphesMEtReader::read() const
{
  const DelphesMEtReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);
  return met_;
}
