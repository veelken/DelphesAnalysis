#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJetReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> DelphesGenJetReader::numInstances_;
std::map<std::string, DelphesGenJetReader *> DelphesGenJetReader::instances_;

DelphesGenJetReader::DelphesGenJetReader()
  : DelphesGenJetReader("GenJet")
{}

DelphesGenJetReader::DelphesGenJetReader(const std::string & branchName_obj)
  : max_nGenJets_(256)
  , branchName_num_(Form("n%s", branchName_obj.data()))
  , branchName_obj_(branchName_obj)
  , pt_(nullptr)
  , eta_(nullptr)
  , phi_(nullptr)
  , mass_(nullptr)
{
  setBranchNames();
}

DelphesGenJetReader::~DelphesGenJetReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    DelphesGenJetReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    delete[] gInstance->pt_;
    delete[] gInstance->eta_;
    delete[] gInstance->phi_;
    delete[] gInstance->mass_;
    instances_[branchName_obj_] = nullptr;
  }
}

void
DelphesGenJetReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    branchName_pt_           = Form("%s_%s", branchName_obj_.data(), "pt");
    branchName_eta_          = Form("%s_%s", branchName_obj_.data(), "eta");
    branchName_phi_          = Form("%s_%s", branchName_obj_.data(), "phi");
    branchName_mass_         = Form("%s_%s", branchName_obj_.data(), "mass");
    instances_[branchName_obj_] = this;
  }
  else
  {
    const DelphesGenJetReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    if(branchName_num_ != gInstance->branchName_num_)
    {
      throw cmsException(this)
        << "Association between configuration parameters 'branchName_num' and 'branchName_obj' must be unique:"
        << " present association 'branchName_num' = " << branchName_num_ << " with 'branchName_obj' = " << branchName_obj_
        << " does not match previous association 'branchName_num' = " << instances_[branchName_obj_]->branchName_num_
        << " with 'branchName_obj' = " << instances_[branchName_obj_]->branchName_obj_ << " !!\n";
    }
  }
  ++numInstances_[branchName_obj_];
}

std::vector<std::string>
DelphesGenJetReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree, max_nGenJets_);
    bai.setBranchAddress(nGenJets_, branchName_num_);
    bai.setBranchAddress(pt_,            branchName_pt_);
    bai.setBranchAddress(eta_,           branchName_eta_);
    bai.setBranchAddress(phi_,           branchName_phi_);
    bai.setBranchAddress(mass_,          branchName_mass_);
    return bai.getBoundBranchNames();
  }
  return {};
}

std::vector<DelphesGenJet>
DelphesGenJetReader::read() const
{
  const DelphesGenJetReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);

  std::vector<DelphesGenJet> genJets;
  const Int_t nGenJets = gInstance->nGenJets_;
  if(nGenJets > max_nGenJets_)
  {
    throw cmsException(this)
      << "Number of generator-level jets stored in Ntuple = " << nGenJets << ", "
         "exceeds max_nGenJets = " << max_nGenJets_ << " !!\n";
  }

  if(nGenJets > 0)
  {
    genJets.reserve(nGenJets);
    for(Int_t idxGenJet = 0; idxGenJet < nGenJets; ++idxGenJet)
    {
      genJets.push_back({
        {
          idxGenJet,
          pt_[idxGenJet],
          eta_[idxGenJet],
          phi_[idxGenJet],
          mass_[idxGenJet],
        },
      });
    } // idxGenJet
  } // nGenJets > 0
  return genJets;
}
