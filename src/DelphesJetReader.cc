#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> DelphesJetReader::numInstances_;
std::map<std::string, DelphesJetReader *> DelphesJetReader::instances_;

DelphesJetReader::DelphesJetReader()
  : DelphesJetReader("Jet")
{}

DelphesJetReader::DelphesJetReader(const std::string & branchName_obj)
  : max_nJets_(256)
  , branchName_num_(Form("n%s", branchName_obj.data()))
  , branchName_obj_(branchName_obj)
  , pt_(nullptr)
  , eta_(nullptr)
  , phi_(nullptr)
  , mass_(nullptr)
  , btag_(nullptr)
  , btagAlgo_(nullptr)
  , btagPhys_(nullptr)
  , flavor_(nullptr)
  , flavorAlgo_(nullptr)
  , flavorPhys_(nullptr)
{
  setBranchNames();
}

DelphesJetReader::~DelphesJetReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    DelphesJetReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    delete[] gInstance->pt_;
    delete[] gInstance->eta_;
    delete[] gInstance->phi_;
    delete[] gInstance->mass_;
    delete[] gInstance->btag_;
    delete[] gInstance->btagAlgo_;
    delete[] gInstance->btagPhys_;
    delete[] gInstance->flavor_;
    delete[] gInstance->flavorAlgo_;
    delete[] gInstance->flavorPhys_;
    instances_[branchName_obj_] = nullptr;
  }
}

void
DelphesJetReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    branchName_pt_         = Form("%s_%s", branchName_obj_.data(), "pt");
    branchName_eta_        = Form("%s_%s", branchName_obj_.data(), "eta");
    branchName_phi_        = Form("%s_%s", branchName_obj_.data(), "phi");
    branchName_mass_       = Form("%s_%s", branchName_obj_.data(), "mass");
    branchName_btag_       = Form("%s_%s", branchName_obj_.data(), "btag");
    branchName_btagAlgo_   = Form("%s_%s", branchName_obj_.data(), "btagAlgo");
    branchName_btagPhys_   = Form("%s_%s", branchName_obj_.data(), "btagPhys");
    branchName_flavor_     = Form("%s_%s", branchName_obj_.data(), "flavor");
    branchName_flavorAlgo_ = Form("%s_%s", branchName_obj_.data(), "flavorAlgo");
    branchName_flavorPhys_ = Form("%s_%s", branchName_obj_.data(), "flavorPhys");
    instances_[branchName_obj_] = this;
  }
  else
  {
    const DelphesJetReader * const gInstance = instances_[branchName_obj_];
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
DelphesJetReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree, max_nJets_);
    bai.setBranchAddress(nJets_, branchName_num_);
    bai.setBranchAddress(pt_,         branchName_pt_);
    bai.setBranchAddress(eta_,        branchName_eta_);
    bai.setBranchAddress(phi_,        branchName_phi_);
    bai.setBranchAddress(mass_,       branchName_mass_);
    bai.setBranchAddress(btag_,       branchName_btag_);
    bai.setBranchAddress(btagAlgo_,   branchName_btagAlgo_);
    bai.setBranchAddress(btagPhys_,   branchName_btagPhys_);
    bai.setBranchAddress(flavor_,     branchName_flavor_);
    bai.setBranchAddress(flavorAlgo_, branchName_flavorAlgo_);
    bai.setBranchAddress(flavorPhys_, branchName_flavorPhys_);
    return bai.getBoundBranchNames();
  }
  return {};
}

std::vector<DelphesJet>
DelphesJetReader::read() const
{
  const DelphesJetReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);

  std::vector<DelphesJet> jets;
  const Int_t nJets = gInstance->nJets_;
  if(nJets > max_nJets_)
  {
    throw cmsException(this)
      << "Number of jets stored in Ntuple = " << nJets << ", "
         "exceeds max_nJets = " << max_nJets_ << " !!\n";
  }

  if(nJets > 0)
  {
    jets.reserve(nJets);
    for(Int_t idxJet = 0; idxJet < nJets; ++idxJet)
    {
      jets.push_back({
        {
          idxJet,
          pt_[idxJet],
          eta_[idxJet],
          phi_[idxJet],
          mass_[idxJet],
        },
        btag_[idxJet],
        btagAlgo_[idxJet],
        btagPhys_[idxJet],
        flavor_[idxJet],
        flavorAlgo_[idxJet],
        flavorPhys_[idxJet],
      });
    } // idxJet
  } // nJets > 0
  return jets;
}

