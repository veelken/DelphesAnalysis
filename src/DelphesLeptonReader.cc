#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> DelphesLeptonReader::numInstances_;
std::map<std::string, DelphesLeptonReader *> DelphesLeptonReader::instances_;

DelphesLeptonReader::DelphesLeptonReader()
  : DelphesLeptonReader("Lepton")
{}

DelphesLeptonReader::DelphesLeptonReader(const std::string & branchName_obj)
  : max_nLeptons_(64)
  , branchName_num_(Form("n%s", branchName_obj.data()))
  , branchName_obj_(branchName_obj)
  , pt_(nullptr)
  , eta_(nullptr)
  , phi_(nullptr)
  , mass_(nullptr)
  , charge_(nullptr)
  , dz_(nullptr)
  , relIso_(nullptr)
  , relIsoRhoCorr_(nullptr)
  , sumPt_(nullptr)
  , sumPtCh_(nullptr)
  , sumPtNeu_(nullptr)
  , sumPtCPU_(nullptr)
{
  setBranchNames();
}

DelphesLeptonReader::~DelphesLeptonReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    DelphesLeptonReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    delete[] gInstance->pt_;
    delete[] gInstance->eta_;
    delete[] gInstance->phi_;
    delete[] gInstance->mass_;
    delete[] gInstance->charge_;
    delete[] gInstance->dz_;
    delete[] gInstance->relIso_;
    delete[] gInstance->relIsoRhoCorr_;
    delete[] gInstance->sumPt_;
    delete[] gInstance->sumPtCh_;
    delete[] gInstance->sumPtNeu_;
    delete[] gInstance->sumPtCPU_;
    instances_[branchName_obj_] = nullptr;
  }
}

void
DelphesLeptonReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    branchName_pt_            = Form("%s_%s", branchName_obj_.data(), "pt");
    branchName_eta_           = Form("%s_%s", branchName_obj_.data(), "eta");
    branchName_phi_           = Form("%s_%s", branchName_obj_.data(), "phi");
    branchName_mass_          = Form("%s_%s", branchName_obj_.data(), "mass");
    branchName_charge_        = Form("%s_%s", branchName_obj_.data(), "charge");
    branchName_dz_            = Form("%s_%s", branchName_obj_.data(), "dz");
    branchName_relIso_        = Form("%s_%s", branchName_obj_.data(), "relIso");
    branchName_relIsoRhoCorr_ = Form("%s_%s", branchName_obj_.data(), "relIsoRhoCorr");
    branchName_sumPt_         = Form("%s_%s", branchName_obj_.data(), "sumPt");
    branchName_sumPtCh_       = Form("%s_%s", branchName_obj_.data(), "sumPtCh");
    branchName_sumPtNeu_      = Form("%s_%s", branchName_obj_.data(), "sumPtNeu");
    branchName_sumPtCPU_      = Form("%s_%s", branchName_obj_.data(), "sumPtCPU");
  }
  else
  {
    const DelphesLeptonReader * const gInstance = instances_[branchName_obj_];
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
DelphesLeptonReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree, max_nLeptons_);
    bai.setBranchAddress(nLeptons_, branchName_num_);
    bai.setBranchAddress(pt_,            branchName_pt_);
    bai.setBranchAddress(eta_,           branchName_eta_);
    bai.setBranchAddress(phi_,           branchName_phi_);
    bai.setBranchAddress(mass_,          branchName_mass_);
    bai.setBranchAddress(charge_,        branchName_charge_);
    bai.setBranchAddress(dz_,            branchName_dz_);
    bai.setBranchAddress(relIso_,        branchName_relIso_);
    bai.setBranchAddress(relIsoRhoCorr_, branchName_relIsoRhoCorr_);
    bai.setBranchAddress(sumPt_,         branchName_sumPt_);
    bai.setBranchAddress(sumPtCh_,       branchName_sumPtCh_);
    bai.setBranchAddress(sumPtNeu_,      branchName_sumPtNeu_);
    bai.setBranchAddress(sumPtCPU_,      branchName_sumPtCPU_);
    return bai.getBoundBranchNames();
  }
  return {};
}

std::vector<DelphesLepton>
DelphesLeptonReader::read() const
{
  const DelphesLeptonReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);

  std::vector<DelphesLepton> leptons;
  const UInt_t nLeptons = gInstance->nLeptons_;
  if(nLeptons > max_nLeptons_)
  {
    throw cmsException(this)
      << "Number of leptons stored in Ntuple = " << nLeptons << ", "
         "exceeds max_nLeptons = " << max_nLeptons_ << " !!\n";
  }

  if(nLeptons > 0)
  {
    leptons.reserve(nLeptons);
    for(UInt_t idxLepton = 0; idxLepton < nLeptons; ++idxLepton)
    {
      leptons.push_back({
        {
          pt_[idxLepton],
          eta_[idxLepton],
          phi_[idxLepton],
          mass_[idxLepton],
        },
        charge_[idxLepton],
        dz_[idxLepton],
        relIso_[idxLepton],
        relIsoRhoCorr_[idxLepton],
        sumPt_[idxLepton],
        sumPtCh_[idxLepton],
        sumPtNeu_[idxLepton],
        sumPtCPU_[idxLepton],
      });
    } // idxLepton
  } // nLeptons > 0
  return leptons;
}
