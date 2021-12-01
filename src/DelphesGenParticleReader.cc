#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticleReader.h"

#include "tthAnalysis/HiggsToTauTau/interface/cmsException.h"             // cmsException()
#include "tthAnalysis/HiggsToTauTau/interface/BranchAddressInitializer.h" // BranchAddressInitializer, TTree, Form()

std::map<std::string, int> DelphesGenParticleReader::numInstances_;
std::map<std::string, DelphesGenParticleReader *> DelphesGenParticleReader::instances_;

DelphesGenParticleReader::DelphesGenParticleReader()
  : DelphesGenParticleReader("GenPart")
{}

DelphesGenParticleReader::DelphesGenParticleReader(const std::string & branchName_obj)
  : max_nGenParticles_(2048)
  , branchName_num_(Form("n%s", branchName_obj.data()))
  , branchName_obj_(branchName_obj)
  , pt_(nullptr)
  , eta_(nullptr)
  , phi_(nullptr)
  , mass_(nullptr)
  , pdgId_(nullptr)
  , status_(nullptr)
  , idxMother_(nullptr)
  , idxDaughter1_(nullptr)
  , idxDaughter2_(nullptr)
{
  setBranchNames();
}

DelphesGenParticleReader::~DelphesGenParticleReader()
{
  --numInstances_[branchName_obj_];
  assert(numInstances_[branchName_obj_] >= 0);
  if(numInstances_[branchName_obj_] == 0)
  {
    DelphesGenParticleReader * const gInstance = instances_[branchName_obj_];
    assert(gInstance);
    delete[] gInstance->pt_;
    delete[] gInstance->eta_;
    delete[] gInstance->phi_;
    delete[] gInstance->mass_;
    delete[] gInstance->pdgId_;
    delete[] gInstance->status_;
    delete[] gInstance->idxMother_;
    delete[] gInstance->idxDaughter1_;
    delete[] gInstance->idxDaughter2_;
    instances_[branchName_obj_] = nullptr;
  }
}

void
DelphesGenParticleReader::setBranchNames()
{
  if(numInstances_[branchName_obj_] == 0)
  {
    branchName_pt_           = Form("%s_%s", branchName_obj_.data(), "pt");
    branchName_eta_          = Form("%s_%s", branchName_obj_.data(), "eta");
    branchName_phi_          = Form("%s_%s", branchName_obj_.data(), "phi");
    branchName_mass_         = Form("%s_%s", branchName_obj_.data(), "mass");
    branchName_pdgId_        = Form("%s_%s", branchName_obj_.data(), "pdgId");
    branchName_status_       = Form("%s_%s", branchName_obj_.data(), "status");
    branchName_idxMother_    = Form("%s_%s", branchName_obj_.data(), "genPartIdxMother");
    branchName_idxDaughter1_ = Form("%s_%s", branchName_obj_.data(), "d1");
    branchName_idxDaughter2_ = Form("%s_%s", branchName_obj_.data(), "d2");
    instances_[branchName_obj_] = this;
  }
  else
  {
    const DelphesGenParticleReader * const gInstance = instances_[branchName_obj_];
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
DelphesGenParticleReader::setBranchAddresses(TTree * tree)
{
  if(instances_[branchName_obj_] == this)
  {
    BranchAddressInitializer bai(tree, max_nGenParticles_);
    bai.setBranchAddress(nGenParticles_, branchName_num_);
    bai.setBranchAddress(pt_,            branchName_pt_);
    bai.setBranchAddress(eta_,           branchName_eta_);
    bai.setBranchAddress(phi_,           branchName_phi_);
    bai.setBranchAddress(mass_,          branchName_mass_);
    bai.setBranchAddress(pdgId_,         branchName_pdgId_);
    bai.setBranchAddress(status_,        branchName_status_);
    bai.setBranchAddress(idxMother_,     branchName_idxMother_);
    bai.setBranchAddress(idxDaughter1_,  branchName_idxDaughter1_);
    bai.setBranchAddress(idxDaughter2_,  branchName_idxDaughter2_);
    return bai.getBoundBranchNames();
  }
  return {};
}

std::vector<DelphesGenParticle>
DelphesGenParticleReader::read() const
{
  const DelphesGenParticleReader * const gInstance = instances_[branchName_obj_];
  assert(gInstance);

  std::vector<DelphesGenParticle> genParticles;
  const Int_t nGenParticles = gInstance->nGenParticles_;
  if(nGenParticles > max_nGenParticles_)
  {
    throw cmsException(this)
      << "Number of generator-level particles stored in Ntuple = " << nGenParticles << ", "
         "exceeds max_nGenParticles = " << max_nGenParticles_ << " !!\n";
  }

  if(nGenParticles > 0)
  {
    genParticles.reserve(nGenParticles);
    for(Int_t idxGenParticle = 0; idxGenParticle < nGenParticles; ++idxGenParticle)
    {
      genParticles.push_back({
        {
          idxGenParticle,
          pt_[idxGenParticle],
          eta_[idxGenParticle],
          phi_[idxGenParticle],
          mass_[idxGenParticle],
        },
        pdgId_[idxGenParticle],
        status_[idxGenParticle],
        idxMother_[idxGenParticle],
        idxDaughter1_[idxGenParticle],
        idxDaughter2_[idxGenParticle],
      });
    } // idxGenParticle
  } // nGenParticles > 0
  return genParticles;
}
