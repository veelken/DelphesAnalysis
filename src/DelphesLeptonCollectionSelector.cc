#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionSelector.h"

DelphesLeptonSelector::DelphesLeptonSelector(Era era, int index, bool debug)
  : min_pt_(7.)
  , max_absEta_(2.5)
  , max_relIsoRhoCorr_(0.10)
{}

void
DelphesLeptonSelector::set_min_pt(float min_pt)
{
  min_pt_ = min_pt;
}

void
DelphesLeptonSelector::set_max_absEta(float max_absEta)
{
  max_absEta_ = max_absEta;
}
 
void
DelphesLeptonSelector::set_max_relIsoRhoCorr(float max_relIsoRhoCorr)
{
  max_relIsoRhoCorr_ = max_relIsoRhoCorr;
}

float
DelphesLeptonSelector::get_min_pt() const
{
  return min_pt_;
}

float
DelphesLeptonSelector::get_max_absEta() const
{
  return max_absEta_;
}
 
float
DelphesLeptonSelector::get_max_relIsoRhoCorr() const
{
  return max_relIsoRhoCorr_;
}

bool
DelphesLeptonSelector::operator()(const DelphesLepton & lepton) const
{
  if(debug_)
  {
    std::cout << get_human_line(this, __func__) << ":\n" << lepton;
  }

  if(lepton.pt() < min_pt_)
  {
    if(debug_)
    {
      std::cout << "FAILS pT = " << lepton.pt() << " >= " << min_pt_ << " cut\n";
    }
    return false;
  }
  if(lepton.absEta() > max_absEta_)
  {
    if(debug_)
    {
      std::cout << "FAILS abs(eta) = " << lepton.absEta() << " <= " << max_absEta_ << " cut\n";
    }
    return false;
  }
  if(lepton.relIsoRhoCorr() > max_relIsoRhoCorr_)
  {
    if(debug_)
    {
      std::cout << "FAILS relIsoRhoCorr = " << lepton.relIsoRhoCorr() << " <= " << max_relIsoRhoCorr_ << " cut\n";
    }
    return false;
  }
  return true;
}
