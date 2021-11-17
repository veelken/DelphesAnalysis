#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionSelector.h"

DelphesJetSelector::DelphesJetSelector(Era era, int index, bool debug)
  : debug_(debug)
  , min_pt_(25.)
  , max_absEta_(2.4)
{}

void
DelphesJetSelector::set_min_pt(float min_pt)
{
  min_pt_ = min_pt;
}

void
DelphesJetSelector::set_max_absEta(float max_absEta)
{
  max_absEta_ = max_absEta;
}

float
DelphesJetSelector::get_min_pt() const
{
  return min_pt_;
}

float
DelphesJetSelector::get_max_absEta() const
{
  return max_absEta_;
}

bool
DelphesJetSelector::operator()(const DelphesJet & jet) const
{
  if(debug_)
  {
    std::cout << get_human_line(this, __func__) << ":\n jet: " << jet << '\n';
  }

  if(jet.pt() < min_pt_)
  {
    if(debug_)
    {
      std::cout << "FAILS pT = " << jet.pt() << " >= " << min_pt_ << " cut\n";
    }
    return false;
  }
  if(max_absEta_ > 0. && jet.absEta() > max_absEta_)
  {
    if(debug_)
    {
      std::cout << "FAILS abs(eta) = " << jet.absEta() << " <= " << max_absEta_ << " cut\n";
    }
    return false;
  }
  return true;
}
