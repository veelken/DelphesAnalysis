#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"

DelphesJet::DelphesJet(const DelphesParticle & particle,
                       Int_t btag,
                       Int_t btagAlgo,
                       Int_t btagPhys,
                       Int_t flavor,
                       Int_t flavorAlgo,
                       Int_t flavorPhys)
  : DelphesParticle(particle)
  , btag_(btag)
  , btagAlgo_(btagAlgo)
  , btagPhys_(btagPhys)
  , flavor_(flavor)
  , flavorAlgo_(flavorAlgo)
  , flavorPhys_(flavorPhys)
{}

DelphesJet::~DelphesJet()
{}

Int_t
DelphesJet::btag() const
{
  return btag_;
}

Int_t
DelphesJet::btagAlgo() const
{
  return btagAlgo_;
}

Int_t
DelphesJet::btagPhys() const
{
  return btagPhys_;
}

Int_t
DelphesJet::flavor() const
{
  return flavor_;
}

Int_t
DelphesJet::flavorAlgo() const
{
  return flavorAlgo_;
}

Int_t
DelphesJet::flavorPhys() const
{
  return flavorPhys_;
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesJet & jet)
{
  stream << static_cast<const DelphesParticle &>(jet) << ",\n"
            " btag = "   << jet.btag()   << ","
            " flavor = " << jet.flavor() << "\n";
  return stream;
}
