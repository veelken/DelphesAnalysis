#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h" 

#include <cmath> // std::fabs()

DelphesParticle::DelphesParticle()
  : DelphesParticle(0., 0., 0., 0.)
{}

DelphesParticle::DelphesParticle(Float_t pt, Float_t eta, Float_t phi, Float_t mass)
  : pt_(pt)
  , eta_(eta)
  , phi_(phi)
  , mass_(mass)
  , absEta_(std::fabs(eta_))
  , p4_{pt_, eta_, phi_, mass_}
{}

DelphesParticle::DelphesParticle(const DelphesParticle::LorentzVector & p4)
  : pt_(p4.pt())
  , eta_(p4.eta())
  , phi_(p4.phi())
  , mass_(p4.mass())
  , absEta_(std::fabs(eta_))
  , p4_(p4)
{}

Float_t
DelphesParticle::pt() const
{
  return pt_;
}

Float_t
DelphesParticle::eta() const
{
  return eta_;
}

Float_t
DelphesParticle::phi() const
{
  return phi_;
}

Float_t
DelphesParticle::mass() const
{
  return mass_;
}

Float_t
DelphesParticle::absEta() const
{
  return absEta_;
}

const DelphesParticle::LorentzVector &
DelphesParticle::p4() const
{
  return p4_;
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesParticle & particle)
{
  stream << " pT = "   << particle.pt()          << ","
            " eta = "  << particle.eta()         << ","
            " phi = "  << particle.phi()         << ","
            " mass = " << particle.mass()        << ","
            " E = "    << particle.p4().energy() << ","
            " |p| = "  << particle.p4().P();
  return stream;
}
