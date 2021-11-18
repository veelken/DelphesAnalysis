#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h"

DelphesGenParticle::DelphesGenParticle(const DelphesParticle & particle,
                                       Int_t pdgId,
                                       Int_t status,
                                       Int_t idxMother,
                                       Int_t idxDaughter1,
                                       Int_t idxDaughter2)
  : DelphesParticle(particle)
  , pdgId_(pdgId)
  , status_(status)
  , idxMother_(idxMother)
  , idxDaughter1_(idxDaughter1)
  , idxDaughter2_(idxDaughter2)
{}

DelphesGenParticle::~DelphesGenParticle()
{}

Int_t 
DelphesGenParticle::pdgId() const
{
  return pdgId_;
}

Int_t 
DelphesGenParticle::status() const
{
  return status_;
}

Int_t 
DelphesGenParticle::idxMother() const
{
  return idxMother_;
}

Int_t 
DelphesGenParticle::idxDaughter1() const
{
  return idxDaughter1_;
}

Int_t 
DelphesGenParticle::idxDaughter2() const
{  
  return idxDaughter2_;
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesGenParticle & genParticle)
{
  stream << "GenParticle #" << genParticle.idx() << ":"
         << " " << static_cast<const DelphesParticle &>(genParticle) << ",\n"
         << " pdgId = "        << genParticle.pdgId()        << ","
         << " status = "       << genParticle.status()       << ","
         << " idxMother = "    << genParticle.idxMother()    << ","
         << " idxDaughter1 = " << genParticle.idxDaughter1() << ","
         << " idxDaughter2 = " << genParticle.idxDaughter2() << "\n";
  return stream;
}
