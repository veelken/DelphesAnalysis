#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"

DelphesLepton::DelphesLepton(const DelphesParticle & particle,
                             Int_t charge,
                             Float_t dz,
                             Float_t relIso,
                             Float_t relIsoRhoCorr,
                             Float_t sumPt,
                             Float_t sumPtCh,
                             Float_t sumPtNeu,
                             Float_t sumPtCPU)
  : DelphesParticle(particle)
  , charge_(charge)
  , dz_(dz)
  , relIso_(relIso)
  , relIsoRhoCorr_(relIsoRhoCorr)
  , sumPt_(sumPt)
  , sumPtCh_(sumPtCh)
  , sumPtNeu_(sumPtNeu)
  , sumPtCPU_(sumPtCPU)
{}

DelphesLepton::~DelphesLepton()
{}

Int_t
DelphesLepton::charge() const
{
  return charge_;
}
  
Float_t
DelphesLepton::dz() const
{
  return dz_;
}

Float_t
DelphesLepton::relIso() const
{
  return relIso_;
}

Float_t 
DelphesLepton::relIsoRhoCorr() const
{
  return relIsoRhoCorr_;
}

Float_t
DelphesLepton::sumPt() const
{
  return sumPt_;
}

Float_t
DelphesLepton::sumPtCh() const
{
  return sumPtCh_;
}

Float_t
DelphesLepton::sumPtNeu() const
{
  return sumPtNeu_;
}

Float_t
DelphesLepton::sumPtCPU() const
{
  return sumPtCPU_;
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesLepton & lepton)
{
  stream << static_cast<const DelphesParticle &>(lepton) << ",\n"
            " dz = "            << lepton.dz()            << ","
            " relIso = "        << lepton.relIso()        << ","
            " relIsoRhoCorr = " << lepton.relIsoRhoCorr() << "\n";
  return stream;
}
