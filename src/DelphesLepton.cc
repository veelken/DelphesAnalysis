#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"

DelphesLepton::DelphesLepton(const DelphesParticle & particle,
                             Type type,
                             Int_t charge,
                             Float_t dz,
                             Float_t relIso,
                             Float_t relIsoRhoCorr,
                             Float_t sumPt,
                             Float_t sumPtCh,
                             Float_t sumPtNeu,
                             Float_t sumPtCPU)
  : DelphesParticle(particle)
  , type_(type)
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

Bool_t
DelphesLepton::is_electron() const
{
  return type_ == kElectron;
}

Bool_t
DelphesLepton::is_muon() const
{
  return type_ == kMuon;
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesLepton & lepton)
{
  if      ( lepton.is_electron() ) stream << "Electron";
  else if ( lepton.is_muon()     ) stream << "Muon";
  else                             stream << "Lepton";
  stream << " #" << lepton.idx() << ":"
         << " " << static_cast<const DelphesParticle &>(lepton) << ",\n"
         << " dz = "            << lepton.dz()            << ","
         << " relIso = "        << lepton.relIso()        << ","
         << " relIsoRhoCorr = " << lepton.relIsoRhoCorr() << "\n";
  return stream;
}
