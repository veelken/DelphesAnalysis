#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h"

#include <iomanip> // std::endl;

DelphesEventInfo::DelphesEventInfo()
{}

DelphesEventInfo::~DelphesEventInfo()
{}

Int_t
DelphesEventInfo::run() const
{
  return run_;
}

Int_t
DelphesEventInfo::luminosityBlock() const
{
  return luminosityBlock_;
}

Long64_t 
DelphesEventInfo::event() const
{
  return event_;
}

Float_t
DelphesEventInfo::scale() const
{
  return scale_;
}

Float_t
DelphesEventInfo::x1() const
{
  return x1_;
}

Float_t
DelphesEventInfo::x2() const
{
  return x2_;
}

Int_t 
DelphesEventInfo::id1() const
{
  return id1_;
}

Int_t
DelphesEventInfo::id2() const
{
  return id2_;
}

Float_t
DelphesEventInfo::alphaQED() const
{
  return alphaQED_;
}

Float_t
DelphesEventInfo::alphaQCD() const
{
  return alphaQCD_;
}

Float_t 
DelphesEventInfo::genweight() const
{
  return genweight_;
}

std::string
DelphesEventInfo::str() const
{
  std::stringstream ss;
  ss << run_ << ':' << luminosityBlock_ << ':' << event_;
  return ss.str();
}

std::ostream &
operator<<(std::ostream & os,
           const DelphesEventInfo & eventInfo)
{
  os << "run = " << eventInfo.run() << ", ls = " << eventInfo.luminosityBlock() << ", event = " << eventInfo.event() << std::endl;
  return os;
}
