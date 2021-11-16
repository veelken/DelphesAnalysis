#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEt.h"

#include <math.h>  // cos, sin
#include <iomanip> // std::endl

DelphesMEt::DelphesMEt()
  : DelphesMEt(0., 0.)
{}

DelphesMEt::DelphesMEt(Float_t pt, Float_t phi)
  : pt_(pt)
  , phi_(phi)
{}

DelphesMEt &
DelphesMEt::operator=(const DelphesMEt & other)
{
  pt_ = other.pt_;
  phi_ = other.phi_;
  return *this;
}

Float_t
DelphesMEt::pt() const
{
  return pt_;
}

Float_t
DelphesMEt::phi() const
{
  return phi_;
}

Float_t
DelphesMEt::px() const
{
  return pt_*cos(phi_);
}

Float_t
DelphesMEt::py() const
{
  return pt_*sin(phi_);
}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesMEt & met)
{
  stream << "pT = " << met.pt() << ", phi = " << met.phi() << std::endl;
  return stream;
}
