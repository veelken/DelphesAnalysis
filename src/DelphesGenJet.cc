#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h"

DelphesGenJet::DelphesGenJet(const DelphesParticle & particle)
  : DelphesParticle(particle)
{}

DelphesGenJet::~DelphesGenJet()
{}

std::ostream &
operator<<(std::ostream & stream,
           const DelphesGenJet & genJet)
{
  stream << "GenJet #" << genJet.idx() << ":"
         << " " << static_cast<const DelphesParticle &>(genJet) << "\n";
  return stream;
}
