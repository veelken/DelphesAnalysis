#ifndef hhAnalysis_DelphesAnalysis_DelphesGenJet_h
#define hhAnalysis_DelphesAnalysis_DelphesGenJet_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h" // DelphesParticle

class DelphesGenJet : public DelphesParticle
{
 public:
  DelphesGenJet() = default;
  DelphesGenJet(const DelphesParticle & particle);

  ~DelphesGenJet();

  friend class DelphesGenJetReader;
};

std::ostream &
operator<<(std::ostream & stream,
           const DelphesGenJet & genJet);

#endif // hhAnalysis_DelphesAnalysis_DelphesGenJet_h
