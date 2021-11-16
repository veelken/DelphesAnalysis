#ifndef hhAnalysis_DelphesAnalysis_DelphesJet_h
#define hhAnalysis_DelphesAnalysis_DelphesJet_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h" // DelphesParticle

class DelphesJet : public DelphesParticle
{
 public:
  DelphesJet() = default;
  DelphesJet(const DelphesParticle & particle,
             Int_t btag,
             Int_t btagAlgo,
             Int_t btagPhys,
             Int_t flavor,
             Int_t flavorAlgo,
             Int_t flavorPhys);

  ~DelphesJet();

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */
  Int_t btag() const;
  Int_t btagAlgo() const;
  Int_t btagPhys() const;
  Int_t flavor() const;
  Int_t flavorAlgo() const;
  Int_t flavorPhys() const;

  friend class DelphesJetReader;

 protected:
  Int_t btag_;
  Int_t btagAlgo_;
  Int_t btagPhys_;
  Int_t flavor_;
  Int_t flavorAlgo_;
  Int_t flavorPhys_;
};

std::ostream &
operator<<(std::ostream & stream,
           const DelphesJet & jet);

#endif // hhAnalysis_DelphesAnalysis_DelphesJet_h

