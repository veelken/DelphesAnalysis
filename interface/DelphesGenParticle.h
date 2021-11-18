#ifndef hhAnalysis_DelphesAnalysis_DelphesGenParticle_h
#define hhAnalysis_DelphesAnalysis_DelphesGenParticle_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h" // DelphesParticle

class DelphesGenParticle : public DelphesParticle
{
 public:
  DelphesGenParticle() = default;
  DelphesGenParticle(const DelphesParticle & particle,
                     Int_t pdgId,
                     Int_t status,
                     Int_t idxMother,
                     Int_t idxDaughter1,
                     Int_t idxDaughter2);

  ~DelphesGenParticle();

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */
  Int_t pdgId() const;
  Int_t status() const;
  Int_t idxMother() const;
  Int_t idxDaughter1() const;
  Int_t idxDaughter2() const;

  friend class DelphesGenParticleReader;

 protected:
  Int_t pdgId_;
  Int_t status_;
  Int_t idxMother_;
  Int_t idxDaughter1_;
  Int_t idxDaughter2_;
};

std::ostream &
operator<<(std::ostream & stream,
           const DelphesGenParticle & genParticle);

#endif // hhAnalysis_DelphesAnalysis_DelphesGenParticle_h
