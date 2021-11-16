#ifndef hhAnalysis_DelphesAnalysis_DelphesParticle_h
#define hhAnalysis_DelphesAnalysis_DelphesParticle_h

#include <Rtypes.h> // Int_t, Float_t
#include <DataFormats/Math/interface/LorentzVector.h> // math::PtEtaPhiMLorentzVector

class DelphesParticle
{
 public:
  typedef math::PtEtaPhiMLorentzVector LorentzVector;

  DelphesParticle();
  DelphesParticle(Float_t pt, Float_t eta, Float_t phi, Float_t mass);
  DelphesParticle(const DelphesParticle::LorentzVector & p4);

  virtual ~DelphesParticle() {}

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */
  Float_t pt() const;
  Float_t eta() const;
  Float_t phi() const;
  Float_t mass() const;
  Float_t absEta() const;

  const LorentzVector & p4() const;

protected:
  Float_t pt_;       ///< pT of the particle
  Float_t eta_;      ///< eta of the particle
  Float_t phi_;      ///< phi of the particle
  Float_t mass_;     ///< mass of the particle

  Float_t absEta_;   ///< |eta| of the particle

  LorentzVector p4_; ///< 4-momentum constructed from the pT, eta, phi and mass
};

std::ostream &
operator<<(std::ostream & stream,
           const DelphesParticle & particle);

#endif // hhAnalysis_DelphesAnalysis_DelphesParticle_h
