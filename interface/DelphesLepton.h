#ifndef hhAnalysis_DelphesAnalysis_DelphesLepton_h
#define hhAnalysis_DelphesAnalysis_DelphesLepton_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h" // DelphesParticle

class DelphesLepton : public DelphesParticle
{
 public:
  DelphesLepton() = default;
  DelphesLepton(const DelphesParticle & particle,
                Int_t charge,
                Float_t dz,
                Float_t relIso,
                Float_t relIsoRhoCorr,
                Float_t sumPt,
                Float_t sumPtCh,
                Float_t sumPtNeu,
                Float_t sumPtCPU);

  ~DelphesLepton();

  /**
   * @brief Funtions to access data-members
   * @return Values of data-members
   */
  Int_t   charge() const;
  Float_t dz() const;
  Float_t relIso() const;
  Float_t relIsoRhoCorr() const;
  Float_t sumPt() const;
  Float_t sumPtCh() const;
  Float_t sumPtNeu() const;
  Float_t sumPtCPU() const;

  friend class DelphesLeptonReader;

 protected:
  Int_t   charge_;
  Float_t dz_;
  Float_t relIso_;
  Float_t relIsoRhoCorr_;
  Float_t sumPt_;
  Float_t sumPtCh_;
  Float_t sumPtNeu_;
  Float_t sumPtCPU_;
};

std::ostream &
operator<<(std::ostream & stream,
           const DelphesLepton & lepton);

#endif // hhAnalysis_DelphesAnalysis_DelphesLepton_h
