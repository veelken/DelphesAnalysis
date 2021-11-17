#ifndef hhAnalysis_DelphesAnalysis_DelphesEventInfo_h
#define hhAnalysis_DelphesAnalysis_DelphesEventInfo_h

#include <Rtypes.h> // Int_t, Long64_t, Float_t

#include <iostream> // std::ostream

class DelphesEventInfo
{
 public:
  DelphesEventInfo();
  ~DelphesEventInfo();

  Int_t run() const;
  Int_t luminosityBlock() const;
  Long64_t event() const;
  Float_t scale() const;
  Float_t x1() const;
  Float_t x2() const;
  Int_t id1() const;
  Int_t id2() const;
  Float_t alphaQED() const;
  Float_t alphaQCD() const;
  Float_t genweight() const;

  std::string str() const;

  friend class DelphesEventInfoReader;

 protected:
  Int_t     run_;             ///< run number
  Int_t     luminosityBlock_; ///< luminosity section
  Long64_t  event_;           ///< event number
  Float_t   scale_;           ///< renormalization scale
  Float_t   x1_;              ///< Bjorken x for first proton
  Float_t   x2_;              ///< Bjorken x for second proton
  Int_t     id1_;             ///< PDG Id of first incoming parton
  Int_t     id2_;             ///< PDG Id of second incoming parton
  Float_t   alphaQED_;        ///< electroweak coupling
  Float_t   alphaQCD_;        ///< strong coupling
  Float_t   genweight_;       ///< generator-level event weight
};

std::ostream &
operator<<(std::ostream& stream,
           const DelphesEventInfo & eventInfo);

#endif // hhAnalysis_DelphesAnalysis_DelphesEventInfo_h
