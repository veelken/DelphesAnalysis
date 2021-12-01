#ifndef hhAnalysis_DelphesAnalysis_DelphesJetCalibration_h
#define hhAnalysis_DelphesAnalysis_DelphesJetCalibration_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h" // DelphesJet

class DelphesJetCalibration
{
 public:
  DelphesJetCalibration();
  ~DelphesJetCalibration();

  DelphesJet operator()(const DelphesJet& jet);
};

#endif // hhAnalysis_DelphesAnalysis_DelphesJetCalibration_h

