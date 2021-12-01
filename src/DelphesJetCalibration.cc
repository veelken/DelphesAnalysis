#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCalibration.h"

DelphesJetCalibration::DelphesJetCalibration()
{}

DelphesJetCalibration::~DelphesJetCalibration()
{}

DelphesJet 
DelphesJetCalibration::operator()(const DelphesJet& jet)
{
  //std::cout << "<DelphesJetCalibration::operator()>:" << std::endl;
  double pt  = jet.pt();
  double eta = jet.eta();
  //std::cout << " pT = " << pt << ", eta = " << eta << std::endl;
  double sf_part1 = 1.;
  if ( abs(eta) <= 5.00 && pt >= 20.00 && pt <= 1000.00 ) {
    sf_part1 *= ((                   abs(eta) <= 1.40) * (pt >= 20.00 && pt <= 1000.00) * 1./(9.87930e-01 - 2.68848e+00*log(pt)/pow(pt, 1.20713e+00)) + 
                 (abs(eta) > 1.40 && abs(eta) <= 2.40) * (pt >= 20.00 && pt <= 1000.00) * 1./(9.84943e-01 - 3.72537e+00*log(pt)/pow(pt, 1.43131e+00)) + 
                 (abs(eta) > 2.40 && abs(eta) <= 3.00) * (pt >= 20.00 && pt <= 1000.00) * 1./(9.89516e-01 - 6.64317e+00*log(pt)/pow(pt, 1.58751e+00)) + 
                 (abs(eta) > 3.00 && abs(eta) <= 5.00) * (pt >= 20.00 && pt <= 1000.00) * 1./(1.00615e+00 - 5.66821e+03*log(pt)/pow(pt, 3.77179e+00)));
  }
  double sf_part2 = 1.;
  if ( jet.btag() >= 1 )
  {
    pt *= sf_part1;
    if ( abs(eta) <= 5.00 && pt >= 20.00 && pt <= 1000.00 ) {
      sf_part2 *= ((                   abs(eta) <= 1.40) * (pt >= 20.00 && pt <= 1000.00) * 1./(1.02022e+00 - 1.46889e+00*log(pt)/pow(pt, 8.77819e-01)) + 
                   (abs(eta) > 1.40 && abs(eta) <= 2.40) * (pt >= 20.00 && pt <= 1000.00) * 1./(1.02747e+00 - 1.23156e+00*log(pt)/pow(pt, 8.33438e-01)) + 
                   (abs(eta) > 2.40 && abs(eta) <= 3.00) * (pt >= 20.00 && pt <= 1000.00) * 1./(1.03500e+00 - 1.04153e+00*log(pt)/pow(pt, 8.01836e-01)) + 
                   (abs(eta) > 3.00 && abs(eta) <= 5.00) * (pt >= 20.00 && pt <= 1000.00) * 1./(1.13454e+00 - 6.36865e-01*log(pt)/pow(pt, 5.80821e-01)));
    }
  }
  double jetPt_calibrated   = sf_part1*sf_part2*jet.pt();
  //std::cout << "jet pT: uncalibrated = " << jet.pt() << ", calibrated = " << jetPt_calibrated << std::endl;
  // CV: do not correct jet mass for energy carried away by neutrinos produced in heavy flavor decays
  double jetMass_calibrated = sf_part1*jet.mass();
  //std::cout << "jet mass: uncalibrated = " << jet.mass() << ", calibrated = " << jetMass_calibrated << std::endl;
  DelphesJet jet_calibrated(
    { jet.idx(), (float)jetPt_calibrated, jet.eta(), jet.phi(), (float)jetMass_calibrated },
    jet.btag(),
    jet.btagAlgo(),
    jet.btagPhys(),
    jet.flavor(),
    jet.flavorAlgo(),
    jet.flavorPhys());
  //std::cout << jet_calibrated;
  return jet_calibrated;
}

