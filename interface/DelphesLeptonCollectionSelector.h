#ifndef hhAnalysis_DelphesAnalysis_DelphesLeptonCollectionSelector_h
#define hhAnalysis_DelphesAnalysis_DelphesLeptonCollectionSelector_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"             // DelphesLepton

#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionSelector.h" // ParticleCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h"       // Era

class DelphesLeptonSelector
{
 public:
  DelphesLeptonSelector(Era era, int index = -1, bool debug = false);
  ~DelphesLeptonSelector() {}

  /**
   * @brief Set cut thresholds
   */
  void set_min_pt(float min_lepton_pt);
  void set_max_absEta(float max_absEta);
  void set_max_relIsoRhoCorr(float max_relIsoRhoCorr);

  /**
   * @brief Get cut thresholds
   */
  float get_min_pt() const;
  float get_max_absEta() const;
  float get_max_relIsoRhoCorr() const;

  /**
   * @brief Check if lepton given as function argument passes selection
   * @return True if lepton passes selection; false otherwise
   */
  bool operator()(const DelphesLepton & lepton) const;
  
 protected:
  bool debug_;

  Float_t min_pt_;            ///< lower cut threshold on reco::GSFElectron pT
  Float_t max_absEta_;        ///< upper cut threshold on absolute value of eta
  Float_t max_relIsoRhoCorr_; ///< upper cut threshold on pileup-corrected relative isolation
};

typedef ParticleCollectionSelector<DelphesLepton, DelphesLeptonSelector> DelphesLeptonCollectionSelector;

#endif // hhAnalysis_DelphesAnalysis_DelphesJetCollectionSelector_h

