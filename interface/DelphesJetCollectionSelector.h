#ifndef hhAnalysis_DelphesAnalysis_DelphesJetCollectionSelector_h
#define hhAnalysis_DelphesAnalysis_DelphesJetCollectionSelector_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                // DelphesJet

#include "tthAnalysis/HiggsToTauTau/interface/ParticleCollectionSelector.h" // ParticleCollectionSelector
#include "tthAnalysis/HiggsToTauTau/interface/analysisAuxFunctions.h"       // Era

class DelphesJetSelector
{
 public:
  explicit DelphesJetSelector(Era era, int index = -1, bool debug = false);
  ~DelphesJetSelector() {}

  /**
   * @brief Set cut thresholds
   */
  void set_min_pt(float min_pt);
  void set_max_absEta(float max_absEta);

  /**
   * @brief Get cut thresholds
   */
  float get_min_pt() const;
  float get_max_absEta() const;

  /**
   * @brief Check if jet given as function argument passes the selection cuts
   * @return True if jet passes selection; false otherwise
   */
  bool
  operator()(const DelphesJet & jet) const;

protected:
  bool debug_;

  Float_t min_pt_;     ///< lower cut threshold on pT
  Float_t max_absEta_; ///< upper cut threshold on absolute value of eta
};

typedef ParticleCollectionSelector<DelphesJet, DelphesJetSelector> DelphesJetCollectionSelector;

#endif // hhAnalysis_DelphesAnalysis_DelphesJetCollectionSelector_h

