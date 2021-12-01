#include "hhAnalysis/DelphesAnalysis/interface/delphesAuxFunctions.h" 

#include "DataFormats/Math/interface/deltaR.h" // deltaR

#include <algorithm>                           // std::sort

std::vector<const DelphesLepton *>
mergeLeptonCollections(const std::vector<const DelphesLepton *> & electrons,
                       const std::vector<const DelphesLepton *> & muons,
                       bool (*sortFunction)(const DelphesLepton *, const DelphesLepton *))
{
  std::vector<const DelphesLepton *> leptons;
  const std::size_t numLeptons = electrons.size() + muons.size();
  if ( numLeptons > 0 )
  {
    leptons.reserve(numLeptons);
    leptons.insert(leptons.end(), electrons.begin(), electrons.end());
    leptons.insert(leptons.end(), muons.begin(), muons.end());
    std::sort(leptons.begin(), leptons.end(), sortFunction);
  }
  return leptons;
}

bool 
isBetterBJet(const DelphesJet * jet1, const DelphesJet * jet2)
{
  // CV: select jets passing b-tag discriminator first,
  //     then select jets by higher pT
  if ( jet1->btag() > jet2->btag() ) return true;
  if ( jet1->btag() < jet2->btag() ) return false;
  return jet1->pt() > jet2->pt();
}

const DelphesGenJet*
get_genJet(const DelphesParticle::LorentzVector& jetP4, const std::vector<DelphesGenJet>& genJets)
{
  const DelphesGenJet* genJet_bestMatch = nullptr;
  double dR_bestMatch = 1.e+3;
  for ( std::vector<DelphesGenJet>::const_iterator genJet = genJets.begin();
        genJet != genJets.end(); ++genJet ) { 
    double dR = deltaR(jetP4, genJet->p4());
    // CV: require that jet directions agree within dR < 0.1 in order to match reconstructed to generator-level jets
    const double dRmax = 0.1;
    if ( dR < dRmax && dR < dR_bestMatch )
    {
      genJet_bestMatch = &(*genJet);
      dR_bestMatch = dR;
    }
  }
  return genJet_bestMatch;
}

const DelphesGenParticle*
get_genBQuarkFromHiggs_or_TopDecay(const DelphesParticle::LorentzVector& jetP4, const std::vector<DelphesGenParticle>& genParticles, int mother_pdgId)
{
  const DelphesGenParticle* genBQuark_bestMatch = nullptr;
  double dR_bestMatch = 1.e+3;
  int numGenParticles = genParticles.size();
  for ( int idxGenParticle = 0; idxGenParticle < numGenParticles; ++idxGenParticle ) {
    const DelphesGenParticle& genParticle = genParticles[idxGenParticle];
    if ( std::abs(genParticle.pdgId()) == 5 ) {
      int idxMother = genParticle.idxMother();
      if ( idxMother >= 0 && idxMother < numGenParticles ) {
        const DelphesGenParticle& mother = genParticles[idxMother];
        if ( std::abs(mother.pdgId()) == mother_pdgId ) {
          double dR = deltaR(jetP4, genParticle.p4());
          // CV: require that jet directions agree within dR < 0.1 in order to match generator-level jets to b-quarks from Higgs boson or top quark decay
          const double dRmax = 0.1;
          if ( dR < dRmax && dR < dR_bestMatch )
          {
            genBQuark_bestMatch = &genParticle;
            dR_bestMatch = dR;
          }
        }
      }
    }
  }
  return genBQuark_bestMatch;
}
