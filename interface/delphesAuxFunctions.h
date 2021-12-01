#ifndef hhAnalysis_DelphesAnalysis_delphesAuxFunctions_h
#define hhAnalysis_DelphesAnalysis_delphesAuxFunctions_h

#include "hhAnalysis/DelphesAnalysis/interface/DelphesParticle.h"                 // DelphesParticle::LorentzVector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"                   // DelphesLepton
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                      // DelphesJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h"              // DelphesGenParticle
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h"                   // DelphesGenJet

#include <vector> // std::vector

template <typename T>
bool
isHigherPt(const T * particle1,
           const T * particle2)
{
  return particle1->pt() > particle2->pt();
}

std::vector<const DelphesLepton *>
mergeLeptonCollections(const std::vector<const DelphesLepton *> & electrons,
                       const std::vector<const DelphesLepton *> & muons,
                       bool (*sortFunction)(const DelphesLepton *, const DelphesLepton *));

bool 
isBetterBJet(const DelphesJet * jet1, const DelphesJet * jet2);

template <typename T>
void 
dumpCollection(const std::vector<T>& particles)
{
  for ( typename std::vector<T>::const_iterator particle = particles.begin();
        particle != particles.end(); ++particle ) {
    std::cout << (*particle);
  }
}

const DelphesGenJet*
get_genJet(const DelphesParticle::LorentzVector& jetP4, const std::vector<DelphesGenJet>& genJets);

const DelphesGenParticle*
get_genBQuarkFromHiggs_or_TopDecay(const DelphesParticle::LorentzVector& jetP4, const std::vector<DelphesGenParticle>& genParticles, int mother_pdgId);

#endif // hhAnalysis_DelphesAnalysis_delphesAuxFunctions_h
