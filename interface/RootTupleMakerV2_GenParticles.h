#ifndef RootTupleMakerV2GenParticles
#define RootTupleMakerV2GenParticles

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class RootTupleMakerV2_GenParticles : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenParticles(const edm::ParameterSet&);
  //
  //Tau Specific functions
  reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>& , int );
  reco::Candidate::LorentzVector getInvisMomentum(const std::vector<const reco::GenParticle*>& , int );
  void  findDaughters(const reco::GenParticle*, std::vector<const reco::GenParticle*>&, int = -1);
  bool  isNeutrino(const reco::GenParticle*);
  int   getGenTauDecayMode(const reco::GenParticle*);
  void  countDecayProducts(const reco::GenParticle*,int&,int&,int&,int&,int&,int&,int&,int&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
};

#endif
