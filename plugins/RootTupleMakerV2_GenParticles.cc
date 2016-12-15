#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_GenParticles.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "TMath.h"

RootTupleMakerV2_GenParticles::RootTupleMakerV2_GenParticles(const edm::ParameterSet& iConfig) :
  genPartInputToken_ (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("InputTag"))),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta"          + suffix );
  produces <std::vector<double> > ( prefix + "Phi"          + suffix );
  produces <std::vector<double> > ( prefix + "P"            + suffix );
  produces <std::vector<double> > ( prefix + "Px"           + suffix );
  produces <std::vector<double> > ( prefix + "Py"           + suffix );
  produces <std::vector<double> > ( prefix + "Pz"           + suffix );
  produces <std::vector<double> > ( prefix + "Pt"           + suffix );
  produces <std::vector<double> > ( prefix + "Energy"       + suffix );
  produces <std::vector<double> > ( prefix + "Mass"         + suffix );
  produces <std::vector<int> >    ( prefix + "PdgId"        + suffix );
  produces <std::vector<double> > ( prefix + "VX"           + suffix );
  produces <std::vector<double> > ( prefix + "VY"           + suffix );
  produces <std::vector<double> > ( prefix + "VZ"           + suffix );
  produces <std::vector<int> >    ( prefix + "NumDaught"    + suffix );
  produces <std::vector<int> >    ( prefix + "Status"       + suffix );
  produces <std::vector<int> >    ( prefix + "MotherIndex"  + suffix );
  produces <std::vector<int> >    ( prefix + "TauDecayMode" + suffix );
  produces <std::vector<double> > ( prefix + "TauVisiblePt" + suffix );
  produces <std::vector<double> > ( prefix + "TauVisibleEta"+ suffix );
  produces <std::vector<double> > ( prefix + "TauVisiblePhi"+ suffix );
  // some genstatusflags
  // see: https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/HepMCCandidate/interface/GenParticle.h
  produces <std::vector<bool> >   ( prefix + "IsPromptFinalState"+ suffix );
  produces <std::vector<bool> >   ( prefix + "IsPromptDecayed"+ suffix );
  produces <std::vector<bool> >   ( prefix + "IsHardProcess"+ suffix );
  produces <std::vector<bool> >   ( prefix + "FromHardProcessFinalState"+ suffix );
  produces <std::vector<bool> >   ( prefix + "FromHardProcessDecayed"+ suffix );
  produces <std::vector<bool> >   ( prefix + "IsLastCopy"+ suffix );
  // Top Pt reweight
  produces <double>               ( prefix + "TopPtWeight"  + suffix );
  // w/z system pt
  produces <double>               ( prefix + "WorZSystemPt"  + suffix );
}

void RootTupleMakerV2_GenParticles::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  p    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  px   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  py   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pz   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  mass  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     pdgId   ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  vx  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vz  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     numDaught  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     status     ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     motherIndex   ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     taudecaymode  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  tauvisiblept  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  tauvisibleeta ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  tauvisiblephi ( new std::vector<double>()  );
  // some genstatusflags
  std::auto_ptr<std::vector<bool> >    isPromptFinalState        ( new std::vector<bool>());
  std::auto_ptr<std::vector<bool> >    isPromptDecayed           ( new std::vector<bool>());
  std::auto_ptr<std::vector<bool> >    isHardProcess             ( new std::vector<bool>());
  std::auto_ptr<std::vector<bool> >    fromHardProcessFinalState ( new std::vector<bool>());
  std::auto_ptr<std::vector<bool> >    fromHardProcessDecayed    ( new std::vector<bool>());
  std::auto_ptr<std::vector<bool> >    isLastCopy                ( new std::vector<bool>());
  std::auto_ptr<double >               topptweight ( new double() );
  *topptweight.get() = 1.0;
  // w/z system Pt
  std::auto_ptr<double >               worzsystempt ( new double() );
  *worzsystempt.get() = -999;
  // not kept in ntuple
  std::auto_ptr<std::vector<int> >     indexInGenPColl( new std::vector<int>()  );
  
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genPartInputToken_, genParticles);

    if( genParticles.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenParticlesInfo") << "Total # GenParticles: " << genParticles->size();

      int nTops = 0;
      float topPt1 = -1;
      float topPt2 = -1;
      std::vector<LorentzVector> hardProcessLeptonLorentzVectors;
      std::vector<float> wOrZGammaPts;
      std::vector<bool> wOrZGammaIsHardProcess;
      std::vector<bool> wOrZGammaIsLastCopy;
      for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
        // keep for topptreweight below
        // see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        if(isTop(&(*it)) && it->status()==62) {
          nTops++;
          if(topPt1<0)
            topPt1 = it->pt();
          else
            topPt2 = it->pt();
        }
        // get the lorentzvectors for leptons
        if(it->isHardProcess()) {
          if(abs(it->pdgId())>=11 && abs(it->pdgId())<=18) // lepton
            hardProcessLeptonLorentzVectors.push_back(it->p4());
        }
        if(abs(it->pdgId())>=22 && abs(it->pdgId())<=24) { // W/Z/gamma
          wOrZGammaPts.push_back(it->pt());
          wOrZGammaIsHardProcess.push_back(it->isHardProcess());
          wOrZGammaIsLastCopy.push_back(it->isLastCopy());
        }

        // don't store any more genparticles when we reach the limit
        // but continue to look for W/Z/gamma and tops (above)
        if(eta->size() >= maxSize)
          continue;

        //don't store LHC protons - pdgId=2212, extremely high eta, phi=0, energy ~= 6500
        if(it->pdgId()==2212 && it->energy()>=6000) continue;
        // actually let's not store any particle with huge eta
        if(fabs(it->eta()) > 10) continue;
        // don't store pythia8 "partons in preparation of hadronization process"
        // details here: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
        if(it->status() >= 71 && it->status()<=79) continue;

        //if(fabs(it->pdgId())==42)
        //  edm::LogInfo("RootTupleMakerV2_GenParticlesInfo") << "LQ found with status: " << it->status();

        // fill in all the vectors
        eta->push_back( it->eta() );
        phi->push_back( it->phi() );
        p->push_back( it->p() );
        px->push_back( it->px() );
        py->push_back( it->py() );
        pz->push_back( it->pz() );
        pt->push_back( it->pt() );
        energy->push_back( it->energy() );
        mass->push_back(it->mass() );
        pdgId->push_back( it->pdgId() );
        vx->push_back( it->vx() );
        vy->push_back( it->vy() );
        vz->push_back( it->vz() );
        numDaught->push_back( it->numberOfDaughters() );
        status->push_back( it->status() );
        indexInGenPColl->push_back(std::distance(genParticles->begin(),it));
        // some genstatus flags
        isPromptFinalState        ->push_back(it->isPromptFinalState());
        isPromptDecayed           ->push_back(it->isPromptDecayed());
        isHardProcess             ->push_back(it->isHardProcess());
        fromHardProcessFinalState ->push_back(it->fromHardProcessFinalState());
        fromHardProcessDecayed    ->push_back(it->fromHardProcessDecayed());
        isLastCopy                ->push_back(it->isLastCopy());

        // >>>>>>  if gen-particle is a tau, check decay mode and fill-in visible momentum parameters:
        int getGenTauDecayMode_ = 0; double tauVisPt  = -999.0; double tauVisEta = -999.0; double tauVisPhi = -999.0;
        if( abs(it->pdgId()) == 15 ){
          getGenTauDecayMode_ =  getGenTauDecayMode( & (*it) );
          if( getGenTauDecayMode_>0 ){//get visible momentum only if decay mode is determined
            std::vector<const reco::GenParticle*> TauDaughters;
            findDaughters( & (*it), TauDaughters, 1);
            LorentzVector tauVis   = getVisMomentum( TauDaughters, 1);
            tauVisPt = tauVis.pt(); tauVisEta = tauVis.eta(); tauVisPhi = tauVis.phi();
            // ------------------------------------------------------------------ DEBUG OUTPUT --------------------------------------------------------------//
            //LorentzVector tauInvis = getInvisMomentum( TauDaughters, 1);
            //std::cout<< "Tau DecayMode/NoOfDaughters: "<<getGenTauDecayMode_<<" / "<<TauDaughters.size()<<std::endl;
            //std::cout<< "          visiblePt/Eta/Phi: "<<tauVisPt<<" / "<<tauVisEta<<" / "<<tauVisPhi <<std::endl;
            //std::cout<< "        invisiblePt/Eta/Phi: "<<tauInvis.pt()<<" / "<<tauInvis.eta()<<" / "<<tauInvis.phi() <<std::endl;
            //std::cout<< "visible+invisiblePt/Eta/Phi: "<<(tauVis+tauInvis).pt()<<" / "<<(tauVis+tauInvis).eta()<<" / "<<(tauVis+tauInvis).phi() <<std::endl;
            //std::cout<< "              genPt/Eta/Phi: "<<it->pt()<<" / "<< it->eta()<<" / "<<it->phi()<<std::endl;
            // ------------------------------------------------------------------ ------------ --------------------------------------------------------------//
          }
        }
        taudecaymode -> push_back( (int)(getGenTauDecayMode_) );
        tauvisiblept -> push_back( tauVisPt  );
        tauvisibleeta-> push_back( tauVisEta );
        tauvisiblephi-> push_back( tauVisPhi );

        int idx = -1;
        for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ) {
          if( it->mother()==&(*mit) ) {
            idx = std::distance(genParticles->begin(),mit);
            break;
          }
        }
        // we have the mother inside the genparticle collection
        // now find the index of mother in the kept/ntuple genparticle collection
        int indexInKeptGenPColl = -1;
        std::vector<int>::iterator genPCollItr = std::find(indexInGenPColl->begin(),indexInGenPColl->end(),idx);
        if(genPCollItr != indexInGenPColl->end())
          indexInKeptGenPColl = std::distance(indexInGenPColl->begin(),genPCollItr);
        motherIndex->push_back( indexInKeptGenPColl );
      } // end of loop over genParticles

      // for topptreweight
      if(nTops==2)
        *topptweight.get() = sqrt(exp(0.156-0.00137*topPt1)*exp(0.156-0.00137*topPt2));

      // for W/Z Pt
      // we always have 2 hard process leptons in DY or W samples; in other samples we may not
      if(hardProcessLeptonLorentzVectors.size()==2) {
        *worzsystempt.get() = (hardProcessLeptonLorentzVectors[0].pt()+hardProcessLeptonLorentzVectors[1].pt());
      }

    } else { // invalid genParticles handle
      edm::LogError("RootTupleMakerV2_GenParticlesError") << "Error! Can't get the genPartInputToken_";
    }
  }

  // put vectors in the event
  iEvent.put( eta,          prefix + "Eta"          + suffix );
  iEvent.put( phi,          prefix + "Phi"          + suffix );
  iEvent.put( p,            prefix + "P"            + suffix );
  iEvent.put( px,           prefix + "Px"           + suffix );
  iEvent.put( py,           prefix + "Py"           + suffix );
  iEvent.put( pz,           prefix + "Pz"           + suffix );
  iEvent.put( pt,           prefix + "Pt"           + suffix );
  iEvent.put( energy,       prefix + "Energy"       + suffix );
  iEvent.put( mass,         prefix + "Mass"       + suffix );
  iEvent.put( pdgId,        prefix + "PdgId"        + suffix );
  iEvent.put( vx,           prefix + "VX"           + suffix );
  iEvent.put( vy,           prefix + "VY"           + suffix );
  iEvent.put( vz,           prefix + "VZ"           + suffix );
  iEvent.put( numDaught,    prefix + "NumDaught"    + suffix );
  iEvent.put( status,       prefix + "Status"       + suffix );
  iEvent.put( motherIndex,  prefix + "MotherIndex"  + suffix );
  iEvent.put( taudecaymode, prefix + "TauDecayMode" + suffix );
  iEvent.put( tauvisiblept, prefix + "TauVisiblePt" + suffix );
  iEvent.put( tauvisibleeta,prefix + "TauVisibleEta"+ suffix );
  iEvent.put( tauvisiblephi,prefix + "TauVisiblePhi"+ suffix );
  // gen status flags
  iEvent.put(isPromptFinalState       ,prefix + "IsPromptFinalState"+ suffix );
  iEvent.put(isPromptDecayed          ,prefix + "IsPromptDecayed"+ suffix );
  iEvent.put(isHardProcess            ,prefix + "IsHardProcess"+ suffix );
  iEvent.put(fromHardProcessFinalState,prefix + "FromHardProcessFinalState"+ suffix );
  iEvent.put(fromHardProcessDecayed   ,prefix + "FromHardProcessDecayed"+ suffix );
  iEvent.put(isLastCopy               ,prefix + "IsLastCopy"+ suffix );
  iEvent.put( topptweight  ,prefix + "TopPtWeight"  + suffix );
  // W/Z system Pt
  iEvent.put( worzsystempt ,prefix + "WorZSystemPt"  + suffix );
}



//------------------------------------------------------------------------------------------------------------------------------------------//
// See: https://hypernews.cern.ch/HyperNews/CMS/get/tauid/301/1/1/1.html                                                                    //
// getVisMomentum(), getInvisMomentum(), findDaughters(), isNeutrino(), countDecayProducts() and getGenTauDecayMode() are reproduced from:  //
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TauAnalysis/CandidateTools/src/candidateAuxFunctions.cc?revision=1.26&view=markup       //
//------------------------------------------------------------------------------------------------------------------------------------------//

reco::Candidate::LorentzVector RootTupleMakerV2_GenParticles::getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status){
  reco::Candidate::LorentzVector p4Vis(0,0,0,0);
  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if ( (status == -1 || (*daughter)->status() == status) && !isNeutrino(*daughter) ) {
      // ------ debug ------
      //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
      //  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
      // ------ ----- ------
      p4Vis += (*daughter)->p4();
    }
  }
  // ------ debug ------
  //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;
  // ------ ----- ------
  return p4Vis;
}

reco::Candidate::LorentzVector RootTupleMakerV2_GenParticles::getInvisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status){
  reco::Candidate::LorentzVector p4Invis(0,0,0,0);
  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if ( (status == -1 || (*daughter)->status() == status) && isNeutrino(*daughter) ) {
      p4Invis += (*daughter)->p4();
    }
  }
  return p4Invis;
}

void RootTupleMakerV2_GenParticles::findDaughters(const reco::GenParticle* mother,  std::vector<const reco::GenParticle*>& daughters, int status){
  unsigned numDaughters = mother->numberOfDaughters();
  for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
    if ( status == -1 || daughter->status() == status ) daughters.push_back(daughter);
    findDaughters(daughter, daughters, status);
  }
}

bool RootTupleMakerV2_GenParticles::isNeutrino(const reco::GenParticle* daughter){
  return ( TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 );
}

bool RootTupleMakerV2_GenParticles::isTop(const reco::GenParticle* particle){
  return ( TMath::Abs(particle->pdgId()) == 6 );
}

void RootTupleMakerV2_GenParticles::countDecayProducts(const reco::GenParticle* genParticle,
						       int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos, 
						       int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons)
{
  int absPdgId = TMath::Abs(genParticle->pdgId());
  int status   = genParticle->status();
  int charge   = genParticle->charge();
  
  if      ( absPdgId == 111 ) ++numPi0s;
  else if ( status   ==   1 ) {
    if      ( absPdgId == 11 ) ++numElectrons;
    else if ( absPdgId == 12 ) ++numElecNeutrinos;
    else if ( absPdgId == 13 ) ++numMuons;
    else if ( absPdgId == 14 ) ++numMuNeutrinos;
    else if ( absPdgId == 15 ) { 
      edm::LogError ("countDecayProducts")
        << "Found tau lepton with status code 1 !!";
      return; 
    }
    else if ( absPdgId == 16 ) return; // no need to count tau neutrinos
    else if ( absPdgId == 22 ) ++numPhotons;
    else if ( charge   !=  0 ) ++numChargedHadrons;
    else                       ++numOtherNeutralHadrons;
  } else {
    unsigned numDaughters = genParticle->numberOfDaughters();
    for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
      const reco::GenParticle* daughter = genParticle->daughterRef(iDaughter).get();
      
      countDecayProducts(daughter, 
			 numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
			 numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);
    }
  }
}

int RootTupleMakerV2_GenParticles::getGenTauDecayMode(const reco::GenParticle* genParticle) 
{
  //--- determine generator level tau decay mode
  //    NOTE: 
  //        (1) function implements logic defined in PhysicsTools/JetMCUtils/src/JetMCTag::genTauDecayMode
  //            for different type of argument 
  //        (2) this implementation should be more robust to handle cases of tau --> tau + gamma radiation

  int numElectrons           = 0;
  int numElecNeutrinos       = 0;
  int numMuons               = 0;
  int numMuNeutrinos         = 0; 
  int numChargedHadrons      = 0;
  int numPi0s                = 0; 
  int numOtherNeutralHadrons = 0;
  int numPhotons             = 0;
  
  countDecayProducts(genParticle,
		     numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
		     numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);
  
  if      ( numElectrons == 1 && numElecNeutrinos == 1 ) return 1;//std::string("electron");
  else if ( numMuons     == 1 && numMuNeutrinos   == 1 ) return 2;//std::string("muon");
  
  switch ( numChargedHadrons ) {
  case 1 : 
    if ( numOtherNeutralHadrons != 0 ) return 11;//std::string("oneProngOther");
    switch ( numPi0s ) {
    case 0:
      return 3;//std::string("oneProng0Pi0");
    case 1:
      return 4;//std::string("oneProng1Pi0");
    case 2:
      return 5;//std::string("oneProng2Pi0");
    default:
      return 6;//std::string("oneProngOther");
    }
  case 3 : 
    if ( numOtherNeutralHadrons != 0 ) return 12;//std::string("threeProngOther");
    switch ( numPi0s ) {
    case 0:
      return 7;//std::string("threeProng0Pi0");
    case 1:
      return 8;//std::string("threeProng1Pi0");
    default:
      return 9;//std::string("threeProngOther");
    }
  default:
    return 13;//std::string("rare");
  }
}

