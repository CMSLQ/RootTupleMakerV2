#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_CaloJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"


RootTupleMakerV2_CaloJets::RootTupleMakerV2_CaloJets(const edm::ParameterSet& iConfig) :
  jetInputToken_ (consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("InputTag"))),
  // inputTagL1Offset(iConfig.getParameter<edm::InputTag>("InputTagL1Offset")),
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
  electronPt (iConfig.getParameter<double>    ("ElectronPt")),
  electronIso (iConfig.getParameter<double>   ("ElectronIso")),
  muonPt (iConfig.getParameter<double>        ("MuonPt")),
  muonIso (iConfig.getParameter<double>       ("MuonIso")),
  jecUncPath(iConfig.getParameter<std::string>("JECUncertainty")),
  readJECuncertainty (iConfig.getParameter<bool>   ("ReadJECuncertainty"))
  //OLD
  //     applyResJEC (iConfig.getParameter<bool>     ("ApplyResidualJEC")),
  //     resJEC (iConfig.getParameter<std::string>   ("ResidualJEC"))
{
  produces <std::vector<float> > ( prefix + "Eta" + suffix );
  produces <std::vector<float> > ( prefix + "Phi" + suffix );
  produces <std::vector<float> > ( prefix + "Pt" + suffix );
  produces <std::vector<float> > ( prefix + "PtRaw" + suffix );
  produces <std::vector<float> > ( prefix + "Energy" + suffix );
  produces <std::vector<float> > ( prefix + "EnergyRaw" + suffix );
  produces <std::vector<float> > ( prefix + "JECUnc" + suffix );
  produces <std::vector<float> > ( prefix + "L2L3ResJEC" + suffix );
  produces <std::vector<float> > ( prefix + "L3AbsJEC" + suffix );
  produces <std::vector<float> > ( prefix + "L2RelJEC" + suffix );
  produces <std::vector<float> > ( prefix + "L1FastJetJEC" + suffix );
  // produces <std::vector<float> > ( prefix + "L1OffsetJEC" + suffix );
  produces <std::vector<int> >    ( prefix + "Overlaps" + suffix );
  produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
  produces <std::vector<float> > ( prefix + "EMF" + suffix );
  produces <std::vector<float> > ( prefix + "resEMF" + suffix );
  produces <std::vector<float> > ( prefix + "HADF" + suffix );
  produces <std::vector<int> >    ( prefix + "n90Hits" + suffix );
  produces <std::vector<float> > ( prefix + "fHPD" + suffix );
  produces <std::vector<float> > ( prefix + "fRBX" + suffix );
  produces <std::vector<float> > ( prefix + "SigmaEta" + suffix );
  produces <std::vector<float> > ( prefix + "SigmaPhi" + suffix );
  produces <std::vector<float> > ( prefix + "TrackCountingHighEffBTag" + suffix );
  produces <std::vector<float> > ( prefix + "TrackCountingHighPurBTag" + suffix );
  produces <std::vector<float> > ( prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
  produces <std::vector<float> > ( prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
  produces <std::vector<float> > ( prefix + "JetProbabilityBTag" + suffix );
  produces <std::vector<float> > ( prefix + "JetBProbabilityBTag" + suffix );
  produces <std::vector<int> >    ( prefix + "PassLooseID" + suffix);
  produces <std::vector<int> >    ( prefix + "PassTightID" + suffix);
}

JetIDSelectionFunctor jetIDLoose( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE );
JetIDSelectionFunctor jetIDTight( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT );

pat::strbitset ret = jetIDLoose.getBitTemplate();

void RootTupleMakerV2_CaloJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<float> >  eta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  phi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt_raw  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  energy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  energy_raw ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  jecUnc_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  l2l3resJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  l3absJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  l2relJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  l1fastjetJEC_vec ( new std::vector<float>()  );
  // std::auto_ptr<std::vector<float> >  l1offsetJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     overlaps ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<float> >  emf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  resEmf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hadf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     n90Hits  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<float> >  fHPD  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  fRBX  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  sigmaEta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  sigmaPhi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackCountingHighEffBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackCountingHighPurBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  simpleSecondaryVertexHighEffBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  simpleSecondaryVertexHighPurBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  jetProbabilityBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  jetBProbabilityBTag  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >  passLooseID  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >  passTightID  ( new std::vector<int>()  );

  //-----------------------------------------------------------------
  //OLD
  //   edm::FileInPath fipUnc(jecUncPath);;
  //   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fipUnc.fullPath());
  //  
  //   JetCorrectorParameters *ResJetCorPar = 0;
  //   FactorizedJetCorrector *JEC = 0;
  //   if(applyResJEC) {
  //     edm::FileInPath fipRes(resJEC);
  //     ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
  //     std::vector<JetCorrectorParameters> vParam;
  //     vParam.push_back(*ResJetCorPar);
  //     JEC = new FactorizedJetCorrector(vParam);
  //   }

  //JEC Uncertainties 
  JetCorrectionUncertainty *jecUnc = 0;
  if(readJECuncertainty)
    {
      //(See https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1075/1.html 
      // and https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2367/1.html)
      // handle the jet corrector parameters collection
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      // get the jet corrector parameters collection from the global tag
      iSetup.get<JetCorrectionsRecord>().get(jecUncPath,JetCorParColl);
      // get the uncertainty parameters from the collection
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      // instantiate the jec uncertainty object
      jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }  

  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByToken(jetInputToken_, jets);

  // edm::Handle<std::vector<pat::Jet> > jetsL1Offset;
  // iEvent.getByLabel(inputTagL1Offset, jetsL1Offset);

  if(jets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_CaloJetsInfo") << "Total # CaloJets: " << jets->size();
    
    for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it ) {
      // exit from loop when you reach the required number of jets
      if(eta->size() >= maxSize)
        break;
      
      ret.set(false);
      int passjetLoose =0;
      if(jetIDLoose( *it, ret )) passjetLoose =1;

      ret.set(false);
      int passjetTight = 0;
      if (jetIDTight( *it, ret)) passjetTight =1;

      if(readJECuncertainty)
	{
	  jecUnc->setJetEta( it->eta() );
	  jecUnc->setJetPt( it->pt() ); // the uncertainty is a function of the corrected pt      
	}

      int ovrlps = 0;
      /* overlaps with good electrons (with different electron IDs) and muons are handled bitwise
         bit 0: eidRobustLoose
         bit 1: eidRobustTight
         bit 2: eidLoose
         bit 3: eidTight
         bit 4: eidRobustHighEnergy
         bit 5: HEEPId
         bit 6: GlobalMuonPromptTight
      */
      const reco::CandidatePtrVector & electrons = it->overlaps("electrons");
      for (size_t i = 0; i < electrons.size(); ++i) {
        // try to cast into pat::Electron
        const pat::Electron *electron = dynamic_cast<const pat::Electron *>(&*electrons[i]);
        if(electron) {
          if( electron->electronID("eidRobustLoose")>0.
              && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso
              && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<0;
          if( electron->electronID("eidRobustTight")>0.
              && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso
              && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<1;
          if( electron->electronID("eidLoose")>0.
              && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso
              && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<2;
          if( electron->electronID("eidTight")>0.
              && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso
              && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<3;
          if( electron->electronID("eidRobustHighEnergy")>0.
              && ((electron->trackIso()+electron->ecalIso()+electron->hcalIso())/electron->pt())<electronIso
              && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<4;
          if( electron->userInt("HEEPId")==0
               && electron->pt()>electronPt ) ovrlps = ovrlps | 1<<5;
        }
      }
      const reco::CandidatePtrVector & muons = it->overlaps("muons");
      for (size_t i = 0; i < muons.size(); ++i) {
        // try to cast into pat::Muon
        const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*muons[i]);
        if(muon) {
          if( muon->muonID("GlobalMuonPromptTight")
              && ((muon->trackIso()+muon->ecalIso()+muon->hcalIso())/muon->pt())<muonIso
              && muon->pt()>muonPt ) ovrlps = ovrlps | 1<<6;
        }
      }
      
      //OLD
      //       float corr = 1.;
      //       if( applyResJEC && iEvent.isRealData() ) {
      //         JEC->setJetEta( it->eta() );
      //         JEC->setJetPt( it->pt() ); // here you put the L2L3 Corrected jet pt
      //         corr = JEC->getCorrection();
      //       }
      //
      //       jecUnc->setJetEta( it->eta() );
      //       jecUnc->setJetPt( it->pt()*corr ); // the uncertainty is a function of the corrected pt

      // Status of JEC
      //std::cout << "CALO: currentJECLevel(): " << it->currentJECLevel() << std::endl;
      //std::cout << "CALO: currentJECSet(): " << it->currentJECSet() << std::endl;
      //-------------------

      // fill in all the vectors

      //OLD
      // pt->push_back( it->pt()*corr );
      // energy->push_back( it->energy()*corr );
      // resJEC_vec->push_back( corr );
      
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      pt_raw->push_back( it->correctedJet("Uncorrected").pt() );
      energy->push_back( it->energy() );
      energy_raw->push_back( it->correctedJet("Uncorrected").energy() );
      l2l3resJEC_vec->push_back( it->pt()/it->correctedJet("L3Absolute").pt() );
      l3absJEC_vec->push_back( it->correctedJet("L3Absolute").pt()/it->correctedJet("L2Relative").pt() );
      l2relJEC_vec->push_back( it->correctedJet("L2Relative").pt()/it->correctedJet("L1FastJet").pt() );
      l1fastjetJEC_vec->push_back( it->correctedJet("L1FastJet").pt()/it->correctedJet("Uncorrected").pt() );
      if(readJECuncertainty)
	jecUnc_vec->push_back( jecUnc->getUncertainty(true) );
      else
	jecUnc_vec->push_back( -999 );
      overlaps->push_back( ovrlps );
      partonFlavour->push_back( it->partonFlavour() );
      emf->push_back( it->emEnergyFraction() );
      resEmf->push_back( it->jetID().restrictedEMF );
      hadf->push_back( it->energyFractionHadronic() );
      n90Hits->push_back( it->jetID().n90Hits );
      fHPD->push_back( it->jetID().fHPD );
      fRBX->push_back( it->jetID().fRBX );
      sigmaEta->push_back( sqrt(it->etaetaMoment()) );
      sigmaPhi->push_back( sqrt(it->phiphiMoment()) );
      trackCountingHighEffBTag->push_back( it->bDiscriminator("trackCountingHighEffBJetTags") );
      trackCountingHighPurBTag->push_back( it->bDiscriminator("trackCountingHighPurBJetTags") );
      simpleSecondaryVertexHighEffBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighEffBJetTags") );
      simpleSecondaryVertexHighPurBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighPurBJetTags") );
      jetProbabilityBTag->push_back( it->bDiscriminator("jetProbabilityBJetTags") );
      jetBProbabilityBTag->push_back( it->bDiscriminator("jetBProbabilityBJetTags") );
      passLooseID->push_back( passjetLoose );
      passTightID->push_back( passjetTight );
    }
  } else {
    edm::LogError("RootTupleMakerV2_CaloJetsError") << "Error! Can't get the jetInputToken_";
  }

  //L1Offset JEC
  // if(jetsL1Offset.isValid())
  //   {
  //     
  //     for( std::vector<pat::Jet>::const_iterator it = jetsL1Offset->begin(); it != jetsL1Offset->end(); ++it )
  // 	{
  // 	  // exit from loop when you reach the required number of jets
  // 	  if(l1offsetJEC_vec->size() >= maxSize)
  // 	    break;
  // 	  
  // 	  l1offsetJEC_vec->push_back( it->correctedJet("L1Offset").pt()/it->correctedJet("Uncorrected").pt() );
  // 	}
  //   }
  // else
  //   {
  //     edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTagL1Offset;
  //   }
  
  //OLD
  delete jecUnc;
  //   delete ResJetCorPar;
  //   delete JEC;

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
  iEvent.put( jecUnc_vec, prefix + "JECUnc" + suffix );
  iEvent.put( l2l3resJEC_vec, prefix + "L2L3ResJEC" + suffix );
  iEvent.put( l3absJEC_vec, prefix + "L3AbsJEC" + suffix );
  iEvent.put( l2relJEC_vec, prefix + "L2RelJEC" + suffix );
  iEvent.put( l1fastjetJEC_vec, prefix + "L1FastJetJEC" + suffix );
  // iEvent.put( l1offsetJEC_vec, prefix + "L1OffsetJEC" + suffix );
  iEvent.put( overlaps, prefix + "Overlaps" + suffix );
  iEvent.put( partonFlavour, prefix + "PartonFlavour" + suffix );
  iEvent.put( emf, prefix + "EMF" + suffix );
  iEvent.put( resEmf, prefix + "resEMF" + suffix );
  iEvent.put( hadf, prefix + "HADF" + suffix );
  iEvent.put( n90Hits, prefix + "n90Hits" + suffix );
  iEvent.put( fHPD, prefix + "fHPD" + suffix );
  iEvent.put( fRBX, prefix + "fRBX" + suffix );
  iEvent.put( sigmaEta, prefix + "SigmaEta" + suffix );
  iEvent.put( sigmaPhi, prefix + "SigmaPhi" + suffix );
  iEvent.put( trackCountingHighEffBTag, prefix + "TrackCountingHighEffBTag" + suffix );
  iEvent.put( trackCountingHighPurBTag, prefix + "TrackCountingHighPurBTag" + suffix );
  iEvent.put( simpleSecondaryVertexHighEffBTag, prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
  iEvent.put( simpleSecondaryVertexHighPurBTag, prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
  iEvent.put( jetProbabilityBTag, prefix + "JetProbabilityBTag" + suffix );
  iEvent.put( jetBProbabilityBTag, prefix + "JetBProbabilityBTag" + suffix );
  iEvent.put( passLooseID, prefix + "PassLooseID" + suffix);
  iEvent.put( passTightID, prefix + "PassTightID" + suffix);
}
