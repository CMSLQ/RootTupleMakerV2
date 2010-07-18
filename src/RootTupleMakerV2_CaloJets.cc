#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_CaloJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


RootTupleMakerV2_CaloJets::RootTupleMakerV2_CaloJets(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    electronPt (iConfig.getParameter<double>    ("ElectronPt")),
    electronIso (iConfig.getParameter<double>   ("ElectronIso")),
    muonPt (iConfig.getParameter<double>        ("MuonPt")),
    muonIso (iConfig.getParameter<double>       ("MuonIso")),
    applyResJEC (iConfig.getParameter<bool>     ("ApplyResidualJEC")),
    resJEC (iConfig.getParameter<std::string>   ("ResidualJEC"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "PtRaw" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "EnergyRaw" + suffix );
  produces <std::vector<int> >    ( prefix + "Overlaps" + suffix );
  produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
  produces <std::vector<double> > ( prefix + "EMF" + suffix );
  produces <std::vector<double> > ( prefix + "resEMF" + suffix );
  produces <std::vector<double> > ( prefix + "HADF" + suffix );
  produces <std::vector<int> >    ( prefix + "n90Hits" + suffix );
  produces <std::vector<double> > ( prefix + "fHPD" + suffix );
  produces <std::vector<double> > ( prefix + "fRBX" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaEta" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaPhi" + suffix );
  produces <std::vector<double> > ( prefix + "TrackCountingHighEffBTag" + suffix );
  produces <std::vector<double> > ( prefix + "SimpleSecondaryVertexBTag" + suffix );
  produces <std::vector<double> > ( prefix + "SoftMuonByPtBTag" + suffix );
  produces <std::vector<double> > ( prefix + "BProbabilityBTag" + suffix );
}

void RootTupleMakerV2_CaloJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt_raw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy_raw ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     overlaps ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  emf  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  resEmf  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadf  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     n90Hits  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  fHPD  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  fRBX  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaPhi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackCountingHighEffBTag  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  simpleSecondaryVertexBTag  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  softMuonByPtBTag  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  bProbabilityBTag  ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(inputTag, jets);

  if(jets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_CaloJetsInfo") << "Total # CaloJets: " << jets->size();

    JetCorrectorParameters *ResJetCorPar = 0;
    FactorizedJetCorrector *JEC = 0;
    if(applyResJEC) {
      edm::FileInPath fipRes(resJEC);
      ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
      std::vector<JetCorrectorParameters> vParam;
      vParam.push_back(*ResJetCorPar);
      JEC = new FactorizedJetCorrector(vParam);
    }

    for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it ) {
      // exit from loop when you reach the required number of jets
      if(eta->size() >= maxSize)
        break;

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

      double corr = 1.;
      if( applyResJEC && iEvent.isRealData() ) {
        JEC->setJetEta( it->eta() );
        JEC->setJetPt( it->pt() ); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      // fill in all the vectors
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt()*corr );
      pt_raw->push_back( it->correctedJet("raw").pt() );
      energy->push_back( it->energy()*corr );
      energy_raw->push_back( it->correctedJet("raw").energy() );
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
      simpleSecondaryVertexBTag->push_back( it->bDiscriminator("simpleSecondaryVertexBJetTag") );
      softMuonByPtBTag->push_back( it->bDiscriminator("softMuonByPtBJetTags") );
      bProbabilityBTag->push_back( it->bDiscriminator("jetBProbabilityBJetTags") );
    }
  } else {
    edm::LogError("RootTupleMakerV2_CaloJetsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
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
  iEvent.put( simpleSecondaryVertexBTag, prefix + "SimpleSecondaryVertexBTag" + suffix );
  iEvent.put( softMuonByPtBTag, prefix + "SoftMuonByPtBTag" + suffix );
  iEvent.put( bProbabilityBTag, prefix + "BProbabilityBTag" + suffix );
}
