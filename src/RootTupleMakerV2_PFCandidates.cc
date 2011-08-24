#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFCandidates.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"

RootTupleMakerV2_PFCandidates::RootTupleMakerV2_PFCandidates(const edm::ParameterSet& iConfig) :
    reducedPFCandidateInputTag(iConfig.getParameter<edm::InputTag>("ReducedPFCandidateInputTag")),
    electronInputTag(iConfig.getParameter<edm::InputTag>("ElectronInputTag")),
    muonInputTag(iConfig.getParameter<edm::InputTag>("MuonInputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    DRmatch (iConfig.getParameter<double>       ("DRmatch"))
{
  produces <std::vector<double> > ( prefix + "EtaLeptLink" + suffix );
  produces <std::vector<double> > ( prefix + "PhiLeptLink" + suffix );
  produces <std::vector<double> > ( prefix + "PtLeptLink" + suffix );
  produces <std::vector<double> > ( prefix + "EnergyLeptLink" + suffix );
  produces <std::vector<int> >    ( prefix + "ChargeLeptLink" + suffix );
}

void RootTupleMakerV2_PFCandidates::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()  );

  //-----------------------------------------------------------------

  // get the PFCandidate Collection
  edm::Handle<edm::View<reco::Candidate> >   reducedpfcands;
  iEvent.getByLabel( reducedPFCandidateInputTag, reducedpfcands );

  // get the GsfElectron Collection
  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(electronInputTag, electrons);
  
  // get the Muon Collection
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(muonInputTag, muons);

  //-------------------------------------------

  if( electrons.isValid() && muons.isValid() )
  { 
             
    for( edm::View<reco::Candidate>::const_iterator icand = reducedpfcands->begin(); icand != reducedpfcands->end(); ++icand )
      {
	
	bool isCloseToRecoEle  = false;
	bool isCloseToRecoMuon = false;
	       
	// look if it close to an electron in a cone
	for( std::vector<pat::Electron>::const_iterator itele = electrons->begin(); itele != electrons->end(); ++itele )
	  {
	    if ( fabs( ROOT::Math::VectorUtil::DeltaR( itele->p4() ,icand->p4() ) ) <= DRmatch ) 
	    {
	      isCloseToRecoEle = true;
	    }
	}//

	// look if it close to a muon in a cone
	for( std::vector<pat::Muon>::const_iterator itmu = muons->begin(); itmu != muons->end(); ++itmu )
	  {
	    if ( fabs( ROOT::Math::VectorUtil::DeltaR( itmu->p4() ,icand->p4() ) ) <= DRmatch ) 
	    {
	      isCloseToRecoMuon = true;
	    }
	}//
	
	if( isCloseToRecoEle || isCloseToRecoMuon )
	  {
	    eta->push_back( icand->eta() );
	    phi->push_back( icand->phi() );
	    pt->push_back( icand->pt() );
	    energy->push_back( icand->energy() );
	    charge->push_back( icand->charge() );	   
	  }

      }//loop over reduced pf candidates
    
  }
  else {
    if( !electrons.isValid() )
      edm::LogError("RootTupleMakerV2_PFCandidatesError") << "Error! Can't get the product " << electronInputTag;
    if( !muons.isValid() )
      edm::LogError("RootTupleMakerV2_PFCandidatesError") << "Error! Can't get the product " << muonInputTag;    
  }
  
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "EtaLeptLink" + suffix );
  iEvent.put( phi, prefix + "PhiLeptLink" + suffix );
  iEvent.put( pt, prefix + "PtLeptLink" + suffix );
  iEvent.put( energy, prefix + "EnergyLeptLink" + suffix );
  iEvent.put( charge, prefix + "ChargeLeptLink" + suffix );
}
