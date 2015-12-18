#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_PFCandidates.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"

#include <set>

RootTupleMakerV2_PFCandidates::RootTupleMakerV2_PFCandidates(const edm::ParameterSet& iConfig) :
    jetInputTag(iConfig.getParameter<edm::InputTag>("JetInputTag")),
    electronInputTag(iConfig.getParameter<edm::InputTag>("ElectronInputTag")),
    muonInputTag(iConfig.getParameter<edm::InputTag>("MuonInputTag")),
    pfcandInputTag(iConfig.getParameter<edm::InputTag>("PFCandInputTag")),
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

  // get the Jet Collection
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(jetInputTag, jets);

  // get the Electron Collection
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByLabel(electronInputTag, electrons);
  
  // get the Muon Collection
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(muonInputTag, muons);
  
  // get the PFCandidate Collection
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByLabel(pfcandInputTag, pfcands);

  //-------------------------------------------

  if( electrons.isValid() && muons.isValid() && pfcands.isValid())
  { 
    // NB: SIC
    // originally the PFCandidates came from the 'pfCandsNotInJet' collection
    // we don't have such a collection in MiniAOD
    // so let's remove the jet constituents by hand
    // loop over the leptons and save their constituents
    //std::set<const reco::CandidatePtr> elecAndMuSourceCands;
    std::vector<reco::CandidatePtr> elecAndMuSourceCands;

    std::vector<const reco::Candidate*> electronCands;
    std::vector<const reco::Candidate*> muonCands;
    std::vector<const reco::Candidate*> jetCands;
    std::vector<const reco::Candidate*> pfCands;
    for(const pat::Electron &el : *electrons) electronCands.push_back(&el);
    for(const pat::Muon &mu : *muons) muonCands.push_back(&mu);
    for(const pat::Jet &jet : *jets) jetCands.push_back(&jet);
    for(const pat::PackedCandidate &pfc : *pfcands) pfCands.push_back(&pfc);


    //first electrons
    for(const reco::Candidate *elec : electronCands)
    {
      for(unsigned int i = 0, n = elec->numberOfSourceCandidatePtrs(); i < n; ++i)
        elecAndMuSourceCands.push_back(elec->sourceCandidatePtr(i));

      // now look for elecAndMuSourceCands inside jets
      for(const reco::Candidate *jet : jetCands)
      {
        std::vector<reco::CandidatePtr> sourceJetCands;
        for(unsigned int i = 0, n = jet->numberOfSourceCandidatePtrs(); i < n; ++i)
        {
          // if jet cand particle is close to electron, see if it is part of the electron
          if(deltaR(*(jet->sourceCandidatePtr(i)),*elec) < DRmatch)
          {
            // if jet constituent is in the electron, remove it from electron constituents set
            std::vector<reco::CandidatePtr>::iterator srcCand = std::find(elecAndMuSourceCands.begin(), elecAndMuSourceCands.end(), jet->sourceCandidatePtr(i));
            if(srcCand != elecAndMuSourceCands.end())
              elecAndMuSourceCands.erase(srcCand);
          }
        }
      }
    }

    //then muons
    for(const reco::Candidate *mu : muonCands)
    {
      for(unsigned int i = 0, n = mu->numberOfSourceCandidatePtrs(); i < n; ++i)
        elecAndMuSourceCands.push_back(mu->sourceCandidatePtr(i));

      // now look for elecAndMuSourceCands inside jets
      for(const reco::Candidate *jet : jetCands)
      {
        std::vector<reco::CandidatePtr> sourceJetCands;
        for(unsigned int i = 0, n = jet->numberOfSourceCandidatePtrs(); i < n; ++i)
        {
          // if jet cand particle is close to electron, see if it is part of the electron
          if(deltaR(*(jet->sourceCandidatePtr(i)),*mu) < DRmatch)
          {
            // if jet constituent is in the electron, remove it from electron constituents set
            std::vector<reco::CandidatePtr>::iterator srcCand = std::find(elecAndMuSourceCands.begin(), elecAndMuSourceCands.end(), jet->sourceCandidatePtr(i));
            if(srcCand != elecAndMuSourceCands.end())
              elecAndMuSourceCands.erase(srcCand);
          }
        }
      }
    }


    for(const reco::Candidate *pfcand : pfCands ){
      bool isCloseToRecoEle  = false;
      bool isCloseToRecoMuon = false;
      //look if it close to an electron in a cone
      for(const pat::Electron &elec : *electrons )
	{
	  if ( fabs( ROOT::Math::VectorUtil::DeltaR( elec.p4() , pfcand->p4() ) ) <= DRmatch ) 
	    {
	      isCloseToRecoEle = true;
	    }
	}//electron loop
      // look if it close to a muon in a cone
      for(const pat::Muon &mu : *muons )
	{
	  if ( fabs( ROOT::Math::VectorUtil::DeltaR( mu.p4() , pfcand->p4() ) ) <= DRmatch ) 
	    {
	      isCloseToRecoMuon = true;
	    }
	}//muon loop
      if( isCloseToRecoEle || isCloseToRecoMuon )
	{
	  eta->push_back( pfcand->eta() );
	  phi->push_back( pfcand->phi() );
	  pt->push_back( pfcand->pt() );
	  energy->push_back( pfcand->energy() );
	  charge->push_back( pfcand->charge() );	   
	}
    }


    //for(const pat::PackedCandidate &pfcand : *packedpfcands )
    //{
    //  bool isCloseToRecoEle  = false;
    //  bool isCloseToRecoMuon = false;

    //  // look if it close to an electron in a cone
    //  for(const pat::Electron &elec : *electrons )
    //  {
    //    if ( fabs( ROOT::Math::VectorUtil::DeltaR( elec.p4() , pfcand.p4() ) ) <= DRmatch ) 
    //    {
    //      isCloseToRecoEle = true;
    //    }
    //  }//

    //  // look if it close to a muon in a cone
    //  for(const pat::Muon &mu : *muons )
    //  {
    //    if ( fabs( ROOT::Math::VectorUtil::DeltaR( mu.p4() , pfcand.p4() ) ) <= DRmatch ) 
    //    {
    //      isCloseToRecoMuon = true;
    //    }
    //  }//

    //  if( isCloseToRecoEle || isCloseToRecoMuon )
    //  {
    //    eta->push_back( pfcand.eta() );
    //    phi->push_back( pfcand.phi() );
    //    pt->push_back( pfcand.pt() );
    //    energy->push_back( pfcand.energy() );
    //    charge->push_back( pfcand.charge() );	   
    //  }

    //}//loop over packed pf candidates

  }
  else {
    if( !electrons.isValid() )
      edm::LogError("RootTupleMakerV2_PFCandidatesError") << "Error! Can't get the product " << electronInputTag;
    if( !muons.isValid() )
      edm::LogError("RootTupleMakerV2_PFCandidatesError") << "Error! Can't get the product " << muonInputTag;    
    if( !pfcands.isValid() )
      edm::LogError("RootTupleMakerV2_PFCandidatesError") << "Error! Can't get the product " << pfcandInputTag;    
  }
  
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "EtaLeptLink" + suffix );
  iEvent.put( phi, prefix + "PhiLeptLink" + suffix );
  iEvent.put( pt, prefix + "PtLeptLink" + suffix );
  iEvent.put( energy, prefix + "EnergyLeptLink" + suffix );
  iEvent.put( charge, prefix + "ChargeLeptLink" + suffix );
}
