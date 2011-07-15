#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Taus.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/Common/interface/Ref.h>
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTau.h"

RootTupleMakerV2_Taus::RootTupleMakerV2_Taus(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt"  + suffix );
  produces <std::vector<double> > ( prefix + "Et"  + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"  + suffix );
  produces <std::vector<int> >    ( prefix + "IsPFTau"  + suffix );
  produces <std::vector<int> >    ( prefix + "IsCaloTau"  + suffix );
  produces <std::vector<int> >    ( prefix + "DecayMode"  + suffix );
  produces <std::vector<double> > ( prefix + "EmFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "Hcal3x3OverPLead"  + suffix );
  produces <std::vector<double> > ( prefix + "HcalMaxOverPLead"  + suffix );
  produces <std::vector<double> > ( prefix + "HcalTotOverPLead"  + suffix );
  produces <std::vector<double> > ( prefix + "IsolationPFChargedHadrCandsPtSum"  + suffix );
  produces <std::vector<double> > ( prefix + "IsolationPFGammaCandsEtSum"  + suffix );
  produces <std::vector<double> > ( prefix + "LeadPFChargedHadrCandsignedSipt"  + suffix );
  produces <std::vector<double> > ( prefix + "EtaLeadCharged"  + suffix );
  produces <std::vector<double> > ( prefix + "PhiLeadCharged"  + suffix );
  produces <std::vector<double> > ( prefix + "PtLeadCharged"  + suffix );
  produces <std::vector<int> >    ( prefix + "AgainstElectronDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "AgainstMuonDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "ByIsolationUsingLeadingPionDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "LeadingPionPtCutDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "LeadingTrackFindingDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "LeadingTrackPtCutDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "TrackIsolationDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "TrackIsolationUsingLeadingPionDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "EcalIsolationDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
  produces <std::vector<int> >    ( prefix + "ByIsolationDiscr"  + suffix );
}

void RootTupleMakerV2_Taus::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  et  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     ispftau  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     iscalotau  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     decaymode  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<double> >  emfraction  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  hcal3x3overplead  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  hcalmaxoverplead  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  hcaltotoverplead  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  isolationpfchargedhadrcandsptsum  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  isolationpfgammacandsetsum  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  leadpfchargedhadrcandsignedsipt ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  etaleadcharged  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  phileadcharged  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  ptleadcharged  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<int> >     againstelectrondiscr  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     againstmuondiscr  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     byisolationusingleadingpiondiscr  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     leadingpionptcutdiscr  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     leadingtrackfindingdiscr  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     leadingtrackptcutdiscr ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     trackisolationdiscr ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     trackisolationusingleadingpiondiscr ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     ecalisolationdiscr ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     ecalisolationusingleadingpiondiscr ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     byisolationdiscr ( new std::vector<int>()   );


  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(inputTag, taus);

  if(taus.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Taus: " << taus->size();

    std::vector<pat::Tau>::const_iterator it     = taus -> begin();
    std::vector<pat::Tau>::const_iterator it_end = taus -> end();
    //
    //
    for (; it != it_end; ++it ) { 
      if ( eta->size() > maxSize ) break;
      //
      // Discriminators are defined in:  
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.27
      if(it->tauID("againstElectron")>0.5){againstelectrondiscr -> push_back ( 1 ) ; }
      if(it->tauID("againstElectron")<0.5){againstelectrondiscr -> push_back ( 0 ) ; }
      if(it->tauID("againstMuon")>0.5){againstmuondiscr -> push_back ( 1 ) ; }
      if(it->tauID("againstMuon")<0.5){againstmuondiscr -> push_back ( 0 ) ; }
      if(it->tauID("byIsolationUsingLeadingPion")>0.5){byisolationusingleadingpiondiscr -> push_back ( 1 ) ; }
      if(it->tauID("byIsolationUsingLeadingPion")<0.5){byisolationusingleadingpiondiscr -> push_back ( 0 ) ; }
      if(it->tauID("leadingPionPtCut")>0.5){leadingpionptcutdiscr -> push_back ( 1 ) ; }
      if(it->tauID("leadingPionPtCut")<0.5){leadingpionptcutdiscr -> push_back ( 0 ) ; }
      if(it->tauID("leadingTrackFinding")>0.5){leadingtrackfindingdiscr -> push_back ( 1 ) ; }
      if(it->tauID("leadingTrackFinding")<0.5){leadingtrackfindingdiscr -> push_back ( 0 ) ; }
      if(it->tauID("leadingTrackPtCut")>0.5){leadingtrackptcutdiscr -> push_back ( 1 ) ; }
      if(it->tauID("leadingTrackPtCut")<0.5){leadingtrackptcutdiscr -> push_back ( 0 ) ; }
      if(it->tauID("trackIsolation")>0.5){trackisolationdiscr -> push_back ( 1 ) ; }
      if(it->tauID("trackIsolation")<0.5){trackisolationdiscr -> push_back ( 0 ) ; }
      if(it->tauID("trackIsolationUsingLeadingPion")>0.5){trackisolationusingleadingpiondiscr -> push_back ( 1 ) ; }
      if(it->tauID("trackIsolationUsingLeadingPion")<0.5){trackisolationusingleadingpiondiscr -> push_back ( 0 ) ; }
      if(it->tauID("ecalIsolation")>0.5){ecalisolationdiscr -> push_back ( 1 ) ; }
      if(it->tauID("ecalIsolation")<0.5){ecalisolationdiscr -> push_back ( 0 ) ; }
      if(it->tauID("ecalIsolationUsingLeadingPion")>0.5){ecalisolationusingleadingpiondiscr -> push_back ( 1 ) ; }
      if(it->tauID("ecalIsolationUsingLeadingPion")<0.5){ecalisolationusingleadingpiondiscr -> push_back ( 0 ) ; }
      if(it->tauID("byIsolation")>0.5){byisolationdiscr -> push_back ( 1 ) ; }
      if(it->tauID("byIsolation")<0.5){byisolationdiscr -> push_back ( 0 ) ; }
      //
      eta -> push_back ( (double)(it -> eta()) ) ;
      phi -> push_back ( (double)(it -> phi()) ) ;
      pt  -> push_back ( (double)(it -> pt() ) ) ;
      et  -> push_back ( (double)(it -> et() ) ) ;
      charge  -> push_back ( (int)(it -> charge() ) ) ;
      if(  it ->isPFTau()  ){ispftau   -> push_back ( 1 ) ;}
      if( !it ->isPFTau()  ){ispftau   -> push_back ( 0 ) ;}
      if(  it ->isCaloTau()){iscalotau -> push_back ( 1 ) ;}
      if( !it ->isCaloTau()){iscalotau -> push_back ( 0 ) ;}
      decaymode                        -> push_back( (double)(it->decayMode())  );
      emfraction                       -> push_back( (double)(it->emFraction()) );
      hcal3x3overplead                 -> push_back( (double)(it->hcal3x3OverPLead()) );
      hcalmaxoverplead                 -> push_back( (double)(it->hcalMaxOverPLead()) );
      hcaltotoverplead                 -> push_back( (double)(it->hcalTotOverPLead()) );
      isolationpfchargedhadrcandsptsum -> push_back( (double)(it->isolationPFChargedHadrCandsPtSum()) );
      isolationpfgammacandsetsum       -> push_back( (double)(it->isolationPFGammaCandsEtSum())       );
      leadpfchargedhadrcandsignedsipt  -> push_back( (double)(it->leadPFChargedHadrCandsignedSipt())  );
      reco::PFCandidateRef leadPFChargedHadrCand_Ref = it->leadPFChargedHadrCand();
      etaleadcharged                   -> push_back( (double)(leadPFChargedHadrCand_Ref->eta()) );
      phileadcharged                   -> push_back( (double)(leadPFChargedHadrCand_Ref->phi()) );
      ptleadcharged                    -> push_back( (double)(leadPFChargedHadrCand_Ref->pt())  );
      //
    }
  } else {
    edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << inputTag;
  }
  
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt , prefix + "Pt"  + suffix );
  iEvent.put( et , prefix + "Et"  + suffix );
  iEvent.put( charge     , prefix + "Charge"      + suffix );
  iEvent.put( ispftau    , prefix + "IsPFTau"     + suffix );
  iEvent.put( iscalotau  , prefix + "IsCaloTau"   + suffix );
  iEvent.put( decaymode  , prefix + "DecayMode"   + suffix );
  iEvent.put( emfraction , prefix + "EmFraction"  + suffix );
  iEvent.put( hcal3x3overplead , prefix + "Hcal3x3OverPLead"  + suffix );
  iEvent.put( hcalmaxoverplead , prefix + "HcalMaxOverPLead"  + suffix );
  iEvent.put( hcaltotoverplead , prefix + "HcalTotOverPLead"  + suffix );
  iEvent.put( isolationpfchargedhadrcandsptsum , prefix + "IsolationPFChargedHadrCandsPtSum"  + suffix );
  iEvent.put( isolationpfgammacandsetsum       , prefix + "IsolationPFGammaCandsEtSum"        + suffix );
  iEvent.put( leadpfchargedhadrcandsignedsipt  , prefix + "LeadPFChargedHadrCandsignedSipt"   + suffix );
  iEvent.put( etaleadcharged, prefix + "EtaLeadCharged" + suffix );
  iEvent.put( phileadcharged, prefix + "PhiLeadCharged" + suffix );
  iEvent.put( ptleadcharged , prefix + "PtLeadCharged"  + suffix );
  iEvent.put( againstelectrondiscr,                prefix + "AgainstElectronDiscr"                + suffix );
  iEvent.put( againstmuondiscr,                    prefix + "AgainstMuonDiscr"                    + suffix );
  iEvent.put( byisolationusingleadingpiondiscr,    prefix + "ByIsolationUsingLeadingPionDiscr"    + suffix );
  iEvent.put( leadingpionptcutdiscr,               prefix + "LeadingPionPtCutDiscr"               + suffix );
  iEvent.put( leadingtrackfindingdiscr,            prefix + "LeadingTrackFindingDiscr"            + suffix );
  iEvent.put( leadingtrackptcutdiscr,              prefix + "LeadingTrackPtCutDiscr"              + suffix );
  iEvent.put( trackisolationdiscr,                 prefix + "TrackIsolationDiscr"                 + suffix );
  iEvent.put( trackisolationusingleadingpiondiscr, prefix + "TrackIsolationUsingLeadingPionDiscr" + suffix );
  iEvent.put( ecalisolationdiscr,                  prefix + "EcalIsolationDiscr"                  + suffix );
  iEvent.put( ecalisolationusingleadingpiondiscr,  prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
  iEvent.put( byisolationdiscr,                    prefix + "ByIsolationDiscr"                    + suffix );

}
