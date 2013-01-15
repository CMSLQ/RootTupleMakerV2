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
  inputTag    (iConfig.getParameter<edm::InputTag>("InputTag")),
  vtxInputTag (iConfig.getParameter<edm::InputTag>("VertexInputTag")),
  prefix      (iConfig.getParameter<std::string>("Prefix")),
  suffix      (iConfig.getParameter<std::string>("Suffix")),
  maxSize     (iConfig.getParameter<unsigned int>("MaxSize")),
  isSCTau     (iConfig.getParameter<bool>("isSCTau")),
  isHPSTau    (iConfig.getParameter<bool>("isHPSTau"))
{
  produces <std::vector<double> > ( prefix + "Eta"                              + suffix );
  produces <std::vector<double> > ( prefix + "Phi"                              + suffix );
  produces <std::vector<double> > ( prefix + "Pt"                               + suffix );
  produces <std::vector<double> > ( prefix + "Et"                               + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"                           + suffix );
  produces <std::vector<int> >    ( prefix + "IsPFTau"                          + suffix );
  produces <std::vector<int> >    ( prefix + "IsCaloTau"                        + suffix );
  produces <std::vector<int> >    ( prefix + "DecayMode"                        + suffix );
  produces <std::vector<double> > ( prefix + "EmFraction"                       + suffix );
  produces <std::vector<double> > ( prefix + "Hcal3x3OverPLead"                 + suffix );
  produces <std::vector<double> > ( prefix + "HcalMaxOverPLead"                 + suffix );
  produces <std::vector<double> > ( prefix + "HcalTotOverPLead"                 + suffix );
  produces <std::vector<double> > ( prefix + "IsolationPFChargedHadrCandsPtSum" + suffix );
  produces <std::vector<double> > ( prefix + "IsolationPFGammaCandsEtSum"       + suffix );
  produces <std::vector<double> > ( prefix + "LeadPFChargedHadrCandsignedSipt"  + suffix );
  produces <std::vector<double> > ( prefix + "EtaLeadCharged"                   + suffix );
  produces <std::vector<double> > ( prefix + "PhiLeadCharged"                   + suffix );
  produces <std::vector<double> > ( prefix + "PtLeadCharged"                    + suffix );
  produces <std::vector<double> > ( prefix + "PhiphiMoment"                     + suffix );
  produces <std::vector<double> > ( prefix + "EtaetaMoment"                     + suffix );
  produces <std::vector<double> > ( prefix + "EtaphiMoment"                     + suffix );
  produces <std::vector<double> > ( prefix + "EcalStripSumEOverPLead"           + suffix ); 
  produces <std::vector<double> > ( prefix + "BremsRecoveryEOverPLead"          + suffix );
  produces <std::vector<double> > ( prefix + "MaximumHCALPFClusterEt"           + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenParticlePt"             + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenParticleEta"            + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenParticlePhi"            + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenJetPt"                  + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenJetEta"                 + suffix ); 
  produces <std::vector<double> > ( prefix + "MatchedGenJetPhi"                 + suffix ); 
  //
  // ShrinkingCone PFTau Specific
  if(isSCTau){
    //shrinkingCone PFTau Discriminators (SCTau)
    produces <std::vector<double> >    ( prefix + "LeadingTrackFindingDiscr"            + suffix );
    produces <std::vector<double> >    ( prefix + "LeadingTrackPtCutDiscr"              + suffix );
    produces <std::vector<double> >    ( prefix + "LeadingPionPtCutDiscr"               + suffix );
    produces <std::vector<double> >    ( prefix + "IsolationDiscr"                      + suffix );
    produces <std::vector<double> >    ( prefix + "TrackIsolationDiscr"                 + suffix );
    produces <std::vector<double> >    ( prefix + "EcalIsolationDiscr"                  + suffix );
    produces <std::vector<double> >    ( prefix + "IsolationUsingLeadingPionDiscr"      + suffix );
    produces <std::vector<double> >    ( prefix + "TrackIsolationUsingLeadingPionDiscr" + suffix );
    produces <std::vector<double> >    ( prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronDiscr"                + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonDiscr"                    + suffix );
    produces <std::vector<double> >    ( prefix + "TaNCDiscr"                           + suffix );
    produces <std::vector<double> >    ( prefix + "TaNCfrOnePercentDiscr"               + suffix );
    produces <std::vector<double> >    ( prefix + "TaNCfrHalfPercentDiscr"              + suffix );
    produces <std::vector<double> >    ( prefix + "TaNCfrQuarterPercentDiscr"           + suffix );
    produces <std::vector<double> >    ( prefix + "TaNCfrTenthPercentDiscr"             + suffix );
  }
  //
  // HPS PFTau Specific
  if(isHPSTau){
    //hps PFTau Discriminators (HPSTau)
    produces <std::vector<double> >    ( prefix + "DecayModeFindingDiscr"  + suffix );
    //
    produces <std::vector<double> >    ( prefix + "AgainstElectronLooseDiscr"        + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMediumDiscr"       + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronTightDiscr"        + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMVADiscr"          + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMVA2rawDiscr"      + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMVA2categoryDiscr" + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronVLooseMVA2Discr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronLooseMVA2Discr"    + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMediumMVA2Discr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronTightMVA2Discr"    + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMVA3rawDiscr"      + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMVA3categoryDiscr" + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronLooseMVA3Discr"    + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronMediumMVA3Discr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronTightMVA3Discr"    + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronVTightMVA3Discr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstElectronDeadECALDiscr"     + suffix );
    //
    produces <std::vector<double> >    ( prefix + "AgainstMuonLooseDiscr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonMediumDiscr"  + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonTightDiscr"   + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonLoose2Discr"  + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonMedium2Discr" + suffix );
    produces <std::vector<double> >    ( prefix + "AgainstMuonTight2Discr"  + suffix );
    //
    produces <std::vector<double> >    ( prefix + "VLooseIsolationDiscr"                           + suffix );
    produces <std::vector<double> >    ( prefix + "LooseIsolationDiscr"                            + suffix );
    produces <std::vector<double> >    ( prefix + "MediumIsolationDiscr"                           + suffix );
    produces <std::vector<double> >    ( prefix + "TightIsolationDiscr"                            + suffix );
    produces <std::vector<double> >    ( prefix + "VLooseIsolationDeltaBetaCorrDiscr"              + suffix );
    produces <std::vector<double> >    ( prefix + "LooseIsolationDeltaBetaCorrDiscr"               + suffix );
    produces <std::vector<double> >    ( prefix + "MediumIsolationDeltaBetaCorrDiscr"              + suffix );
    produces <std::vector<double> >    ( prefix + "TightIsolationDeltaBetaCorrDiscr"               + suffix );
    produces <std::vector<double> >    ( prefix + "VLooseCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    produces <std::vector<double> >    ( prefix + "LooseCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    produces <std::vector<double> >    ( prefix + "MediumCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    produces <std::vector<double> >    ( prefix + "TightCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    produces <std::vector<double> >    ( prefix + "CombinedIsolationDeltaBetaCorr3HitsDiscr"       + suffix );
    produces <std::vector<double> >    ( prefix + "LooseCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    produces <std::vector<double> >    ( prefix + "MediumCombinedIsolationDeltaBetaCorr3HitsDiscr" + suffix );
    produces <std::vector<double> >    ( prefix + "TightCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    produces <std::vector<double> >    ( prefix + "IsolationMVArawDiscr"                           + suffix );
    produces <std::vector<double> >    ( prefix + "LooseIsolationMVADiscr"                         + suffix );
    produces <std::vector<double> >    ( prefix + "MediumIsolationMVADiscr"                        + suffix );
    produces <std::vector<double> >    ( prefix + "TightIsolationMVADiscr"                         + suffix );
    produces <std::vector<double> >    ( prefix + "LooseIsolationMVA2Discr"                        + suffix );
    produces <std::vector<double> >    ( prefix + "MediumIsolationMVA2Discr"                       + suffix );
    produces <std::vector<double> >    ( prefix + "TightIsolationMVA2Discr"                        + suffix );
    //
    // HPSTau Signal PFCandidates Info
    produces <std::vector<double> >    ( prefix + "SignalPFChargedHadrCandsPt"    + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFChargedHadrCandsEta"   + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFChargedHadrCandsPhi"   + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFChargedHadrCandsCount" + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFNeutrHadrCandsPt"      + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFNeutrHadrCandsEta"     + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFNeutrHadrCandsPhi"     + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFNeutrHadrCandsCount"   + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFGammaCandsPt"          + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFGammaCandsEta"         + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFGammaCandsPhi"         + suffix );
    produces <std::vector<double> >    ( prefix + "SignalPFGammaCandsCount"       + suffix );
    //
    // HPSTau Vertex Info
    produces <std::vector<int> >    ( prefix + "VtxIndex"      + suffix );
    produces <std::vector<double> > ( prefix + "VtxDistXY"     + suffix );
    produces <std::vector<double> > ( prefix + "VtxDistZ"      + suffix );
    produces <std::vector<double> > ( prefix + "LeadVtxDistXY" + suffix );
    produces <std::vector<double> > ( prefix + "LeadVtxDistZ"  + suffix );
    //
    // --------------------------------------------------------------------------------------- //
    // HPS Tau Optional Isolation information
    //produces <std::vector<double> >    ( prefix + "IsolationPFChargedHadrCandsPt"    + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFChargedHadrCandsEta"   + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFChargedHadrCandsPhi"   + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFChargedHadrCandsCount" + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFNeutrHadrCandsPt"      + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFNeutrHadrCandsEta"     + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFNeutrHadrCandsPhi"     + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFNeutrHadrCandsCount"   + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFGammaCandsPt"          + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFGammaCandsEta"         + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFGammaCandsPhi"         + suffix );
    //produces <std::vector<double> >    ( prefix + "IsolationPFGammaCandsCount"       + suffix );
    // --------------------------------------------------------------------------------------- //
    //
  }
  //
}

void RootTupleMakerV2_Taus::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //
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
  std::auto_ptr<std::vector<double> >  phiphimoment  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  etaetamoment  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  etaphimoment  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  ecalstripsumeoverplead  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  bremsrecoveryeoverplead  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  maximumhcalpfclusteret  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenparticlept ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenparticleeta ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenparticlephi ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenjetpt ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenjeteta ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedgenjetphi ( new std::vector<double>()   );
  //
  //shrinkingCone PFTau Discriminators (SCTau)
  std::auto_ptr<std::vector<double> >  leadingtrackfindingdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  leadingtrackptcutdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  leadingpionptcutdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  isolationdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  trackisolationdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  ecalisolationdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  isolationusingleadingpiondiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  trackisolationusingleadingpiondiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  ecalisolationusingleadingpiondiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectrondiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuondiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tancdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tancfronepercentdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tancfrhalfpercentdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tancfrquarterpercentdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tancfrtenthpercentdiscr  ( new std::vector<double>()   );
  //
  //hps PFTau Discriminators (HPSTau)
  std::auto_ptr<std::vector<double> >  decaymodefindingdiscr  ( new std::vector<double>()   );
  //
  std::auto_ptr<std::vector<double> >  againstelectronloosediscr        ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmediumdiscr       ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectrontightdiscr        ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmvadiscr          ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmva2rawdiscr      ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmva2categorydiscr ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronvloosemva2discr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronloosemva2discr    ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmediummva2discr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectrontightmva2discr    ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmva3rawdiscr      ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmva3categorydiscr ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronloosemva3discr    ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronmediummva3discr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectrontightmva3discr    ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectronvtightmva3discr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstelectrondeadecaldiscr     ( new std::vector<double>()   );
  //
  std::auto_ptr<std::vector<double> >  againstmuonloosediscr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuonmediumdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuontightdiscr   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuonloose2discr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuonmedium2discr ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  againstmuontight2discr  ( new std::vector<double>()   );
  //
  std::auto_ptr<std::vector<double> >  vlooseisolationdiscr                           ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  looseisolationdiscr                            ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumisolationdiscr                           ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightisolationdiscr                            ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  vlooseisolationdeltabetacorrdiscr              ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  looseisolationdeltabetacorrdiscr               ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumisolationdeltabetacorrdiscr              ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightisolationdeltabetacorrdiscr               ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  vloosecombinedisolationdeltabetacorrdiscr      ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  loosecombinedisolationdeltabetacorrdiscr       ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumcombinedisolationdeltabetacorrdiscr      ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightcombinedisolationdeltabetacorrdiscr       ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  combinedisolationdeltabetacorr3hitsdiscr       ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  loosecombinedisolationdeltabetacorr3hitsdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumcombinedisolationdeltabetacorr3hitsdiscr ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightcombinedisolationdeltabetacorr3hitsdiscr  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  isolationmvarawdiscr                           ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  looseisolationmvadiscr                         ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumisolationmvadiscr                        ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightisolationmvadiscr                         ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  looseisolationmva2discr                        ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  mediumisolationmva2discr                       ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  tightisolationmva2discr                        ( new std::vector<double>()   );
  //
  //Signal Particles (HPSTau)  
  std::auto_ptr<std::vector<double> >  signalpfchargedhadrcandspt    ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfchargedhadrcandseta   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfchargedhadrcandsphi   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfchargedhadrcandscount ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfneutrhadrcandspt      ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfneutrhadrcandseta     ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfneutrhadrcandsphi     ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfneutrhadrcandscount   ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfgammacandspt          ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfgammacandseta         ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfgammacandsphi         ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  signalpfgammacandscount       ( new std::vector<double>()   );
  //                                                                                                                                                                          
  // HPSTau Vertex Info 
  std::auto_ptr<std::vector<int> >     vtxIndex   ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<double> >  vtxDistXY  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtxDistZ   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtx0DistXY ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtx0DistZ  ( new std::vector<double>()  );
  //
  // --------------------------------------------------------------------------------------- //
  // HPS Tau Optional Isolation information
  //std::auto_ptr<std::vector<double> >  isolationpfchargedhadrcandspt    ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfchargedhadrcandseta   ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfchargedhadrcandsphi   ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfchargedhadrcandscount ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfneutrhadrcandspt      ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfneutrhadrcandseta     ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfneutrhadrcandsphi     ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfneutrhadrcandscount   ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfgammacandspt          ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfgammacandseta         ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfgammacandsphi         ( new std::vector<double>()   );
  //std::auto_ptr<std::vector<double> >  isolationpfgammacandscount       ( new std::vector<double>()   );
  // --------------------------------------------------------------------------------------- //
  //  
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);
  //
  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(inputTag, taus);
  //
  if(taus.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Taus: " << taus->size();

    std::vector<pat::Tau>::const_iterator it     = taus -> begin();
    std::vector<pat::Tau>::const_iterator it_end = taus -> end();
    //
    for (; it != it_end; ++it ) { 
      if ( eta->size() > maxSize ) break;
      //
      // SCTau Discriminators are at: (cvs up -r 1.53 PhysicsTools/PatAlgos/python/tools/tauTools.py)
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/tools/tauTools.py?revision=1.53&view=markup
      // Out of the Box, Clean SC Taus are given by switchToPFTau:
      // tauID("leadingTrackFinding") > 0.5 
      // tauID("leadingPionPtCut") > 0.5 
      // tauID("byIsolationUsingLeadingPion") > 0.5
      // tauID("againstMuon") > 0.5 
      // tauID("againstElectron") > 0.5'
      // (signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3)'
      //
      if(isSCTau){
	againstelectrondiscr                -> push_back( it->tauID("againstElectron" )               );//applied in switchToPFTau
	againstmuondiscr                    -> push_back( it->tauID("againstMuon")                    );//applied in switchToPFTau 
	ecalisolationdiscr                  -> push_back( it->tauID("ecalIsolation")                  );
	ecalisolationusingleadingpiondiscr  -> push_back( it->tauID("byIsolationUsingLeadingPion")    );
	isolationdiscr                      -> push_back( it->tauID("byIsolation")                    );
	isolationusingleadingpiondiscr      -> push_back( it->tauID("byIsolationUsingLeadingPion")    );//applied in switchToPFTau 
	leadingtrackfindingdiscr            -> push_back( it->tauID("leadingTrackFinding")            );//applied in switchToPFTau
	leadingtrackptcutdiscr              -> push_back( it->tauID("leadingTrackPtCut")              );
	leadingpionptcutdiscr               -> push_back( it->tauID("leadingPionPtCut")               );//applied in switchToPFTau 
	trackisolationdiscr                 -> push_back( it->tauID("trackIsolation")                 );
	trackisolationusingleadingpiondiscr -> push_back( it->tauID("trackIsolationUsingLeadingPion") );
	//
	tancdiscr                 -> push_back( it->tauID("byTaNC")                 );
	tancfronepercentdiscr     -> push_back( it->tauID("byTaNCfrOnePercent")     );
	tancfrhalfpercentdiscr    -> push_back( it->tauID("byTaNCfrHalfPercent")    );
	tancfrquarterpercentdiscr -> push_back( it->tauID("byTaNCfrQuarterPercent") );
	tancfrtenthpercentdiscr   -> push_back( it->tauID("byTaNCfrTenthPercent")   );
      }
      //
      // --  PAT V08-09-51 -- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X#V08_09_51  --
      //
      // HPSTau Discriminators are at: 
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.31.6.4&view=markup
      //
      // Clean HPS Taus are defined at:
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/cleaningLayer1/tauCleaner_cfi.py?revision=1.10&view=markup
      //
      // Out of the Box, Clean HPS Taus are given by cleanPatTaus: 
      // tauID("decayModeFinding") > 0.5                            ------ KEPT
      // tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5       ------ REMOVED
      // tauID("againstMuonMedium") > 0.5                           ------ REMOVED
      // tauID("againstElectronMedium") > 0.5                       ------ REMOVED
      // pt > 20.0                                                  ------ CHANGED to pt>15.0
      // abs(eta) < 2.3                                             ------ CHANGED to abs(eta)<2.5
      //
      if(isHPSTau){
	decaymodefindingdiscr      -> push_back( it->tauID("decayModeFinding")      );//applied in cleanPatTaus, KEPT in cmssw_cfg.py
	//
	againstelectronloosediscr        -> push_back( it->tauID("againstElectronLoose")        );
	againstelectronmediumdiscr       -> push_back( it->tauID("againstElectronMedium")       );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	againstelectrontightdiscr        -> push_back( it->tauID("againstElectronTight")        );
	againstelectronmvadiscr          -> push_back( it->tauID("againstElectronMVA")          );
	againstelectronmva2rawdiscr      -> push_back( it->tauID("againstElectronMVA2raw")      );
	againstelectronmva2categorydiscr -> push_back( it->tauID("againstElectronMVA2category") );
	againstelectronvloosemva2discr   -> push_back( it->tauID("againstElectronVLooseMVA2")   );
	againstelectronloosemva2discr    -> push_back( it->tauID("againstElectronLooseMVA2")    );
	againstelectronmediummva2discr   -> push_back( it->tauID("againstElectronMediumMVA2")   );
	againstelectrontightmva2discr    -> push_back( it->tauID("againstElectronTightMVA2")    );
	againstelectronmva3rawdiscr      -> push_back( it->tauID("againstElectronMVA3raw")      );
	againstelectronmva3categorydiscr -> push_back( it->tauID("againstElectronMVA3category") );
	againstelectronloosemva3discr    -> push_back( it->tauID("againstElectronLooseMVA3")    );
	againstelectronmediummva3discr   -> push_back( it->tauID("againstElectronMediumMVA3")   );
	againstelectrontightmva3discr    -> push_back( it->tauID("againstElectronTightMVA3")    );
	againstelectronvtightmva3discr   -> push_back( it->tauID("againstElectronVTightMVA3")   );
	againstelectrondeadecaldiscr     -> push_back( it->tauID("againstElectronDeadECAL")     );
	//
	againstmuonloosediscr   -> push_back( it->tauID("againstMuonLoose")   );
	againstmuonmediumdiscr  -> push_back( it->tauID("againstMuonMedium")  );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	againstmuontightdiscr   -> push_back( it->tauID("againstMuonTight")   );
	againstmuonloose2discr  -> push_back( it->tauID("againstMuonLoose2")  );
	againstmuonmedium2discr -> push_back( it->tauID("againstMuonMedium2") );
	againstmuontight2discr  -> push_back( it->tauID("againstMuonTight2")  );
	//
	vlooseisolationdiscr              -> push_back( it->tauID("byVLooseIsolation") );
	looseisolationdiscr               -> push_back( it->tauID("byLooseIsolation")  );
	mediumisolationdiscr              -> push_back( it->tauID("byMediumIsolation") );
	tightisolationdiscr               -> push_back( it->tauID("byTightIsolation" ) );
	vlooseisolationdeltabetacorrdiscr -> push_back( it->tauID("byVLooseIsolationDeltaBetaCorr") );
	looseisolationdeltabetacorrdiscr  -> push_back( it->tauID("byLooseIsolationDeltaBetaCorr")  );
	mediumisolationdeltabetacorrdiscr -> push_back( it->tauID("byMediumIsolationDeltaBetaCorr") );
	tightisolationdeltabetacorrdiscr  -> push_back( it->tauID("byTightIsolationDeltaBetaCorr")  );
	vloosecombinedisolationdeltabetacorrdiscr -> push_back( it->tauID("byVLooseCombinedIsolationDeltaBetaCorr") );
	loosecombinedisolationdeltabetacorrdiscr  -> push_back( it->tauID("byLooseCombinedIsolationDeltaBetaCorr" ) );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	mediumcombinedisolationdeltabetacorrdiscr -> push_back( it->tauID("byMediumCombinedIsolationDeltaBetaCorr") );
	tightcombinedisolationdeltabetacorrdiscr  -> push_back( it->tauID("byTightCombinedIsolationDeltaBetaCorr" ) );
	combinedisolationdeltabetacorr3hitsdiscr       -> push_back( it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits" )   );
	loosecombinedisolationdeltabetacorr3hitsdiscr  -> push_back( it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")  );
	mediumcombinedisolationdeltabetacorr3hitsdiscr -> push_back( it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") );
	tightcombinedisolationdeltabetacorr3hitsdiscr  -> push_back( it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")  );
	isolationmvarawdiscr     -> push_back( it->tauID("byIsolationMVAraw")     );
	looseisolationmvadiscr   -> push_back( it->tauID("byLooseIsolationMVA")   );
	mediumisolationmvadiscr  -> push_back( it->tauID("byMediumIsolationMVA")  );
	tightisolationmvadiscr   -> push_back( it->tauID("byTightIsolationMVA")   );
	looseisolationmva2discr  -> push_back( it->tauID("byLooseIsolationMVA2")  );
	mediumisolationmva2discr -> push_back( it->tauID("byMediumIsolationMVA2") );
	tightisolationmva2discr  -> push_back( it->tauID("byTightIsolationMVA2")  );
	//
      }
      //
      //
      eta -> push_back ( (double)(it -> eta()) ) ;
      phi -> push_back ( (double)(it -> phi()) ) ;
      pt  -> push_back ( (double)(it -> pt() ) ) ;
      et  -> push_back ( (double)(it -> et() ) ) ;
      charge  -> push_back ( (int)(it -> charge() ) ) ;
      if(  it ->isPFTau()  ){ispftau   -> push_back ( 1 ) ;}// all should be PF Tau
      if( !it ->isPFTau()  ){ispftau   -> push_back ( 0 ) ;}// all should be PF Tau
      if(  it ->isCaloTau()){iscalotau -> push_back ( 1 ) ;}// all should be PF Tau
      if( !it ->isCaloTau()){iscalotau -> push_back ( 0 ) ;}// all should be PF Tau
      decaymode                        -> push_back( (double)(it->decayMode())  );
      emfraction                       -> push_back( (double)(it->emFraction()) );
      hcal3x3overplead                 -> push_back( (double)(it->hcal3x3OverPLead()) );
      hcalmaxoverplead                 -> push_back( (double)(it->hcalMaxOverPLead()) );
      hcaltotoverplead                 -> push_back( (double)(it->hcalTotOverPLead()) );
      isolationpfchargedhadrcandsptsum -> push_back( (double)(it->isolationPFChargedHadrCandsPtSum()) );
      isolationpfgammacandsetsum       -> push_back( (double)(it->isolationPFGammaCandsEtSum())       );
      leadpfchargedhadrcandsignedsipt  -> push_back( (double)(it->leadPFChargedHadrCandsignedSipt())  );
      reco::PFCandidateRef leadPFChargedHadrCand_Ref = it->leadPFChargedHadrCand();
      if(leadPFChargedHadrCand_Ref.isNonnull()){// this check is needed in case hpsTau fails decayModeFinding.
	etaleadcharged                   -> push_back( (double)(leadPFChargedHadrCand_Ref->eta()) );
	phileadcharged                   -> push_back( (double)(leadPFChargedHadrCand_Ref->phi()) );
	ptleadcharged                    -> push_back( (double)(leadPFChargedHadrCand_Ref->pt())  );
      }
      phiphimoment                     -> push_back( (double)(it->phiphiMoment()) );
      etaetamoment                     -> push_back( (double)(it->etaetaMoment()) );
      etaphimoment                     -> push_back( (double)(it->etaphiMoment()) );
      ecalstripsumeoverplead           -> push_back( (double)(it->ecalStripSumEOverPLead()) );
      bremsrecoveryeoverplead          -> push_back( (double)(it->bremsRecoveryEOverPLead()) );
      maximumhcalpfclusteret           -> push_back( (double)(it->maximumHCALPFClusterEt()) );
      //
      // ----------- Vertex association ----------- //
      double minVtxDist3D = 9999.;
      int    vtxIndex_    = -1;
      double vtxDistXY_   = -9999.;
      double vtxDistZ_    = -9999.;      
      double vtx0DistXY_  = -9999.;
      double vtx0DistZ_   = -9999.;

      if( primaryVertices.isValid() ) {
	edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Primary Vertices: " << primaryVertices->size();

	int i_vertex = 0;

        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {
	  double distX  = (it->vertex()).x()-(v_it->position()).x();
	  double distY  = (it->vertex()).y()-(v_it->position()).y();
	  double distZ  = (it->vertex()).z()-(v_it->position()).z();
	  double distXY = sqrt(pow(distX,2) + pow(distY,2));
	  double dist3D = sqrt(pow(distXY,2) + pow(distZ,2));
	  
	  if ( i_vertex == 0 ) {  //leading vertex, by default sorted by sum(pt^2)
	    vtx0DistXY_ = distXY;
	    vtx0DistZ_  = distZ ;
	  }

          if( dist3D<minVtxDist3D ) {
            minVtxDist3D = dist3D;
            vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
            vtxDistXY_   = distXY;
            vtxDistZ_    = distZ;
          }

	  i_vertex++;
        }
      } else {
	edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << vtxInputTag;
      }
      vtxIndex                 -> push_back( vtxIndex_  );
      vtxDistXY                -> push_back( vtxDistXY_ );
      vtxDistZ                 -> push_back( vtxDistZ_  );
      vtx0DistXY               -> push_back( vtx0DistXY_);
      vtx0DistZ                -> push_back( vtx0DistZ_ );
      //
      // ----------- Gen-Reco Matching ----------- //
      double genparPt =-999.; double genjetPt =-999.;
      double genparEta=-999.; double genjetEta=-999.;
      double genparPhi=-999.; double genjetPhi=-999.;
      if ( !iEvent.isRealData() ) {
	for(uint igen = 0 ; igen < it->genParticleRefs().size() ; ++igen ){//genParticleRefs().size() is either 0 or 1
	  genparPt=it->genParticle(igen)->pt();
	  genparEta=it->genParticle(igen)->eta();
	  genparPhi=it->genParticle(igen)->phi();
	}
	if( it->genJet() ){
	  genjetPt=it->genJet()->pt();
	  genjetEta=it->genJet()->eta();
	  genjetPhi=it->genJet()->phi();
	}
      }
      matchedgenparticlept   -> push_back ( (double)(genparPt)  );
      matchedgenparticleeta  -> push_back ( (double)(genparEta) );
      matchedgenparticlephi  -> push_back ( (double)(genparPhi) );
      matchedgenjetpt        -> push_back ( (double)(genjetPt)  );
      matchedgenjeteta       -> push_back ( (double)(genjetEta) );
      matchedgenjetphi       -> push_back ( (double)(genjetPhi) );
      //
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //  --  User Isolation and isoDeposit Methods -- Feb 2012
      // 
      // User Isolation definitions and "keys" are at:
      //  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/PatCandidates/interface/Isolation.h?revision=1.7&view=markup
      //  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.27&view=markup
      //   it -> userIsolation(pat::PfChargedHadronIso);
      //   it -> userIsolation(pat::PfNeutralHadronIso);
      //   it -> userIsolation(pat::PfGammaIso);
      //
      // IsoDeposit is defined at:
      //  http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_2_5/doc/html/d4/d0a/classreco_1_1IsoDeposit.html
      //   (it -> isoDeposit(pat::PfChargedHadronIso)) -> candEnergy();
      //   (it -> isoDeposit(pat::PfChargedHadronIso)) -> depositWithin(0.5);
      //
      //
      // In comparison:
      // "pat::tau -> isolationPFCands()" is taken directly from the reco::PFTau object.
      // The parameters (cone size, pt threshold) are defined at:
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoTauTag/RecoTau/python/PFRecoTauQualityCuts_cfi.py?revision=1.4&view=markup
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoTauTag/RecoTau/python/PFRecoTauProducer_cfi.py?revision=1.23&view=markup
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/IsolationAlgos/plugins/PFTauExtractor.cc?revision=1.1&view=markup
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      if(isHPSTau){
	reco::PFCandidateRefVector signalPFChargedHadrCands_RefVector = it->signalPFChargedHadrCands();
	double signalPFChargedHadrCands_RefVector_size=0;
	if(signalPFChargedHadrCands_RefVector.isNonnull()){
	  signalPFChargedHadrCands_RefVector_size=(double)(signalPFChargedHadrCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalChargedHad     = signalPFChargedHadrCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalChargedHad_end = signalPFChargedHadrCands_RefVector.end();
	  for (; itSignalChargedHad != itSignalChargedHad_end; ++itSignalChargedHad ) {
	    signalpfchargedhadrcandspt  -> push_back((double)(*itSignalChargedHad)->pt());
	    signalpfchargedhadrcandseta -> push_back((double)(*itSignalChargedHad)->eta());
	    signalpfchargedhadrcandsphi -> push_back((double)(*itSignalChargedHad)->phi());
	  }
	}
	signalpfchargedhadrcandscount -> push_back(signalPFChargedHadrCands_RefVector_size);
	//
	reco::PFCandidateRefVector signalPFNeutrHadrCands_RefVector = it->signalPFNeutrHadrCands();
	double signalPFNeutrHadrCands_RefVector_size=0;
	if(signalPFNeutrHadrCands_RefVector.isNonnull()){
  	  signalPFNeutrHadrCands_RefVector_size=(double)(signalPFNeutrHadrCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalNeutrHadr     = signalPFNeutrHadrCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalNeutrHadr_end = signalPFNeutrHadrCands_RefVector.end();
	  for (; itSignalNeutrHadr != itSignalNeutrHadr_end; ++itSignalNeutrHadr ) {
	    signalpfneutrhadrcandspt  -> push_back((double)(*itSignalNeutrHadr)->pt());	
	    signalpfneutrhadrcandseta -> push_back((double)(*itSignalNeutrHadr)->eta());	
	    signalpfneutrhadrcandsphi -> push_back((double)(*itSignalNeutrHadr)->phi());	
	  }
	}
	signalpfneutrhadrcandscount -> push_back(signalPFNeutrHadrCands_RefVector_size);
	//
	reco::PFCandidateRefVector signalPFGammaCands_RefVector = it->signalPFGammaCands();
	double signalPFGammaCands_RefVector_size=0;
	if(signalPFGammaCands_RefVector.isNonnull()){
	  signalPFGammaCands_RefVector_size=(double)(signalPFGammaCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalGamma     = signalPFGammaCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalGamma_end = signalPFGammaCands_RefVector.end();
	  for (; itSignalGamma != itSignalGamma_end; ++itSignalGamma ) {
	    signalpfgammacandspt  -> push_back((double)(*itSignalGamma)->pt());
	    signalpfgammacandseta -> push_back((double)(*itSignalGamma)->eta());
	    signalpfgammacandsphi -> push_back((double)(*itSignalGamma)->phi());
	  }
	}
	signalpfgammacandscount -> push_back(signalPFGammaCands_RefVector_size);
	//	
	// --------------------------------------------------------------------------------------- //
	// HPS Tau Optional Isolation information 
	//reco::PFCandidateRefVector isoPFChargedHadrCands_RefVector = it->isolationPFChargedHadrCands();
	//isolationpfchargedhadrcandscount ->push_back((double)(isoPFChargedHadrCands_RefVector.size()));
	//if(isoPFChargedHadrCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoChargedHad     = isoPFChargedHadrCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoChargedHad_end = isoPFChargedHadrCands_RefVector.end();
	//  for (; itIsoChargedHad != itIsoChargedHad_end; ++itIsoChargedHad ) {
	//    isolationpfchargedhadrcandspt  -> push_back((double)(*itIsoChargedHad)->pt());
	//    isolationpfchargedhadrcandseta -> push_back((double)(*itIsoChargedHad)->eta());
	//    isolationpfchargedhadrcandsphi -> push_back((double)(*itIsoChargedHad)->phi());
	//  }
	//}
	//reco::PFCandidateRefVector isoPFNeutrHadrCands_RefVector = it->isolationPFNeutrHadrCands();
	//isolationpfneutrhadrcandscount -> push_back((double)(isoPFNeutrHadrCands_RefVector.size()));
	//if(isoPFNeutrHadrCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoNeutrHadr     = isoPFNeutrHadrCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoNeutrHadr_end = isoPFNeutrHadrCands_RefVector.end();
	//  for (; itIsoNeutrHadr != itIsoNeutrHadr_end; ++itIsoNeutrHadr ) {
	//    isolationpfneutrhadrcandspt  -> push_back((double)(*itIsoNeutrHadr)->pt());
	//    isolationpfneutrhadrcandseta -> push_back((double)(*itIsoNeutrHadr)->eta());
	//    isolationpfneutrhadrcandsphi -> push_back((double)(*itIsoNeutrHadr)->phi());
	//  }
	//}
	//reco::PFCandidateRefVector isoPFGammaCands_RefVector = it->isolationPFGammaCands();
	//isolationpfgammacandscount -> push_back((double)(isoPFGammaCands_RefVector.size()));
	//if(isoPFGammaCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoGamma     = isoPFGammaCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoGamma_end = isoPFGammaCands_RefVector.end();
	//  for (; itIsoGamma != itIsoGamma_end; ++itIsoGamma ) {
	//    isolationpfgammacandspt  -> push_back((double)(*itIsoGamma)->pt());
	//    isolationpfgammacandseta -> push_back((double)(*itIsoGamma)->eta());
	//    isolationpfgammacandsphi -> push_back((double)(*itIsoGamma)->phi());
	//  }
	//}
	// --------------------------------------------------------------------------------------- //
	//
      }
      //
    }
  } else {
    edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << inputTag;
  }
  iEvent.put( eta,                              prefix + "Eta"                               + suffix );
  iEvent.put( phi,                              prefix + "Phi"                               + suffix );
  iEvent.put( pt,                               prefix + "Pt"                                + suffix );
  iEvent.put( et,                               prefix + "Et"                                + suffix );
  iEvent.put( charge,                           prefix + "Charge"                            + suffix );
  iEvent.put( ispftau,                          prefix + "IsPFTau"                           + suffix );
  iEvent.put( iscalotau,                        prefix + "IsCaloTau"                         + suffix );
  iEvent.put( decaymode,                        prefix + "DecayMode"                         + suffix );
  iEvent.put( emfraction,                       prefix + "EmFraction"                        + suffix );
  iEvent.put( hcal3x3overplead,                 prefix + "Hcal3x3OverPLead"                  + suffix );
  iEvent.put( hcalmaxoverplead,                 prefix + "HcalMaxOverPLead"                  + suffix );
  iEvent.put( hcaltotoverplead,                 prefix + "HcalTotOverPLead"                  + suffix );
  iEvent.put( isolationpfchargedhadrcandsptsum, prefix + "IsolationPFChargedHadrCandsPtSum"  + suffix );
  iEvent.put( isolationpfgammacandsetsum,       prefix + "IsolationPFGammaCandsEtSum"        + suffix );
  iEvent.put( leadpfchargedhadrcandsignedsipt,  prefix + "LeadPFChargedHadrCandsignedSipt"   + suffix );
  iEvent.put( etaleadcharged,                   prefix + "EtaLeadCharged"                    + suffix );
  iEvent.put( phileadcharged,                   prefix + "PhiLeadCharged"                    + suffix );
  iEvent.put( ptleadcharged,                    prefix + "PtLeadCharged"                     + suffix );
  iEvent.put( phiphimoment,                     prefix + "PhiphiMoment"                      + suffix );
  iEvent.put( etaetamoment,                     prefix + "EtaetaMoment"                      + suffix );
  iEvent.put( etaphimoment,                     prefix + "EtaphiMoment"                      + suffix );
  iEvent.put( ecalstripsumeoverplead,           prefix + "EcalStripSumEOverPLead"            + suffix );
  iEvent.put( bremsrecoveryeoverplead,          prefix + "BremsRecoveryEOverPLead"           + suffix );
  iEvent.put( maximumhcalpfclusteret,           prefix + "MaximumHCALPFClusterEt"            + suffix );
  iEvent.put( matchedgenparticlept,             prefix + "MatchedGenParticlePt"              + suffix );
  iEvent.put( matchedgenparticleeta,            prefix + "MatchedGenParticleEta"             + suffix );
  iEvent.put( matchedgenparticlephi,            prefix + "MatchedGenParticlePhi"             + suffix );
  iEvent.put( matchedgenjetpt,                  prefix + "MatchedGenJetPt"                   + suffix );
  iEvent.put( matchedgenjeteta,                 prefix + "MatchedGenJetEta"                  + suffix );
  iEvent.put( matchedgenjetphi,                 prefix + "MatchedGenJetPhi"                  + suffix );
  //
  if(isSCTau){
    iEvent.put( leadingtrackfindingdiscr,               prefix + "LeadingTrackFindingDiscr"            + suffix );
    iEvent.put( leadingtrackptcutdiscr,                 prefix + "LeadingTrackPtCutDiscr"              + suffix );
    iEvent.put( leadingpionptcutdiscr,                  prefix + "LeadingPionPtCutDiscr"               + suffix );
    iEvent.put( isolationdiscr,                         prefix + "IsolationDiscr"                      + suffix );
    iEvent.put( trackisolationdiscr,                    prefix + "TrackIsolationDiscr"                 + suffix );
    iEvent.put( ecalisolationdiscr,                     prefix + "EcalIsolationDiscr"                  + suffix );
    iEvent.put( isolationusingleadingpiondiscr,         prefix + "IsolationUsingLeadingPionDiscr"      + suffix );
    iEvent.put( trackisolationusingleadingpiondiscr,    prefix + "TrackIsolationUsingLeadingPionDiscr" + suffix );
    iEvent.put( ecalisolationusingleadingpiondiscr,     prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
    iEvent.put( againstelectrondiscr,                   prefix + "AgainstElectronDiscr"                + suffix );
    iEvent.put( againstmuondiscr,                       prefix + "AgainstMuonDiscr"                    + suffix );
    iEvent.put( tancdiscr,                              prefix + "TaNCDiscr"                           + suffix );
    iEvent.put( tancfronepercentdiscr,                  prefix + "TaNCfrOnePercentDiscr"               + suffix );
    iEvent.put( tancfrhalfpercentdiscr,                 prefix + "TaNCfrHalfPercentDiscr"              + suffix );
    iEvent.put( tancfrquarterpercentdiscr,              prefix + "TaNCfrQuarterPercentDiscr"           + suffix );
    iEvent.put( tancfrtenthpercentdiscr,                prefix + "TaNCfrTenthPercentDiscr"             + suffix );
  }
  //
  if(isHPSTau){
    iEvent.put( decaymodefindingdiscr,                     prefix + "DecayModeFindingDiscr"                       + suffix );
    //
    iEvent.put( againstelectronloosediscr,                 prefix + "AgainstElectronLooseDiscr"                   + suffix );
    iEvent.put( againstelectronmediumdiscr,                prefix + "AgainstElectronMediumDiscr"                  + suffix );
    iEvent.put( againstelectrontightdiscr,                 prefix + "AgainstElectronTightDiscr"                   + suffix );
    iEvent.put( againstelectronmvadiscr,                   prefix + "AgainstElectronMVADiscr"                     + suffix );
    iEvent.put( againstelectronmva2rawdiscr,               prefix + "AgainstElectronMVA2rawDiscr"                 + suffix );
    iEvent.put( againstelectronmva2categorydiscr,          prefix + "AgainstElectronMVA2categoryDiscr"            + suffix );
    iEvent.put( againstelectronvloosemva2discr,            prefix + "AgainstElectronVLooseMVA2Discr"              + suffix );
    iEvent.put( againstelectronloosemva2discr,             prefix + "AgainstElectronLooseMVA2Discr"               + suffix );
    iEvent.put( againstelectronmediummva2discr,            prefix + "AgainstElectronMediumMVA2Discr"              + suffix );
    iEvent.put( againstelectrontightmva2discr,             prefix + "AgainstElectronTightMVA2Discr"               + suffix );
    iEvent.put( againstelectronmva3rawdiscr,               prefix + "AgainstElectronMVA3rawDiscr"                 + suffix );
    iEvent.put( againstelectronmva3categorydiscr,          prefix + "AgainstElectronMVA3categoryDiscr"            + suffix );
    iEvent.put( againstelectronloosemva3discr,             prefix + "AgainstElectronLooseMVA3Discr"               + suffix );
    iEvent.put( againstelectronmediummva3discr,            prefix + "AgainstElectronMediumMVA3Discr"              + suffix );
    iEvent.put( againstelectrontightmva3discr,             prefix + "AgainstElectronTightMVA3Discr"               + suffix );
    iEvent.put( againstelectronvtightmva3discr,            prefix + "AgainstElectronVTightMVA3Discr"              + suffix );
    iEvent.put( againstelectrondeadecaldiscr,              prefix + "AgainstElectronDeadECALDiscr"                + suffix );
    //
    iEvent.put( againstmuonloosediscr,                     prefix + "AgainstMuonLooseDiscr"                       + suffix );
    iEvent.put( againstmuonmediumdiscr,                    prefix + "AgainstMuonMediumDiscr"                      + suffix );
    iEvent.put( againstmuontightdiscr,                     prefix + "AgainstMuonTightDiscr"                       + suffix );
    iEvent.put( againstmuonloose2discr,                    prefix + "AgainstMuonLoose2Discr"                      + suffix );
    iEvent.put( againstmuonmedium2discr,                   prefix + "AgainstMuonMedium2Discr"                     + suffix );
    iEvent.put( againstmuontight2discr,                    prefix + "AgainstMuonTight2Discr"                      + suffix );

    iEvent.put( vlooseisolationdiscr,                           prefix + "VLooseIsolationDiscr"                           + suffix );
    iEvent.put( looseisolationdiscr,                            prefix + "LooseIsolationDiscr"                            + suffix );
    iEvent.put( mediumisolationdiscr,                           prefix + "MediumIsolationDiscr"                           + suffix );
    iEvent.put( tightisolationdiscr,                            prefix + "TightIsolationDiscr"                            + suffix );
    iEvent.put( vlooseisolationdeltabetacorrdiscr,              prefix + "VLooseIsolationDeltaBetaCorrDiscr"              + suffix );
    iEvent.put( looseisolationdeltabetacorrdiscr,               prefix + "LooseIsolationDeltaBetaCorrDiscr"               + suffix );
    iEvent.put( mediumisolationdeltabetacorrdiscr,              prefix + "MediumIsolationDeltaBetaCorrDiscr"              + suffix );
    iEvent.put( tightisolationdeltabetacorrdiscr,               prefix + "TightIsolationDeltaBetaCorrDiscr"               + suffix );
    iEvent.put( vloosecombinedisolationdeltabetacorrdiscr,      prefix + "VLooseCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    iEvent.put( loosecombinedisolationdeltabetacorrdiscr,       prefix + "LooseCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    iEvent.put( mediumcombinedisolationdeltabetacorrdiscr,      prefix + "MediumCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    iEvent.put( tightcombinedisolationdeltabetacorrdiscr,       prefix + "TightCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    iEvent.put( combinedisolationdeltabetacorr3hitsdiscr,       prefix + "CombinedIsolationDeltaBetaCorr3HitsDiscr"       + suffix );
    iEvent.put( loosecombinedisolationdeltabetacorr3hitsdiscr,  prefix + "LooseCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    iEvent.put( mediumcombinedisolationdeltabetacorr3hitsdiscr, prefix + "MediumCombinedIsolationDeltaBetaCorr3HitsDiscr" + suffix );
    iEvent.put( tightcombinedisolationdeltabetacorr3hitsdiscr,  prefix + "TightCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    iEvent.put( isolationmvarawdiscr,                           prefix + "IsolationMVArawDiscr"                           + suffix );
    iEvent.put( looseisolationmvadiscr,                         prefix + "LooseIsolationMVADiscr"                         + suffix );
    iEvent.put( mediumisolationmvadiscr,                        prefix + "MediumIsolationMVADiscr"                        + suffix );
    iEvent.put( tightisolationmvadiscr,                         prefix + "TightIsolationMVADiscr"                         + suffix );
    iEvent.put( looseisolationmva2discr,                        prefix + "LooseIsolationMVA2Discr"                        + suffix );
    iEvent.put( mediumisolationmva2discr,                       prefix + "MediumIsolationMVA2Discr"                       + suffix );
    iEvent.put( tightisolationmva2discr,                        prefix + "TightIsolationMVA2Discr"                        + suffix );
    //
    iEvent.put( vtxIndex,    prefix + "VtxIndex"       + suffix );
    iEvent.put( vtxDistXY,   prefix + "VtxDistXY"      + suffix );
    iEvent.put( vtxDistZ,    prefix + "VtxDistZ"       + suffix );
    iEvent.put( vtx0DistXY,  prefix + "LeadVtxDistXY"  + suffix );
    iEvent.put( vtx0DistZ,   prefix + "LeadVtxDistZ"   + suffix );
    //
    iEvent.put( signalpfchargedhadrcandspt,   prefix +  "SignalPFChargedHadrCandsPt"   + suffix );
    iEvent.put( signalpfchargedhadrcandseta,  prefix +  "SignalPFChargedHadrCandsEta"  + suffix );
    iEvent.put( signalpfchargedhadrcandsphi,  prefix +  "SignalPFChargedHadrCandsPhi"  + suffix );
    iEvent.put( signalpfchargedhadrcandscount,prefix +  "SignalPFChargedHadrCandsCount"+ suffix );
    iEvent.put( signalpfneutrhadrcandspt,     prefix +  "SignalPFNeutrHadrCandsPt"     + suffix );
    iEvent.put( signalpfneutrhadrcandseta,    prefix +  "SignalPFNeutrHadrCandsEta"    + suffix );
    iEvent.put( signalpfneutrhadrcandsphi,    prefix +  "SignalPFNeutrHadrCandsPhi"    + suffix );
    iEvent.put( signalpfneutrhadrcandscount,  prefix +  "SignalPFNeutrHadrCandsCount"  + suffix );
    iEvent.put( signalpfgammacandspt,         prefix +  "SignalPFGammaCandsPt"         + suffix );
    iEvent.put( signalpfgammacandseta,        prefix +  "SignalPFGammaCandsEta"        + suffix );
    iEvent.put( signalpfgammacandsphi,        prefix +  "SignalPFGammaCandsPhi"        + suffix );
    iEvent.put( signalpfgammacandscount,      prefix +  "SignalPFGammaCandsCount"      + suffix );
    //
    // --------------------------------------------------------------------------------------- //
    // HPS Tau Optional Isolation information 
    //iEvent.put( isolationpfchargedhadrcandspt,   prefix +  "IsolationPFChargedHadrCandsPt"   + suffix );
    //iEvent.put( isolationpfchargedhadrcandseta,  prefix +  "IsolationPFChargedHadrCandsEta"  + suffix );
    //iEvent.put( isolationpfchargedhadrcandsphi,  prefix +  "IsolationPFChargedHadrCandsPhi"  + suffix );
    //iEvent.put( isolationpfchargedhadrcandscount,prefix +  "IsolationPFChargedHadrCandsCount"+ suffix );
    //iEvent.put( isolationpfneutrhadrcandspt,     prefix +  "IsolationPFNeutrHadrCandsPt"     + suffix );
    //iEvent.put( isolationpfneutrhadrcandseta,    prefix +  "IsolationPFNeutrHadrCandsEta"    + suffix );
    //iEvent.put( isolationpfneutrhadrcandsphi,    prefix +  "IsolationPFNeutrHadrCandsPhi"    + suffix );
    //iEvent.put( isolationpfneutrhadrcandscount,  prefix +  "IsolationPFNeutrHadrCandsCount"  + suffix );
    //iEvent.put( isolationpfgammacandspt,         prefix +  "IsolationPFGammaCandsPt"         + suffix );
    //iEvent.put( isolationpfgammacandseta,        prefix +  "IsolationPFGammaCandsEta"        + suffix );
    //iEvent.put( isolationpfgammacandsphi,        prefix +  "IsolationPFGammaCandsPhi"        + suffix );
    //iEvent.put( isolationpfgammacandscount,      prefix +  "IsolationPFGammaCandsCount"      + suffix );
    // --------------------------------------------------------------------------------------- //
  }
}
