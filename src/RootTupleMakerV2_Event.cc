#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Event.h"
#include "FWCore/Framework/interface/Event.h"

RootTupleMakerV2_Event::RootTupleMakerV2_Event(const edm::ParameterSet& iConfig) :
  //fastJetForIsolationInputTag(iConfig.getParameter<edm::InputTag>("FastJetForIsolationInputTag")),//not used in 2012
  fastJetForJECInputTag(iConfig.getParameter<edm::InputTag>("FastJetForJECInputTag")),
  fastJetForHEEPInputTag(iConfig.getParameter<edm::InputTag>("FastJetForHEEPInputTag")),
  fastJetForJECCCPUInputTag(iConfig.getParameter<edm::InputTag>("FastJetForJECCCPUInputTag")),
  fastJetForJECCNInputTag(iConfig.getParameter<edm::InputTag>("FastJetForJECCNInputTag")),
  fastJetForJECCNTInputTag(iConfig.getParameter<edm::InputTag>("FastJetForJECCNTInputTag"))
{
  produces <unsigned int> ( "run"    );
  produces <unsigned int> ( "event"  );
  produces <unsigned int> ( "bunch"  );
  produces <unsigned int> ( "ls"     );
  produces <unsigned int> ( "orbit"  );
  produces <double>       ( "time"   );
  produces <bool>         ( "isData" );
  //produces <double>     ( "rhoIso" );//not used in 2012 
  produces <double>       ( "rhoJets"     );
  produces <double>       ( "rhoForHEEP"  );
  produces <double>       ( "rhoJetsCCPU" );
  produces <double>       ( "rhoJetsCN"   );
  produces <double>       ( "rhoJetsCNT"  );
}

void RootTupleMakerV2_Event::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int >  run   ( new unsigned int(iEvent.id().run()        ) );
  std::auto_ptr<unsigned int >  event ( new unsigned int(iEvent.id().event()      ) );
  std::auto_ptr<unsigned int >  ls    ( new unsigned int(iEvent.luminosityBlock() ) );

  double sec  = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;

  std::auto_ptr<unsigned int >  bunch ( new unsigned int(iEvent.bunchCrossing()   ) );
  std::auto_ptr<unsigned int >  orbit ( new unsigned int(iEvent.orbitNumber()     ) );
  std::auto_ptr<double >        time  ( new double(sec+usec/conv));

  std::auto_ptr<bool >          isdata  ( new bool(iEvent.isRealData()));

  // edm::Handle<double> rhoH;
  // iEvent.getByLabel(fastJetForIsolationInputTag,rhoH);
  // std::auto_ptr<double >        rhoIso  ( new double( *rhoH.product() ) );

  edm::Handle<double> rhoHJets;
  iEvent.getByLabel(fastJetForJECInputTag,rhoHJets);
  std::auto_ptr<double >        rhoJets     ( new double( *rhoHJets.product() )     );

  edm::Handle<double> rhoHForHEEP;
  iEvent.getByLabel(fastJetForHEEPInputTag, rhoHForHEEP);
  std::auto_ptr<double >        rhoForHEEP  ( new double( *rhoHForHEEP.product() )     );

  edm::Handle<double> rhoHJetsCCPU;
  iEvent.getByLabel(fastJetForJECCCPUInputTag,rhoHJetsCCPU);
  std::auto_ptr<double >        rhoJetsCCPU ( new double( *rhoHJetsCCPU.product() ) );

  edm::Handle<double> rhoHJetsCN;
  iEvent.getByLabel(fastJetForJECCNInputTag,rhoHJetsCN);
  std::auto_ptr<double >        rhoJetsCN   ( new double( *rhoHJetsCN.product() )   );

  edm::Handle<double> rhoHJetsCNT;
  iEvent.getByLabel(fastJetForJECCNTInputTag,rhoHJetsCNT);
  std::auto_ptr<double >        rhoJetsCNT  ( new double( *rhoHJetsCNT.product() )  );



  //-----------------------------------------------------------------
  iEvent.put( run,   "run"   );
  iEvent.put( event, "event" );
  iEvent.put( ls   , "ls"    );
  iEvent.put( bunch, "bunch" );
  iEvent.put( orbit, "orbit" );
  iEvent.put( time,  "time"  );
  iEvent.put( isdata,"isData");
  //iEvent.put( rhoIso,   "rhoIso"   );
  iEvent.put( rhoJets,     "rhoJets"     );
  iEvent.put( rhoForHEEP,  "rhoForHEEP"  );
  iEvent.put( rhoJetsCCPU, "rhoJetsCCPU" );
  iEvent.put( rhoJetsCN,   "rhoJetsCN"   );
  iEvent.put( rhoJetsCNT,  "rhoJetsCNT"  );
}
