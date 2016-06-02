#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Event.h"
#include "FWCore/Framework/interface/Event.h"

RootTupleMakerV2_Event::RootTupleMakerV2_Event(const edm::ParameterSet& iConfig) :
  globalTagLabel_(iConfig.getParameter<std::string>("globalTag")),
  globalTagToken_(consumes<std::string>(iConfig.getParameter<std::string>("globalTag"))),
  fixedGridRhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoAllInputTag"))),
  fixedGridRhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetAllInputTag"))),
  fixedGridRhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetAllCaloInputTag"))),
  fixedGridRhoFastjetCentralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralInputTag"))),
  fixedGridRhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralCaloInputTag"))),
  fixedGridRhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralNeutralInputTag"))),
  fixedGridRhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralChargedPileUpInputTag")))
{
  //produces <std::string>  ("globalTag");
  //produces <std::string>  ("cmsswRelease");
  produces <unsigned int> ( "run"    );
  produces <unsigned int> ( "event"  );
  produces <unsigned int> ( "bunch"  );
  produces <unsigned int> ( "ls"     );
  produces <unsigned int> ( "orbit"  );
  produces <double>       ( "time"   );
  produces <bool>         ( "isData" );
  produces <double>       ( "fixedGridRhoAll"                         );
  produces <double>       ( "fixedGridRhoFastjetAll"              );
  produces <double>       ( "fixedGridRhoFastjetAllCalo"              );
  produces <double>       ( "fixedGridRhoFastjetCentral"          );
  produces <double>       ( "fixedGridRhoFastjetCentralCalo"          );
  produces <double>       ( "fixedGridRhoFastjetCentralNeutral"       );
  produces <double>       ( "fixedGridRhoFastjetCentralChargedPileUp" );
}

void RootTupleMakerV2_Event::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //std::string releaseVersion = (iEvent.processHistory().rbegin())->releaseVersion();
  //std::string releaseVersion = "";
  //std::auto_ptr<std::string>  globaltag   ( new std::string(globalTagLabel_) );
  //std::auto_ptr<std::string>  cmsswrelease   ( new std::string(releaseVersion) );
  //std::auto_ptr<std::string>  cmsswrelease   ( new std::string((iEvent.processHistory().rbegin())->releaseVersion()) );
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

  edm::Handle<double> fixedGridRhoAllH;
  iEvent.getByToken(fixedGridRhoAllToken_,fixedGridRhoAllH);
  std::auto_ptr<double> fixedGridRhoAll (new double(*fixedGridRhoAllH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetAllH;
  iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAllH);
  std::auto_ptr<double> fixedGridRhoFastjetAll (new double(*fixedGridRhoFastjetAllH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetAllCaloH;
  iEvent.getByToken(fixedGridRhoFastjetAllCaloToken_,fixedGridRhoFastjetAllCaloH);
  std::auto_ptr<double> fixedGridRhoFastjetAllCalo (new double(*fixedGridRhoFastjetAllCaloH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetCentralH;
  iEvent.getByToken(fixedGridRhoFastjetCentralToken_,fixedGridRhoFastjetCentralH);
  std::auto_ptr<double> fixedGridRhoFastjetCentral (new double(*fixedGridRhoFastjetCentralH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetCentralCaloH;
  iEvent.getByToken(fixedGridRhoFastjetCentralCaloToken_,fixedGridRhoFastjetCentralCaloH);
  std::auto_ptr<double> fixedGridRhoFastjetCentralCalo (new double(*fixedGridRhoFastjetCentralCaloH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetCentralNeutralH;
  iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken_,fixedGridRhoFastjetCentralNeutralH);
  std::auto_ptr<double> fixedGridRhoFastjetCentralNeutral (new double(*fixedGridRhoFastjetCentralNeutralH.product() ) );

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUpH;
  iEvent.getByToken(fixedGridRhoFastjetCentralChargedPileUpToken_,fixedGridRhoFastjetCentralChargedPileUpH);
  std::auto_ptr<double> fixedGridRhoFastjetCentralChargedPileUp (new double(*fixedGridRhoFastjetCentralChargedPileUpH.product() ) );




  //-----------------------------------------------------------------
  //iEvent.put( globaltag, "globalTag");
  //iEvent.put( cmsswrelease, "cmsswRelease");
  iEvent.put( run,   "run"   );
  iEvent.put( event, "event" );
  iEvent.put( ls   , "ls"    );
  iEvent.put( bunch, "bunch" );
  iEvent.put( orbit, "orbit" );
  iEvent.put( time,  "time"  );
  iEvent.put( isdata,"isData");
  iEvent.put( fixedGridRhoAll,                         "fixedGridRhoAll"                         );
  iEvent.put( fixedGridRhoFastjetAll,                  "fixedGridRhoFastjetAll"                  );
  iEvent.put( fixedGridRhoFastjetAllCalo,              "fixedGridRhoFastjetAllCalo"              );
  iEvent.put( fixedGridRhoFastjetCentral,              "fixedGridRhoFastjetCentral"              );
  iEvent.put( fixedGridRhoFastjetCentralCalo,          "fixedGridRhoFastjetCentralCalo"          );
  iEvent.put( fixedGridRhoFastjetCentralNeutral,       "fixedGridRhoFastjetCentralNeutral"       );
  iEvent.put( fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp" );
}
