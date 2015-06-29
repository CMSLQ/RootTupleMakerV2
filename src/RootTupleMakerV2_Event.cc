#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Event.h"
#include "FWCore/Framework/interface/Event.h"

RootTupleMakerV2_Event::RootTupleMakerV2_Event(const edm::ParameterSet& iConfig) :
  fixedGridRhoAllInputTag(iConfig.getParameter<edm::InputTag>("FixedGridRhoAllInputTag")),
  fixedGridRhoFastjetAllCaloInputTag(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetAllCaloInputTag")),
  fixedGridRhoFastjetCentralCaloInputTag(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralCaloInputTag")),
  fixedGridRhoFastjetCentralChargedPileUpInputTag(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralChargedPileUpInputTag")),
  fixedGridRhoFastjetCentralNeutralInputTag(iConfig.getParameter<edm::InputTag>("FixedGridRhoFastjetCentralNeutralInputTag")),
  fillRhos (iConfig.getUntrackedParameter<bool>("FillRhos",true))
{
  produces <unsigned int> ( "run"    );
  produces <unsigned int> ( "event"  );
  produces <unsigned int> ( "bunch"  );
  produces <unsigned int> ( "ls"     );
  produces <unsigned int> ( "orbit"  );
  produces <double>       ( "time"   );
  produces <bool>         ( "isData" );
  produces <double>       ( "fixedGridRhoAll"                         );
  produces <double>       ( "fixedGridRhoFastjetAllCalo"              );
  produces <double>       ( "fixedGridRhoFastjetCentralCalo"          );
  produces <double>       ( "fixedGridRhoFastjetCentralChargedPileUp" );
  produces <double>       ( "fixedGridRhoFastjetCentralNeutral"       );
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

  std::auto_ptr<double> fixedGridRhoAll (new double(-999));
  std::auto_ptr<double> fixedGridRhoFastjetAllCalo (new double(-999));
  std::auto_ptr<double> fixedGridRhoFastjetCentralCalo (new double(-999));
  std::auto_ptr<double> fixedGridRhoFastjetCentralChargedPileUp (new double(-999));
  std::auto_ptr<double> fixedGridRhoFastjetCentralNeutral (new double(-999));
  if(fillRhos)
  {
    edm::Handle<double> fixedGridRhoAllH;
    iEvent.getByLabel(fixedGridRhoAllInputTag,fixedGridRhoAllH);
    fixedGridRhoAll.reset(new double(*fixedGridRhoAllH.product()));

    edm::Handle<double> fixedGridRhoFastjetAllCaloH;
    iEvent.getByLabel(fixedGridRhoFastjetAllCaloInputTag,fixedGridRhoFastjetAllCaloH);
    fixedGridRhoFastjetAllCalo.reset(new double(*fixedGridRhoFastjetAllCaloH.product() ) );

    edm::Handle<double> fixedGridRhoFastjetCentralCaloH;
    iEvent.getByLabel(fixedGridRhoFastjetCentralCaloInputTag,fixedGridRhoFastjetCentralCaloH);
    fixedGridRhoFastjetCentralCalo.reset(new double(*fixedGridRhoFastjetCentralCaloH.product() ) );

    edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUpH;
    iEvent.getByLabel(fixedGridRhoFastjetCentralChargedPileUpInputTag,fixedGridRhoFastjetCentralChargedPileUpH);
    fixedGridRhoFastjetCentralChargedPileUp.reset(new double(*fixedGridRhoFastjetCentralChargedPileUpH.product() ) );

    edm::Handle<double> fixedGridRhoFastjetCentralNeutralH;
    iEvent.getByLabel(fixedGridRhoFastjetCentralNeutralInputTag,fixedGridRhoFastjetCentralNeutralH);
    fixedGridRhoFastjetCentralNeutral.reset(new double(*fixedGridRhoFastjetCentralNeutralH.product() ) );
  }




  //-----------------------------------------------------------------
  iEvent.put( run,   "run"   );
  iEvent.put( event, "event" );
  iEvent.put( ls   , "ls"    );
  iEvent.put( bunch, "bunch" );
  iEvent.put( orbit, "orbit" );
  iEvent.put( time,  "time"  );
  iEvent.put( isdata,"isData");
  iEvent.put( fixedGridRhoAll,                         "fixedGridRhoAll"                         );
  iEvent.put( fixedGridRhoFastjetAllCalo,              "fixedGridRhoFastjetAllCalo"              );
  iEvent.put( fixedGridRhoFastjetCentralCalo,          "fixedGridRhoFastjetCentralCalo"          );
  iEvent.put( fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp" );
  iEvent.put( fixedGridRhoFastjetCentralNeutral,       "fixedGridRhoFastjetCentralNeutral"       );
}
