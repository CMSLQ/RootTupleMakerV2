#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


RootTupleMakerV2_GenEventInfo::RootTupleMakerV2_GenEventInfo(const edm::ParameterSet& iConfig) :
    genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
    storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
    pdfWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFWeightsInputTag")),
    pileupInfoSrc(iConfig.getParameter<edm::InputTag>("pileupInfo"))
{
  produces <unsigned int> ( "ProcessID" );
  produces <double>       ( "PtHat" );
  produces <std::vector<double> > ( "PDFWeights" );
  produces <std::vector<int> > ( "PileUpInteractions");
  produces <std::vector<int> > ( "PileUpOriginBX" ) ;
}

void RootTupleMakerV2_GenEventInfo::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int >         processID   ( new unsigned int() );
  std::auto_ptr<double >               ptHat ( new double() );
  std::auto_ptr<std::vector<double> >  pdfWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int >  >   Number_interactions  ( new std::vector<int>() );
  std::auto_ptr<std::vector<int >  >   OriginBX( new std::vector<int>() );

  *processID.get() = 0;
  *ptHat.get() = 0.;

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(genEvtInfoInputTag, genEvtInfoProduct);

    if( genEvtInfoProduct.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << genEvtInfoInputTag;

      *processID.get() = genEvtInfoProduct->signalProcessID();
      *ptHat.get() = ( genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : 0. );

    } else {
      edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << genEvtInfoInputTag;
    }
    // PDF Weights Part
    if( storePDFWeights ) {
      edm::Handle<std::vector<double> > pdfWeightsHandle;
      iEvent.getByLabel(pdfWeightsInputTag, pdfWeightsHandle);

      if( pdfWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << pdfWeightsInputTag;

        *pdfWeights.get() = *pdfWeightsHandle;

      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << pdfWeightsInputTag;
      }
    }
    // PileupSummary Part
    edm::Handle<std::vector<PileupSummaryInfo> >  puInfo;
    iEvent.getByLabel(pileupInfoSrc, puInfo);
    
    if(puInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = puInfo->begin(); it != puInfo->end(); ++it ) {
	Number_interactions -> push_back ( it->getPU_NumInteractions() );
	OriginBX -> push_back ( it -> getBunchCrossing());
      }
    }
    else {
      edm::LogError("RootTupleMakerV2_PileUpError") << "Error! Can't get the product " << pileupInfoSrc;
    }
  }

  //-----------------------------------------------------------------
  iEvent.put( processID, "ProcessID" );
  iEvent.put( ptHat, "PtHat" );
  iEvent.put( pdfWeights, "PDFWeights" );
  iEvent.put( Number_interactions,   "PileUpInteractions"   );
  iEvent.put( OriginBX,   "PileUpOriginBX" );
}
