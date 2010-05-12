#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


RootTupleMakerV2_GenEventInfo::RootTupleMakerV2_GenEventInfo(const edm::ParameterSet& iConfig) :
    genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
    storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
    pdfWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFWeightsInputTag"))
{
  produces <unsigned int> ( "ProcessID" );
  produces <double>       ( "PtHat" );
  produces <std::vector<double> > ( "PDFWeights" );
}

void RootTupleMakerV2_GenEventInfo::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int >         processID   ( new unsigned int() );
  std::auto_ptr<double >               ptHat ( new double() );
  std::auto_ptr<std::vector<double> >  pdfWeights  ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(genEvtInfoInputTag, genEvtInfoProduct);

    if( genEvtInfoProduct.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << genEvtInfoInputTag;

      *processID.get() = genEvtInfoProduct->signalProcessID();
      *ptHat.get() = ( genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : 0.0 );

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
  }

  //-----------------------------------------------------------------
  iEvent.put( processID, "ProcessID" );
  iEvent.put( ptHat, "PtHat" );
  iEvent.put( pdfWeights, "PDFWeights" );
}
