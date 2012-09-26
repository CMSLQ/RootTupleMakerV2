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
    pdfCTEQWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFCTEQWeightsInputTag")),
    pdfMSTWWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFMSTWWeightsInputTag")),
    pdfNNPDFWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFNNPDFWeightsInputTag")),
    pileupInfoSrc(iConfig.getParameter<edm::InputTag>("pileupInfo"))
{
  produces <unsigned int> ( "ProcessID" );
  produces <double>       ( "PtHat" );
  produces <std::vector<double> > ( "PDFCTEQWeights" );
  produces <std::vector<double> > ( "PDFMSTWWeights" );
  produces <std::vector<double> > ( "PDFNNPDFWeights" );
  produces <std::vector<int> > ( "PileUpInteractions");
  produces <std::vector<int> > ( "PileUpOriginBX" ) ;
  produces <std::vector<float> > ( "PileUpInteractionsTrue" );
  produces <double>       ( "Weight" );
//   produces <double>       ( "AlphaQCD" );
//   produces <double>       ( "AlphaQED" );
}

void RootTupleMakerV2_GenEventInfo::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int >         processID   ( new unsigned int() );
  std::auto_ptr<double >               ptHat ( new double() );
  std::auto_ptr<std::vector<double> >  pdfCTEQWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfMSTWWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfNNPDFWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int >  >   Number_interactions  ( new std::vector<int>() );
  std::auto_ptr<std::vector<float> >   trueNumberInteractions ( new std::vector<float>() );
  std::auto_ptr<std::vector<int >  >   OriginBX( new std::vector<int>() );
  std::auto_ptr<double >               weight ( new double() );
//   std::auto_ptr<double >               alphaQCD ( new double() );
//   std::auto_ptr<double >               alphaQED ( new double() );

  *processID.get() = 0;
  *ptHat.get() = 0.;
  *weight.get() = 0.;
//   *alphaQCD.get() = 0.;
//   *alphaQED.get() = 0.;

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(genEvtInfoInputTag, genEvtInfoProduct);

    if( genEvtInfoProduct.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << genEvtInfoInputTag;

      *processID.get() = genEvtInfoProduct->signalProcessID();
      *ptHat.get() = ( genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : 0. );
      *weight.get() = genEvtInfoProduct->weight();
//       *alphaQCD.get() = genEvtInfoProduct->alphaQCD();
//       *alphaQED.get() = genEvtInfoProduct->alphaQED();

    } else {
      edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << genEvtInfoInputTag;
    }
    // PDF Weights Part
    if( storePDFWeights ) {

      edm::Handle<std::vector<double> > pdfCTEQWeightsHandle;
      edm::Handle<std::vector<double> > pdfMSTWWeightsHandle;
      edm::Handle<std::vector<double> > pdfNNPDFWeightsHandle;

      iEvent.getByLabel(pdfCTEQWeightsInputTag, pdfCTEQWeightsHandle);
      iEvent.getByLabel(pdfMSTWWeightsInputTag, pdfMSTWWeightsHandle);
      iEvent.getByLabel(pdfNNPDFWeightsInputTag, pdfNNPDFWeightsHandle);

      if( pdfCTEQWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << pdfCTEQWeightsInputTag;
        *pdfCTEQWeights.get() = *pdfCTEQWeightsHandle;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << pdfCTEQWeightsInputTag;
      }

      if( pdfMSTWWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << pdfMSTWWeightsInputTag;
        *pdfMSTWWeights.get() = *pdfMSTWWeightsHandle;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << pdfMSTWWeightsInputTag;
      }

      if( pdfNNPDFWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained " << pdfNNPDFWeightsInputTag;

        *pdfNNPDFWeights.get() = *pdfNNPDFWeightsHandle;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the product " << pdfNNPDFWeightsInputTag;
      }

    }
    // PileupSummary Part
    edm::Handle<std::vector<PileupSummaryInfo> >  puInfo;
    iEvent.getByLabel(pileupInfoSrc, puInfo);
    
    if(puInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = puInfo->begin(); it != puInfo->end(); ++it ) {
	trueNumberInteractions -> push_back ( it -> getTrueNumInteractions() );
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
  iEvent.put( pdfCTEQWeights, "PDFCTEQWeights" );
  iEvent.put( pdfMSTWWeights, "PDFMSTWWeights" );
  iEvent.put( pdfNNPDFWeights, "PDFNNPDFWeights" );
  iEvent.put( Number_interactions,   "PileUpInteractions"   );
  iEvent.put( trueNumberInteractions, "PileUpInteractionsTrue" );
  iEvent.put( OriginBX,   "PileUpOriginBX" );
  iEvent.put( weight, "Weight" );
//   iEvent.put( alphaQCD, "alphaQCD" );
//   iEvent.put( alphaQED, "alphaQED" );
}
