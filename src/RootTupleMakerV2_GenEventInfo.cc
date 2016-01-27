#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include <iostream>
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include <string>


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
  produces <std::vector<double> > ( "ScaleWeights" );
  produces <double> ( "amcNLOWeight" );
  produces <std::vector<int> > ( "PileUpInteractions");
  produces <std::vector<int> > ( "PileUpOriginBX" ) ;
  produces <std::vector<float> > ( "PileUpInteractionsTrue" );
  produces <double>       ( "Weight" );
//   produces <double>       ( "AlphaQCD" );
//   produces <double>       ( "AlphaQED" );
}

void RootTupleMakerV2_GenEventInfo::
endRun(edm::Run const& iRun, edm::EventSetup const&) {
  
  //Logging the LHE event info to check, as per https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights
  //The assumption will be made that it will be in this format:
  // <weightgroup combine="envelope" type="Central scale variation">
  //   <weight id="1"> mur=1 muf=1 </weight>
  //   <weight id="2"> mur=1 muf=2 </weight>
  //   <weight id="3"> mur=1 muf=0.5 </weight>
  //   <weight id="4"> mur=2 muf=1 </weight>
  //   <weight id="5"> mur=2 muf=2 </weight>
  //   <weight id="6"> mur=2 muf=0.5 </weight>
  //   <weight id="7"> mur=0.5 muf=1 </weight>
  //   <weight id="8"> mur=0.5 muf=2 </weight>
  //   <weight id="9"> mur=0.5 muf=0.5 </weight>
  // </weightgroup>
  // Thus meaning the 0'th weight is the central weight - this was checked for
  // /TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1/MINIAODSIM
  // /DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM 
  // /WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM

  edm::Handle<LHERunInfoProduct> run; 
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator; 
  iRun.getByLabel( "externalLHEProducer", run );

  if (run.isValid()){
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());

    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      if ( (iter->tag().compare("initrwgt") ==0) ){
	edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") <<"LHE File tag: " << iter->tag() << std::endl;
	std::vector<std::string> lines = iter->lines();
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
	  edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << lines.at(iLine);
	}
      }
    }
  }
}
void RootTupleMakerV2_GenEventInfo::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<unsigned int >         processID   ( new unsigned int() );
  std::auto_ptr<double >               ptHat ( new double() );
  std::auto_ptr<std::vector<double> >  pdfCTEQWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfMSTWWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfNNPDFWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scaleWeights  ( new std::vector<double>()  );
  std::auto_ptr<double>  amcNLOweight  ( new double()  );
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
	/*//for debugging to see what the weight variations are
	std::vector<double> weights = (*pdfCTEQWeightsHandle);
	std::cout << "Event weight for central PDF:" << weights[0] << std::endl;
	unsigned int nmembers = weights.size();
	for (unsigned int j=1; j<nmembers; j+=2) {
	  std::cout << "Event weight for PDF variation +" << (j+1)/2 << ": " << weights[j] << std::endl;
	  std::cout << "Event weight for PDF variation -" << (j+1)/2 << ": " << weights[j+1] << std::endl;
	}
	*/
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
    ///Weights from LHE part
    edm::Handle<LHEEventProduct> EvtHandle;
    iEvent.getByLabel( "externalLHEProducer" , EvtHandle );

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel(genEvtInfoInputTag, genEvtInfo);

    //Non-madgraph samples may not have this information, if so, skip
    if (EvtHandle.isValid() && genEvtInfo.isValid()){
      //Powheg samples have valid EvtHandle but seg fault when trying to access weights, so skip if the weights vector is empty
      if (EvtHandle->weights().size()>0){

      	double theWeight = genEvtInfo->weight();
      	double thisWeight = -1.;
	
      	//This follows the suggestion here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights
      	//WARNING: This assumes the first 9 weights are renormalization/factorization-related  This seems to be true for all Madgraph samples
      	for (unsigned int i=0; i <= 8; i++) {
      		thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
		scaleWeights->push_back(thisWeight);
      	}
	//This is for MG PDF weights.  PYTHIA doesn't have externalLHEProducer, need to get them in a different way
      	for (unsigned int i=9; i <= 109; i++) {
	  thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
	  pdfNNPDFWeights->push_back(thisWeight);
	}
      	for (unsigned int i=315; i <= 365; i++) {
	  thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
	  pdfMSTWWeights->push_back(thisWeight);
	}
      	for (unsigned int i=392; i <= 444; i++) {
	  thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
	  pdfCTEQWeights->push_back(thisWeight);
	}
	
      	EvtHandle->weights()[0].wgt < 0 ? *amcNLOweight.get()=-1. : *amcNLOweight.get()=1.;
      
      	/*
		unsigned int num_whichWeight = EvtHandle->weights().size();
		for (unsigned int iWeight = 0; iWeight < num_whichWeight; iWeight++) {
		amcNLOWeights->push_back( EvtHandle->weights()[iWeight].wgt/EvtHandle->originalXWGTUP() ); 
		if(iWeight==0){
		std::cout << "           weightLHE[" << iWeight << "] = " << EvtHandle->weights()[iWeight].wgt << std::endl;
		std::cout << "        newweightLHE[" << iWeight << "] = " << EvtHandle->weights()[iWeight].wgt/EvtHandle->originalXWGTUP() << std::endl;
		std::cout << " weight*newweightLHE[" << iWeight << "] = " << theWeight*EvtHandle->weights()[iWeight].wgt/EvtHandle->originalXWGTUP() << std::endl;
		int i = iWeight;
		if(iWeight<=8)std::cout << " scaleWeight[" << iWeight << "] = " << theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()) << std::endl<<std::endl<<std::endl;
		}
		}
      	*/
      }
    }
  }
  
  
  //-----------------------------------------------------------------
  iEvent.put( processID, "ProcessID" );
  iEvent.put( ptHat, "PtHat" );
  iEvent.put( pdfCTEQWeights, "PDFCTEQWeights" );
  iEvent.put( pdfMSTWWeights, "PDFMSTWWeights" );
  iEvent.put( pdfNNPDFWeights, "PDFNNPDFWeights" );
  iEvent.put( scaleWeights, "ScaleWeights" );
  iEvent.put( amcNLOweight, "amcNLOWeight" );
  iEvent.put( Number_interactions,   "PileUpInteractions"   );
  iEvent.put( trueNumberInteractions, "PileUpInteractionsTrue" );
  iEvent.put( OriginBX,   "PileUpOriginBX" );
  iEvent.put( weight, "Weight" );
//   iEvent.put( alphaQCD, "alphaQCD" );
//   iEvent.put( alphaQED, "alphaQED" );
}
