#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_GenEventInfo.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <iostream>
#include <string>


RootTupleMakerV2_GenEventInfo::RootTupleMakerV2_GenEventInfo(const edm::ParameterSet& iConfig) :
  genEvtInfoInputToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag"))),
  storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
  pdfCTEQWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFCTEQWeightsInputTag"))),
  pdfMMTHWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFMMTHWeightsInputTag"))),
  pdfNNPDFWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFNNPDFWeightsInputTag"))),
  pdfPDF4LHCWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFPDF4LHCWeightsInputTag"))),
  pileupInfoSrcToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
  LHERunInfoToken_(consumes<LHERunInfoProduct, edm::InRun >(iConfig.getParameter<edm::InputTag>("LHERunInfoProductInputTag"))),
  LHEEventProductToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInputTag")))
{
  produces <unsigned int> ( "ProcessID" );
  produces <double>       ( "PtHat" );
  produces <std::vector<double> > ( "PDFCTEQWeights" );
  produces <std::vector<double> > ( "PDFMMTHWeights" );
  produces <std::vector<double> > ( "PDFNNPDFWeights" );
  produces <std::vector<double> > ( "PDFPDF4LHCWeights" );
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
  iRun.getByToken(LHERunInfoToken_, run );

  if (run.isValid()){
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());

    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      if ( (iter->tag().compare("initrwgt") ==0) ){
	edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") <<"LHE File tag: " << iter->tag() << std::endl;
	std::vector<std::string> lines = iter->lines();
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
	  //edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << lines.at(iLine);
	  std::cout << lines.at(iLine);
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
  std::auto_ptr<std::vector<double> >  pdfMMTHWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfNNPDFWeights  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pdfPDF4LHCWeights  ( new std::vector<double>()  );
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
    iEvent.getByToken(genEvtInfoInputToken_, genEvtInfoProduct);

    if( genEvtInfoProduct.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained genEvtInfoInputToken_";

      *processID.get() = genEvtInfoProduct->signalProcessID();
      *ptHat.get() = ( genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : 0. );
      *weight.get() = genEvtInfoProduct->weight();
//       *alphaQCD.get() = genEvtInfoProduct->alphaQCD();
//       *alphaQED.get() = genEvtInfoProduct->alphaQED();

    } else {
      edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the genEvtInfoInputToken_";
    }
    // PDF Weights Part
    if( storePDFWeights ) {

      edm::Handle<std::vector<double> > pdfCTEQWeightsHandle;
      edm::Handle<std::vector<double> > pdfMMTHWeightsHandle;
      edm::Handle<std::vector<double> > pdfNNPDFWeightsHandle;
      edm::Handle<std::vector<double> > pdfPDF4LHCWeightsHandle;

      iEvent.getByToken(pdfCTEQWeightsInputToken_, pdfCTEQWeightsHandle);
      iEvent.getByToken(pdfMMTHWeightsInputToken_, pdfMMTHWeightsHandle);
      iEvent.getByToken(pdfNNPDFWeightsInputToken_, pdfNNPDFWeightsHandle);
      iEvent.getByToken(pdfPDF4LHCWeightsInputToken_, pdfPDF4LHCWeightsHandle);

      if( pdfCTEQWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained pdfCTEQWeightsInputToken_";
        //*pdfCTEQWeights.get() = *pdfCTEQWeightsHandle;
	//Instead of the above line, have to divide by central value - if 0, set to 1
	std::vector<double> weights = (*pdfCTEQWeightsHandle);
	std::vector<double> reWeightsCTEQ;
	unsigned int nmembers = weights.size();	
	for (unsigned int i=0; i<nmembers; i++) {
	  if(weights[0]!=0)reWeightsCTEQ.push_back(weights[i]/weights[0]);
	  else reWeightsCTEQ.push_back(1.0);
	}
        *pdfCTEQWeights.get() = reWeightsCTEQ;
	//for debugging to see what the weight variations are
	/*
	for (unsigned int j=0; j<nmembers; j+=2) {
	  if(j==0){std::cout << "Event weight for PDF central value" << j << ": " << weights[j] << std::endl;j++;}
	  else{
	    std::cout << "Event weight for PDF variation +" << (j+1)/2 << ": " << weights[j] << std::endl;
	    std::cout << "Event weight for PDF variation -" << (j+1)/2 << ": " << weights[j+1] << std::endl;
	  }
	}
	*/
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the pdfCTEQWeightsInputToken_";
      }

      if( pdfMMTHWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained pdfMMTHWeightsInputToken_";
        //*pdfMMTHWeights.get() = *pdfMMTHWeightsHandle;
	//Instead of the above line, have to divide by central value - if 0, set to 1
	std::vector<double> weights = (*pdfMMTHWeightsHandle);
	std::vector<double> reWeightsMMTH;
	unsigned int nmembers = weights.size();	
	for (unsigned int i=0; i<nmembers; i++) {
	  if(weights[0]!=0)reWeightsMMTH.push_back(weights[i]/weights[0]);
	  else reWeightsMMTH.push_back(1.0);
	}
        *pdfMMTHWeights.get() = reWeightsMMTH;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the pdfMMTHWeightsInputToken_";
      }

      if( pdfNNPDFWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained pdfNNPDFWeightsInputToken_";

        //*pdfNNPDFWeights.get() = *pdfNNPDFWeightsHandle;
	//Instead of the above line, have to divide by central value - if 0, set to 1
	std::vector<double> weights = (*pdfNNPDFWeightsHandle);
	std::vector<double> reWeightsNNPDF;
	unsigned int nmembers = weights.size();	
	for (unsigned int i=0; i<nmembers; i++) {
	  if(weights[0]!=0)reWeightsNNPDF.push_back(weights[i]/weights[0]);
	  else reWeightsNNPDF.push_back(1.0);
	}
        *pdfNNPDFWeights.get() = reWeightsNNPDF;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the pdfNNPDFWeightsInputToken_";
      }
     
      if( pdfPDF4LHCWeightsHandle.isValid() ) {
        edm::LogInfo("RootTupleMakerV2_GenEventInfoInfo") << "Successfully obtained pdfPDF4LHCWeightsInputToken_";

        //*pdfPDF4LHCWeights.get() = *pdfPDF4LHCWeightsHandle;
	//Instead of the above line, have to divide by central value - if 0, set to 1
	std::vector<double> weights = (*pdfPDF4LHCWeightsHandle);
	std::vector<double> reWeightsPDF4LHC;
	unsigned int nmembers = weights.size();	
	for (unsigned int i=0; i<nmembers; i++) {
	  if(weights[0]!=0)reWeightsPDF4LHC.push_back(weights[i]/weights[0]);
	  else reWeightsPDF4LHC.push_back(1.0);
	}
        *pdfPDF4LHCWeights.get() = reWeightsPDF4LHC;
      } else {
        edm::LogError("RootTupleMakerV2_GenEventInfoError") << "Error! Can't get the pdfPDF4LHCWeightsInputToken_";
      }
      
    }
    // PileupSummary Part
    edm::Handle<std::vector<PileupSummaryInfo> >  puInfo;
    iEvent.getByToken(pileupInfoSrcToken_, puInfo);
    
    if(puInfo.isValid()) {
      for( std::vector<PileupSummaryInfo>::const_iterator it = puInfo->begin(); it != puInfo->end(); ++it ) {
	trueNumberInteractions -> push_back ( it -> getTrueNumInteractions() );
	Number_interactions -> push_back ( it->getPU_NumInteractions() );
	OriginBX -> push_back ( it -> getBunchCrossing());
      }
    }
    else {
      edm::LogError("RootTupleMakerV2_PileUpError") << "Error! Can't get the pileupInfoSrcToken_";
    }
    ///Weights from LHE part
    edm::Handle<LHEEventProduct> EvtHandle;
    iEvent.getByToken(LHEEventProductToken_, EvtHandle );

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(genEvtInfoInputToken_, genEvtInfo);

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
	  pdfMMTHWeights->push_back(thisWeight);
	}
      	for (unsigned int i=392; i <= 444; i++) {
	  thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
	  pdfCTEQWeights->push_back(thisWeight);
	}
      	for (unsigned int i=392; i <= 444; i++) {//fixme todo: MG samples don't currently have PDF4LHC - will they be added?
	  thisWeight = theWeight * (EvtHandle->weights()[i].wgt/EvtHandle->originalXWGTUP()); 
	  pdfPDF4LHCWeights->push_back(thisWeight);
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
  iEvent.put( pdfMMTHWeights, "PDFMMTHWeights" );
  iEvent.put( pdfNNPDFWeights, "PDFNNPDFWeights" );
  iEvent.put( pdfPDF4LHCWeights, "PDFPDF4LHCWeights" );
  iEvent.put( scaleWeights, "ScaleWeights" );
  iEvent.put( amcNLOweight, "amcNLOWeight" );
  iEvent.put( Number_interactions,   "PileUpInteractions"   );
  iEvent.put( trueNumberInteractions, "PileUpInteractionsTrue" );
  iEvent.put( OriginBX,   "PileUpOriginBX" );
  iEvent.put( weight, "Weight" );
//   iEvent.put( alphaQCD, "alphaQCD" );
//   iEvent.put( alphaQED, "alphaQED" );
}
