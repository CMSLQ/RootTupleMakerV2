#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

RootTupleMakerV2_MET::RootTupleMakerV2_MET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    metToken(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("InputTag"))),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET")),
    store_MET_significance (iConfig.getParameter<bool>  ("StoreMETSignificance")),
    uncertainty (iConfig.getParameter<std::string>  ("Uncertainty")),
    corLevel       (iConfig.getParameter<std::string>  ("CorrectionLevel"))
{
  produces <std::vector<float> > ( prefix + "MET" + suffix );
  produces <std::vector<float> > ( prefix + "METPhi" + suffix );
  produces <std::vector<float> > ( prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    produces <std::vector<float> > ( prefix + "METUncorr" + suffix );
    produces <std::vector<float> > ( prefix + "METPhiUncorr" + suffix );
    produces <std::vector<float> > ( prefix + "SumETUncorr" + suffix );
  }
  if ( store_MET_significance ) {
    produces <std::vector<float> > ( prefix + "METSig" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDXX" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDXY" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDYX" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDYY" + suffix );
  }
}

void RootTupleMakerV2_MET::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::unique_ptr<std::vector<float> >  met  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  metphi  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  sumet  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  metuncorr  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  metphiuncorr  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  sumetuncorr  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  metsig  ( new std::vector<float>()  );  
  std::unique_ptr<std::vector<float> >  metsigmatrixdxx  ( new std::vector<float>()  );  
  std::unique_ptr<std::vector<float> >  metsigmatrixdxy  ( new std::vector<float>()  );  
  std::unique_ptr<std::vector<float> >  metsigmatrixdyx  ( new std::vector<float>()  );  
  std::unique_ptr<std::vector<float> >  metsigmatrixdyy  ( new std::vector<float>()  );  

  //-----------------------------------------------------------------
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken, mets);

  if(mets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_METInfo") << "Total # METs: " << mets->size();

    for( pat::METCollection::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

      pat::MET::METUncertainty shift = pat::MET::NoShift;
      pat::MET::METCorrectionLevel level = pat::MET::Type1;

      if(     uncertainty=="NoShift")          shift=pat::MET::NoShift;
      else if(uncertainty=="JetResUp")         shift=pat::MET::JetResUp;
      else if(uncertainty=="JetResDown")       shift=pat::MET::JetResDown;
      else if(uncertainty=="JetEnUp")          shift=pat::MET::JetEnUp;
      else if(uncertainty=="JetEnDown")        shift=pat::MET::JetEnDown;
      else if(uncertainty=="MuonEnUp")         shift=pat::MET::MuonEnUp;
      else if(uncertainty=="MuonEnDown")       shift=pat::MET::MuonEnDown;
      else if(uncertainty=="ElectronEnUp")     shift=pat::MET::ElectronEnUp;
      else if(uncertainty=="ElectronEnDown")   shift=pat::MET::ElectronEnDown;
      else if(uncertainty=="TauEnUp")          shift=pat::MET::TauEnUp;
      else if(uncertainty=="TauEnDown")        shift=pat::MET::TauEnDown;
      else if(uncertainty=="UnclusteredEnUp")  shift=pat::MET::UnclusteredEnUp;
      else if(uncertainty=="UnclusteredEnDown")shift=pat::MET::UnclusteredEnDown;
      else if(uncertainty=="PhotonEnUp")       shift=pat::MET::PhotonEnUp;
      else if(uncertainty=="PhotonEnDown")     shift=pat::MET::PhotonEnDown;
      else if(uncertainty=="JetResUpSmear")    shift=pat::MET::JetResUpSmear;
      else if(uncertainty=="JetResDownSmear")  shift=pat::MET::JetResDownSmear;
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET uncertainty label: " << uncertainty;

      bool applyCorrection = true;
      if(corLevel=="Raw")                level = pat::MET::Raw;
      else if(corLevel=="Type1")         level = pat::MET::Type1;
      else if(corLevel=="Type01")        level = pat::MET::Type01;
      else if(corLevel=="TypeXY")        level = pat::MET::TypeXY;
      else if(corLevel=="Type1XY")       level = pat::MET::Type1XY;
      else if(corLevel=="Type01XY")      level = pat::MET::Type01XY;//dont wan't to use 01XY because it overcorrects
      else if(corLevel=="Type1Smear")    level = pat::MET::Type1Smear;
      else if(corLevel=="Type01Smear")   level = pat::MET::Type01Smear;
      else if(corLevel=="Type1SmearXY")  level = pat::MET::Type1SmearXY;
      else if(corLevel=="Type01SmearXY") level = pat::MET::Type01SmearXY;
      else if(corLevel=="RawCalo")       level = pat::MET::RawCalo;
      // this is our adhoc way of not running the shiftedPt, shiftedPhi, shiftedSumEt on the recorrected MET collections
      //  (which we don't want to do because they aren't filled, so accessing them causes a seg fault)
      else if(corLevel=="NoCorrection")  applyCorrection = false;
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET correction level label: " << corLevel;

      if(applyCorrection)
	{
	  met->push_back( it->shiftedPt(shift,level) );
	  metphi->push_back( it->shiftedPhi(shift,level) );
	  sumet->push_back( it->shiftedSumEt(shift,level) );
	}
      else
	{
	  met->push_back( it->pt() );
	  metphi->push_back( it->phi() );
	  sumet->push_back( it->sumEt() );
	}
      
      if ( store_uncorrected_MET ) {
        metuncorr->push_back( it->uncorPt() );
        metphiuncorr->push_back( it->uncorPhi() );
        sumetuncorr->push_back( it->uncorSumEt() );
      }

      if ( store_MET_significance ) {
        //-- See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance#Known_Issues
	float sigmaX2= it->getSignificanceMatrix()(0,0);
	float sigmaY2= it->getSignificanceMatrix()(1,1);
	float significance = -1;
	if(sigmaX2<1.e10 && sigmaY2<1.e10) 
	  significance = it->metSignificance();
	//--

	metsig->push_back( significance );
	metsigmatrixdxx->push_back( it->getSignificanceMatrix()(0,0) );
	metsigmatrixdxy->push_back( it->getSignificanceMatrix()(0,1) );
	metsigmatrixdyx->push_back( it->getSignificanceMatrix()(1,0) );
	metsigmatrixdyy->push_back( it->getSignificanceMatrix()(1,1) );
	//See DataFormats/METReco/src/MET.cc
      }

    }
  } else {
    edm::LogError("RootTupleMakerV2_METError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put(std::move(met), prefix + "MET" + suffix );
  iEvent.put(std::move(metphi), prefix + "METPhi" + suffix );
  iEvent.put(std::move(sumet), prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    iEvent.put(std::move(metuncorr), prefix + "METUncorr" + suffix );
    iEvent.put(std::move(metphiuncorr), prefix + "METPhiUncorr" + suffix );
    iEvent.put(std::move(sumetuncorr), prefix + "SumETUncorr" + suffix );
  }
  if ( store_MET_significance ) {
    iEvent.put(std::move(metsig), prefix + "METSig" + suffix );
    iEvent.put(std::move(metsigmatrixdxx), prefix + "METSigMatrixDXX" + suffix );
    iEvent.put(std::move(metsigmatrixdxy), prefix + "METSigMatrixDXY" + suffix );
    iEvent.put(std::move(metsigmatrixdyx), prefix + "METSigMatrixDYX" + suffix );
    iEvent.put(std::move(metsigmatrixdyy), prefix + "METSigMatrixDYY" + suffix );
  }
}
