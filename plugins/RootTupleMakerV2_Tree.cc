#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Tree.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "DataFormats/Provenance/interface/SelectedProducts.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include "boost/foreach.hpp"
#include <TBranch.h>
#include <TLorentzVector.h>
	
RootTupleMakerV2_Tree::~RootTupleMakerV2_Tree() {}

void RootTupleMakerV2_Tree::
analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  BOOST_FOREACH( BranchConnector* connector, connectors)
    connector->connect(iEvent);
  tree->Fill();
}

template <class T>
void RootTupleMakerV2_Tree::TypedBranchConnector<T>::
connect(const edm::Event& iEvent) {
  edm::Handle<T> handle_;
  //iEvent.getByLabel(ml, pin, handle_);
  iEvent.getByToken(token_,handle_);
  object_ = *handle_;
}

template <class T> 
RootTupleMakerV2_Tree::TypedBranchConnector<T>::
TypedBranchConnector(edm::BranchDescription const* desc,
                     std::string t,
                     TTree * tree,
                     edm::EDGetTokenT<T> token)
  :  ml( desc->moduleLabel() ),
     pin( desc->productInstanceName() ),
     token_(token)
{
  object_ptr_ = &object_;
  std::string s=pin+t;
  if(t!="")  { tree->Branch(pin.c_str(),  object_ptr_, s.c_str() );}  //raw type
  else       { tree->Branch(pin.c_str(), &object_ptr_            );}  //vector<type>
}

void RootTupleMakerV2_Tree::beginJob() {
  //usesResource("TFileService");
  tree = fs->make<TTree>("tree", "");

  typedef std::map<std::string,       bool> mapStringBool;
  typedef std::map<std::string,        int> mapStringInt;
  typedef std::map<std::string,std::string> mapStringString;
  typedef std::map<std::string,std::vector<float> > mapStringDoubles;
  typedef std::vector<std::vector<float> > vectorVectorFloats;
  typedef std::vector<std::vector<int> > vectorVectorInts;
  typedef std::vector<std::vector<bool> > vectorVectorBools;
  typedef std::vector<std::vector<std::string> > vectorVectorStrings;

  std::map<std::string, LEAFTYPE> leafmap;
  leafmap["bool"]      = BOOL;       leafmap["bools"]     = BOOL_V;
  leafmap["short int"] = SHORT;      leafmap["shorts"]    = SHORT_V;
  leafmap["ushort int"]= U_SHORT;    leafmap["ushorts"]   = U_SHORT_V;
  leafmap["int"]       = INT;        leafmap["ints"]      = INT_V;
  leafmap["uint"]      = U_INT;      leafmap["uints"]     = U_INT_V;
  leafmap["float"]     = FLOAT;      leafmap["floats"]    = FLOAT_V;
  leafmap["double"]    = DOUBLE;     leafmap["doubles"]   = DOUBLE_V;
  leafmap["lint"]      = LONG;       leafmap["longs"]     = LONG_V;
  leafmap["ulint"]     = U_LONG;     leafmap["ulongs"]    = U_LONG_V;
  leafmap["StringStringstdmap"] = STRING_STRING_M;
  leafmap["Stringboolstdmap"  ] = STRING_BOOL_M;
  leafmap["Stringintstdmap"   ] = STRING_INT_M;
  leafmap["Stringfloatsstdmap"] = STRING_FLOAT_V_M;
  leafmap["floatss"] = FLOAT_V_V;
  leafmap["intss"] = INT_V_V;
  leafmap["boolss"] = BOOL_V_V;
  leafmap["Stringss"] = STRING_V_V;
  // leafmap[""] = LORENTZ_V_V;
  
  //
  leafmap["String"]     = STRING;     leafmap["Strings"]    = STRING_V;

  edm::Service<edm::ConstProductRegistry> reg;
  edm::SelectedProducts allBranches = reg->allBranchDescriptions();
  //DEBUG
  //edm::LogWarning("RootTupleMakerV2_Tree") << cfg;
  //DEBUG
  edm::ProductSelectorRules productSelectorRules_(cfg, "outputCommands", "RootTupleMakerV2_Tree");
  edm::ProductSelector productSelector_;
  productSelector_.initialize(productSelectorRules_, allBranches);

  //DEBUG
  edm::LogWarning("RootTupleMakerV2_Tree") << productSelector_;
  //DEBUG
  std::set<std::string> branchnames;

  BOOST_FOREACH( const edm::SelectedProducts::value_type& selection, allBranches) {
    //DEBUG
    edm::LogWarning("RootTupleMakerV2_Tree") << selection->productInstanceName();
    //DEBUG
    if(productSelector_.selected(*selection)) {

      //Check for duplicate branch names
      if (branchnames.find( selection->productInstanceName()) != branchnames.end() ) {
        throw edm::Exception(edm::errors::Configuration)
          << "More than one branch named: "
          << selection->productInstanceName() << std::endl
          << "Exception thrown from RootTupleMakerV2_Tree::beginJob" << std::endl;
      }
      else {
        branchnames.insert( selection->productInstanceName() );
      }
      

      //Create RootTupleMakerV2_Tree branch

      std::string tokenName=selection->moduleLabel()+"_"+selection->productInstanceName();
      
      switch(leafmap.find( selection->friendlyClassName() )->second) {
      case BOOL            : { connectors.push_back( new TypedBranchConnector                      <bool>         (selection, "/O", tree, tokenMap.find(tokenName)->second) ); break; }
      case BOOL_V          : { connectors.push_back( new TypedBranchConnector<std::vector          <bool> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case INT             : { connectors.push_back( new TypedBranchConnector                       <int>         (selection, "/I", tree, tokenMap.find(tokenName)->second) ); break; }
      case INT_V           : { connectors.push_back( new TypedBranchConnector<std::vector           <int> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_INT           : { connectors.push_back( new TypedBranchConnector              <unsigned int>         (selection, "/i", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_INT_V         : { connectors.push_back( new TypedBranchConnector<std::vector  <unsigned int> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case SHORT           : { connectors.push_back( new TypedBranchConnector                     <short>         (selection, "/S", tree, tokenMap.find(tokenName)->second) ); break; }
      case SHORT_V         : { connectors.push_back( new TypedBranchConnector<std::vector         <short> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_SHORT         : { connectors.push_back( new TypedBranchConnector            <unsigned short>         (selection, "/s", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_SHORT_V       : { connectors.push_back( new TypedBranchConnector<std::vector<unsigned short> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case FLOAT           : { connectors.push_back( new TypedBranchConnector                     <float>         (selection, "/F", tree, tokenMap.find(tokenName)->second) ); break; }
      case FLOAT_V         : { connectors.push_back( new TypedBranchConnector<std::vector         <float> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case DOUBLE          : { connectors.push_back( new TypedBranchConnector                    <double>         (selection, "/D", tree, tokenMap.find(tokenName)->second) ); break; }
      case DOUBLE_V        : { connectors.push_back( new TypedBranchConnector<std::vector        <double> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case LONG            : { connectors.push_back( new TypedBranchConnector                      <long>         (selection, "/L", tree, tokenMap.find(tokenName)->second) ); break; }
      case LONG_V          : { connectors.push_back( new TypedBranchConnector<std::vector          <long> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_LONG          : { connectors.push_back( new TypedBranchConnector             <unsigned long>         (selection, "/l", tree, tokenMap.find(tokenName)->second) ); break; }
      case U_LONG_V        : { connectors.push_back( new TypedBranchConnector<std::vector <unsigned long> >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
	//                     {       }
      case STRING          : { connectors.push_back( new TypedBranchConnector             <std::string  >         (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case STRING_V        : { connectors.push_back( new TypedBranchConnector<std::vector <std::string  > >       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
	{       }
      case STRING_INT_M    : { connectors.push_back( new TypedBranchConnector<mapStringInt>        (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case STRING_BOOL_M   : { connectors.push_back( new TypedBranchConnector<mapStringBool>       (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case STRING_STRING_M : { connectors.push_back( new TypedBranchConnector<mapStringString>     (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case STRING_FLOAT_V_M: { connectors.push_back( new TypedBranchConnector<mapStringDoubles>    (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case FLOAT_V_V       : { connectors.push_back( new TypedBranchConnector<vectorVectorFloats>  (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case INT_V_V         : { connectors.push_back( new TypedBranchConnector<vectorVectorInts>    (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case BOOL_V_V        : { connectors.push_back( new TypedBranchConnector<vectorVectorBools>   (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }
      case STRING_V_V      : { connectors.push_back( new TypedBranchConnector<vectorVectorStrings> (selection,   "", tree, tokenMap.find(tokenName)->second) ); break; }

      default:
        {
          std::string leafstring = "";
          typedef std::pair<std::string, LEAFTYPE> pair_t;
          BOOST_FOREACH( const pair_t& leaf, leafmap)
            leafstring+= "\t" + leaf.first + "\n";

          throw edm::Exception(edm::errors::Configuration)
            << "class RootTupleMakerV2_Tree does not handle leaves of type " << selection->className() << " like\n"
            <<   selection->friendlyClassName()   << "_"
            <<   selection->moduleLabel()         << "_"
            <<   selection->productInstanceName() << "_"
            <<   selection->processName()         << std::endl
            << "Valid leaf types are (friendlyClassName):\n"
            <<   leafstring
            << "Exception thrown from RootTupleMakerV2_Tree::beginJob\n";
        }
      }
    }
  }}


/*
  void
  RootTupleMakerV2_Tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  }
*/



RootTupleMakerV2_Tree::RootTupleMakerV2_Tree(const edm::ParameterSet& pset) : 
  cfg(pset)
{
  tokenMap["rootTupleElectrons_ElectronCutFlowHashesEGammaIDHEEP"]=rootTupleElectrons_ElectronCutFlowHashesEGammaIDHEEP_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowHashesEGammaIDLoose"]=rootTupleElectrons_ElectronCutFlowHashesEGammaIDLoose_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowHashesEGammaIDMedium"]=rootTupleElectrons_ElectronCutFlowHashesEGammaIDMedium_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowHashesEGammaIDTight"]=rootTupleElectrons_ElectronCutFlowHashesEGammaIDTight_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowHashesEGammaIDVeto"]=rootTupleElectrons_ElectronCutFlowHashesEGammaIDVeto_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowNamesEGammaIDHEEP"]=rootTupleElectrons_ElectronCutFlowNamesEGammaIDHEEP_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowNamesEGammaIDLoose"]=rootTupleElectrons_ElectronCutFlowNamesEGammaIDLoose_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowNamesEGammaIDMedium"]=rootTupleElectrons_ElectronCutFlowNamesEGammaIDMedium_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowNamesEGammaIDTight"]=rootTupleElectrons_ElectronCutFlowNamesEGammaIDTight_Token_;
  tokenMap["rootTupleElectrons_ElectronCutFlowNamesEGammaIDVeto"]=rootTupleElectrons_ElectronCutFlowNamesEGammaIDVeto_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjCollectionName"]=rootTupleTriggerObjects_HLTriggerObjCollectionName_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjFilterNames"]=rootTupleTriggerObjects_HLTriggerObjFilterNames_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjPathNames"]=rootTupleTriggerObjects_HLTriggerObjPathNames_Token_;
  tokenMap["rootTupleEventSelection_isBPTX0"]=rootTupleEventSelection_isBPTX0_Token_;
  tokenMap["rootTupleEventSelection_isBSCBeamHalo"]=rootTupleEventSelection_isBSCBeamHalo_Token_;
  tokenMap["rootTupleEventSelection_isBSCMinBias"]=rootTupleEventSelection_isBSCMinBias_Token_;
  tokenMap["rootTupleEventSelection_isPhysDeclared"]=rootTupleEventSelection_isPhysDeclared_Token_;
  tokenMap["rootTupleEventSelection_isPrimaryVertex"]=rootTupleEventSelection_isPrimaryVertex_Token_;
  tokenMap["rootTupleEventSelection_isTrackingFailure"]=rootTupleEventSelection_isTrackingFailure_Token_;
  tokenMap["rootTupleEventSelection_passBadEESupercrystalFilter"]=rootTupleEventSelection_passBadEESupercrystalFilter_Token_;
  tokenMap["rootTupleEventSelection_passBeamHaloFilterTight"]=rootTupleEventSelection_passBeamHaloFilterTight_Token_;
  tokenMap["rootTupleEventSelection_passEcalDeadCellTriggerPrimitiveFilter"]=rootTupleEventSelection_passEcalDeadCellTriggerPrimitiveFilter_Token_;
  tokenMap["rootTupleEventSelection_passEcalLaserCorrFilter"]=rootTupleEventSelection_passEcalLaserCorrFilter_Token_;
  tokenMap["rootTupleEventSelection_passHBHENoiseFilter"]=rootTupleEventSelection_passHBHENoiseFilter_Token_;
  tokenMap["rootTupleEventSelection_passHBHENoiseIsoFilter"]=rootTupleEventSelection_passHBHENoiseIsoFilter_Token_;
  tokenMap["rootTupleEventSelection_passHcalLaserEventFilter"]=rootTupleEventSelection_passHcalLaserEventFilter_Token_;
  tokenMap["rootTupleEventSelection_passLogErrorTooManyClusters"]=rootTupleEventSelection_passLogErrorTooManyClusters_Token_;
  tokenMap["rootTupleEventSelection_passManyStripClus53X"]=rootTupleEventSelection_passManyStripClus53X_Token_;
  tokenMap["rootTupleEventSelection_passTooManyStripClus53X"]=rootTupleEventSelection_passTooManyStripClus53X_Token_;
  tokenMap["rootTupleEventSelection_passTrackingFailureFilter"]=rootTupleEventSelection_passTrackingFailureFilter_Token_;
  tokenMap["rootTupleEvent_isData"]=rootTupleEvent_isData_Token_;
  tokenMap["rootTupleMuons_hasVeryForwardPFMuon"]=rootTupleMuons_hasVeryForwardPFMuon_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJethasJetWithBadUncAK4CHS"]=rootTuplePFJetsAK4CHS_PFJethasJetWithBadUncAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJethasJetWithBadUncAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJethasJetWithBadUncAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJethasJetWithBadUncAK5CHS"]=rootTuplePFJetsAK5CHS_PFJethasJetWithBadUncAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJethasJetWithBadUncAK5"]=rootTuplePFJetsAK5_PFJethasJetWithBadUncAK5_Token_;
  tokenMap["rootTupleElectrons_ElectronGsfCtfCharge"]=rootTupleElectrons_ElectronGsfCtfCharge_Token_;
  tokenMap["rootTupleElectrons_ElectronGsfCtfScPixCharge"]=rootTupleElectrons_ElectronGsfCtfScPixCharge_Token_;
  tokenMap["rootTupleElectrons_ElectronGsfScPixCharge"]=rootTupleElectrons_ElectronGsfScPixCharge_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTDoubleEleMatched"]=rootTupleElectrons_ElectronHLTDoubleEleMatched_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTEleJetJetMatched"]=rootTupleElectrons_ElectronHLTEleJetJetMatched_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleMatched"]=rootTupleElectrons_ElectronHLTSingleEleMatched_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleWP85Matched"]=rootTupleElectrons_ElectronHLTSingleEleWP85Matched_Token_;
  tokenMap["rootTupleElectrons_ElectronHasEcalDrivenSeed"]=rootTupleElectrons_ElectronHasEcalDrivenSeed_Token_;
  tokenMap["rootTupleElectrons_ElectronHasMatchedConvPhot"]=rootTupleElectrons_ElectronHasMatchedConvPhot_Token_;
  tokenMap["rootTupleElectrons_ElectronHasTrackerDrivenSeed"]=rootTupleElectrons_ElectronHasTrackerDrivenSeed_Token_;
  tokenMap["rootTupleElectrons_ElectronIsEB"]=rootTupleElectrons_ElectronIsEB_Token_;
  tokenMap["rootTupleElectrons_ElectronIsEE"]=rootTupleElectrons_ElectronIsEE_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDEoP"]=rootTupleElectrons_ElectronPassEGammaIDEoP_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDLoose"]=rootTupleElectrons_ElectronPassEGammaIDLoose_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDMedium"]=rootTupleElectrons_ElectronPassEGammaIDMedium_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDTight"]=rootTupleElectrons_ElectronPassEGammaIDTight_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDTrigTight"]=rootTupleElectrons_ElectronPassEGammaIDTrigTight_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDTrigWP70"]=rootTupleElectrons_ElectronPassEGammaIDTrigWP70_Token_;
  tokenMap["rootTupleElectrons_ElectronPassEGammaIDVeto"]=rootTupleElectrons_ElectronPassEGammaIDVeto_Token_;
  tokenMap["rootTupleElectrons_ElectronPassHEEPID"]=rootTupleElectrons_ElectronPassHEEPID_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleIsoMuonMatched"]=rootTupleMuons_MuonHLTSingleIsoMuonMatched_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleMuonMatched"]=rootTupleMuons_MuonHLTSingleMuonMatched_Token_;
  tokenMap["rootTupleVertex_VertexIsFake"]=rootTupleVertex_VertexIsFake_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjPassedPathL3Filter"]=rootTupleTriggerObjects_HLTriggerObjPassedPathL3Filter_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjPassedPathLastFilter"]=rootTupleTriggerObjects_HLTriggerObjPassedPathLastFilter_Token_;
  tokenMap["rootTupleEvent_fixedGridRhoAll"]=rootTupleEvent_fixedGridRhoAll_Token_;
  tokenMap["rootTupleEvent_fixedGridRhoFastjetAllCalo"]=rootTupleEvent_fixedGridRhoFastjetAllCalo_Token_;
  tokenMap["rootTupleEvent_fixedGridRhoFastjetCentralCalo"]=rootTupleEvent_fixedGridRhoFastjetCentralCalo_Token_;
  tokenMap["rootTupleEvent_fixedGridRhoFastjetCentralChargedPileUp"]=rootTupleEvent_fixedGridRhoFastjetCentralChargedPileUp_Token_;
  tokenMap["rootTupleEvent_fixedGridRhoFastjetCentralNeutral"]=rootTupleEvent_fixedGridRhoFastjetCentralNeutral_Token_;
  tokenMap["rootTupleEvent_time"]=rootTupleEvent_time_Token_;
  tokenMap["rootTupleGenEventInfo_PtHat"]=rootTupleGenEventInfo_PtHat_Token_;
  tokenMap["rootTupleGenEventInfo_Weight"]=rootTupleGenEventInfo_Weight_Token_;
  tokenMap["rootTupleGenEventInfo_amcNLOWeight"]=rootTupleGenEventInfo_amcNLOWeight_Token_;
  tokenMap["rootTupleElectrons_ElectronBeamSpotDXYError"]=rootTupleElectrons_ElectronBeamSpotDXYError_Token_;
  tokenMap["rootTupleElectrons_ElectronBeamSpotDXY"]=rootTupleElectrons_ElectronBeamSpotDXY_Token_;
  tokenMap["rootTupleElectrons_ElectronCaloEnergy"]=rootTupleElectrons_ElectronCaloEnergy_Token_;
  tokenMap["rootTupleElectrons_ElectronDCotTheta"]=rootTupleElectrons_ElectronDCotTheta_Token_;
  tokenMap["rootTupleElectrons_ElectronDeltaEtaTrkSC"]=rootTupleElectrons_ElectronDeltaEtaTrkSC_Token_;
  tokenMap["rootTupleElectrons_ElectronDeltaEtaTrkSeedSC"]=rootTupleElectrons_ElectronDeltaEtaTrkSeedSC_Token_;
  tokenMap["rootTupleElectrons_ElectronDeltaPhiTrkSC"]=rootTupleElectrons_ElectronDeltaPhiTrkSC_Token_;
  tokenMap["rootTupleElectrons_ElectronDist"]=rootTupleElectrons_ElectronDist_Token_;
  tokenMap["rootTupleElectrons_ElectronE1x5OverE5x5"]=rootTupleElectrons_ElectronE1x5OverE5x5_Token_;
  tokenMap["rootTupleElectrons_ElectronE2x5OverE5x5"]=rootTupleElectrons_ElectronE2x5OverE5x5_Token_;
  tokenMap["rootTupleElectrons_ElectronESuperClusterOverP"]=rootTupleElectrons_ElectronESuperClusterOverP_Token_;
  tokenMap["rootTupleElectrons_ElectronEcalEnergy"]=rootTupleElectrons_ElectronEcalEnergy_Token_;
  tokenMap["rootTupleElectrons_ElectronEcalIsoDR03"]=rootTupleElectrons_ElectronEcalIsoDR03_Token_;
  tokenMap["rootTupleElectrons_ElectronEcalIsoPAT"]=rootTupleElectrons_ElectronEcalIsoPAT_Token_;
  tokenMap["rootTupleElectrons_ElectronEnergy"]=rootTupleElectrons_ElectronEnergy_Token_;
  tokenMap["rootTupleElectrons_ElectronEta"]=rootTupleElectrons_ElectronEta_Token_;
  tokenMap["rootTupleElectrons_ElectronFbrem"]=rootTupleElectrons_ElectronFbrem_Token_;
  tokenMap["rootTupleElectrons_ElectronFull5x5E1x5OverE5x5"]=rootTupleElectrons_ElectronFull5x5E1x5OverE5x5_Token_;
  tokenMap["rootTupleElectrons_ElectronFull5x5E2x5OverE5x5"]=rootTupleElectrons_ElectronFull5x5E2x5OverE5x5_Token_;
  tokenMap["rootTupleElectrons_ElectronFull5x5SigmaIEtaIEta"]=rootTupleElectrons_ElectronFull5x5SigmaIEtaIEta_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTDoubleEleMatchEta"]=rootTupleElectrons_ElectronHLTDoubleEleMatchEta_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTDoubleEleMatchPhi"]=rootTupleElectrons_ElectronHLTDoubleEleMatchPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTDoubleEleMatchPt"]=rootTupleElectrons_ElectronHLTDoubleEleMatchPt_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTEleJetJetMatchEta"]=rootTupleElectrons_ElectronHLTEleJetJetMatchEta_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTEleJetJetMatchPhi"]=rootTupleElectrons_ElectronHLTEleJetJetMatchPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTEleJetJetMatchPt"]=rootTupleElectrons_ElectronHLTEleJetJetMatchPt_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleMatchEta"]=rootTupleElectrons_ElectronHLTSingleEleMatchEta_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleMatchPhi"]=rootTupleElectrons_ElectronHLTSingleEleMatchPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleMatchPt"]=rootTupleElectrons_ElectronHLTSingleEleMatchPt_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleWP85MatchEta"]=rootTupleElectrons_ElectronHLTSingleEleWP85MatchEta_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleWP85MatchPhi"]=rootTupleElectrons_ElectronHLTSingleEleWP85MatchPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronHLTSingleEleWP85MatchPt"]=rootTupleElectrons_ElectronHLTSingleEleWP85MatchPt_Token_;
  tokenMap["rootTupleElectrons_ElectronHcalIsoD1DR03"]=rootTupleElectrons_ElectronHcalIsoD1DR03_Token_;
  tokenMap["rootTupleElectrons_ElectronHcalIsoD2DR03"]=rootTupleElectrons_ElectronHcalIsoD2DR03_Token_;
  tokenMap["rootTupleElectrons_ElectronHcalIsoDR03FullCone"]=rootTupleElectrons_ElectronHcalIsoDR03FullCone_Token_;
  tokenMap["rootTupleElectrons_ElectronHcalIsoDR03"]=rootTupleElectrons_ElectronHcalIsoDR03_Token_;
  tokenMap["rootTupleElectrons_ElectronHcalIsoPAT"]=rootTupleElectrons_ElectronHcalIsoPAT_Token_;
  tokenMap["rootTupleElectrons_ElectronHoE"]=rootTupleElectrons_ElectronHoE_Token_;
  tokenMap["rootTupleElectrons_ElectronLeadVtxDistXY"]=rootTupleElectrons_ElectronLeadVtxDistXY_Token_;
  tokenMap["rootTupleElectrons_ElectronLeadVtxDistZ"]=rootTupleElectrons_ElectronLeadVtxDistZ_Token_;
  tokenMap["rootTupleElectrons_ElectronMatchedGenParticleEta"]=rootTupleElectrons_ElectronMatchedGenParticleEta_Token_;
  tokenMap["rootTupleElectrons_ElectronMatchedGenParticlePhi"]=rootTupleElectrons_ElectronMatchedGenParticlePhi_Token_;
  tokenMap["rootTupleElectrons_ElectronMatchedGenParticlePt"]=rootTupleElectrons_ElectronMatchedGenParticlePt_Token_;
  tokenMap["rootTupleElectrons_ElectronPFChargedHadronIso03"]=rootTupleElectrons_ElectronPFChargedHadronIso03_Token_;
  tokenMap["rootTupleElectrons_ElectronPFChargedHadronIso04"]=rootTupleElectrons_ElectronPFChargedHadronIso04_Token_;
  tokenMap["rootTupleElectrons_ElectronPFNeutralHadronIso03"]=rootTupleElectrons_ElectronPFNeutralHadronIso03_Token_;
  tokenMap["rootTupleElectrons_ElectronPFNeutralHadronIso04"]=rootTupleElectrons_ElectronPFNeutralHadronIso04_Token_;
  tokenMap["rootTupleElectrons_ElectronPFPUIso03"]=rootTupleElectrons_ElectronPFPUIso03_Token_;
  tokenMap["rootTupleElectrons_ElectronPFPhotonIso03"]=rootTupleElectrons_ElectronPFPhotonIso03_Token_;
  tokenMap["rootTupleElectrons_ElectronPFPhotonIso04"]=rootTupleElectrons_ElectronPFPhotonIso04_Token_;
  tokenMap["rootTupleElectrons_ElectronPhi"]=rootTupleElectrons_ElectronPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronPrimaryVertexDXYError"]=rootTupleElectrons_ElectronPrimaryVertexDXYError_Token_;
  tokenMap["rootTupleElectrons_ElectronPrimaryVertexDXY"]=rootTupleElectrons_ElectronPrimaryVertexDXY_Token_;
  tokenMap["rootTupleElectrons_ElectronPtHeep"]=rootTupleElectrons_ElectronPtHeep_Token_;
  tokenMap["rootTupleElectrons_ElectronPt"]=rootTupleElectrons_ElectronPt_Token_;
  tokenMap["rootTupleElectrons_ElectronR9"]=rootTupleElectrons_ElectronR9_Token_;
  tokenMap["rootTupleElectrons_ElectronRelIsoPAT"]=rootTupleElectrons_ElectronRelIsoPAT_Token_;
  tokenMap["rootTupleElectrons_ElectronSCEnergy"]=rootTupleElectrons_ElectronSCEnergy_Token_;
  tokenMap["rootTupleElectrons_ElectronSCEta"]=rootTupleElectrons_ElectronSCEta_Token_;
  tokenMap["rootTupleElectrons_ElectronSCPhi"]=rootTupleElectrons_ElectronSCPhi_Token_;
  tokenMap["rootTupleElectrons_ElectronSCPt"]=rootTupleElectrons_ElectronSCPt_Token_;
  tokenMap["rootTupleElectrons_ElectronSCRawEnergy"]=rootTupleElectrons_ElectronSCRawEnergy_Token_;
  tokenMap["rootTupleElectrons_ElectronSigmaEtaEta"]=rootTupleElectrons_ElectronSigmaEtaEta_Token_;
  tokenMap["rootTupleElectrons_ElectronSigmaIEtaIEta"]=rootTupleElectrons_ElectronSigmaIEtaIEta_Token_;
  tokenMap["rootTupleElectrons_ElectronTrackPt"]=rootTupleElectrons_ElectronTrackPt_Token_;
  tokenMap["rootTupleElectrons_ElectronTrackValidFractionOfHits"]=rootTupleElectrons_ElectronTrackValidFractionOfHits_Token_;
  tokenMap["rootTupleElectrons_ElectronTrackVx"]=rootTupleElectrons_ElectronTrackVx_Token_;
  tokenMap["rootTupleElectrons_ElectronTrackVy"]=rootTupleElectrons_ElectronTrackVy_Token_;
  tokenMap["rootTupleElectrons_ElectronTrackVz"]=rootTupleElectrons_ElectronTrackVz_Token_;
  tokenMap["rootTupleElectrons_ElectronTrkIsoDR03"]=rootTupleElectrons_ElectronTrkIsoDR03_Token_;
  tokenMap["rootTupleElectrons_ElectronTrkIsoPAT"]=rootTupleElectrons_ElectronTrkIsoPAT_Token_;
  tokenMap["rootTupleElectrons_ElectronVtxDistXY"]=rootTupleElectrons_ElectronVtxDistXY_Token_;
  tokenMap["rootTupleElectrons_ElectronVtxDistZ"]=rootTupleElectrons_ElectronVtxDistZ_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronEnergy"]=rootTupleGenElectronsFromWs_GenWElectronEnergy_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronEta"]=rootTupleGenElectronsFromWs_GenWElectronEta_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronMass"]=rootTupleGenElectronsFromWs_GenWElectronMass_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronP"]=rootTupleGenElectronsFromWs_GenWElectronP_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPhi"]=rootTupleGenElectronsFromWs_GenWElectronPhi_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPt"]=rootTupleGenElectronsFromWs_GenWElectronPt_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPx"]=rootTupleGenElectronsFromWs_GenWElectronPx_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPy"]=rootTupleGenElectronsFromWs_GenWElectronPy_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPz"]=rootTupleGenElectronsFromWs_GenWElectronPz_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronTauVisibleEta"]=rootTupleGenElectronsFromWs_GenWElectronTauVisibleEta_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronTauVisiblePhi"]=rootTupleGenElectronsFromWs_GenWElectronTauVisiblePhi_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronTauVisiblePt"]=rootTupleGenElectronsFromWs_GenWElectronTauVisiblePt_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronVX"]=rootTupleGenElectronsFromWs_GenWElectronVX_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronVY"]=rootTupleGenElectronsFromWs_GenWElectronVY_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronVZ"]=rootTupleGenElectronsFromWs_GenWElectronVZ_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronEnergy"]=rootTupleGenElectronsFromZs_GenZElectronEnergy_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronEta"]=rootTupleGenElectronsFromZs_GenZElectronEta_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronMass"]=rootTupleGenElectronsFromZs_GenZElectronMass_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronP"]=rootTupleGenElectronsFromZs_GenZElectronP_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPhi"]=rootTupleGenElectronsFromZs_GenZElectronPhi_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPt"]=rootTupleGenElectronsFromZs_GenZElectronPt_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPx"]=rootTupleGenElectronsFromZs_GenZElectronPx_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPy"]=rootTupleGenElectronsFromZs_GenZElectronPy_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPz"]=rootTupleGenElectronsFromZs_GenZElectronPz_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronTauVisibleEta"]=rootTupleGenElectronsFromZs_GenZElectronTauVisibleEta_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronTauVisiblePhi"]=rootTupleGenElectronsFromZs_GenZElectronTauVisiblePhi_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronTauVisiblePt"]=rootTupleGenElectronsFromZs_GenZElectronTauVisiblePt_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronVX"]=rootTupleGenElectronsFromZs_GenZElectronVX_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronVY"]=rootTupleGenElectronsFromZs_GenZElectronVY_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronVZ"]=rootTupleGenElectronsFromZs_GenZElectronVZ_Token_;
  tokenMap["rootTupleGenEventInfo_PDFCTEQWeights"]=rootTupleGenEventInfo_PDFCTEQWeights_Token_;
  tokenMap["rootTupleGenEventInfo_PDFMSTWWeights"]=rootTupleGenEventInfo_PDFMSTWWeights_Token_;
  tokenMap["rootTupleGenEventInfo_PDFNNPDFWeights"]=rootTupleGenEventInfo_PDFNNPDFWeights_Token_;
  tokenMap["rootTupleGenEventInfo_ScaleWeights"]=rootTupleGenEventInfo_ScaleWeights_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetEMFAK4"]=rootTupleGenJetsAK4_GenJetEMFAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetEnergyAK4"]=rootTupleGenJetsAK4_GenJetEnergyAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetEtaAK4"]=rootTupleGenJetsAK4_GenJetEtaAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetHADFAK4"]=rootTupleGenJetsAK4_GenJetHADFAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetPAK4"]=rootTupleGenJetsAK4_GenJetPAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetPhiAK4"]=rootTupleGenJetsAK4_GenJetPhiAK4_Token_;
  tokenMap["rootTupleGenJetsAK4_GenJetPtAK4"]=rootTupleGenJetsAK4_GenJetPtAK4_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetEMFAK5"]=rootTupleGenJetsAK5_GenJetEMFAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetEnergyAK5"]=rootTupleGenJetsAK5_GenJetEnergyAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetEtaAK5"]=rootTupleGenJetsAK5_GenJetEtaAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetHADFAK5"]=rootTupleGenJetsAK5_GenJetHADFAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetPAK5"]=rootTupleGenJetsAK5_GenJetPAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetPhiAK5"]=rootTupleGenJetsAK5_GenJetPhiAK5_Token_;
  tokenMap["rootTupleGenJetsAK5_GenJetPtAK5"]=rootTupleGenJetsAK5_GenJetPtAK5_Token_;
  tokenMap["rootTupleGenMETTrue_GenMETPhiTrue"]=rootTupleGenMETTrue_GenMETPhiTrue_Token_;
  tokenMap["rootTupleGenMETTrue_GenMETTrue"]=rootTupleGenMETTrue_GenMETTrue_Token_;
  tokenMap["rootTupleGenMETTrue_GenSumETTrue"]=rootTupleGenMETTrue_GenSumETTrue_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuEnergy"]=rootTupleGenMuonsFromWs_GenWMuEnergy_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuEta"]=rootTupleGenMuonsFromWs_GenWMuEta_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuMass"]=rootTupleGenMuonsFromWs_GenWMuMass_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuP"]=rootTupleGenMuonsFromWs_GenWMuP_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPhi"]=rootTupleGenMuonsFromWs_GenWMuPhi_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPt"]=rootTupleGenMuonsFromWs_GenWMuPt_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPx"]=rootTupleGenMuonsFromWs_GenWMuPx_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPy"]=rootTupleGenMuonsFromWs_GenWMuPy_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPz"]=rootTupleGenMuonsFromWs_GenWMuPz_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuTauVisibleEta"]=rootTupleGenMuonsFromWs_GenWMuTauVisibleEta_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuTauVisiblePhi"]=rootTupleGenMuonsFromWs_GenWMuTauVisiblePhi_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuTauVisiblePt"]=rootTupleGenMuonsFromWs_GenWMuTauVisiblePt_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuVX"]=rootTupleGenMuonsFromWs_GenWMuVX_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuVY"]=rootTupleGenMuonsFromWs_GenWMuVY_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuVZ"]=rootTupleGenMuonsFromWs_GenWMuVZ_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuEnergy"]=rootTupleGenMuonsFromZs_GenZMuEnergy_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuEta"]=rootTupleGenMuonsFromZs_GenZMuEta_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuMass"]=rootTupleGenMuonsFromZs_GenZMuMass_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuP"]=rootTupleGenMuonsFromZs_GenZMuP_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPhi"]=rootTupleGenMuonsFromZs_GenZMuPhi_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPt"]=rootTupleGenMuonsFromZs_GenZMuPt_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPx"]=rootTupleGenMuonsFromZs_GenZMuPx_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPy"]=rootTupleGenMuonsFromZs_GenZMuPy_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPz"]=rootTupleGenMuonsFromZs_GenZMuPz_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuTauVisibleEta"]=rootTupleGenMuonsFromZs_GenZMuTauVisibleEta_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuTauVisiblePhi"]=rootTupleGenMuonsFromZs_GenZMuTauVisiblePhi_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuTauVisiblePt"]=rootTupleGenMuonsFromZs_GenZMuTauVisiblePt_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuVX"]=rootTupleGenMuonsFromZs_GenZMuVX_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuVY"]=rootTupleGenMuonsFromZs_GenZMuVY_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuVZ"]=rootTupleGenMuonsFromZs_GenZMuVZ_Token_;
  tokenMap["rootTupleGenParticles_GenParticleEnergy"]=rootTupleGenParticles_GenParticleEnergy_Token_;
  tokenMap["rootTupleGenParticles_GenParticleEta"]=rootTupleGenParticles_GenParticleEta_Token_;
  tokenMap["rootTupleGenParticles_GenParticleMass"]=rootTupleGenParticles_GenParticleMass_Token_;
  tokenMap["rootTupleGenParticles_GenParticleP"]=rootTupleGenParticles_GenParticleP_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePhi"]=rootTupleGenParticles_GenParticlePhi_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePt"]=rootTupleGenParticles_GenParticlePt_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePx"]=rootTupleGenParticles_GenParticlePx_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePy"]=rootTupleGenParticles_GenParticlePy_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePz"]=rootTupleGenParticles_GenParticlePz_Token_;
  tokenMap["rootTupleGenParticles_GenParticleTauVisibleEta"]=rootTupleGenParticles_GenParticleTauVisibleEta_Token_;
  tokenMap["rootTupleGenParticles_GenParticleTauVisiblePhi"]=rootTupleGenParticles_GenParticleTauVisiblePhi_Token_;
  tokenMap["rootTupleGenParticles_GenParticleTauVisiblePt"]=rootTupleGenParticles_GenParticleTauVisiblePt_Token_;
  tokenMap["rootTupleGenParticles_GenParticleVX"]=rootTupleGenParticles_GenParticleVX_Token_;
  tokenMap["rootTupleGenParticles_GenParticleVY"]=rootTupleGenParticles_GenParticleVY_Token_;
  tokenMap["rootTupleGenParticles_GenParticleVZ"]=rootTupleGenParticles_GenParticleVZ_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauEnergy"]=rootTupleGenTausFromWs_GenWTauEnergy_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauEta"]=rootTupleGenTausFromWs_GenWTauEta_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauMass"]=rootTupleGenTausFromWs_GenWTauMass_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauP"]=rootTupleGenTausFromWs_GenWTauP_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPhi"]=rootTupleGenTausFromWs_GenWTauPhi_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPt"]=rootTupleGenTausFromWs_GenWTauPt_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPx"]=rootTupleGenTausFromWs_GenWTauPx_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPy"]=rootTupleGenTausFromWs_GenWTauPy_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPz"]=rootTupleGenTausFromWs_GenWTauPz_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauTauVisibleEta"]=rootTupleGenTausFromWs_GenWTauTauVisibleEta_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauTauVisiblePhi"]=rootTupleGenTausFromWs_GenWTauTauVisiblePhi_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauTauVisiblePt"]=rootTupleGenTausFromWs_GenWTauTauVisiblePt_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauVX"]=rootTupleGenTausFromWs_GenWTauVX_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauVY"]=rootTupleGenTausFromWs_GenWTauVY_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauVZ"]=rootTupleGenTausFromWs_GenWTauVZ_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauEnergy"]=rootTupleGenTausFromZs_GenZTauEnergy_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauEta"]=rootTupleGenTausFromZs_GenZTauEta_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauMass"]=rootTupleGenTausFromZs_GenZTauMass_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauP"]=rootTupleGenTausFromZs_GenZTauP_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPhi"]=rootTupleGenTausFromZs_GenZTauPhi_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPt"]=rootTupleGenTausFromZs_GenZTauPt_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPx"]=rootTupleGenTausFromZs_GenZTauPx_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPy"]=rootTupleGenTausFromZs_GenZTauPy_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPz"]=rootTupleGenTausFromZs_GenZTauPz_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauTauVisibleEta"]=rootTupleGenTausFromZs_GenZTauTauVisibleEta_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauTauVisiblePhi"]=rootTupleGenTausFromZs_GenZTauTauVisiblePhi_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauTauVisiblePt"]=rootTupleGenTausFromZs_GenZTauTauVisiblePt_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauVX"]=rootTupleGenTausFromZs_GenZTauVX_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauVY"]=rootTupleGenTausFromZs_GenZTauVY_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauVZ"]=rootTupleGenTausFromZs_GenZTauVZ_Token_;
  tokenMap["rootTupleMuons_MuonBackToBackCompatibility"]=rootTupleMuons_MuonBackToBackCompatibility_Token_;
  tokenMap["rootTupleMuons_MuonBeamSpotDXYError"]=rootTupleMuons_MuonBeamSpotDXYError_Token_;
  tokenMap["rootTupleMuons_MuonBeamSpotDXY"]=rootTupleMuons_MuonBeamSpotDXY_Token_;
  tokenMap["rootTupleMuons_MuonBestTrackVtxDistXY"]=rootTupleMuons_MuonBestTrackVtxDistXY_Token_;
  tokenMap["rootTupleMuons_MuonBestTrackVtxDistZ"]=rootTupleMuons_MuonBestTrackVtxDistZ_Token_;
  tokenMap["rootTupleMuons_MuonCocktailEtaError"]=rootTupleMuons_MuonCocktailEtaError_Token_;
  tokenMap["rootTupleMuons_MuonCocktailEta"]=rootTupleMuons_MuonCocktailEta_Token_;
  tokenMap["rootTupleMuons_MuonCocktailGlobalChi2"]=rootTupleMuons_MuonCocktailGlobalChi2_Token_;
  tokenMap["rootTupleMuons_MuonCocktailP"]=rootTupleMuons_MuonCocktailP_Token_;
  tokenMap["rootTupleMuons_MuonCocktailPhiError"]=rootTupleMuons_MuonCocktailPhiError_Token_;
  tokenMap["rootTupleMuons_MuonCocktailPhi"]=rootTupleMuons_MuonCocktailPhi_Token_;
  tokenMap["rootTupleMuons_MuonCocktailPtError"]=rootTupleMuons_MuonCocktailPtError_Token_;
  tokenMap["rootTupleMuons_MuonCocktailPt"]=rootTupleMuons_MuonCocktailPt_Token_;
  tokenMap["rootTupleMuons_MuonCocktailQOverPError"]=rootTupleMuons_MuonCocktailQOverPError_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkD0Error"]=rootTupleMuons_MuonCocktailTrkD0Error_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkD0"]=rootTupleMuons_MuonCocktailTrkD0_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkDzError"]=rootTupleMuons_MuonCocktailTrkDzError_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkDz"]=rootTupleMuons_MuonCocktailTrkDz_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkValidFractionOfHits"]=rootTupleMuons_MuonCocktailTrkValidFractionOfHits_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkVtxDXY"]=rootTupleMuons_MuonCocktailTrkVtxDXY_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkVtxDZ"]=rootTupleMuons_MuonCocktailTrkVtxDZ_Token_;
  tokenMap["rootTupleMuons_MuonCosmicCompatibility"]=rootTupleMuons_MuonCosmicCompatibility_Token_;
  tokenMap["rootTupleMuons_MuonEcalIso"]=rootTupleMuons_MuonEcalIso_Token_;
  tokenMap["rootTupleMuons_MuonEcalVetoIso"]=rootTupleMuons_MuonEcalVetoIso_Token_;
  tokenMap["rootTupleMuons_MuonEnergy"]=rootTupleMuons_MuonEnergy_Token_;
  tokenMap["rootTupleMuons_MuonEtaError"]=rootTupleMuons_MuonEtaError_Token_;
  tokenMap["rootTupleMuons_MuonEta"]=rootTupleMuons_MuonEta_Token_;
  tokenMap["rootTupleMuons_MuonGlobalChi2"]=rootTupleMuons_MuonGlobalChi2_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleIsoMuonMatchEta"]=rootTupleMuons_MuonHLTSingleIsoMuonMatchEta_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleIsoMuonMatchPhi"]=rootTupleMuons_MuonHLTSingleIsoMuonMatchPhi_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleIsoMuonMatchPt"]=rootTupleMuons_MuonHLTSingleIsoMuonMatchPt_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleMuonMatchEta"]=rootTupleMuons_MuonHLTSingleMuonMatchEta_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleMuonMatchPhi"]=rootTupleMuons_MuonHLTSingleMuonMatchPhi_Token_;
  tokenMap["rootTupleMuons_MuonHLTSingleMuonMatchPt"]=rootTupleMuons_MuonHLTSingleMuonMatchPt_Token_;
  tokenMap["rootTupleMuons_MuonHOIso"]=rootTupleMuons_MuonHOIso_Token_;
  tokenMap["rootTupleMuons_MuonHcalIso"]=rootTupleMuons_MuonHcalIso_Token_;
  tokenMap["rootTupleMuons_MuonHcalVetoIso"]=rootTupleMuons_MuonHcalVetoIso_Token_;
  tokenMap["rootTupleMuons_MuonMatchedGenParticleEta"]=rootTupleMuons_MuonMatchedGenParticleEta_Token_;
  tokenMap["rootTupleMuons_MuonMatchedGenParticlePhi"]=rootTupleMuons_MuonMatchedGenParticlePhi_Token_;
  tokenMap["rootTupleMuons_MuonMatchedGenParticlePt"]=rootTupleMuons_MuonMatchedGenParticlePt_Token_;
  tokenMap["rootTupleMuons_MuonOverlapCompatibility"]=rootTupleMuons_MuonOverlapCompatibility_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03ChargedHadron"]=rootTupleMuons_MuonPFIsoR03ChargedHadron_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03ChargedParticle"]=rootTupleMuons_MuonPFIsoR03ChargedParticle_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03NeutralHadronHT"]=rootTupleMuons_MuonPFIsoR03NeutralHadronHT_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03NeutralHadron"]=rootTupleMuons_MuonPFIsoR03NeutralHadron_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03PU"]=rootTupleMuons_MuonPFIsoR03PU_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03PhotonHT"]=rootTupleMuons_MuonPFIsoR03PhotonHT_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR03Photon"]=rootTupleMuons_MuonPFIsoR03Photon_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04ChargedHadron"]=rootTupleMuons_MuonPFIsoR04ChargedHadron_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04ChargedParticle"]=rootTupleMuons_MuonPFIsoR04ChargedParticle_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04NeutralHadronHT"]=rootTupleMuons_MuonPFIsoR04NeutralHadronHT_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04NeutralHadron"]=rootTupleMuons_MuonPFIsoR04NeutralHadron_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04PU"]=rootTupleMuons_MuonPFIsoR04PU_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04PhotonHT"]=rootTupleMuons_MuonPFIsoR04PhotonHT_Token_;
  tokenMap["rootTupleMuons_MuonPFIsoR04Photon"]=rootTupleMuons_MuonPFIsoR04Photon_Token_;
  tokenMap["rootTupleMuons_MuonP"]=rootTupleMuons_MuonP_Token_;
  tokenMap["rootTupleMuons_MuonPhiError"]=rootTupleMuons_MuonPhiError_Token_;
  tokenMap["rootTupleMuons_MuonPhi"]=rootTupleMuons_MuonPhi_Token_;
  tokenMap["rootTupleMuons_MuonPrimaryVertexDXYError"]=rootTupleMuons_MuonPrimaryVertexDXYError_Token_;
  tokenMap["rootTupleMuons_MuonPrimaryVertexDXY"]=rootTupleMuons_MuonPrimaryVertexDXY_Token_;
  tokenMap["rootTupleMuons_MuonPtError"]=rootTupleMuons_MuonPtError_Token_;
  tokenMap["rootTupleMuons_MuonPt"]=rootTupleMuons_MuonPt_Token_;
  tokenMap["rootTupleMuons_MuonQOverPError"]=rootTupleMuons_MuonQOverPError_Token_;
  tokenMap["rootTupleMuons_MuonTimeCompatibility"]=rootTupleMuons_MuonTimeCompatibility_Token_;
  tokenMap["rootTupleMuons_MuonTrackChi2"]=rootTupleMuons_MuonTrackChi2_Token_;
  tokenMap["rootTupleMuons_MuonTrackerIsoSumPT"]=rootTupleMuons_MuonTrackerIsoSumPT_Token_;
  tokenMap["rootTupleMuons_MuonTrkD0Error"]=rootTupleMuons_MuonTrkD0Error_Token_;
  tokenMap["rootTupleMuons_MuonTrkD0"]=rootTupleMuons_MuonTrkD0_Token_;
  tokenMap["rootTupleMuons_MuonTrkDzError"]=rootTupleMuons_MuonTrkDzError_Token_;
  tokenMap["rootTupleMuons_MuonTrkDz"]=rootTupleMuons_MuonTrkDz_Token_;
  tokenMap["rootTupleMuons_MuonTrkEtaError"]=rootTupleMuons_MuonTrkEtaError_Token_;
  tokenMap["rootTupleMuons_MuonTrkEta"]=rootTupleMuons_MuonTrkEta_Token_;
  tokenMap["rootTupleMuons_MuonTrkIso"]=rootTupleMuons_MuonTrkIso_Token_;
  tokenMap["rootTupleMuons_MuonTrkPhiError"]=rootTupleMuons_MuonTrkPhiError_Token_;
  tokenMap["rootTupleMuons_MuonTrkPhi"]=rootTupleMuons_MuonTrkPhi_Token_;
  tokenMap["rootTupleMuons_MuonTrkPtError"]=rootTupleMuons_MuonTrkPtError_Token_;
  tokenMap["rootTupleMuons_MuonTrkPt"]=rootTupleMuons_MuonTrkPt_Token_;
  tokenMap["rootTupleMuons_MuonTrkValidFractionOfHits"]=rootTupleMuons_MuonTrkValidFractionOfHits_Token_;
  tokenMap["rootTupleMuons_MuonTrkVx"]=rootTupleMuons_MuonTrkVx_Token_;
  tokenMap["rootTupleMuons_MuonTrkVy"]=rootTupleMuons_MuonTrkVy_Token_;
  tokenMap["rootTupleMuons_MuonTrkVz"]=rootTupleMuons_MuonTrkVz_Token_;
  tokenMap["rootTupleMuons_MuonVtxDistXY"]=rootTupleMuons_MuonVtxDistXY_Token_;
  tokenMap["rootTupleMuons_MuonVtxDistZ"]=rootTupleMuons_MuonVtxDistZ_Token_;
  tokenMap["rootTuplePFCandidates_PFCandEnergyLeptLink"]=rootTuplePFCandidates_PFCandEnergyLeptLink_Token_;
  tokenMap["rootTuplePFCandidates_PFCandEtaLeptLink"]=rootTuplePFCandidates_PFCandEtaLeptLink_Token_;
  tokenMap["rootTuplePFCandidates_PFCandPhiLeptLink"]=rootTuplePFCandidates_PFCandPhiLeptLink_Token_;
  tokenMap["rootTuplePFCandidates_PFCandPtLeptLink"]=rootTuplePFCandidates_PFCandPtLeptLink_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBestVertexTrackAssociationFactorAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBestVertexTrackAssociationFactorAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBetaAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBetaAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBetaClassicAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBetaClassicAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBetaStarAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBetaStarAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBetaStarClassicAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBetaStarClassicAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetChargedEmEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetChargedEmEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetChargedHadronEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetChargedHadronEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetChargedMuEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetChargedMuEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertexWeighted3DSeparationAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertexWeighted3DSeparationAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertexWeightedXYSeparationAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertexWeightedXYSeparationAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertexWeightedZSeparationAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertexWeightedZSeparationAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetCombinedInclusiveSecondaryVertexBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetCombinedInclusiveSecondaryVertexBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetCombinedMVABTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetCombinedMVABTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetCombinedSecondaryVertexBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetCombinedSecondaryVertexBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetElectronEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetElectronEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetEnergyAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetEnergyAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetEnergyRawAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetEnergyRawAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetEtaAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetEtaAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetHFEMEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetHFEMEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetHFHadronEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetHFHadronEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetJECUncAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetJECUncAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetJetBProbabilityBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetJetBProbabilityBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetJetProbabilityBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetJetProbabilityBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetL1FastJetJECAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetL1FastJetJECAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetL2L3ResJECAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetL2L3ResJECAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetL2RelJECAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetL2RelJECAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetL3AbsJECAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetL3AbsJECAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetMuonEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetMuonEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetNeutralEmEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetNeutralEmEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetNeutralHadronEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetNeutralHadronEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPhiAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPhiAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPhotonEnergyFractionAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPhotonEnergyFractionAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPileupMVAAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPileupMVAAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPtAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPtAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPtRawAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPtRawAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetScaledDownEnergyAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetScaledDownEnergyAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetScaledDownPtAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetScaledDownPtAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetScaledUpEnergyAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetScaledUpEnergyAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetScaledUpPtAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetScaledUpPtAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSimpleSecondaryVertexHighEffBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSimpleSecondaryVertexHighEffBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSimpleSecondaryVertexHighPurBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSimpleSecondaryVertexHighPurBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSmearedDownEnergyAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSmearedDownEnergyAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSmearedDownPtAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSmearedDownPtAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSmearedUpEnergyAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSmearedUpEnergyAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetSmearedUpPtAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetSmearedUpPtAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetTrackCountingHighEffBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetTrackCountingHighEffBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetTrackCountingHighPurBTagAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetTrackCountingHighPurBTagAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBestVertexTrackAssociationFactorAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBestVertexTrackAssociationFactorAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBetaAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBetaAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBetaClassicAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBetaClassicAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBetaStarAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBetaStarAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBetaStarClassicAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBetaStarClassicAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetChargedEmEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetChargedEmEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetChargedHadronEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetChargedHadronEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetChargedMuEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetChargedMuEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeighted3DSeparationAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeighted3DSeparationAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeightedXYSeparationAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeightedXYSeparationAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeightedZSeparationAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertexWeightedZSeparationAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetCombinedInclusiveSecondaryVertexBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetCombinedInclusiveSecondaryVertexBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetCombinedMVABTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetCombinedMVABTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetCombinedSecondaryVertexBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetCombinedSecondaryVertexBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetElectronEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetElectronEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetEnergyAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetEnergyAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetEnergyRawAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetEnergyRawAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetEtaAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetEtaAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetHFEMEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetHFEMEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetHFHadronEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetHFHadronEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetJECUncAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetJECUncAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetJetBProbabilityBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetJetBProbabilityBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetJetProbabilityBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetJetProbabilityBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetL1FastJetJECAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetL1FastJetJECAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetL2L3ResJECAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetL2L3ResJECAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetL2RelJECAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetL2RelJECAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetL3AbsJECAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetL3AbsJECAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetMuonEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetMuonEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetNeutralEmEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetNeutralEmEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetNeutralHadronEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetNeutralHadronEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPhiAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPhiAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPhotonEnergyFractionAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPhotonEnergyFractionAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPtAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPtAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPtRawAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPtRawAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetScaledDownEnergyAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetScaledDownEnergyAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetScaledDownPtAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetScaledDownPtAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetScaledUpEnergyAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetScaledUpEnergyAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetScaledUpPtAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetScaledUpPtAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSimpleSecondaryVertexHighEffBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSimpleSecondaryVertexHighEffBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSimpleSecondaryVertexHighPurBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSimpleSecondaryVertexHighPurBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSmearedDownEnergyAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSmearedDownEnergyAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSmearedDownPtAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSmearedDownPtAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSmearedUpEnergyAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSmearedUpEnergyAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetSmearedUpPtAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetSmearedUpPtAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetTrackCountingHighEffBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetTrackCountingHighEffBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetTrackCountingHighPurBTagAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetTrackCountingHighPurBTagAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBestVertexTrackAssociationFactorAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBestVertexTrackAssociationFactorAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBetaAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBetaAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBetaClassicAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBetaClassicAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBetaStarAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBetaStarAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBetaStarClassicAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBetaStarClassicAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetChargedEmEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetChargedEmEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetChargedHadronEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetChargedHadronEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetChargedMuEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetChargedMuEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertexWeighted3DSeparationAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertexWeighted3DSeparationAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertexWeightedXYSeparationAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertexWeightedXYSeparationAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertexWeightedZSeparationAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertexWeightedZSeparationAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetCombinedInclusiveSecondaryVertexBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetCombinedInclusiveSecondaryVertexBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetCombinedMVABTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetCombinedMVABTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetCombinedSecondaryVertexBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetCombinedSecondaryVertexBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetElectronEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetElectronEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetEnergyAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetEnergyAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetEnergyRawAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetEnergyRawAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetEtaAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetEtaAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetHFEMEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetHFEMEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetHFHadronEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetHFHadronEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetJECUncAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetJECUncAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetJetBProbabilityBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetJetBProbabilityBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetJetProbabilityBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetJetProbabilityBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetL1FastJetJECAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetL1FastJetJECAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetL2L3ResJECAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetL2L3ResJECAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetL2RelJECAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetL2RelJECAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetL3AbsJECAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetL3AbsJECAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetMuonEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetMuonEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetNeutralEmEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetNeutralEmEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetNeutralHadronEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetNeutralHadronEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPhiAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPhiAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPhotonEnergyFractionAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPhotonEnergyFractionAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPileupMVAAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPileupMVAAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPtAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPtAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPtRawAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPtRawAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetScaledDownEnergyAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetScaledDownEnergyAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetScaledDownPtAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetScaledDownPtAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetScaledUpEnergyAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetScaledUpEnergyAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetScaledUpPtAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetScaledUpPtAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSimpleSecondaryVertexHighEffBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSimpleSecondaryVertexHighEffBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSimpleSecondaryVertexHighPurBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSimpleSecondaryVertexHighPurBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSmearedDownEnergyAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSmearedDownEnergyAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSmearedDownPtAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSmearedDownPtAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSmearedUpEnergyAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSmearedUpEnergyAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetSmearedUpPtAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetSmearedUpPtAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetTrackCountingHighEffBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetTrackCountingHighEffBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetTrackCountingHighPurBTagAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetTrackCountingHighPurBTagAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBestVertexTrackAssociationFactorAK5"]=rootTuplePFJetsAK5_PFJetBestVertexTrackAssociationFactorAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBetaAK5"]=rootTuplePFJetsAK5_PFJetBetaAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBetaClassicAK5"]=rootTuplePFJetsAK5_PFJetBetaClassicAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBetaStarAK5"]=rootTuplePFJetsAK5_PFJetBetaStarAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBetaStarClassicAK5"]=rootTuplePFJetsAK5_PFJetBetaStarClassicAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetChargedEmEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetChargedEmEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetChargedHadronEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetChargedHadronEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetChargedMuEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetChargedMuEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertexWeighted3DSeparationAK5"]=rootTuplePFJetsAK5_PFJetClosestVertexWeighted3DSeparationAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertexWeightedXYSeparationAK5"]=rootTuplePFJetsAK5_PFJetClosestVertexWeightedXYSeparationAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertexWeightedZSeparationAK5"]=rootTuplePFJetsAK5_PFJetClosestVertexWeightedZSeparationAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetCombinedInclusiveSecondaryVertexBTagAK5"]=rootTuplePFJetsAK5_PFJetCombinedInclusiveSecondaryVertexBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetCombinedMVABTagAK5"]=rootTuplePFJetsAK5_PFJetCombinedMVABTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetCombinedSecondaryVertexBTagAK5"]=rootTuplePFJetsAK5_PFJetCombinedSecondaryVertexBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetElectronEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetElectronEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetEnergyAK5"]=rootTuplePFJetsAK5_PFJetEnergyAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetEnergyRawAK5"]=rootTuplePFJetsAK5_PFJetEnergyRawAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetEtaAK5"]=rootTuplePFJetsAK5_PFJetEtaAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetHFEMEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetHFEMEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetHFHadronEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetHFHadronEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetJECUncAK5"]=rootTuplePFJetsAK5_PFJetJECUncAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetJetBProbabilityBTagAK5"]=rootTuplePFJetsAK5_PFJetJetBProbabilityBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetJetProbabilityBTagAK5"]=rootTuplePFJetsAK5_PFJetJetProbabilityBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetL1FastJetJECAK5"]=rootTuplePFJetsAK5_PFJetL1FastJetJECAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetL2L3ResJECAK5"]=rootTuplePFJetsAK5_PFJetL2L3ResJECAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetL2RelJECAK5"]=rootTuplePFJetsAK5_PFJetL2RelJECAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetL3AbsJECAK5"]=rootTuplePFJetsAK5_PFJetL3AbsJECAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetMuonEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetMuonEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetNeutralEmEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetNeutralEmEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetNeutralHadronEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetNeutralHadronEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPhiAK5"]=rootTuplePFJetsAK5_PFJetPhiAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPhotonEnergyFractionAK5"]=rootTuplePFJetsAK5_PFJetPhotonEnergyFractionAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPileupMVAAK5"]=rootTuplePFJetsAK5_PFJetPileupMVAAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPtAK5"]=rootTuplePFJetsAK5_PFJetPtAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPtRawAK5"]=rootTuplePFJetsAK5_PFJetPtRawAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetScaledDownEnergyAK5"]=rootTuplePFJetsAK5_PFJetScaledDownEnergyAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetScaledDownPtAK5"]=rootTuplePFJetsAK5_PFJetScaledDownPtAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetScaledUpEnergyAK5"]=rootTuplePFJetsAK5_PFJetScaledUpEnergyAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetScaledUpPtAK5"]=rootTuplePFJetsAK5_PFJetScaledUpPtAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSimpleSecondaryVertexHighEffBTagAK5"]=rootTuplePFJetsAK5_PFJetSimpleSecondaryVertexHighEffBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSimpleSecondaryVertexHighPurBTagAK5"]=rootTuplePFJetsAK5_PFJetSimpleSecondaryVertexHighPurBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSmearedDownEnergyAK5"]=rootTuplePFJetsAK5_PFJetSmearedDownEnergyAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSmearedDownPtAK5"]=rootTuplePFJetsAK5_PFJetSmearedDownPtAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSmearedUpEnergyAK5"]=rootTuplePFJetsAK5_PFJetSmearedUpEnergyAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetSmearedUpPtAK5"]=rootTuplePFJetsAK5_PFJetSmearedUpPtAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetTrackCountingHighEffBTagAK5"]=rootTuplePFJetsAK5_PFJetTrackCountingHighEffBTagAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetTrackCountingHighPurBTagAK5"]=rootTuplePFJetsAK5_PFJetTrackCountingHighPurBTagAK5_Token_;
  tokenMap["rootTuplePFMETPuppiType1Cor_PFMETPhiPuppiType1Cor"]=rootTuplePFMETPuppiType1Cor_PFMETPhiPuppiType1Cor_Token_;
  tokenMap["rootTuplePFMETPuppiType1Cor_PFMETPuppiType1Cor"]=rootTuplePFMETPuppiType1Cor_PFMETPuppiType1Cor_Token_;
  tokenMap["rootTuplePFMETPuppiType1Cor_PFSumETPuppiType1Cor"]=rootTuplePFMETPuppiType1Cor_PFSumETPuppiType1Cor_Token_;
  tokenMap["rootTuplePFMETPuppi_PFMETPhiPuppi"]=rootTuplePFMETPuppi_PFMETPhiPuppi_Token_;
  tokenMap["rootTuplePFMETPuppi_PFMETPuppi"]=rootTuplePFMETPuppi_PFMETPuppi_Token_;
  tokenMap["rootTuplePFMETPuppi_PFSumETPuppi"]=rootTuplePFMETPuppi_PFSumETPuppi_Token_;
  tokenMap["rootTuplePFMETType01Cor_PFMETPhiType01Cor"]=rootTuplePFMETType01Cor_PFMETPhiType01Cor_Token_;
  tokenMap["rootTuplePFMETType01Cor_PFMETType01Cor"]=rootTuplePFMETType01Cor_PFMETType01Cor_Token_;
  tokenMap["rootTuplePFMETType01Cor_PFSumETType01Cor"]=rootTuplePFMETType01Cor_PFSumETType01Cor_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnDown_PFMETPhiType1CorElectronEnDown"]=rootTuplePFMETType01XYCorElectronEnDown_PFMETPhiType1CorElectronEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnDown_PFMETType1CorElectronEnDown"]=rootTuplePFMETType01XYCorElectronEnDown_PFMETType1CorElectronEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnDown_PFSumETType1CorElectronEnDown"]=rootTuplePFMETType01XYCorElectronEnDown_PFSumETType1CorElectronEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnUp_PFMETPhiType1CorElectronEnUp"]=rootTuplePFMETType01XYCorElectronEnUp_PFMETPhiType1CorElectronEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnUp_PFMETType1CorElectronEnUp"]=rootTuplePFMETType01XYCorElectronEnUp_PFMETType1CorElectronEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorElectronEnUp_PFSumETType1CorElectronEnUp"]=rootTuplePFMETType01XYCorElectronEnUp_PFSumETType1CorElectronEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnDown_PFMETPhiType1CorJetEnDown"]=rootTuplePFMETType01XYCorJetEnDown_PFMETPhiType1CorJetEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnDown_PFMETType1CorJetEnDown"]=rootTuplePFMETType01XYCorJetEnDown_PFMETType1CorJetEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnDown_PFSumETType1CorJetEnDown"]=rootTuplePFMETType01XYCorJetEnDown_PFSumETType1CorJetEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnUp_PFMETPhiType1CorJetEnUp"]=rootTuplePFMETType01XYCorJetEnUp_PFMETPhiType1CorJetEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnUp_PFMETType1CorJetEnUp"]=rootTuplePFMETType01XYCorJetEnUp_PFMETType1CorJetEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetEnUp_PFSumETType1CorJetEnUp"]=rootTuplePFMETType01XYCorJetEnUp_PFSumETType1CorJetEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResDown_PFMETPhiType1CorJetResDown"]=rootTuplePFMETType01XYCorJetResDown_PFMETPhiType1CorJetResDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResDown_PFMETType1CorJetResDown"]=rootTuplePFMETType01XYCorJetResDown_PFMETType1CorJetResDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResDown_PFSumETType1CorJetResDown"]=rootTuplePFMETType01XYCorJetResDown_PFSumETType1CorJetResDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResUp_PFMETPhiType1CorJetResUp"]=rootTuplePFMETType01XYCorJetResUp_PFMETPhiType1CorJetResUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResUp_PFMETType1CorJetResUp"]=rootTuplePFMETType01XYCorJetResUp_PFMETType1CorJetResUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorJetResUp_PFSumETType1CorJetResUp"]=rootTuplePFMETType01XYCorJetResUp_PFSumETType1CorJetResUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnDown_PFMETPhiType1CorMuonEnDown"]=rootTuplePFMETType01XYCorMuonEnDown_PFMETPhiType1CorMuonEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnDown_PFMETType1CorMuonEnDown"]=rootTuplePFMETType01XYCorMuonEnDown_PFMETType1CorMuonEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnDown_PFSumETType1CorMuonEnDown"]=rootTuplePFMETType01XYCorMuonEnDown_PFSumETType1CorMuonEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnUp_PFMETPhiType1CorMuonEnUp"]=rootTuplePFMETType01XYCorMuonEnUp_PFMETPhiType1CorMuonEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnUp_PFMETType1CorMuonEnUp"]=rootTuplePFMETType01XYCorMuonEnUp_PFMETType1CorMuonEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorMuonEnUp_PFSumETType1CorMuonEnUp"]=rootTuplePFMETType01XYCorMuonEnUp_PFSumETType1CorMuonEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnDown_PFMETPhiType1CorTauEnDown"]=rootTuplePFMETType01XYCorTauEnDown_PFMETPhiType1CorTauEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnDown_PFMETType1CorTauEnDown"]=rootTuplePFMETType01XYCorTauEnDown_PFMETType1CorTauEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnDown_PFSumETType1CorTauEnDown"]=rootTuplePFMETType01XYCorTauEnDown_PFSumETType1CorTauEnDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnUp_PFMETPhiType1CorTauEnUp"]=rootTuplePFMETType01XYCorTauEnUp_PFMETPhiType1CorTauEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnUp_PFMETType1CorTauEnUp"]=rootTuplePFMETType01XYCorTauEnUp_PFMETType1CorTauEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorTauEnUp_PFSumETType1CorTauEnUp"]=rootTuplePFMETType01XYCorTauEnUp_PFSumETType1CorTauEnUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredDown_PFMETPhiType1CorUnclusteredDown"]=rootTuplePFMETType01XYCorUnclusteredDown_PFMETPhiType1CorUnclusteredDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredDown_PFMETType1CorUnclusteredDown"]=rootTuplePFMETType01XYCorUnclusteredDown_PFMETType1CorUnclusteredDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredDown_PFSumETType1CorUnclusteredDown"]=rootTuplePFMETType01XYCorUnclusteredDown_PFSumETType1CorUnclusteredDown_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredUp_PFMETPhiType1CorUnclusteredUp"]=rootTuplePFMETType01XYCorUnclusteredUp_PFMETPhiType1CorUnclusteredUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredUp_PFMETType1CorUnclusteredUp"]=rootTuplePFMETType01XYCorUnclusteredUp_PFMETType1CorUnclusteredUp_Token_;
  tokenMap["rootTuplePFMETType01XYCorUnclusteredUp_PFSumETType1CorUnclusteredUp"]=rootTuplePFMETType01XYCorUnclusteredUp_PFSumETType1CorUnclusteredUp_Token_;
  tokenMap["rootTuplePFMETType01XYCor_PFMETPhiType01XYCor"]=rootTuplePFMETType01XYCor_PFMETPhiType01XYCor_Token_;
  tokenMap["rootTuplePFMETType01XYCor_PFMETType01XYCor"]=rootTuplePFMETType01XYCor_PFMETType01XYCor_Token_;
  tokenMap["rootTuplePFMETType01XYCor_PFSumETType01XYCor"]=rootTuplePFMETType01XYCor_PFSumETType01XYCor_Token_;
  tokenMap["rootTuplePFMETType1Cor_PFMETPhiType1Cor"]=rootTuplePFMETType1Cor_PFMETPhiType1Cor_Token_;
  tokenMap["rootTuplePFMETType1Cor_PFMETType1Cor"]=rootTuplePFMETType1Cor_PFMETType1Cor_Token_;
  tokenMap["rootTuplePFMETType1Cor_PFSumETType1Cor"]=rootTuplePFMETType1Cor_PFSumETType1Cor_Token_;
  tokenMap["rootTuplePFMET_PFMETPhi"]=rootTuplePFMET_PFMETPhi_Token_;
  tokenMap["rootTuplePFMET_PFMET"]=rootTuplePFMET_PFMET_Token_;
  tokenMap["rootTuplePFMET_PFSumET"]=rootTuplePFMET_PFSumET_Token_;
  tokenMap["rootTupleVertex_VertexChi2"]=rootTupleVertex_VertexChi2_Token_;
  tokenMap["rootTupleVertex_VertexNDF"]=rootTupleVertex_VertexNDF_Token_;
  tokenMap["rootTupleVertex_VertexRho"]=rootTupleVertex_VertexRho_Token_;
  tokenMap["rootTupleVertex_VertexXErr"]=rootTupleVertex_VertexXErr_Token_;
  tokenMap["rootTupleVertex_VertexX"]=rootTupleVertex_VertexX_Token_;
  tokenMap["rootTupleVertex_VertexYErr"]=rootTupleVertex_VertexYErr_Token_;
  tokenMap["rootTupleVertex_VertexY"]=rootTupleVertex_VertexY_Token_;
  tokenMap["rootTupleVertex_VertexZErr"]=rootTupleVertex_VertexZErr_Token_;
  tokenMap["rootTupleVertex_VertexZ"]=rootTupleVertex_VertexZ_Token_;
  tokenMap["rootTupleElectrons_ElectronRhoIsoHEEP"]=rootTupleElectrons_ElectronRhoIsoHEEP_Token_;
  tokenMap["rootTupleGenEventInfo_PileUpInteractionsTrue"]=rootTupleGenEventInfo_PileUpInteractionsTrue_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjEta"]=rootTupleTriggerObjects_HLTriggerObjEta_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjPhi"]=rootTupleTriggerObjects_HLTriggerObjPhi_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjPt"]=rootTupleTriggerObjects_HLTriggerObjPt_Token_;
  tokenMap["rootTupleElectrons_ElectronCharge"]=rootTupleElectrons_ElectronCharge_Token_;
  tokenMap["rootTupleElectrons_ElectronClassif"]=rootTupleElectrons_ElectronClassif_Token_;
  tokenMap["rootTupleElectrons_ElectronMissingHitsEG"]=rootTupleElectrons_ElectronMissingHitsEG_Token_;
  tokenMap["rootTupleElectrons_ElectronMissingHits"]=rootTupleElectrons_ElectronMissingHits_Token_;
  tokenMap["rootTupleElectrons_ElectronNumberOfBrems"]=rootTupleElectrons_ElectronNumberOfBrems_Token_;
  tokenMap["rootTupleElectrons_ElectronOverlaps"]=rootTupleElectrons_ElectronOverlaps_Token_;
  tokenMap["rootTupleElectrons_ElectronPassId"]=rootTupleElectrons_ElectronPassId_Token_;
  tokenMap["rootTupleElectrons_ElectronPassIsoPAT"]=rootTupleElectrons_ElectronPassIsoPAT_Token_;
  tokenMap["rootTupleElectrons_ElectronVtxIndex"]=rootTupleElectrons_ElectronVtxIndex_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronMotherIndex"]=rootTupleGenElectronsFromWs_GenWElectronMotherIndex_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronNumDaught"]=rootTupleGenElectronsFromWs_GenWElectronNumDaught_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronPdgId"]=rootTupleGenElectronsFromWs_GenWElectronPdgId_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronStatus"]=rootTupleGenElectronsFromWs_GenWElectronStatus_Token_;
  tokenMap["rootTupleGenElectronsFromWs_GenWElectronTauDecayMode"]=rootTupleGenElectronsFromWs_GenWElectronTauDecayMode_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronMotherIndex"]=rootTupleGenElectronsFromZs_GenZElectronMotherIndex_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronNumDaught"]=rootTupleGenElectronsFromZs_GenZElectronNumDaught_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronPdgId"]=rootTupleGenElectronsFromZs_GenZElectronPdgId_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronStatus"]=rootTupleGenElectronsFromZs_GenZElectronStatus_Token_;
  tokenMap["rootTupleGenElectronsFromZs_GenZElectronTauDecayMode"]=rootTupleGenElectronsFromZs_GenZElectronTauDecayMode_Token_;
  tokenMap["rootTupleGenEventInfo_PileUpInteractions"]=rootTupleGenEventInfo_PileUpInteractions_Token_;
  tokenMap["rootTupleGenEventInfo_PileUpOriginBX"]=rootTupleGenEventInfo_PileUpOriginBX_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuMotherIndex"]=rootTupleGenMuonsFromWs_GenWMuMotherIndex_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuNumDaught"]=rootTupleGenMuonsFromWs_GenWMuNumDaught_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuPdgId"]=rootTupleGenMuonsFromWs_GenWMuPdgId_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuStatus"]=rootTupleGenMuonsFromWs_GenWMuStatus_Token_;
  tokenMap["rootTupleGenMuonsFromWs_GenWMuTauDecayMode"]=rootTupleGenMuonsFromWs_GenWMuTauDecayMode_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuMotherIndex"]=rootTupleGenMuonsFromZs_GenZMuMotherIndex_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuNumDaught"]=rootTupleGenMuonsFromZs_GenZMuNumDaught_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuPdgId"]=rootTupleGenMuonsFromZs_GenZMuPdgId_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuStatus"]=rootTupleGenMuonsFromZs_GenZMuStatus_Token_;
  tokenMap["rootTupleGenMuonsFromZs_GenZMuTauDecayMode"]=rootTupleGenMuonsFromZs_GenZMuTauDecayMode_Token_;
  tokenMap["rootTupleGenParticles_GenParticleMotherIndex"]=rootTupleGenParticles_GenParticleMotherIndex_Token_;
  tokenMap["rootTupleGenParticles_GenParticleNumDaught"]=rootTupleGenParticles_GenParticleNumDaught_Token_;
  tokenMap["rootTupleGenParticles_GenParticlePdgId"]=rootTupleGenParticles_GenParticlePdgId_Token_;
  tokenMap["rootTupleGenParticles_GenParticleStatus"]=rootTupleGenParticles_GenParticleStatus_Token_;
  tokenMap["rootTupleGenParticles_GenParticleTauDecayMode"]=rootTupleGenParticles_GenParticleTauDecayMode_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauMotherIndex"]=rootTupleGenTausFromWs_GenWTauMotherIndex_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauNumDaught"]=rootTupleGenTausFromWs_GenWTauNumDaught_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauPdgId"]=rootTupleGenTausFromWs_GenWTauPdgId_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauStatus"]=rootTupleGenTausFromWs_GenWTauStatus_Token_;
  tokenMap["rootTupleGenTausFromWs_GenWTauTauDecayMode"]=rootTupleGenTausFromWs_GenWTauTauDecayMode_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauMotherIndex"]=rootTupleGenTausFromZs_GenZTauMotherIndex_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauNumDaught"]=rootTupleGenTausFromZs_GenZTauNumDaught_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauPdgId"]=rootTupleGenTausFromZs_GenZTauPdgId_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauStatus"]=rootTupleGenTausFromZs_GenZTauStatus_Token_;
  tokenMap["rootTupleGenTausFromZs_GenZTauTauDecayMode"]=rootTupleGenTausFromZs_GenZTauTauDecayMode_Token_;
  tokenMap["rootTupleMuons_MuonBestTrackVtxIndex"]=rootTupleMuons_MuonBestTrackVtxIndex_Token_;
  tokenMap["rootTupleMuons_MuonCharge"]=rootTupleMuons_MuonCharge_Token_;
  tokenMap["rootTupleMuons_MuonCocktailCharge"]=rootTupleMuons_MuonCocktailCharge_Token_;
  tokenMap["rootTupleMuons_MuonCocktailRefitID"]=rootTupleMuons_MuonCocktailRefitID_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkHits"]=rootTupleMuons_MuonCocktailTrkHits_Token_;
  tokenMap["rootTupleMuons_MuonCocktailTrkVtxIndex"]=rootTupleMuons_MuonCocktailTrkVtxIndex_Token_;
  tokenMap["rootTupleMuons_MuonGlobalTrkValidHits"]=rootTupleMuons_MuonGlobalTrkValidHits_Token_;
  tokenMap["rootTupleMuons_MuonIsGlobal"]=rootTupleMuons_MuonIsGlobal_Token_;
  tokenMap["rootTupleMuons_MuonIsPF"]=rootTupleMuons_MuonIsPF_Token_;
  tokenMap["rootTupleMuons_MuonIsTracker"]=rootTupleMuons_MuonIsTracker_Token_;
  tokenMap["rootTupleMuons_MuonPassID"]=rootTupleMuons_MuonPassID_Token_;
  tokenMap["rootTupleMuons_MuonPixelHits"]=rootTupleMuons_MuonPixelHits_Token_;
  tokenMap["rootTupleMuons_MuonSegmentMatches"]=rootTupleMuons_MuonSegmentMatches_Token_;
  tokenMap["rootTupleMuons_MuonStationMatches"]=rootTupleMuons_MuonStationMatches_Token_;
  tokenMap["rootTupleMuons_MuonTrackLayersWithMeasurement"]=rootTupleMuons_MuonTrackLayersWithMeasurement_Token_;
  tokenMap["rootTupleMuons_MuonTrkHitsTrackerOnly"]=rootTupleMuons_MuonTrkHitsTrackerOnly_Token_;
  tokenMap["rootTupleMuons_MuonTrkHits"]=rootTupleMuons_MuonTrkHits_Token_;
  tokenMap["rootTupleMuons_MuonTrkPixelHits"]=rootTupleMuons_MuonTrkPixelHits_Token_;
  tokenMap["rootTupleMuons_MuonVtxIndex"]=rootTupleMuons_MuonVtxIndex_Token_;
  tokenMap["rootTuplePFCandidates_PFCandChargeLeptLink"]=rootTuplePFCandidates_PFCandChargeLeptLink_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetBestVertexTrackAssociationIndexAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetBestVertexTrackAssociationIndexAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetChargedHadronMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetChargedHadronMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetChargedMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetChargedMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertex3DIndexAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertex3DIndexAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertexXYIndexAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertexXYIndexAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetClosestVertexZIndexAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetClosestVertexZIndexAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetElectronMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetElectronMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetHFEMMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetHFEMMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetHFHadronMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetHFHadronMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetMuonMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetMuonMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetNConstituentsAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetNConstituentsAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetNeutralHadronMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetNeutralHadronMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetNeutralMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetNeutralMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPartonFlavourAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPartonFlavourAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPassLooseIDAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPassLooseIDAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPassTightIDAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPassTightIDAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4CHS_PFJetPhotonMultiplicityAK4CHS"]=rootTuplePFJetsAK4CHS_PFJetPhotonMultiplicityAK4CHS_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetBestVertexTrackAssociationIndexAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetBestVertexTrackAssociationIndexAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetChargedHadronMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetChargedHadronMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetChargedMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetChargedMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertex3DIndexAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertex3DIndexAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertexXYIndexAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertexXYIndexAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetClosestVertexZIndexAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetClosestVertexZIndexAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetElectronMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetElectronMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetHFEMMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetHFEMMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetHFHadronMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetHFHadronMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetMuonMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetMuonMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetNConstituentsAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetNConstituentsAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetNeutralHadronMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetNeutralHadronMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetNeutralMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetNeutralMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPartonFlavourAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPartonFlavourAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPassLooseIDAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPassLooseIDAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPassTightIDAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPassTightIDAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK4Puppi_PFJetPhotonMultiplicityAK4Puppi"]=rootTuplePFJetsAK4Puppi_PFJetPhotonMultiplicityAK4Puppi_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetBestVertexTrackAssociationIndexAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetBestVertexTrackAssociationIndexAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetChargedHadronMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetChargedHadronMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetChargedMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetChargedMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertex3DIndexAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertex3DIndexAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertexXYIndexAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertexXYIndexAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetClosestVertexZIndexAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetClosestVertexZIndexAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetElectronMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetElectronMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetHFEMMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetHFEMMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetHFHadronMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetHFHadronMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetMuonMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetMuonMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetNConstituentsAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetNConstituentsAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetNeutralHadronMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetNeutralHadronMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetNeutralMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetNeutralMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPartonFlavourAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPartonFlavourAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPassLooseIDAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPassLooseIDAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPassTightIDAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPassTightIDAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5CHS_PFJetPhotonMultiplicityAK5CHS"]=rootTuplePFJetsAK5CHS_PFJetPhotonMultiplicityAK5CHS_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetBestVertexTrackAssociationIndexAK5"]=rootTuplePFJetsAK5_PFJetBestVertexTrackAssociationIndexAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetChargedHadronMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetChargedHadronMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetChargedMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetChargedMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertex3DIndexAK5"]=rootTuplePFJetsAK5_PFJetClosestVertex3DIndexAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertexXYIndexAK5"]=rootTuplePFJetsAK5_PFJetClosestVertexXYIndexAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetClosestVertexZIndexAK5"]=rootTuplePFJetsAK5_PFJetClosestVertexZIndexAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetElectronMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetElectronMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetHFEMMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetHFEMMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetHFHadronMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetHFHadronMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetMuonMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetMuonMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetNConstituentsAK5"]=rootTuplePFJetsAK5_PFJetNConstituentsAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetNeutralHadronMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetNeutralHadronMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetNeutralMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetNeutralMultiplicityAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPartonFlavourAK5"]=rootTuplePFJetsAK5_PFJetPartonFlavourAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPassLooseIDAK5"]=rootTuplePFJetsAK5_PFJetPassLooseIDAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPassTightIDAK5"]=rootTuplePFJetsAK5_PFJetPassTightIDAK5_Token_;
  tokenMap["rootTuplePFJetsAK5_PFJetPhotonMultiplicityAK5"]=rootTuplePFJetsAK5_PFJetPhotonMultiplicityAK5_Token_;
  tokenMap["rootTupleVertex_VertexNTracksW05"]=rootTupleVertex_VertexNTracksW05_Token_;
  tokenMap["rootTupleVertex_VertexNTracks"]=rootTupleVertex_VertexNTracks_Token_;
  tokenMap["rootTupleTriggerObjects_HLTriggerObjTypeIds"]=rootTupleTriggerObjects_HLTriggerObjTypeIds_Token_;
  tokenMap["rootTupleEvent_bunch"]=rootTupleEvent_bunch_Token_;
  tokenMap["rootTupleEvent_event"]=rootTupleEvent_event_Token_;
  tokenMap["rootTupleEvent_ls"]=rootTupleEvent_ls_Token_;
  tokenMap["rootTupleEvent_orbit"]=rootTupleEvent_orbit_Token_;
  tokenMap["rootTupleEvent_run"]=rootTupleEvent_run_Token_;
  tokenMap["rootTupleGenEventInfo_ProcessID"]=rootTupleGenEventInfo_ProcessID_Token_;
}
