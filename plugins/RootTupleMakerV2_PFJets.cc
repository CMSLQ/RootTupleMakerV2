
#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_PFJets.h"


RootTupleMakerV2_PFJets::RootTupleMakerV2_PFJets(const edm::ParameterSet& iConfig) :
  jetInputToken_ (consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("InputTag"))),
  inputTag (iConfig.getParameter<edm::InputTag>("InputTag")),//fixme using old way to find PUPPI, is this a problem?
  prefix  (iConfig.getParameter<std::string>  ("Prefix")),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  mvaPileupIDname  (iConfig.getParameter<std::string>  ("MVAPileupIDName")),
  maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
  jecUncPath (iConfig.getParameter<std::string>("JECUncertainty")),
  readJECuncertainty (iConfig.getParameter<bool>   ("ReadJECuncertainty")),
  readJERuncertainty (iConfig.getParameter<bool>   ("ReadJERuncertainty")),
  //
  vtxInputToken_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexInputTag"))),
  jerUncPath (iConfig.getParameter<std::string>("JERUncertainty")),
  jer_from_gt (iConfig.getParameter<bool>      ("ReadJERFromGT")),
  jer_resolutions_file (iConfig.getParameter<edm::FileInPath>  ("JERResolutionsFile")),
  jer_scale_factors_file (iConfig.getParameter<edm::FileInPath>("JERScaleFactorsFile")),
  rhoToken (consumes<double>(iConfig.getParameter<edm::InputTag>("RhoCollection")))
{
  if(upperCase(inputTag.label()).find(upperCase("PUPPI")) != std::string::npos)
    isPuppiJetColl = true;
  else
    isPuppiJetColl = false;
  produces <bool>                 ( prefix + "hasJetWithBadUnc" + suffix );
  produces <std::vector<float> >  ( prefix + "Eta" + suffix );
  produces <std::vector<float> >  ( prefix + "Phi" + suffix );
  produces <std::vector<float> >  ( prefix + "Pt" + suffix );
  produces <std::vector<float> >  ( prefix + "Mt" + suffix );
  produces <std::vector<float> >  ( prefix + "SmearedUpPt" + suffix );
  produces <std::vector<float> >  ( prefix + "SmearedDownPt" + suffix );
  produces <std::vector<float> >  ( prefix + "ScaledUpPt" + suffix );
  produces <std::vector<float> >  ( prefix + "ScaledDownPt" + suffix );
  produces <std::vector<float> >  ( prefix + "PtRaw" + suffix );
  produces <std::vector<float> >  ( prefix + "Energy" + suffix );
  produces <std::vector<float> >  ( prefix + "SmearedUpEnergy" + suffix );
  produces <std::vector<float> >  ( prefix + "SmearedDownEnergy" + suffix );
  produces <std::vector<float> >  ( prefix + "ScaledUpEnergy" + suffix );
  produces <std::vector<float> >  ( prefix + "ScaledDownEnergy" + suffix );
  produces <std::vector<float> >  ( prefix + "EnergyRaw" + suffix );
  produces <std::vector<float> >  ( prefix + "JECUnc" + suffix );
  produces <std::vector<float> >  ( prefix + "L2L3ResJEC" + suffix );
  produces <std::vector<float> >  ( prefix + "L3AbsJEC" + suffix );
  produces <std::vector<float> >  ( prefix + "L2RelJEC" + suffix );
  produces <std::vector<float> >  ( prefix + "L1FastJetJEC" + suffix );
  produces <std::vector<float> >  ( prefix + "JERRes" + suffix );
  produces <std::vector<float> >  ( prefix + "JERResSF" + suffix );
  produces <std::vector<float> >  ( prefix + "JERResSFUp" + suffix );
  produces <std::vector<float> >  ( prefix + "JERResSFDown" + suffix );
  produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
  produces <std::vector<int> >    ( prefix + "HadronFlavour" + suffix );
  produces <std::vector<float> >  ( prefix + "ChargedEmEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "ChargedHadronEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "ChargedMuEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "ElectronEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "MuonEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "NeutralEmEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "NeutralHadronEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "PhotonEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "HFHadronEnergyFraction"  + suffix );
  produces <std::vector<float> >  ( prefix + "HFEMEnergyFraction"  + suffix );
  produces <std::vector<int> >    ( prefix + "ChargedHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "ChargedMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "ElectronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "MuonMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "NeutralHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "NeutralMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "PhotonMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "HFHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "HFEMMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "NConstituents"  + suffix );
  //produces <std::vector<float> >  ( prefix + "TrackCountingHighEffBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "TrackCountingHighPurBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "JetProbabilityBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "JetBProbabilityBTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "CombinedSecondaryVertexBTag" + suffix );    
  //produces <std::vector<float> >  ( prefix + "CombinedSecondaryVertexMVABTag" + suffix ); 
  //produces <std::vector<float> >  ( prefix + "SoftElectronByPtBTag" + suffix );           
  //produces <std::vector<float> >  ( prefix + "SoftElectronByIP3dBTag" + suffix );         
  //produces <std::vector<float> >  ( prefix + "SoftMuonBTag" + suffix );                   
  //produces <std::vector<float> >  ( prefix + "SoftMuonByPtBTag" + suffix );               
  //produces <std::vector<float> >  ( prefix + "SoftMuonByIP3dBTag" + suffix );             
  produces <std::vector<float> >  ( prefix + "CombinedInclusiveSecondaryVertexBTag" + suffix );
  produces <std::vector<float> >  ( prefix + "CombinedMVABTag" + suffix );
  //produces <std::vector<float> >  ( prefix + "CombinedCvsLJetTags" + suffix );
  //produces <std::vector<float> >  ( prefix + "CombinedCvsBJetTags" + suffix );
  produces <std::vector<int> >    ( prefix + "PassLooseID" + suffix);
  produces <std::vector<int> >    ( prefix + "PassTightID" + suffix);
  // for non-PUPPI jets, we get the MVA pileup ID discriminator
  if(!isPuppiJetColl){
    produces <std::vector<float> >   ( prefix + "PileupMVA" + suffix);
    produces <std::vector<int> >     ( prefix + "PileupMVApassesLoose"  + suffix);
    produces <std::vector<int> >     ( prefix + "PileupMVApassesMedium" + suffix);
    produces <std::vector<int> >     ( prefix + "PileupMVApassesTight"  + suffix);
    produces <std::vector<float> >   (prefix + "VtxNtracks"  + suffix);
    produces <std::vector<float> >   (prefix + "VtxMass"  + suffix);
    produces <std::vector<float> >   (prefix + "VtxPx"  + suffix);
    produces <std::vector<float> >   (prefix + "VtxPy"  + suffix);
    produces <std::vector<float> >   (prefix + "Vtx3DVal"  + suffix);
    produces <std::vector<float> >   (prefix + "Vtx3DSig"  + suffix);
    produces <std::vector<float> >   (prefix + "LeadTrackPt"  + suffix);
    produces <std::vector<float> >   (prefix + "SoftLepPt"  + suffix);
    produces <std::vector<float> >   (prefix + "SoftLepRatio"  + suffix);
    produces <std::vector<float> >   (prefix + "SoftLepDr"  + suffix);
  }
  produces <std::vector<float> >  ( prefix + "BestVertexTrackAssociationFactor" + suffix );
  produces <std::vector<int> >    ( prefix + "BestVertexTrackAssociationIndex" + suffix);
  produces <std::vector<float> >  ( prefix + "ClosestVertexWeighted3DSeparation" + suffix );
  produces <std::vector<float> >  ( prefix + "ClosestVertexWeightedXYSeparation" + suffix );
  produces <std::vector<float> >  ( prefix + "ClosestVertexWeightedZSeparation" + suffix );
  produces <std::vector<int> >    ( prefix + "ClosestVertex3DIndex" + suffix);
  produces <std::vector<int> >    ( prefix + "ClosestVertexXYIndex" + suffix);
  produces <std::vector<int> >    ( prefix + "ClosestVertexZIndex" + suffix);
  produces <std::vector<float> >  ( prefix + "Beta" + suffix ) ;
  produces <std::vector<float> >  ( prefix + "BetaStar" + suffix ) ;
  produces <std::vector<float> >  ( prefix + "BetaClassic" + suffix ) ;
  produces <std::vector<float> >  ( prefix + "BetaStarClassic" + suffix ) ;
}


PFJetIDSelectionFunctor pfjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
PFJetIDSelectionFunctor pfjetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );

pat::strbitset retpf = pfjetIDLoose.getBitTemplate();


std::string RootTupleMakerV2_PFJets::upperCase(std::string input)
{
  for (std::string::iterator it = input.begin(); it != input.end(); ++ it)
    *it = toupper(*it);
  return input;
}


void RootTupleMakerV2_PFJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::auto_ptr<bool>                  hasJetWithBadUnc ( new bool() );
  std::auto_ptr<std::vector<float> >   eta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   phi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   pt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   mt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   ptSmearedUp  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   ptSmearedDown  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   ptScaledUp  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   ptScaledDown  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energySmearedUp  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energySmearedDown  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energyScaledUp  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energyScaledDown  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   pt_raw  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   energy_raw ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   jecUnc_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   l2l3resJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   l3absJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   l2relJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   l1fastjetJEC_vec ( new std::vector<float>()  );
  // std::auto_ptr<std::vector<float> >   l1offsetJEC_vec ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   jerRes_vec ( new std::vector<float>() );
  std::auto_ptr<std::vector<float> >   jerResSF_vec ( new std::vector<float>() );
  std::auto_ptr<std::vector<float> >   jerResSFup_vec ( new std::vector<float>() );
  std::auto_ptr<std::vector<float> >   jerResSFdown_vec ( new std::vector<float>() );
  std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     hadronFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<float> >   chargedEmEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   chargedHadronEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   chargedMuEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   electronEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   muonEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   neutralEmEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   neutralHadronEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   photonEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   hfHadronEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   hfEMEnergyFraction  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedHadronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     electronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     muonMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralHadronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     photonMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     hfHadronMultiplicity ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     hfEMMultiplicity ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     nConstituents  ( new std::vector<int>()  ) ;
  //std::auto_ptr<std::vector<float> >   trackCountingHighEffBTag  ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   trackCountingHighPurBTag  ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   simpleSecondaryVertexHighEffBTag  ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   simpleSecondaryVertexHighPurBTag  ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   jetProbabilityBTag  ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   jetBProbabilityBTag  ( new std::vector<float>()  );

  //std::auto_ptr<std::vector<float> >   combinedSecondaryVertexBTag          ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   combinedSecondaryVertexMVABTag       ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   softElectronByPtBTag                 ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   softElectronByIP3dBTag               ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   softMuonBTag                         ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   softMuonByPtBTag                     ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   softMuonByIP3dBTag                   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   combinedInclusiveSecondaryVertexBTag ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >   combinedMVABTag                      ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   combinedCvsLJetTag                   ( new std::vector<float>()  );
  //std::auto_ptr<std::vector<float> >   combinedCvsBJetTag                   ( new std::vector<float>()  );
	
  std::auto_ptr<std::vector<int> >   passLooseID  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >   passTightID  ( new std::vector<int>()  );

  std::auto_ptr<std::vector<float> > jetpileup_mva  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >   jetpileup_mva_passesLoose  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >   jetpileup_mva_passesMedium ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >   jetpileup_mva_passesTight  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<float> > vtxNtracks  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > vtxMass  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > vtxPx  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > vtxPy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > vtx3DVal  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > vtx3DSig  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > leadTrackPt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > softLepPt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > softLepRatio  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> > softLepDr  ( new std::vector<float>()  );


  std::auto_ptr <std::vector<float> >   bestVertexTrackAssociationFactor  ( new std::vector<float>()  );
  std::auto_ptr <std::vector<int> >     bestVertexTrackAssociationIndex   ( new std::vector<int>()  );
  std::auto_ptr <std::vector<float> >   closestVertexWeighted3DSeparation  ( new std::vector<float>()  );
  std::auto_ptr <std::vector<float> >   closestVertexWeightedXYSeparation  ( new std::vector<float>()  );
  std::auto_ptr <std::vector<float> >   closestVertexWeightedZSeparation  ( new std::vector<float>()  );
  std::auto_ptr <std::vector<int> >     closestVertex3DIndex            ( new std::vector<int>()  );
  std::auto_ptr <std::vector<int> >     closestVertexXYIndex           ( new std::vector<int>()  );
  std::auto_ptr <std::vector<int> >     closestVertexZIndex            ( new std::vector<int>()  );

  std::auto_ptr <std::vector<float > > betaStar        ( new std::vector<float>());
  std::auto_ptr <std::vector<float > > betaStarClassic ( new std::vector<float>());
  std::auto_ptr <std::vector<float > > beta            ( new std::vector<float>());
  std::auto_ptr <std::vector<float > > betaClassic     ( new std::vector<float>());
	
  //-----------------------------------------------------------------

  //JEC Uncertainties
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
  *hasJetWithBadUnc.get() = false;
  JetCorrectionUncertainty *jecUnc = 0;
  if(readJECuncertainty)
    {
      // handle the jet corrector parameters collection
      edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
      // get the jet corrector parameters collection from the global tag
      iSetup.get<JetCorrectionsRecord>().get(jecUncPath,JetCorParColl);
      // get the uncertainty parameters from the collection
      JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
      // instantiate the jec uncertainty object
      jecUnc = new JetCorrectionUncertainty(JetCorPar);
    }
  //JER Uncertainties
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
  // Access jet resolution and scale factor from the condition database
  // or from text files
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if (readJERuncertainty)
  {
    // Two differents way to create a class instance
    if (jer_from_gt)
    {
      // First way, using the get() static method
      resolution = JME::JetResolution::get(iSetup, jerUncPath);
      res_sf = JME::JetResolutionScaleFactor::get(iSetup, jecUncPath); //same label as jecUncPath
    }
    else
    {
      // Second way, using the constructor
      resolution = JME::JetResolution(jer_resolutions_file.fullPath());
      res_sf = JME::JetResolutionScaleFactor(jer_scale_factors_file.fullPath());
      edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Reading scale factors from: " << jer_scale_factors_file;
    }// get the factors from the global tag or from test file
  }
  // need rho for JER
  edm::Handle<double> rho;
  iEvent.getByToken(rhoToken, rho);
  
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByToken(jetInputToken_, jets);
	
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vtxInputToken_,primaryVertices);

  if(jets.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # PFJets: " << jets->size();

      int ijet = -1;
      for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it )
	{
	  ijet++;
	  // Only look at jets with pt>=17 GeV
	  if( it->pt()<17 ) 
	    continue;

	  // exit from loop when you reach the required number of jets
	  if(eta->size() >= maxSize)
	    break;

	  retpf.set(false);
	  int passjetLoose =0;
	  if(pfjetIDLoose( *it, retpf )) passjetLoose =1;

	  retpf.set(false);
	  int passjetTight = 0;
	  if (pfjetIDTight( *it, retpf)) passjetTight =1;


	  if(readJECuncertainty)
	    {
	      jecUnc->setJetEta( it->eta() );
	      // the uncertainty is a function of the corrected pt
	      jecUnc->setJetPt( it->pt() );
	    }

	  // Status of JEC
	  //std::cout << "PF: currentJECLevel(): " << it->currentJECLevel() << std::endl;
	  //std::cout << "PF: currentJECSet(): " << it->currentJECSet() << std::endl;
	  //-------------------
	  
	  // Vertex association

	  int bestVtxIndex3Ddist = -1;
	  int bestVtxIndexXYdist = -1;
	  int bestVtxIndexZdist = -1;

	  int bestVtxIndexSharedTracks = -1;

	  float minVtxDist3D = 999999.;
	  float minVtxDistXY = -99999.;
	  float minVtxDistZ  = -99999.;
	  float maxTrackAssocRatio = -9999.;

	  // Loop on primary Vertices and jets and perform associations 

	  reco::VertexCollection::const_iterator lead_vertex = primaryVertices->end();
	  bool found_lead_vertex = false;

	  if(primaryVertices.isValid())
	    {
	      edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

	      // Main Vertex Loop
	      for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it )
		{

		  float sumweights = 0.0;
		  float dist3Dweighted = 0.0;
		  float distXYweighted = 0.0;
		  float distZweighted = 0.0;
		  float assocsumpttracks = 0.0;
		  float trackassociationratio = 0.000001;

		  if ( v_it -> isFake() || v_it -> ndof() < 4 ) continue;
		  if ( !found_lead_vertex ) { 
		    lead_vertex = v_it;
		    found_lead_vertex = true;
		  }

		  // Loop on tracks in jet, calculate PT weighted 3D distance to vertex and PT weighted shared track ratio
		  const reco::TrackRefVector &jtracks=it->associatedTracks();
		  for(reco::TrackRefVector::const_iterator jtIt=jtracks.begin(); jtIt!=jtracks.end(); ++jtIt)
		    {
		      if( jtIt->isNull() ) continue;
		      const reco::Track *jtrack=jtIt->get();
		      float trackptweight = jtrack->pt();
		      sumweights += trackptweight;

		      // Weighted Distance Calculation
		      float distXY= jtrack->dxy(v_it->position());
		      float distZ = jtrack->dz(v_it->position());
		      dist3Dweighted = trackptweight*(sqrt(pow(distXY,2) + pow(distZ,2)));
		      distXYweighted = trackptweight*distXY;
		      distZweighted = trackptweight*distZ;


		      // Loop on vertex tracks, find PT weighted shared tracks. 
		      for(reco::Vertex::trackRef_iterator vtIt=v_it->tracks_begin(); vtIt!=v_it->tracks_end(); ++vtIt)
			{
			  if( vtIt->isNull() ) continue;
			  const reco::Track *vtrack=vtIt->get();
			  if(vtrack!=jtrack) continue;
			  assocsumpttracks+=jtrack->pt();
			  break;
			}

		      trackassociationratio = assocsumpttracks/sumweights;

		    }

		  // Divide distances by sum of weights. 
		  dist3Dweighted = dist3Dweighted / sumweights;
		  distXYweighted = distXYweighted / sumweights;
		  distZweighted  = distZweighted  / sumweights;	

		  // Find vertex with minimum weighted distance. 
		  if( dist3Dweighted < minVtxDist3D )
		    {
		      minVtxDist3D = dist3Dweighted;
		      bestVtxIndex3Ddist = int(std::distance(primaryVertices->begin(),v_it));

		    }

		  if( distXYweighted < minVtxDistXY )
		    {
		      minVtxDistXY = distXYweighted;
		      bestVtxIndexXYdist = int(std::distance(primaryVertices->begin(),v_it));
		    }

		  if( distZweighted < minVtxDistZ )
		    {
		      minVtxDistZ = distZweighted;
		      bestVtxIndexZdist = int(std::distance(primaryVertices->begin(),v_it));
		    }

		  // Find vertex with minimum weighted distance. 
		  if( trackassociationratio > maxTrackAssocRatio )
		    {
		      maxTrackAssocRatio = trackassociationratio ;
		      bestVtxIndexSharedTracks = int(std::distance(primaryVertices->begin(),v_it));
		    }						

		  //std::cout<<dist3Dweighted<<"  "<<distXYweighted<<"  "<<distZweighted<<"  "<<trackassociationratio<<"  "<<int(std::distance(primaryVertices->begin(),v_it))<<std::endl;


		}
	      //std::cout<<"---------------------"<<std::endl;
	    }	
	  else
	    {
	      edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product primary vertices";
	    }

	  // Get the constituents of the PFJet
	  // On MiniAOD, we must use daughters instead

	  //std::vector <reco::PFCandidatePtr> constituents  = it -> getPFConstituents();
	  //std::vector <reco::PFCandidatePtr>::iterator i_constituent   = constituents.begin();
	  //std::vector <reco::PFCandidatePtr>::iterator end_constituent = constituents.end();
	  int numberOfDaughters = it->numberOfDaughters();
	  float sum_track_pt = 0.;

	  float jetBetaStar        = 0.0 ;
	  float jetBetaStarClassic = 0.0 ;
	  float jetBeta            = 0.0 ;
	  float jetBetaClassic     = 0.0 ;

	  if ( found_lead_vertex ) { 
	    for(int dau = 0; dau < numberOfDaughters; ++dau) {
	      edm::Ptr<pat::PackedCandidate> constituent(it->daughterPtr(dau));

	      try { 
		// constituent->pseudoTrack().pt() in MiniAOD (best we can have) is the same as constituent pt
		//   so no need to make the pseudoTrack here
		float track_pt = constituent->pt();
		sum_track_pt += track_pt;

		// If it's used in the fit of the primary vertex or associated to it, take it
		// See: https://hypernews.cern.ch/HyperNews/CMS/get/csa14/85/1/1/1.html
		// Note that for CHS jets, all daughters will be at least PV Loose (i.e., not from a "pileup" vertex)
		bool track_from_lead_vertex = (constituent->fromPV() == pat::PackedCandidate::PVUsedInFit) ||
		  (constituent->fromPV() == pat::PackedCandidate::PVTight);

		bool track_from_other_vertex = false;

		float dZ0 = fabs(constituent->dz()); // constituent dz is only kept at float precision (MiniAOD trick)
		float dZ = dZ0; 

		// look at other vertices in PV collection, not primary ( pv[0] )
		for(reco::VertexCollection::const_iterator v_it=primaryVertices->begin()+1; v_it!=primaryVertices->end(); ++v_it )
		  {
		    if( v_it -> isFake() || v_it -> ndof() < 4 ) continue;
		    // is track from a non-PV non-pileup vertex in the PV collection?
		    if(!track_from_other_vertex)
		      track_from_other_vertex = constituent->fromPV() == pat::PackedCandidate::PVLoose;
		    dZ = std::min(dZ, constituent->dz(v_it -> position()));
		  }

		if      (  track_from_lead_vertex && !track_from_other_vertex ) jetBetaClassic     += track_pt;
		else if ( !track_from_lead_vertex &&  track_from_other_vertex ) jetBetaStarClassic += track_pt;

		if      ( dZ0 < 0.2 ) jetBeta     += track_pt;
		else if ( dZ  < 0.2 ) jetBetaStar += track_pt;

	      }

	      catch (cms::Exception & e) { std::cout << e << std::endl; } 

	    }
	  }

	  if ( sum_track_pt != 0. ) { 
	    jetBetaStar        /= sum_track_pt ;
	    jetBetaStarClassic /= sum_track_pt ;
	    jetBeta            /= sum_track_pt ;
	    jetBetaClassic     /= sum_track_pt ;
	  } 
	  else { 
	    assert ( jetBetaStar        == 0.0 );
	    assert ( jetBetaStarClassic == 0.0 );
	    assert ( jetBeta            == 0.0 );
	    assert ( jetBetaClassic     == 0.0 );
	  }

	  betaStar        -> push_back ( jetBetaStar        ) ;
	  betaStarClassic -> push_back ( jetBetaStarClassic ) ;
	  beta            -> push_back ( jetBeta            ) ;
	  betaClassic     -> push_back ( jetBetaClassic     ) ;

	  bestVertexTrackAssociationFactor ->push_back( maxTrackAssocRatio );
	  bestVertexTrackAssociationIndex ->push_back( bestVtxIndexSharedTracks );
	  closestVertexWeighted3DSeparation ->push_back( minVtxDist3D);
	  closestVertexWeightedXYSeparation ->push_back( minVtxDistXY );
	  closestVertexWeightedZSeparation ->push_back( minVtxDistZ);
	  closestVertex3DIndex ->push_back( bestVtxIndex3Ddist);
	  closestVertexXYIndex ->push_back( bestVtxIndexXYdist);
	  closestVertexZIndex ->push_back( bestVtxIndexZdist);

	  eta->push_back( it->eta() );
	  phi->push_back( it->phi() );
	  pt->push_back( it->pt() );
	  mt->push_back( it->mt() );
	  
	  // JER
	  if(readJERuncertainty)
	    {
	      // JER
	      // create a JetParameters object
	      // You *must* know in advance which variables are needed for getting the resolution. For the moment, only pt and eta are needed, but this will
	      // probably change in the future when PU dependency is added. Please keep an eye on the twiki page:
	      //     https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
	      JME::JetParameters jetParameters;
	      jetParameters.setJetPt(it->pt());
	      jetParameters.setJetEta(it->eta());
	      jetParameters.setRho(*rho);
	      // get resolution
	      float jetRes = resolution.getResolution(jetParameters);
	      // We do the same thing to access the scale factor
	      float jetResScaleFactor = res_sf.getScaleFactor(jetParameters);
	      // Access up and down variation of the scale factor
	      float jetResScaleFactorUp = res_sf.getScaleFactor(jetParameters, Variation::UP);
	      float jetResScaleFactorDown = res_sf.getScaleFactor(jetParameters, Variation::DOWN);
	      
	      //// JER smearing: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
	      // requires matching to gen jet, so just store the info needed to do the smearing here
	      //if ( !iEvent.isRealData() ) { 
	      //  // make the collections
	      //}
	      //else { 
	      //  ptSmearedUp       -> push_back ( it -> pt()     );
	      //  energySmearedUp   -> push_back ( it -> energy() );
	      //  ptSmearedDown     -> push_back ( it -> pt()     );
	      //  energySmearedDown -> push_back ( it -> energy() );
	      //  ptScaledUp        -> push_back ( it -> pt()     );
	      //  energyScaledUp    -> push_back ( it -> energy() );
	      //  ptScaledDown      -> push_back ( it -> pt()     );
	      //  energyScaledDown  -> push_back ( it -> energy() );
	      //}
	      jerRes_vec->push_back(jetRes);
	      jerResSF_vec->push_back(jetResScaleFactor);
	      jerResSFup_vec->push_back(jetResScaleFactorUp);
	      jerResSFdown_vec->push_back(jetResScaleFactorDown);
	    }
	  else
	    {
	      jerRes_vec->push_back(-999);
	      jerResSF_vec->push_back(-999);
	      jerResSFup_vec->push_back(-999);
	      jerResSFdown_vec->push_back(-999);
	    }
	  
	  pt_raw->push_back( it->correctedJet("Uncorrected").pt() );
	  energy->push_back( it->energy() );
	  energy_raw->push_back( it->correctedJet("Uncorrected").energy() );
	  //JEC
	  //// test
	  //std::cout << "inside module with label: " << inputTag.label() << std::endl;
	  //for(std::vector<std::string>::const_iterator jecItr = it->availableJECSets().begin();
	  //    jecItr != it->availableJECSets().end(); ++jecItr)
	  //{
	  //  std::cout << "Found available JEC: " << *jecItr << "; corrections available: ";
	  //  for(std::vector<std::string>::const_iterator jecLItr = it->availableJECLevels().begin();
	  //      jecLItr != it->availableJECLevels().end(); ++jecLItr)
	  //    std::cout << *jecLItr << ",";
	  //  std::cout << std::endl;
	  //}
	  //// end test
	  if(iEvent.isRealData())
	    l2l3resJEC_vec->push_back( it->correctedJet("L2L3Residual").pt()/it->correctedJet("L3Absolute").pt() );//FIXME do we really want these as ratios to each other? Does L3Absolute e.g. contain L2Relative already? If not they are useless...
	  else
	    l2l3resJEC_vec->push_back( 1.0 );
	  l3absJEC_vec->push_back( it->correctedJet("L3Absolute").pt()/it->correctedJet("L2Relative").pt() );
	  // puppi appears to have no L1FastJet correction available for 2012D PromptReco-v3
	  if(!isPuppiJetColl)
	    {
	      l2relJEC_vec->push_back( it->correctedJet("L2Relative").pt()/it->correctedJet("L1FastJet").pt() );
	      l1fastjetJEC_vec->push_back( it->correctedJet("L1FastJet").pt()/it->correctedJet("Uncorrected").pt() );
	    }
	  if(readJECuncertainty){ 
	    float uncertainty = -999.;
	    try { 
	      uncertainty = jecUnc->getUncertainty(true);
	    } 
	    catch ( cms::Exception & e ) { 
	      edm::LogWarning("RootTupleMakerV2_PFJetsError") << "Warning! For PFJet with eta = " << it -> eta() << " caught JEC unc exception: " << e;
	      uncertainty = -999.;
	      *hasJetWithBadUnc.get() = true;
	    }
	    jecUnc_vec->push_back( uncertainty );
	  }
	  else {
	    jecUnc_vec->push_back( -999 );
	  }

	  partonFlavour->push_back( it->partonFlavour() );
	  hadronFlavour->push_back( it->hadronFlavour() );
	  chargedEmEnergyFraction->push_back( it->chargedEmEnergyFraction() );
	  chargedHadronEnergyFraction->push_back( it->chargedHadronEnergyFraction() );
	  // same as : it->chargedHadronEnergy() / it->correctedJet("Uncorrected").energy()
	  chargedMuEnergyFraction->push_back( it->chargedMuEnergyFraction() );
	  electronEnergyFraction->push_back( it->electronEnergy() / it->correctedJet("Uncorrected").energy() );
	  // 'const class pat::Jet' has no member named 'electronEnergyFraction'
	  muonEnergyFraction->push_back( it->muonEnergyFraction() );
	  neutralEmEnergyFraction->push_back( it->neutralEmEnergyFraction() );
	  neutralHadronEnergyFraction->push_back( it->neutralHadronEnergyFraction() );
	  photonEnergyFraction->push_back( it->photonEnergyFraction() );
	  hfHadronEnergyFraction->push_back( it->HFHadronEnergyFraction() );
	  hfEMEnergyFraction->push_back( it->HFEMEnergyFraction() );
	  chargedHadronMultiplicity->push_back( it->chargedHadronMultiplicity() );
	  chargedMultiplicity->push_back( it->chargedMultiplicity() );
	  electronMultiplicity->push_back( it->electronMultiplicity() );
	  muonMultiplicity->push_back( it->muonMultiplicity() );
	  neutralHadronMultiplicity->push_back( it->neutralHadronMultiplicity() );
	  neutralMultiplicity->push_back( it->neutralMultiplicity() );
	  photonMultiplicity->push_back( it->photonMultiplicity() );
	  hfHadronMultiplicity->push_back( it->HFHadronMultiplicity() );
	  hfEMMultiplicity->push_back( it->HFEMMultiplicity() );
	  nConstituents->push_back( it->numberOfDaughters() ); // same as it->getPFConstituents().size()
	  //trackCountingHighEffBTag            ->push_back( it->bDiscriminator("pfTrackCountingHighEffBJetTags") );
	  //trackCountingHighPurBTag            ->push_back( it->bDiscriminator("pfTrackCountingHighPurBJetTags") );
	  //simpleSecondaryVertexHighEffBTag    ->push_back( it->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags") );
	  //simpleSecondaryVertexHighPurBTag    ->push_back( it->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags") );
	  //jetProbabilityBTag                  ->push_back( it->bDiscriminator("pfJetProbabilityBJetTags") );
	  //jetBProbabilityBTag                 ->push_back( it->bDiscriminator("pfJetBProbabilityBJetTags") );
	  //combinedSecondaryVertexBTag         ->push_back( it->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags") );
	  combinedInclusiveSecondaryVertexBTag->push_back( it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
	  combinedMVABTag                     ->push_back( it->bDiscriminator("pfCombinedMVAV2BJetTags") );
	  //combinedCvsLJetTag                  ->push_back( it->bDiscriminator("pfCombinedCvsLJetTags") );
	  //combinedCvsBJetTag                  ->push_back( it->bDiscriminator("pfCombinedCvsBJetTags") );
	  //combinedSecondaryVertexMVABTag      ->push_back( it->bDiscriminator("combinedSecondaryVertexMVABJetTags") );
	  //softElectronByPtBTag                ->push_back( it->bDiscriminator("softElectronByPtBJetTags") );                
	  //softElectronByIP3dBTag              ->push_back( it->bDiscriminator("softElectronByIP3dBJetTags") );
	  //softMuonBTag                        ->push_back( it->bDiscriminator("softMuonBJetTags") );
	  //softMuonByPtBTag                    ->push_back( it->bDiscriminator("softMuonByPtBJetTags") );                
	  //softMuonByIP3dBTag                  ->push_back( it->bDiscriminator("softMuonByIP3dBJetTags") );
	  passLooseID->push_back( passjetLoose );
	  passTightID->push_back( passjetTight );
	  //// JET Pileup1 MVA
	  if(!isPuppiJetColl){
	    jetpileup_mva->push_back(it->userFloat(mvaPileupIDname));
	    std::string idName = "pileupJetId:fullId";
	    jetpileup_mva_passesLoose ->push_back(bool(it->userInt(idName) & (1 << 2)));
	    jetpileup_mva_passesMedium->push_back(bool(it->userInt(idName) & (1 << 1)));
	    jetpileup_mva_passesTight ->push_back(bool(it->userInt(idName) & (1 << 0)));

	    vtxNtracks->push_back(it->userFloat("vtxNtracks"));
	    vtxMass->push_back(it->userFloat("vtxMass"));
	    vtxPx->push_back(it->userFloat("vtxPx"));
	    vtxPy->push_back(it->userFloat("vtxPy"));
	    vtx3DVal->push_back(it->userFloat("vtx3DVal"));
	    vtx3DSig->push_back(it->userFloat("vtx3DSig"));


	    float leadTrackPt_ = 0, softLepPt_ = 0, softLepRatio_ = 0, softLepDr_ = 0;
	    for ( unsigned k = 0; k < it->numberOfSourceCandidatePtrs(); ++k ) {
	      reco::CandidatePtr pfJetConstituent = it->sourceCandidatePtr(k);
	      
	      const reco::Candidate* kcand = pfJetConstituent.get();
	      const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>( kcand );
	      //if ( !lPack ) throw cms::Exception( "NoPackedConstituent" ) << " For jet " << i << " failed to get constituent " << k << std::endl;
	      float candPt = kcand->pt();
	      float candDr = ROOT::Math::VectorUtil::DeltaR(kcand->p4(),it->p4());
	      
	      if(lPack->charge() != 0 && candPt > leadTrackPt_) leadTrackPt_ = candPt;
	      
	      if(abs(lPack->pdgId()) == 11 || abs(lPack->pdgId()) == 13) {
		if(candPt > softLepPt_){
		  softLepPt_ = candPt;
		  softLepRatio_ = candPt/it->pt();
		  softLepDr_ = candDr;
		}
	      }
	    }
             
	    leadTrackPt->push_back(leadTrackPt_);
	    softLepPt->push_back(softLepPt_);
	    softLepRatio->push_back(softLepRatio_);
	    softLepDr->push_back(softLepDr_);



	    //std::cout<<"Passes loose: "<<bool(it->userInt(idName) & (1 << 2))<< " medium: "<<bool(it->userInt(idName) & (1 << 1))<< " tight: "<<bool(it->userInt(idName) & (1 << 0))<<std::endl;
	  }

	  // 			//////////////////////////////////////////////////////////////////// 
	  // 			if( fabs(it->eta()) > 3) 
	  // 			  {
	  //  			    float SUM = it->chargedEmEnergyFraction() + it->chargedHadronEnergyFraction() + it->neutralEmEnergyFraction() + it->neutralHadronEnergyFraction() + it->chargedMuEnergyFraction() + it->HFHadronEnergyFraction() + it->HFEMEnergyFraction() ; 

	  //  			    std::cout << "eta,chargedEmEnergy,chargedHadronEnergy,neutralEmEnergy,neutralHadronEnergy,chargedMuEnergy,HFHadronEnergy,HFEMEnergy,SUM: "  
	  //  				      << it->eta() << " , "
	  //  				      << it->chargedEmEnergyFraction() << " , " 
	  //  				      << it->chargedHadronEnergyFraction() << " , " 
	  //  				      << it->neutralEmEnergyFraction() << " , "  
	  //  				      << it->neutralHadronEnergyFraction() << " , "
	  //  				      << it->chargedMuEnergyFraction() << " , "
	  //  				      << it->HFHadronEnergyFraction() << " , "
	  //  				      << it->HFEMEnergyFraction() << " , "
	  //  				      << SUM
	  //  				      << std::endl; 
	  // 			  }
	  // 			////////////////////////////////////////////////////////////////////

	  // TODO: Gen matching?

	  // Trigger matching
	  // FIXME? matcher and embedder are not run during cmsRun for some reason
	  //TEST
	  //std::cout << "size of trigger matches: " << it->triggerObjectMatches().size() << std::endl;
	  //TEST
	  //// The Ele+Jet+Jet path
	  //const pat::TriggerObjectStandAloneCollection matchesEleJetJet = it->triggerObjectMatchesByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*");
	  //if(matchesEleJetJet.size() > 0)
	  //{
	  //  // do stuff, like fill ntuple vars
	  //}

	} // end loop over jets
    }
  else
    {
      edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the pf jets";
    }

  //L1Offset JEC
  // if(jetsL1Offset.isValid())
  // {
  // 
  //   for( std::vector<pat::Jet>::const_iterator it = jetsL1Offset->begin(); it != jetsL1Offset->end(); ++it )
  //     {
  //	// Only look at jets with pt>=20 GeV
  //	if( it->pt()<20 )
  //	  continue;
  //       // exit from loop when you reach the required number of jets
  //       if(l1offsetJEC_vec->size() >= maxSize)
  // 	break;
  // 
  //       l1offsetJEC_vec->push_back( it->correctedJet("L1Offset").pt()/it->correctedJet("Uncorrected").pt() );
  //     }
  // }
  // else
  // {
  // 	edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTagL1Offset;
  // }

  delete jecUnc;

  //-----------------------------------------------------------------
  // put vectors in the event
	
	
  iEvent.put( hasJetWithBadUnc,prefix + "hasJetWithBadUnc" + suffix );
  iEvent.put( bestVertexTrackAssociationFactor,prefix + "BestVertexTrackAssociationFactor" + suffix );
  iEvent.put( bestVertexTrackAssociationIndex,prefix + "BestVertexTrackAssociationIndex" + suffix);
  iEvent.put( closestVertexWeighted3DSeparation,prefix + "ClosestVertexWeighted3DSeparation" + suffix );
  iEvent.put( closestVertexWeightedXYSeparation,prefix + "ClosestVertexWeightedXYSeparation" + suffix );
  iEvent.put( closestVertexWeightedZSeparation,prefix + "ClosestVertexWeightedZSeparation" + suffix );
  iEvent.put( closestVertex3DIndex,prefix + "ClosestVertex3DIndex" + suffix);
  iEvent.put( closestVertexXYIndex,prefix + "ClosestVertexXYIndex" + suffix);
  iEvent.put( closestVertexZIndex,prefix + "ClosestVertexZIndex" + suffix);
	
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( mt, prefix + "Mt" + suffix );
  iEvent.put( ptSmearedUp, prefix + "SmearedUpPt" + suffix );
  iEvent.put( ptSmearedDown, prefix + "SmearedDownPt" + suffix );
  iEvent.put( ptScaledUp, prefix + "ScaledUpPt" + suffix );
  iEvent.put( ptScaledDown, prefix + "ScaledDownPt" + suffix );
  iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( energySmearedUp, prefix + "SmearedUpEnergy" + suffix );
  iEvent.put( energySmearedDown, prefix + "SmearedDownEnergy" + suffix );
  iEvent.put( energyScaledUp, prefix + "ScaledUpEnergy" + suffix );
  iEvent.put( energyScaledDown, prefix + "ScaledDownEnergy" + suffix );
  iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
  iEvent.put( jecUnc_vec, prefix + "JECUnc" + suffix );
  iEvent.put( l2l3resJEC_vec, prefix + "L2L3ResJEC" + suffix );
  iEvent.put( l3absJEC_vec, prefix + "L3AbsJEC" + suffix );
  iEvent.put( l2relJEC_vec, prefix + "L2RelJEC" + suffix );
  iEvent.put( l1fastjetJEC_vec, prefix + "L1FastJetJEC" + suffix );
  iEvent.put( jerRes_vec, prefix + "JERRes" + suffix );
  iEvent.put( jerResSF_vec, prefix + "JERResSF" + suffix );
  iEvent.put( jerResSFup_vec, prefix + "JERResSFUp" + suffix );
  iEvent.put( jerResSFdown_vec, prefix + "JERResSFDown" + suffix );
  iEvent.put( partonFlavour, prefix + "PartonFlavour" + suffix );
  iEvent.put( hadronFlavour, prefix + "HadronFlavour" + suffix );
  iEvent.put( chargedEmEnergyFraction,  prefix + "ChargedEmEnergyFraction"  + suffix );
  iEvent.put( chargedHadronEnergyFraction,  prefix + "ChargedHadronEnergyFraction"  + suffix );
  iEvent.put( chargedMuEnergyFraction,  prefix + "ChargedMuEnergyFraction"  + suffix );
  iEvent.put( electronEnergyFraction,  prefix + "ElectronEnergyFraction"  + suffix );
  iEvent.put( muonEnergyFraction,  prefix + "MuonEnergyFraction"  + suffix );
  iEvent.put( neutralEmEnergyFraction,  prefix + "NeutralEmEnergyFraction"  + suffix );
  iEvent.put( neutralHadronEnergyFraction,  prefix + "NeutralHadronEnergyFraction"  + suffix );
  iEvent.put( photonEnergyFraction,  prefix + "PhotonEnergyFraction"  + suffix );
  iEvent.put( hfHadronEnergyFraction,  prefix + "HFHadronEnergyFraction"  + suffix );
  iEvent.put( hfEMEnergyFraction,  prefix + "HFEMEnergyFraction"  + suffix );
  iEvent.put( chargedHadronMultiplicity,  prefix + "ChargedHadronMultiplicity"  + suffix );
  iEvent.put( chargedMultiplicity,  prefix + "ChargedMultiplicity"  + suffix );
  iEvent.put( electronMultiplicity,  prefix + "ElectronMultiplicity"  + suffix );
  iEvent.put( muonMultiplicity,  prefix + "MuonMultiplicity"  + suffix );
  iEvent.put( neutralHadronMultiplicity,  prefix + "NeutralHadronMultiplicity"  + suffix );
  iEvent.put( neutralMultiplicity,  prefix + "NeutralMultiplicity"  + suffix );
  iEvent.put( photonMultiplicity,  prefix + "PhotonMultiplicity"  + suffix );
  iEvent.put( hfHadronMultiplicity,  prefix + "HFHadronMultiplicity"  + suffix );
  iEvent.put( hfEMMultiplicity,  prefix + "HFEMMultiplicity"  + suffix );
  iEvent.put( nConstituents,  prefix + "NConstituents"  + suffix );
  //iEvent.put( trackCountingHighEffBTag, prefix + "TrackCountingHighEffBTag" + suffix );
  //iEvent.put( trackCountingHighPurBTag, prefix + "TrackCountingHighPurBTag" + suffix );
  //iEvent.put( simpleSecondaryVertexHighEffBTag, prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
  //iEvent.put( simpleSecondaryVertexHighPurBTag, prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
  //iEvent.put( jetProbabilityBTag, prefix + "JetProbabilityBTag" + suffix );
  //iEvent.put( jetBProbabilityBTag, prefix + "JetBProbabilityBTag" + suffix );
  //iEvent.put( combinedSecondaryVertexBTag         ,prefix + "CombinedSecondaryVertexBTag" + suffix );    
  //iEvent.put( combinedSecondaryVertexMVABTag      ,prefix + "CombinedSecondaryVertexMVABTag" + suffix ); 
  //iEvent.put( softElectronByPtBTag                ,prefix + "SoftElectronByPtBTag" + suffix );           
  //iEvent.put( softElectronByIP3dBTag              ,prefix + "SoftElectronByIP3dBTag" + suffix );         
  //iEvent.put( softMuonBTag                        ,prefix + "SoftMuonBTag" + suffix );                   
  //iEvent.put( softMuonByPtBTag                    ,prefix + "SoftMuonByPtBTag" + suffix );               
  //iEvent.put( softMuonByIP3dBTag                  ,prefix + "SoftMuonByIP3dBTag" + suffix );            
  iEvent.put( combinedInclusiveSecondaryVertexBTag,prefix + "CombinedInclusiveSecondaryVertexBTag" + suffix );
  iEvent.put( combinedMVABTag                     ,prefix + "CombinedMVABTag" + suffix ) ;
  //iEvent.put( combinedCvsLJetTag                  ,prefix + "CombinedCvsLJetTags" + suffix ) ;
  //iEvent.put( combinedCvsBJetTag                  ,prefix + "CombinedCvsBJetTags" + suffix ) ;
  iEvent.put( passLooseID, prefix + "PassLooseID" + suffix);
  iEvent.put( passTightID, prefix + "PassTightID" + suffix);
	
  if(!isPuppiJetColl){
    iEvent.put( jetpileup_mva, prefix + "PileupMVA" + suffix);
    iEvent.put( jetpileup_mva_passesLoose,  prefix + "PileupMVApassesLoose"  + suffix);
    iEvent.put( jetpileup_mva_passesMedium, prefix + "PileupMVApassesMedium" + suffix);
    iEvent.put( jetpileup_mva_passesTight,  prefix + "PileupMVApassesTight"  + suffix);
    iEvent.put( vtxNtracks,  prefix + "VtxNtracks"  + suffix);
    iEvent.put( vtxMass,  prefix + "VtxMass"  + suffix);
    iEvent.put( vtxPx,  prefix + "VtxPx"  + suffix);
    iEvent.put( vtxPy,  prefix + "VtxPy"  + suffix);
    iEvent.put( vtx3DVal,  prefix + "Vtx3DVal"  + suffix);
    iEvent.put( vtx3DSig,  prefix + "Vtx3DSig"  + suffix);
    iEvent.put( leadTrackPt,  prefix + "LeadTrackPt"  + suffix);
    iEvent.put( softLepPt,  prefix + "SoftLepPt"  + suffix);
    iEvent.put( softLepRatio,  prefix + "SoftLepRatio"  + suffix);
    iEvent.put( softLepDr,  prefix + "SoftLepDr"  + suffix);
  }
  iEvent.put(betaStar       , prefix + "BetaStar"        + suffix ) ;
  iEvent.put(betaStarClassic, prefix + "BetaStarClassic" + suffix ) ;
  iEvent.put(beta           , prefix + "Beta"            + suffix ) ;
  iEvent.put(betaClassic    , prefix + "BetaClassic"     + suffix ) ;




}
