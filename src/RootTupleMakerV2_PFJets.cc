#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

RootTupleMakerV2_PFJets::RootTupleMakerV2_PFJets(const edm::ParameterSet& iConfig) :
inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
inputTagL1Offset(iConfig.getParameter<edm::InputTag>("InputTagL1Offset")),
prefix  (iConfig.getParameter<std::string>  ("Prefix")),
suffix  (iConfig.getParameter<std::string>  ("Suffix")),
maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
jecUncPath(iConfig.getParameter<std::string>("JECUncertainty")),
readJECuncertainty (iConfig.getParameter<bool>   ("ReadJECuncertainty")),
vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag"))

//OLD
//    applyResJEC (iConfig.getParameter<bool>     ("ApplyResidualJEC")),
//    resJEC (iConfig.getParameter<std::string>   ("ResidualJEC"))
{
        produces <bool>                 ( "hasJetWithBadUnc" );
	produces <std::vector<double> > ( prefix + "Eta" + suffix );
	produces <std::vector<double> > ( prefix + "Phi" + suffix );
	produces <std::vector<double> > ( prefix + "Pt" + suffix );
	produces <std::vector<double> > ( prefix + "PtRaw" + suffix );
	produces <std::vector<double> > ( prefix + "Energy" + suffix );
	produces <std::vector<double> > ( prefix + "EnergyRaw" + suffix );
	produces <std::vector<double> > ( prefix + "JECUnc" + suffix );
	produces <std::vector<double> > ( prefix + "L2L3ResJEC" + suffix );
	produces <std::vector<double> > ( prefix + "L3AbsJEC" + suffix );
	produces <std::vector<double> > ( prefix + "L2RelJEC" + suffix );
	produces <std::vector<double> > ( prefix + "L1FastJetJEC" + suffix );
	produces <std::vector<double> > ( prefix + "L1OffsetJEC" + suffix );
	produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
	produces <std::vector<double> > ( prefix + "ChargedEmEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "ChargedHadronEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "ChargedMuEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "ElectronEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "MuonEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "NeutralEmEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "NeutralHadronEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "PhotonEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "HFHadronEnergyFraction"  + suffix );
	produces <std::vector<double> > ( prefix + "HFEMEnergyFraction"  + suffix );
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
	produces <std::vector<double> > ( prefix + "TrackCountingHighEffBTag" + suffix );
	produces <std::vector<double> > ( prefix + "TrackCountingHighPurBTag" + suffix );
	produces <std::vector<double> > ( prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
	produces <std::vector<double> > ( prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
	produces <std::vector<double> > ( prefix + "JetProbabilityBTag" + suffix );
	produces <std::vector<double> > ( prefix + "JetBProbabilityBTag" + suffix );
	produces <std::vector<double> > ( prefix + "CombinedSecondaryVertexBTag" + suffix );    
	produces <std::vector<double> > ( prefix + "CombinedSecondaryVertexMVABTag" + suffix ); 
	produces <std::vector<double> > ( prefix + "SoftElectronByPtBTag" + suffix );           
	produces <std::vector<double> > ( prefix + "SoftElectronByIP3dBTag" + suffix );         
	produces <std::vector<double> > ( prefix + "SoftMuonBTag" + suffix );                   
	produces <std::vector<double> > ( prefix + "SoftMuonByPtBTag" + suffix );               
	produces <std::vector<double> > ( prefix + "SoftMuonByIP3dBTag" + suffix );             
	produces <std::vector<double> > ( prefix + "CombinedInclusiveSecondaryVertexBTag" + suffix );
	produces <std::vector<double> > ( prefix + "CombinedMVABTag" + suffix );
	produces <std::vector<int> >    ( prefix + "PassLooseID" + suffix);
	produces <std::vector<int> >    ( prefix + "PassTightID" + suffix);
	produces <std::vector<double> > ( prefix + "BestVertexTrackAssociationFactor" + suffix );
	produces <std::vector<int> >    ( prefix + "BestVertexTrackAssociationIndex" + suffix);
	produces <std::vector<double> > ( prefix + "ClosestVertexWeighted3DSeparation" + suffix );
	produces <std::vector<double> > ( prefix + "ClosestVertexWeightedXYSeparation" + suffix );
	produces <std::vector<double> > ( prefix + "ClosestVertexWeightedZSeparation" + suffix );
	produces <std::vector<int> >    ( prefix + "ClosestVertex3DIndex" + suffix);
	produces <std::vector<int> >    ( prefix + "ClosestVertexXYIndex" + suffix);
	produces <std::vector<int> >    ( prefix + "ClosestVertexZIndex" + suffix);
	produces <std::vector<double> > ( prefix + "Beta" + suffix ) ;
	produces <std::vector<double> > ( prefix + "BetaStar" + suffix ) ;
	produces <std::vector<double> > ( prefix + "BetaClassic" + suffix ) ;
	produces <std::vector<double> > ( prefix + "BetaStarClassic" + suffix ) ;
}


PFJetIDSelectionFunctor pfjetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
PFJetIDSelectionFunctor pfjetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );

pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

void RootTupleMakerV2_PFJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

        std::auto_ptr<bool>                  hasJetWithBadUnc ( new bool() );
	std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  pt_raw  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  energy_raw ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  jecUnc_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  l2l3resJEC_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  l3absJEC_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  l2relJEC_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  l1fastjetJEC_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  l1offsetJEC_vec ( new std::vector<double>()  );
	std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
	std::auto_ptr<std::vector<double> >  chargedEmEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  chargedHadronEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  chargedMuEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  electronEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  muonEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  neutralEmEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  neutralHadronEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  photonEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  hfHadronEnergyFraction  ( new std::vector<double>()  ) ;
	std::auto_ptr<std::vector<double> >  hfEMEnergyFraction  ( new std::vector<double>()  ) ;
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
	std::auto_ptr<std::vector<double> >  trackCountingHighEffBTag  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  trackCountingHighPurBTag  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  simpleSecondaryVertexHighEffBTag  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  simpleSecondaryVertexHighPurBTag  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  jetProbabilityBTag  ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  jetBProbabilityBTag  ( new std::vector<double>()  );

	std::auto_ptr<std::vector<double> >  combinedSecondaryVertexBTag          ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  combinedSecondaryVertexMVABTag       ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  softElectronByPtBTag                 ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  softElectronByIP3dBTag               ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  softMuonBTag                         ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  softMuonByPtBTag                     ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  softMuonByIP3dBTag                   ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  combinedInclusiveSecondaryVertexBTag ( new std::vector<double>()  );
	std::auto_ptr<std::vector<double> >  combinedMVABTag                      ( new std::vector<double>()  );
	
	std::auto_ptr<std::vector<int> >  passLooseID  ( new std::vector<int>()  );
	std::auto_ptr<std::vector<int> >  passTightID  ( new std::vector<int>()  );
	std::auto_ptr <std::vector<double> >  bestVertexTrackAssociationFactor  ( new std::vector<double>()  );
	std::auto_ptr <std::vector<int> >     bestVertexTrackAssociationIndex   ( new std::vector<int>()  );
	std::auto_ptr <std::vector<double> >  closestVertexWeighted3DSeparation  ( new std::vector<double>()  );
	std::auto_ptr <std::vector<double> >  closestVertexWeightedXYSeparation  ( new std::vector<double>()  );
	std::auto_ptr <std::vector<double> >  closestVertexWeightedZSeparation  ( new std::vector<double>()  );
	std::auto_ptr <std::vector<int> >     closestVertex3DIndex            ( new std::vector<int>()  );
	std::auto_ptr <std::vector<int> >     closestVertexXYIndex           ( new std::vector<int>()  );
	std::auto_ptr <std::vector<int> >     closestVertexZIndex            ( new std::vector<int>()  );

	std::auto_ptr <std::vector<double > > betaStar        ( new std::vector<double>());
	std::auto_ptr <std::vector<double > > betaStarClassic ( new std::vector<double>());
	std::auto_ptr <std::vector<double > > beta            ( new std::vector<double>());
	std::auto_ptr <std::vector<double > > betaClassic     ( new std::vector<double>());
	
	//-----------------------------------------------------------------

	// OLD
	//   edm::FileInPath fipUnc(jecUncPath);;
	//   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fipUnc.fullPath());
	//
	//   JetCorrectorParameters *ResJetCorPar = 0;
	//   FactorizedJetCorrector *JEC = 0;
	//   if(applyResJEC) {
	//     edm::FileInPath fipRes(resJEC);
	//     ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
	//     std::vector<JetCorrectorParameters> vParam;
	//     vParam.push_back(*ResJetCorPar);
	//     JEC = new FactorizedJetCorrector(vParam);
	//   }

	//JEC Uncertainties

	*hasJetWithBadUnc.get() = false;
	JetCorrectionUncertainty *jecUnc = 0;
	if(readJECuncertainty)
	{
		//(See https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1075/1.html
		// and https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2367/1.html)
		// handle the jet corrector parameters collection
		edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
		// get the jet corrector parameters collection from the global tag
		iSetup.get<JetCorrectionsRecord>().get(jecUncPath,JetCorParColl);
		// get the uncertainty parameters from the collection
		JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
		// instantiate the jec uncertainty object
		jecUnc = new JetCorrectionUncertainty(JetCorPar);
	}

	edm::Handle<std::vector<pat::Jet> > jets;
	iEvent.getByLabel(inputTag, jets);

	edm::Handle<std::vector<pat::Jet> > jetsL1Offset;
	iEvent.getByLabel(inputTagL1Offset, jetsL1Offset);
	
	edm::Handle<reco::VertexCollection> primaryVertices;  // DB
	iEvent.getByLabel(vtxInputTag,primaryVertices);       // DB

	if(jets.isValid())
	{
		edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # PFJets: " << jets->size();

		for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it )
		{
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

			// OLD
			//       double corr = 1.;
			//       if( applyResJEC && iEvent.isRealData() ) {
			//         JEC->setJetEta( it->eta() );
			//         JEC->setJetPt( it->pt() ); // here you put the L2L3 Corrected jet pt
			//         corr = JEC->getCorrection();
			//       }

			// Status of JEC
			//std::cout << "PF: currentJECLevel(): " << it->currentJECLevel() << std::endl;
			//std::cout << "PF: currentJECSet(): " << it->currentJECSet() << std::endl;
			//-------------------

			// fill in all the vectors

			// OLD
			// pt->push_back( it->pt()*corr );
			// energy->push_back( it->energy()*corr );
			// resJEC_vec->push_back( corr );



			// Vertex association

			int bestVtxIndex3Ddist = -1;
			int bestVtxIndexXYdist = -1;
			int bestVtxIndexZdist = -1;

			int bestVtxIndexSharedTracks = -1;
			
			double minVtxDist3D = 999999.;
			double minVtxDistXY = -99999.;
			double minVtxDistZ  = -99999.;
			double maxTrackAssocRatio = -9999.;
			
			
			// Loop on primary Vertices and jets and perform associations 
			
			reco::VertexCollection::const_iterator lead_vertex = primaryVertices->end();
			bool found_lead_vertex = false;

			if(primaryVertices.isValid())
			{
				edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # Primary Vertices: " << primaryVertices->size();
				
				// Main Vertex Loop
				for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it )
				{

					double sumweights = 0.0;
					double dist3Dweighted = 0.0;
					double distXYweighted = 0.0;
					double distZweighted = 0.0;
					double assocsumpttracks = 0.0;
					double trackassociationratio = 0.000001;
	
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
						double trackptweight = jtrack->pt();
						sumweights += trackptweight;

						// Weighted Distance Calculation
						double distXY= jtrack->dxy(v_it->position());
						double distZ = jtrack->dz(v_it->position());
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
				edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << vtxInputTag;
			}

			// Get the constituents of the PFJet

			std::vector <reco::PFCandidatePtr> constituents  = it -> getPFConstituents();
			std::vector <reco::PFCandidatePtr>::iterator i_constituent   = constituents.begin();
			std::vector <reco::PFCandidatePtr>::iterator end_constituent = constituents.end();
			double sum_track_pt = 0.;
	
			double jetBetaStar        = 0.0 ;
			double jetBetaStarClassic = 0.0 ;
			double jetBeta            = 0.0 ;
			double jetBetaClassic     = 0.0 ;

			if ( found_lead_vertex ) { 
			  for (; i_constituent != end_constituent; ++i_constituent ) { 
			    reco::PFCandidatePtr & constituent = *i_constituent;
			    if ( ! constituent -> trackRef().isNonnull()   ) continue;
			    if ( ! constituent -> trackRef().isAvailable() ) continue;
			    
			    try { 
			      double track_pt = constituent -> trackRef() -> pt();
			      sum_track_pt += track_pt;
			      
			      bool track_from_lead_vertex = find( lead_vertex->tracks_begin(), 
								  lead_vertex->tracks_end()  , 
								  reco::TrackBaseRef( constituent -> trackRef())) != lead_vertex -> tracks_end();
			      
			      bool track_from_other_vertex = false;
			      
			      double dZ0 = fabs(constituent ->trackRef()->dz(lead_vertex->position()));
			      double dZ = dZ0; 
			      
			      for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ){
				if( v_it -> isFake() || v_it -> ndof() < 4 ) continue;
				bool is_lead_vertex  = (v_it -> position() - lead_vertex -> position()).r() < 0.02;
				if( ! is_lead_vertex && ! track_from_other_vertex ) {
				  track_from_other_vertex = find( v_it -> tracks_begin(), 
								  v_it -> tracks_end  (), 
								  reco::TrackBaseRef(constituent -> trackRef())) != v_it -> tracks_end(); 
				}
				dZ = std::min(dZ,fabs(constituent->trackRef()->dz( v_it -> position())));
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
			pt_raw->push_back( it->correctedJet("Uncorrected").pt() );
			energy->push_back( it->energy() );
			energy_raw->push_back( it->correctedJet("Uncorrected").energy() );
			l2l3resJEC_vec->push_back( it->pt()/it->correctedJet("L3Absolute").pt() );
			l3absJEC_vec->push_back( it->correctedJet("L3Absolute").pt()/it->correctedJet("L2Relative").pt() );
			l2relJEC_vec->push_back( it->correctedJet("L2Relative").pt()/it->correctedJet("L1FastJet").pt() );
			l1fastjetJEC_vec->push_back( it->correctedJet("L1FastJet").pt()/it->correctedJet("Uncorrected").pt() );
			if(readJECuncertainty){ 
			  double uncertainty = -999.;
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
			trackCountingHighEffBTag->push_back( it->bDiscriminator("trackCountingHighEffBJetTags") );
			trackCountingHighPurBTag->push_back( it->bDiscriminator("trackCountingHighPurBJetTags") );
			simpleSecondaryVertexHighEffBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighEffBJetTags") );
			simpleSecondaryVertexHighPurBTag->push_back( it->bDiscriminator("simpleSecondaryVertexHighPurBJetTags") );
			jetProbabilityBTag->push_back( it->bDiscriminator("jetProbabilityBJetTags") );
			jetBProbabilityBTag->push_back( it->bDiscriminator("jetBProbabilityBJetTags") );
			combinedSecondaryVertexBTag         ->push_back( it->bDiscriminator("combinedSecondaryVertexBJetTags"         ));
			combinedSecondaryVertexMVABTag      ->push_back( it->bDiscriminator("combinedSecondaryVertexMVABJetTags"      ));
			softElectronByPtBTag                ->push_back( it->bDiscriminator("softElectronByPtBJetTags"                ));                
			softElectronByIP3dBTag              ->push_back( it->bDiscriminator("softElectronByIP3dBJetTags"              ));
			softMuonBTag                        ->push_back( it->bDiscriminator("softMuonBJetTags"                        ));
			softMuonByPtBTag                    ->push_back( it->bDiscriminator("softMuonByPtBJetTags"                    ));                
			softMuonByIP3dBTag                  ->push_back( it->bDiscriminator("softMuonByIP3dBJetTags"                  ));
			combinedInclusiveSecondaryVertexBTag->push_back( it->bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));
			combinedMVABTag                     ->push_back( it->bDiscriminator("combinedMVABJetTags"                     ));
			passLooseID->push_back( passjetLoose );
			passTightID->push_back( passjetTight );
			

// 			//////////////////////////////////////////////////////////////////// 
// 			if( fabs(it->eta()) > 3) 
// 			  {
//  			    double SUM = it->chargedEmEnergyFraction() + it->chargedHadronEnergyFraction() + it->neutralEmEnergyFraction() + it->neutralHadronEnergyFraction() + it->chargedMuEnergyFraction() + it->HFHadronEnergyFraction() + it->HFEMEnergyFraction() ; 
			    
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

		}
	}
	else
	{
		edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTag;
	}

	//L1Offset JEC
	if(jetsL1Offset.isValid())
	{

	  for( std::vector<pat::Jet>::const_iterator it = jetsL1Offset->begin(); it != jetsL1Offset->end(); ++it )
	    {
	      // exit from loop when you reach the required number of jets
	      if(l1offsetJEC_vec->size() >= maxSize)
		break;

	      l1offsetJEC_vec->push_back( it->correctedJet("L1Offset").pt()/it->correctedJet("Uncorrected").pt() );
	    }
	}
	else
	{
		edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTagL1Offset;
	}

	//OLD
	delete jecUnc;
	//   delete ResJetCorPar;
	//   delete JEC;

	//-----------------------------------------------------------------
	// put vectors in the event
	
	
	iEvent.put( hasJetWithBadUnc, "hasJetWithBadUnc" );
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
	iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
	iEvent.put( energy, prefix + "Energy" + suffix );
	iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
	iEvent.put( jecUnc_vec, prefix + "JECUnc" + suffix );
	iEvent.put( l2l3resJEC_vec, prefix + "L2L3ResJEC" + suffix );
	iEvent.put( l3absJEC_vec, prefix + "L3AbsJEC" + suffix );
	iEvent.put( l2relJEC_vec, prefix + "L2RelJEC" + suffix );
	iEvent.put( l1fastjetJEC_vec, prefix + "L1FastJetJEC" + suffix );
	iEvent.put( l1offsetJEC_vec, prefix + "L1OffsetJEC" + suffix );
	iEvent.put( partonFlavour, prefix + "PartonFlavour" + suffix );
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
	iEvent.put( trackCountingHighEffBTag, prefix + "TrackCountingHighEffBTag" + suffix );
	iEvent.put( trackCountingHighPurBTag, prefix + "TrackCountingHighPurBTag" + suffix );
	iEvent.put( simpleSecondaryVertexHighEffBTag, prefix + "SimpleSecondaryVertexHighEffBTag" + suffix );
	iEvent.put( simpleSecondaryVertexHighPurBTag, prefix + "SimpleSecondaryVertexHighPurBTag" + suffix );
	iEvent.put( jetProbabilityBTag, prefix + "JetProbabilityBTag" + suffix );
	iEvent.put( jetBProbabilityBTag, prefix + "JetBProbabilityBTag" + suffix );
	iEvent.put( combinedSecondaryVertexBTag         ,prefix + "CombinedSecondaryVertexBTag" + suffix );    
	iEvent.put( combinedSecondaryVertexMVABTag      ,prefix + "CombinedSecondaryVertexMVABTag" + suffix ); 
	iEvent.put( softElectronByPtBTag                ,prefix + "SoftElectronByPtBTag" + suffix );           
	iEvent.put( softElectronByIP3dBTag              ,prefix + "SoftElectronByIP3dBTag" + suffix );         
	iEvent.put( softMuonBTag                        ,prefix + "SoftMuonBTag" + suffix );                   
	iEvent.put( softMuonByPtBTag                    ,prefix + "SoftMuonByPtBTag" + suffix );               
	iEvent.put( softMuonByIP3dBTag                  ,prefix + "SoftMuonByIP3dBTag" + suffix );            
	iEvent.put( combinedInclusiveSecondaryVertexBTag,prefix + "CombinedInclusiveSecondaryVertexBTag" + suffix );
	iEvent.put( combinedMVABTag                     ,prefix + "CombinedMVABTag" + suffix ) ;
	iEvent.put( passLooseID, prefix + "PassLooseID" + suffix);
	iEvent.put( passTightID, prefix + "PassTightID" + suffix);
       
	iEvent.put(betaStar       , prefix + "BetaStar"        + suffix ) ;
	iEvent.put(betaStarClassic, prefix + "BetaStarClassic" + suffix ) ;
	iEvent.put(beta           , prefix + "Beta"            + suffix ) ;
	iEvent.put(betaClassic    , prefix + "BetaClassic"     + suffix ) ;




}
