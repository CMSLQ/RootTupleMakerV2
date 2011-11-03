#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Tree.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Event.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_EventSelection.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_CaloJets.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFJets.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Electrons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Taus.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Muons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Trigger.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_TriggerObjects.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Vertex.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenEventInfo.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenParticles.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenJets.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenMET.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Photons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFCandidates.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PhysicsDSTStream2011.h"

DEFINE_FWK_MODULE(RootTupleMakerV2_Tree);
DEFINE_FWK_MODULE(RootTupleMakerV2_Event);
DEFINE_FWK_MODULE(RootTupleMakerV2_EventSelection);
DEFINE_FWK_MODULE(RootTupleMakerV2_CaloJets);
DEFINE_FWK_MODULE(RootTupleMakerV2_PFJets);
DEFINE_FWK_MODULE(RootTupleMakerV2_Electrons);
DEFINE_FWK_MODULE(RootTupleMakerV2_Taus);
DEFINE_FWK_MODULE(RootTupleMakerV2_MET);
DEFINE_FWK_MODULE(RootTupleMakerV2_Muons);
DEFINE_FWK_MODULE(RootTupleMakerV2_Trigger);
DEFINE_FWK_MODULE(RootTupleMakerV2_TriggerObjects);
DEFINE_FWK_MODULE(RootTupleMakerV2_Vertex);
DEFINE_FWK_MODULE(RootTupleMakerV2_GenEventInfo);
DEFINE_FWK_MODULE(RootTupleMakerV2_GenParticles);
DEFINE_FWK_MODULE(RootTupleMakerV2_GenJets);
DEFINE_FWK_MODULE(RootTupleMakerV2_GenMET);
DEFINE_FWK_MODULE(RootTupleMakerV2_Photons);
DEFINE_FWK_MODULE(RootTupleMakerV2_PFCandidates);
DEFINE_FWK_MODULE(RootTupleMakerV2_PhysicsDSTStream2011);
