#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();

#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Tree.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Event.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_EventSelection.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_CaloJets.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFJets.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Electrons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Muons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_SuperClusters.h"
#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Trigger.h"

DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_Tree);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_Event);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_EventSelection);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_CaloJets);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_PFJets);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_Electrons);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_MET);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_Muons);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_SuperClusters);
DEFINE_ANOTHER_FWK_MODULE(RootTupleMakerV2_Trigger);
