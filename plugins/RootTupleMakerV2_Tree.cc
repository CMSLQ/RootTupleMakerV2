#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Tree.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "DataFormats/Provenance/interface/SelectedProducts.h"
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <map>
#include "boost/foreach.hpp"
#include <TBranch.h>
#include <TLorentzVector.h>

RootTupleMakerV2_Tree::~RootTupleMakerV2_Tree() {}

void RootTupleMakerV2_Tree::
analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  for( BranchConnector* connector : connectors)
    connector->connect(iEvent);
  tree->Fill();
}

template <class T>
void RootTupleMakerV2_Tree::TypedBranchConnector<T>::
connect(const edm::Event& iEvent) {
  edm::Handle<T> handle_;
  iEvent.getByLabel(ml, pin, handle_);
  //iEvent.getByToken(token_,handle_);
  object_ = *handle_;
}

template <class T> 
RootTupleMakerV2_Tree::TypedBranchConnector<T>::
TypedBranchConnector(edm::BranchDescription const* desc,
                     std::string t,
                     TTree * tree)
  :  ml( desc->moduleLabel() ),
     pin( desc->productInstanceName() )
{
  object_ptr_ = &object_;
  std::string s=pin+t;
  if(t!="")  { tree->Branch(pin.c_str(),  object_ptr_, s.c_str() );}  //raw type
  else       { tree->Branch(pin.c_str(), &object_ptr_            );}  //vector<type>
}

void RootTupleMakerV2_Tree::beginJob() {}

void RootTupleMakerV2_Tree::endJob() {}

//void RootTupleMakerV2_Tree::
//beginJob() {
RootTupleMakerV2_Tree::
RootTupleMakerV2_Tree(const edm::ParameterSet& cfg) {
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
  //auto allBranches = reg->allBranchDescriptions();
  //DEBUG
  //cout << "RootTupleMakerV2_Tree  " << cfg << endl;
  //DEBUG
  edm::ProductSelectorRules productSelectorRules_(cfg, "outputCommands", "RootTupleMakerV2_Tree");
  edm::ProductSelector productSelector_;
  productSelector_.initialize(productSelectorRules_, allBranches);

  ////DEBUG
  //for(auto const& selection : allBranches)
  //  edm::LogWarning("RootTupleMakerV2_Tree") << "IN allBranches: FOUND PRODUCT INSTANCE NAME: " << selection->productInstanceName();
  ////DEBUG
  ////DEBUG
  //edm::LogWarning("RootTupleMakerV2_Tree") << "ProductSelector looks like:" << productSelector_ << " END OF PRODUCT SELECTOR";
  ////DEBUG
  ////DEBUG
  //for(auto const& selection : allBranches)
  //  edm::LogWarning("RootTupleMakerV2_Tree") << "IN allBranches: FOUND PRODUCT INSTANCE NAME: " << selection->productInstanceName();
  ////DEBUG

  std::set<std::string> branchnames;

  for(auto const& selection : allBranches) {
    //DEBUG
    std::cout << "RootTupleMakerV2_Tree:  FOUND PRODUCT INSTANCE NAME: " << selection->productInstanceName() << std::endl;
    //DEBUG
    if(productSelector_.selected(*selection)) {

      //Check for duplicate branch names
      if (branchnames.find( selection->productInstanceName()) != branchnames.end() ) {
        throw edm::Exception(edm::errors::Configuration)
          << "More than one branch named: "
          << selection->productInstanceName() << std::endl
          << "Exception thrown from RootTupleMakerV2_Tree::RootTupleMakerV2_Tree" << std::endl;
      }
      else {
        branchnames.insert( selection->productInstanceName() );
      }
      

      ////Create RootTupleMakerV2_Tree branch
      //switch(leafmap.find( selection->friendlyClassName() )->second) {
      //case BOOL            : { edm::EDGetTokenT<bool> token = consumes<bool>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName())); connectors.push_back( new TypedBranchConnector                      <bool>         (selection, "/O", tree, token) ); break; }
      //case BOOL_V          : { edm::EDGetTokenT<std::vector<bool> > token = consumes<std::vector<bool> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector          <bool> >       (selection,   "", tree, token) ); break; }
      //case INT             : { edm::EDGetTokenT<int> token = consumes<int>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector                       <int>         (selection, "/I", tree, token) ); break; }
      //case INT_V           : { edm::EDGetTokenT<std::vector<int> >token = consumes<std::vector<int> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector           <int> >       (selection,   "", tree, token) ); break; }
      //case U_INT           : { edm::EDGetTokenT<unsigned int> token = consumes<unsigned int>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector              <unsigned int>         (selection, "/i", tree, token) ); break; }
      //case U_INT_V         : { edm::EDGetTokenT<std::vector<unsigned int> > token = consumes<std::vector<unsigned int> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector  <unsigned int> >       (selection,   "", tree, token) ); break; }
      //case SHORT           : { edm::EDGetTokenT<short> token = consumes<short>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector                     <short>         (selection, "/S", tree, token) ); break; }
      //case SHORT_V         : { edm::EDGetTokenT<std::vector<short> > token = consumes<std::vector<short> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector         <short> >       (selection,   "", tree, token) ); break; }
      //case U_SHORT         : { edm::EDGetTokenT<unsigned short> token = consumes<unsigned short>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector            <unsigned short>         (selection, "/s", tree, token) ); break; }
      //case U_SHORT_V       : { edm::EDGetTokenT<std::vector<unsigned short> > token = consumes<std::vector<unsigned short> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector<unsigned short> >       (selection,   "", tree, token) ); break; }
      //case FLOAT           : { edm::EDGetTokenT<float> token = consumes<float>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector                     <float>         (selection, "/F", tree, token) ); break; }
      //case FLOAT_V         : { edm::EDGetTokenT<std::vector<float> > token = consumes<std::vector<float> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector         <float> >       (selection,   "", tree, token) ); break; }
      //case DOUBLE          : { edm::EDGetTokenT<double> token = consumes<double>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector                    <double>         (selection, "/D", tree, token) ); break; }
      //case DOUBLE_V        : { edm::EDGetTokenT<std::vector<double> > token = consumes<std::vector<double> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector        <double> >       (selection,   "", tree, token) ); break; }
      //case LONG            : { edm::EDGetTokenT<long> token = consumes<long>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector                      <long>         (selection, "/L", tree, token) ); break; }
      //case LONG_V          : { edm::EDGetTokenT<std::vector<long> > token = consumes<std::vector<long> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector          <long> >       (selection,   "", tree, token) ); break; }
      //case U_LONG          : { edm::EDGetTokenT<unsigned long> token = consumes<unsigned long>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector             <unsigned long>         (selection, "/l", tree, token) ); break; }
      //case U_LONG_V        : { edm::EDGetTokenT<std::vector<unsigned long> > token = consumes<std::vector<unsigned long> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector <unsigned long> >       (selection,   "", tree, token) ); break; }
      ////                     {       }
      //case STRING          : { edm::EDGetTokenT<std::string> token = consumes<std::string>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector             <std::string  >         (selection,   "", tree, token) ); break; }
      //case STRING_V        : { edm::EDGetTokenT<std::vector<std::string> > token = consumes<std::vector<std::string> >(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<std::vector <std::string  > >       (selection,   "", tree, token) ); break; }
	    //                       {       }
      //case STRING_INT_M    : { edm::EDGetTokenT<mapStringInt> token = consumes<mapStringInt>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<mapStringInt>        (selection,   "", tree, token) ); break; }
      //case STRING_BOOL_M   : { edm::EDGetTokenT<mapStringBool> token = consumes<mapStringBool>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<mapStringBool>       (selection,   "", tree, token) ); break; }
      //case STRING_STRING_M : { edm::EDGetTokenT<mapStringString> token = consumes<mapStringString>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<mapStringString>     (selection,   "", tree, token) ); break; }
      //case STRING_FLOAT_V_M: { edm::EDGetTokenT<mapStringDoubles> token = consumes<mapStringDoubles>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<mapStringDoubles>    (selection,   "", tree, token) ); break; }
      //case FLOAT_V_V       : { edm::EDGetTokenT<vectorVectorFloats> token = consumes<vectorVectorFloats>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<vectorVectorFloats>  (selection,   "", tree, token) ); break; }
      //case INT_V_V         : { edm::EDGetTokenT<vectorVectorInts> token = consumes<vectorVectorInts>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<vectorVectorInts>    (selection,   "", tree, token) ); break; }
      //case BOOL_V_V        : { edm::EDGetTokenT<vectorVectorBools> token = consumes<vectorVectorBools>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<vectorVectorBools>   (selection,   "", tree, token) ); break; }
      //case STRING_V_V      : { edm::EDGetTokenT<vectorVectorStrings> token = consumes<vectorVectorStrings>(edm::InputTag(selection->moduleLabel(),selection->productInstanceName()));  connectors.push_back( new TypedBranchConnector<vectorVectorStrings> (selection,   "", tree, token) ); break; }

      //default:
      //  {
      //    std::string leafstring = "";
      //    typedef std::pair<std::string, LEAFTYPE> pair_t;
      //    BOOST_FOREACH( const pair_t& leaf, leafmap)
      //      leafstring+= "\t" + leaf.first + "\n";

      //    throw edm::Exception(edm::errors::Configuration)
      //      << "class RootTupleMakerV2_Tree does not handle leaves of type " << selection->className() << " like\n"
      //      <<   selection->friendlyClassName()   << "_"
      //      <<   selection->moduleLabel()         << "_"
      //      <<   selection->productInstanceName() << "_"
      //      <<   selection->processName()         << std::endl
      //      << "Valid leaf types are (friendlyClassName):\n"
      //      <<   leafstring
      //      << "Exception thrown from RootTupleMakerV2_Tree::beginJob\n";
      //  }
      //}
    }
  }
}

void
RootTupleMakerV2_Tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
