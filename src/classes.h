#include <map>
#include <string>
#include <vector>
#include "DataFormats/Common/interface/Wrapper.h"

namespace {

  struct dictionary {
    std::map<std::string, std::vector<float > > dummy;
    edm::Wrapper< std::map<std::string, std::vector<float> > > dummy1;

    std::vector<std::vector<float > > dummy2;
    edm::Wrapper< std::vector<std::vector<float> > > dummy3;

    std::vector<std::vector<int > > dummy4;
    edm::Wrapper< std::vector<std::vector<int> > > dummy5;

  };

}
