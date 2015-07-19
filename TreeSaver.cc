#include "bToDDoubleStar/TreeSaver.h"
#include <string>

//the static variables

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  TTree* TreeSaver::pDataTree;

  bool TreeSaver::initialized=false;
  vector<void*> TreeSaver::treeData;

  vector<float> TreeSaver::dataF;

  vector<int> TreeSaver::dataI;
  vector<string> TreeSaver::fieldNamesF;
  vector<string> TreeSaver::fieldNamesI;

#if defined(BELLE_NAMESPACE)
}
#endif
