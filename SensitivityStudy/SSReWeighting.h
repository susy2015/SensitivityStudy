#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include "TChain.h"

#include "SusyAnaTools/Tools/NTupleReader.h"

#include "ConstantsSnippet.h"

struct SSSampleInfo
{
  std::string Tag;
  double weight;
  TChain *chain;
};

class SSSampleWeight
{
 public:
  std::vector<SSSampleInfo> SSSampleInfos;
  void SSSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, double kf, const TString &inputFileList );
 private:
  bool FillChain(TChain *chain, const TString &inputFileList, std::string tag);
};
