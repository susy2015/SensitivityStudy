#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
//define lumi in pb-1
//2015 data 13TeV
//#define LUMI 2262.946
//2016 ICHEP Exp
//#define LUMI 8000
//2016 End of year Exp
#define LUMI 36000

//Fill chain from txt file
bool FillChain(TChain *chain, const TString &inputFileList, std::string tag)
{
  std::ifstream infile( inputFileList, std::ifstream::in );
  std::string buffer;

  if(!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  std::cout << "TreeUtilities : FillChain " << tag << std::endl;
  while(1)
  {
    buffer.clear();
    infile >> buffer;

    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    //std::cout << (buffer.find(tag) != std::string::npos) << std::endl;
    if (buffer.find(tag) != std::string::npos) 
    {
      //std::cout << tag << " found!" << std::endl;
      chain->Add(buffer.c_str());
    }
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return true;
}

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
};

void SSSampleWeight::SSSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, double kf, const TString &inputFileList)
{
  SSSampleInfo oneInfo;

  oneInfo.Tag = tag;
  oneInfo.weight = xsec*lumi*kf/nevents;
  //weight is one if we are reading data
  //if( tag.find("HTMHT") != std::string::npos ) oneInfo.weight = 1;
  //negative weight for the sample other than QCD and HTMHT
  //if( !(tag.find("QCD") != std::string::npos) ) oneInfo.weight = -xsec*lumi/nevents;
  //if( tag.find("HTMHT") != std::string::npos ) oneInfo.weight = 1;
  oneInfo.chain= new TChain("stopTreeMaker/SSTree");
  //oneInfo.chain= new TChain("stopTreeMaker/AUX");
  if(!FillChain(oneInfo.chain, inputFileList, oneInfo.Tag))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }
  SSSampleInfos.push_back(oneInfo);
  oneInfo = {};
}
