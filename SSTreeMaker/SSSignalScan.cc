#include <fstream>
#include <map>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdlib.h>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/searchBins.h"

#include "SSSignalScan.h"

int main(int argc, char* argv[])
{
  if (argc < 1)
  {
    std::cerr <<"Please give 1 argument " << "inputFileName " << std::endl;
    std::cerr <<"Valid configurations are: " << std::endl;
    std::cerr <<"./SSSignalScan root://cmseos.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Feb_2016_Ntp_v6X_forMoriond/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Spring15_74X_Feb_2016_Ntp_v6p0_forMoriond_TTJets_SingleLeptFromT/160213_195624/0000/stopFlatNtuples_2.root" << std::endl;
    return -1;
  }
  std::string input_str(argv[1]);
  std::string trim;
  //std::string outputpath = "/eos/uscms/store/group/lpcsusyhad/hua/Skimmed_2015Nov15/";

  std::string output_str;
  //here is a little bit tricky when dealing with the slash... need to improve
  //for all the data samples and ttbar leptonic MC samples
  std::string tag = input_str.substr(find_Nth(input_str,10,"/") + 1,find_Nth(input_str,11,"/")-find_Nth(input_str,10,"/")-1);
  std::size_t idpos = input_str.find("stopFlatNtuples");
  std::string fileid = input_str.substr (idpos);

  output_str = "SSSignalScan_" + tag + "_" + fileid;
  std::cout << "Output File Name: " << output_str << std::endl;

  TChain *originalTree = new TChain("stopTreeMaker/AUX");
  originalTree->Add(input_str.c_str());
  //originalTree->SetBranchStatus("*", 1);
   
  std::shared_ptr<topTagger::type3TopTagger>type3Ptr(nullptr);
  NTupleReader *tr=0;
  //initialize the type3Ptr defined in the customize.h
  AnaFunctions::prepareForNtupleReader();
  tr = new NTupleReader(originalTree, AnaConsts::activatedBranchNames);
  const std::string spec = "lostlept";
  BaselineVessel *myBaselineVessel = 0;
  myBaselineVessel = new BaselineVessel(*tr, spec, "fastsim");
  myBaselineVessel->toptaggerCfgFile = "Example_TopTagger.cfg";
  //The passBaseline is registered here
  tr->registerFunction(*myBaselineVessel);

  SignalScanHistgram mySignalScanHistgram;
  mySignalScanHistgram.BookHistgram( (output_str).c_str() );

  while(tr->getNextEvent())
  {
    bool passSSTrimAndSlim = tr->getVar<bool>("passBaseline"+spec);
    double SusyMotherMass = tr->getVar<double>("SusyMotherMass");
    double SusyLSPMass    = tr->getVar<double>("SusyLSPMass");
    (mySignalScanHistgram.h_totEvt_xSusyMotherMass_ySusyLSPMass)->Fill(SusyMotherMass,SusyLSPMass,1);
  }

  (mySignalScanHistgram.oFile)->Write();
  (mySignalScanHistgram.oFile)->Close();

  if (originalTree) delete originalTree;

  std::string d = "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15";
  std::system(("xrdcp " + output_str + " " + d).c_str());
  std::system(("rm " + output_str).c_str());

  return 0;
}
