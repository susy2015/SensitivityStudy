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

#include "SSTrimAndSlimCommon.h"

int main(int argc, char* argv[])
{
  if (argc < 1)
  {
    std::cerr <<"Please give 1 argument " << "inputFileName " << std::endl;
    std::cerr <<"Valid configurations are: " << std::endl;
    std::cerr <<"./SSTrimAndSlim_Signal root://cmseos.fnal.gov//store/user/lpcsusyhad/Spring15_74X_Feb_2016_Ntp_v6X_forMoriond/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/Spring15_74X_Feb_2016_Ntp_v6p0_forMoriond_TTJets_SingleLeptFromT/160213_195624/0000/stopFlatNtuples_2.root" << std::endl;
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

  output_str = "SSTrimAndSlimmed_" + tag + "_" + fileid;
  std::cout << "Output File Name: " << output_str << std::endl;

  TChain *originalTree = new TChain("stopTreeMaker/AUX");
  originalTree->Add(input_str.c_str());
  //originalTree->SetBranchStatus("*", 1);
   
  TFile* output = new TFile((output_str).c_str(), "RECREATE");
  TDirectory *mydict = output->mkdir("stopTreeMaker");
  mydict->cd();
  TTree* selectedTree = new TTree("SSTree","SSTree");
  //TTree* selectedTree = originalTree->CloneTree(0);
  //search bin variables
  Double_t met,mt2; Int_t ntopjets, nbotjets;
  selectedTree->Branch("met",&met,"met/D");
  selectedTree->Branch("mt2",&mt2,"mt2/D");
  selectedTree->Branch("nTop",&ntopjets,"nTop/I");
  selectedTree->Branch("nBot",&nbotjets,"nBot/I");
  //muon and electron variables, and bool of passLeptVeto
  Int_t nmus,nels; Bool_t passLeptVeto;
  selectedTree->Branch("nMuons"    ,&nmus,"nMuons/I"    );
  selectedTree->Branch("nElectrons",&nels,"nElectrons/I");
  selectedTree->Branch("passLeptVeto",&passLeptVeto,"passLeptVeto/O");
  //AUX variables maybe useful for research
  Int_t njets30,njets50; Double_t ht,htTops;
  selectedTree->Branch("nJets30",&njets30,"nJets30/I");
  selectedTree->Branch("nJets50",&njets50,"nJets50/I");
  selectedTree->Branch("ht",&ht,"ht/D");
  selectedTree->Branch("htTops",&htTops,"htTops/D");

  //LSP and mediator mass for SUS Signal
  Double_t SusyMotherMass, SusyLSPMass;
  selectedTree->Branch("SusyMotherMass",&SusyMotherMass,"SusyMotherMass/D");
  selectedTree->Branch("SusyLSPMass"   ,&SusyLSPMass   ,"SusyLSPMass/D");

  std::shared_ptr<topTagger::type3TopTagger>type3Ptr(nullptr);
  NTupleReader *tr=0;
  //initialize the type3Ptr defined in the customize.h
  AnaFunctions::prepareForNtupleReader();
  tr = new NTupleReader(originalTree, AnaConsts::activatedBranchNames);
  const std::string spec = "lostlept";
  BaselineVessel *myBaselineVessel = 0;
  myBaselineVessel = new BaselineVessel(*tr, spec, "fastsim");
  if( !useNewTagger ){ myBaselineVessel->SetupTopTagger(false, "Legacy_TopTagger.cfg" ); }
  else
  {
    if( useLegacycfg ){ myBaselineVessel->SetupTopTagger(true, "Legacy_TopTagger.cfg" ); }
    else{ myBaselineVessel->SetupTopTagger(true, "TopTagger.cfg" ); }
  }
  //The passBaseline is registered here
  tr->registerFunction(*myBaselineVessel);

  while(tr->getNextEvent())
  {
    /*
    bool passLeptVeto = tr.getVar<bool>("passLeptVeto"+spec);
    bool passnJets = tr.getVar<bool>("passnJets"+spec);
    bool passMET = tr.getVar<bool>("passMET"+spec);
    bool passHT = tr.getVar<bool>("passHT"+spec);
    bool passMT2 = tr.getVar<bool>("passMT2"+spec);
    bool passTagger = tr.getVar<bool>("passTagger"+spec);
    bool passBJets = tr.getVar<bool>("passBJets"+spec);
    bool passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter"+spec);
    bool passQCDHighMETFilter = tr.getVar<bool>("passQCDHighMETFilter"+spec);
    bool passdPhis = tr.getVar<bool>("passdPhis"+spec);
    */
    bool passSSTrimAndSlim = tr->getVar<bool>("passBaseline"+spec);
    /*
    passSSTrimAndSlim = ( met > 200)
                && passnJets
                && passHT
                && passMT2
                //&& passTagger
                && passBJets
                && passNoiseEventFilter;
    */
    if(passSSTrimAndSlim)
    {
      //searchbin variables
      met = tr->getVar<double>("met");
      mt2 = tr->getVar<double>("best_had_brJet_MT2"+spec);       
      ntopjets = tr->getVar<int>("nTopCandSortedCnt"+spec);
      nbotjets = tr->getVar<int>("cntCSVS"+spec);
      //Muon and Electron variables
      nmus = tr->getVar<int>("nMuons_CUT"+spec);
      nels = tr->getVar<int>("nElectrons_CUT"+spec);
      passLeptVeto = tr->getVar<bool>("passLeptVeto"+spec);

      //AUX variables
      njets30 = tr->getVar<int>("cntNJetsPt30Eta24"+spec);
      njets50 = tr->getVar<int>("cntNJetsPt50Eta24"+spec);
      ht = tr->getVar<double>("HT"+spec);
      htTops = GetHTTops( tr->getVec<TLorentzVector>("vTops"+spec) );
      //double mht = tr->getVar<double>("mht"); 
 
      SusyMotherMass = tr->getVar<double>("SusyMotherMass");
      SusyLSPMass    = tr->getVar<double>("SusyLSPMass");
     
      selectedTree->Fill();
    }
    else continue;
  }
  selectedTree->Write();
  output->Write(); 
  output->Close();

  if (originalTree) delete originalTree;

  std::string d = "root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15";
  std::system(("xrdcp " + output_str + " " + d).c_str());
  std::system(("rm " + output_str).c_str());

  return 0;
}
