#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <vector>

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/xSec.h"

#include "TStopwatch.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TPad.h"
#include "TStyle.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "SS.h"

//SearchBins mySearchBins("SB_v1_2017");
SearchBins mySearchBins("SB_59_2016");

//const double Scale = 591.5/2153.736;
//                             NJets =      1        2       3        4        5        6       7      >=8
const double NJetRweightingFactor[8] = {1.02845,1.08559,1.06879,0.922173,0.871796,0.99674,0.993756,0.539612};//2016 ICHEP, v8 MC and v9 data, 12.9 fb-1
//const double NJetRweightingFactor[8] = {0.926542,1.03995,0.919711,0.723581,0.869969,0.95682,0.584418,0.874059};//2016 Moriond

void LoopSSCS( SSSampleWeight& mySSSampleWeight )
{
  //clock to monitor the run time
  size_t t0 = clock();
  std::vector<SSSampleInfo>::iterator iter_SSSampleInfos;

  double mucs[NSB] = {0}, elcs[NSB]={0}, mucs_NMC[NSB] = {0}, elcs_NMC[NSB]={0};
  std::cout << "Let's do sensitivity study: " << std::endl;
  
  for(iter_SSSampleInfos = mySSSampleWeight.SSSampleInfos.begin(); iter_SSSampleInfos != mySSSampleWeight.SSSampleInfos.end(); iter_SSSampleInfos++)
  {    
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_SSSampleInfos).chain);

    double thisweight = (*iter_SSSampleInfos).weight;
    std::cout <<"Sample Type: "<< (*iter_SSSampleInfos).Tag << "; Weight: " << thisweight << std::endl;

    while(tr.getNextEvent())
    {
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      //searchbin variables
      int ntopjets = tr.getVar<int>("nTop");
      int nbotjets = tr.getVar<int>("nBot");
      double mt2 = tr.getVar<double>("mt2");
      double met = tr.getVar<double>("met");
      double ht = tr.getVar<double>("ht");

      //Get electron and muon for LL study
      int nElectrons = tr.getVar<int>("nElectrons");
      int nMuons = tr.getVar<int>("nMuons");
      int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met );
      //int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met, ht );
      //if( searchbin_id > 61 ) std::cout << "Get it!" << std::endl;
      if(searchbin_id<0) continue;

      if( ((*iter_SSSampleInfos).Tag).find("TTJets") != std::string::npos )
      {
        double TTJetsSF = 0.8;
        //Counting the muon control sample
        if (nElectrons == 0 && nMuons == 1)
        {
          if(searchbin_id>=0){ mucs[searchbin_id]+=(thisweight*TTJetsSF); mucs_NMC[searchbin_id]++; }
        }
        if (nElectrons == 1 && nMuons == 0)
        {
          if(searchbin_id>=0){ elcs[searchbin_id]+=(thisweight*TTJetsSF); elcs_NMC[searchbin_id]++; }
        }
      }
    }//end of inner loop
  }//end of Samples loop

  std::ofstream DC_sb_LL_Header;
  DC_sb_LL_Header.open ( "DC_sb_LL_Header.h" );
  DC_sb_LL_Header << "  const double head_DC_sb_MC_LL_cs[" << NSB << "] = ";
  for( int i_cal = 0 ; i_cal < NSB ; i_cal++ )
  {
    if( i_cal == 0 ) { DC_sb_LL_Header << "{"; }
    DC_sb_LL_Header << mucs[i_cal];
    if( i_cal != NSB-1 ) { DC_sb_LL_Header << ","; }
    else{ DC_sb_LL_Header << "};"; }
  }
  DC_sb_LL_Header << std::endl;
  DC_sb_LL_Header << "  const double head_DC_sb_NMC_LL_cs[" << NSB << "] = ";
  for( int i_cal = 0 ; i_cal < NSB ; i_cal++ )
  {
    if( i_cal == 0 ) { DC_sb_LL_Header << "{"; }
    DC_sb_LL_Header << mucs_NMC[i_cal];
    if( i_cal != NSB-1 ) { DC_sb_LL_Header << ","; }
    else{ DC_sb_LL_Header << "};"; }
  }

  DC_sb_LL_Header.close();
  return ;
}

void LoopDSB( SSSampleWeight& mySSSampleWeight , SSSampleWeight& mySSSampleWeightSignal)
{

  // TFile *pre_trim_file = pre_trim_file = new TFile("signalScan_SMS-T2tt_forHua.root");
  //TFile *pre_trim_file = pre_trim_file = new TFile("signalScan_SMS-T1tttt_forHua.root");
 
  //TH1D * h1_thisSig_totEntries = (TH1D*) pre_trim_file->Get("totEntries_1300_1000");
  //double thisSig_totEntries = h1_thisSig_totEntries->GetBinContent(1);
  //std::cout << "totEntries_1000_50 = " << thisSig_totEntries << std::endl;


  //clock to monitor the run time
  size_t t0 = clock();
  std::vector<SSSampleInfo>::iterator iter_SSSampleInfos;
  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file

  SSDataCard mySSDataCard;
  //std::cout << "Let's do sensitivity study: " << std::endl;
  const int nmetbin=17;
  const int nmt2bin=16;
  double nBGEvents[nmetbin][nmt2bin]; // met,mt2
  double nSEvents[nmetbin][nmt2bin];
  double nSEvents2[nmetbin][nmt2bin];
  double nSEvents3[nmetbin][nmt2bin];
  double nSEvents4[nmetbin][nmt2bin];
  double nmuCSEvents[nmetbin][nmt2bin];

  const int nHTbin=34;
  double nBGEventsMETHT[nmetbin][nHTbin]; // met,mt2
  double nSEventsMETHT[nmetbin][nHTbin];
  double nSEvents2METHT[nmetbin][nHTbin];
  double nSEvents3METHT[nmetbin][nHTbin];
  double nmuCSEventsMETHT[nmetbin][nHTbin];
  double nBGEventsMETHTtop[nmetbin][nHTbin]; // met,mt2
  double nSEventsMETHTtop[nmetbin][nHTbin];
  double nmuCSEventsMETHTtop[nmetbin][nHTbin];

  double nBGEventsMT2HT[nmt2bin][nHTbin];
  double nmuCSEventsMT2HT[nmt2bin][nHTbin];
  double nSEventsMT2HT[nmt2bin][nHTbin];


  for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
  {
    for (int metbinc=0;metbinc<nmetbin;++metbinc)
    {
      nBGEvents[metbinc][mt2binc]=0.0;
      nSEvents[metbinc][mt2binc]=0.0;
      nSEvents2[metbinc][mt2binc]=0.0;
      nSEvents3[metbinc][mt2binc]=0.0;
      nSEvents4[metbinc][mt2binc]=0.0;
      nmuCSEvents[metbinc][mt2binc]=0.0;
    }
    for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
    {
      nBGEventsMT2HT[mt2binc][HTbinc]=0.0;
      nmuCSEventsMT2HT[mt2binc][HTbinc]=0.0;
      nSEventsMT2HT[mt2binc][HTbinc]=0.0;
    }
  }

  for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
  {
    for (int metbinc=0;metbinc<nmetbin;++metbinc)
    {
      nBGEventsMETHT[metbinc][HTbinc]=0.0;
      nSEventsMETHT[metbinc][HTbinc]=0.0;
      nSEvents2METHT[metbinc][HTbinc]=0.0;
      nSEvents3METHT[metbinc][HTbinc]=0.0;
      nmuCSEventsMETHT[metbinc][HTbinc]=0.0;
      nBGEventsMETHTtop[metbinc][HTbinc]=0.0;
      nSEventsMETHTtop[metbinc][HTbinc]=0.0;
      nmuCSEventsMETHTtop[metbinc][HTbinc]=0.0;
    }
  }

  const int ntopcut=2;
  const int nbcut=2;
  //const double minmetcut=250.0;
  const double minmt2cut=200.0;
  const double minHTcut=300.0;
  const double minHTtopcut=0.0;

  //begin to loop over SM events
  for(iter_SSSampleInfos = mySSSampleWeight.SSSampleInfos.begin(); iter_SSSampleInfos != mySSSampleWeight.SSSampleInfos.end(); iter_SSSampleInfos++)
  {    
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_SSSampleInfos).chain);

    double thisweight = (*iter_SSSampleInfos).weight;
    std::cout <<"Sample Type: "<< (*iter_SSSampleInfos).Tag << "; Weight: " << thisweight << std::endl;


    while(tr.getNextEvent())
    {
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      //searchbin variables
      int ntopjets = tr.getVar<int>("nTop");
      int nbotjets = tr.getVar<int>("nBot");
      double mt2 = tr.getVar<double>("mt2");
      double met = tr.getVar<double>("met");
      double HT = tr.getVar<double>("ht");
      double HTtop = tr.getVar<double>("htTops");
      //closure plots variables
      //int njets30 = tr.getVar<int>("nJets30");
      //int njets50 = tr.getVar<int>("nJets50");
      //double ht = tr.getVar<double>("ht");
      int nMuons = tr.getVar<int>("nMuons");
      int nElectrons = tr.getVar<int>("nElectrons");

      bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
      //if(!passLeptVeto) continue;

      //std::cout << "ntopjets = " << ntopjets << std::endl;
      //if (passLeptVeto && ntopjets==1 && nbotjets==1)
      //if (nElectrons==0 && ntopjets==ntopcut && nbotjets==nbcut)
      //if (nElectrons==0 && ntopjets>=ntopcut && nbotjets==nbcut)
      //if (nElectrons==0 && met>250.0 && ntopjets>=ntopcut && nbotjets>=nbcut)
      //if (nElectrons==0 && met>250.0 && ntopjets>=ntopcut && nbotjets==nbcut)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut && met<500.0 && mt2<650.0)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut && met<600.0)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets>=nbcut && met<450.0 && HT<1000.0)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut && met<450.0 && mt2<400.0)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut && mt2<550.0)
      //if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut && met<550.0 && mt2<400.0)
	//if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets>=nbcut && met<450.0)
      if (nElectrons==0 && met>250.0 && ntopjets==ntopcut && nbotjets==nbcut)
      {
	double metcutc=250.0;
	double mt2cutc=200.0;
	double HTcutc=300.0;
	double HTtopcutc=0.0;

	for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
	{
	  mt2cutc=minmt2cut+mt2binc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
	  {
	    if (met>metcutc+metbinc*50.0 && mt2>mt2cutc)
	    {
	      if (passLeptVeto) nBGEvents[metbinc][mt2binc]+=thisweight;
	      if (nMuons==1) nmuCSEvents[metbinc][mt2binc]+=thisweight;
	    }
	  }

	  for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	  {
	    if (HT>minHTcut+HTbinc*50.0 && mt2>mt2cutc)
	    {
	      if (passLeptVeto) nBGEventsMT2HT[mt2binc][HTbinc]+=thisweight;
	      if (nMuons==1) nmuCSEventsMT2HT[mt2binc][HTbinc]+=thisweight;
	    }
	  }

	}

	for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	{
	  HTcutc=minHTcut+HTbinc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
	  {
	    if (met>metcutc+metbinc*50.0 && HT>HTcutc)
	    //if (met<metcutc+metbinc*50.0 && HT<HTcutc)
	    //if (met>metcutc+metbinc*50.0 && HT<HTcutc)
	    //if (met<metcutc+metbinc*50.0 && HT>HTcutc)
	    {
	      if (passLeptVeto) nBGEventsMETHT[metbinc][HTbinc]+=thisweight;
	      if (nMuons==1) nmuCSEventsMETHT[metbinc][HTbinc]+=thisweight;
	    }
	  }
	}

	for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	{
	  HTtopcutc=minHTtopcut+HTbinc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
	  {
	    if (met>metcutc+metbinc*50.0 && HTtop>HTtopcutc)
	    {
	      if (passLeptVeto) nBGEventsMETHTtop[metbinc][HTbinc]+=thisweight;
	      if (nMuons==1) nmuCSEventsMETHTtop[metbinc][HTbinc]+=thisweight;
	    }
	  }
	}


      }
    }//end of inner loop
  }//end of Samples loop


  //begin to loop over signal events
  std::vector<SSSampleInfo>::iterator iter_SSSampleInfosSignal;
  for(iter_SSSampleInfosSignal = mySSSampleWeightSignal.SSSampleInfos.begin(); iter_SSSampleInfosSignal != mySSSampleWeightSignal.SSSampleInfos.end(); iter_SSSampleInfosSignal++)
  {    
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_SSSampleInfosSignal).chain);

    double thisweight = (*iter_SSSampleInfosSignal).weight;
    std::cout <<"Sample Type: "<< (*iter_SSSampleInfosSignal).Tag << "; Weight: " << thisweight << std::endl;


//  for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
//  {
//    for (int metbinc=0;metbinc<nmetbin;++metbinc)
//    {
//      std::cout << "nSEvents[" << metbinc << "][" << mt2binc << "] = " << nSEvents[mt2binc][metbinc] << std::endl;
//    }
//  }


    while(tr.getNextEvent())
    {
      if(tr.getEvtNum()%200000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      //searchbin variables
      int ntopjets = tr.getVar<int>("nTop");
      int nbotjets = tr.getVar<int>("nBot");
      double mt2 = tr.getVar<double>("mt2");
      double met = tr.getVar<double>("met");
      double StopMass = tr.getVar<double>("SusyMotherMass");
      double LSPMass = tr.getVar<double>("SusyLSPMass");
      double HT = tr.getVar<double>("ht");
      double HTtop = tr.getVar<double>("htTops");
     //closure plots variables
      //int njets30 = tr.getVar<int>("nJets30");
      //int njets50 = tr.getVar<int>("nJets50");
      //double ht = tr.getVar<double>("ht");
      int nMuons = tr.getVar<int>("nMuons");
      int nElectrons = tr.getVar<int>("nElectrons");
      bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
      //if(!passLeptVeto) continue;

      //std::cout << "ntopjets = " << ntopjets << std::endl;
      //if (StopMass==1900.0 && LSPMass==50.0 && nMuons==0 && nElectrons==0 && ntopjets==ntopcut && nbotjets==nbcut)
      //if (StopMass==1300.0 && LSPMass==1000.0 && met>250.0 && passLeptVeto && ntopjets>=ntopcut && nbotjets==nbcut)
	//if (StopMass==1900.0 && LSPMass==100.0 && met>250.0 && passLeptVeto && ntopjets>=ntopcut && nbotjets>=nbcut)
      //if (StopMass==1300.0 && LSPMass==1000.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets>=nbcut)
      //if (StopMass==1300.0 && LSPMass==1000.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets>=nbcut && met<450.0 && HT<1000.0)
      //if (StopMass==450.0 && LSPMass==300.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut)
      //if (StopMass==1000.0 && LSPMass==50.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut && mt2<550.0)
      //if (StopMass==1900.0 && LSPMass==100.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut && met<550.0 && mt2<400.0)
      //if (StopMass==700.0 && LSPMass==450.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut && met<500.0 && mt2<650.0)
      //if (StopMass==1000.0 && LSPMass==50.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut && met<600.0)
	//if (StopMass==1900.0 && LSPMass==100.0 && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets>=nbcut && met<450.0)
      //if (StopMass==1900.0 && LSPMass==100.0 && met>250.0 && passLeptVeto && ntopjets>=ntopcut && nbotjets>=nbcut)
      bool isbench=false;
      double corrWeight=thisweight;
      //if (StopMass==1900.0 && LSPMass==100.0) isbench=true;
      //if (StopMass==1900.0 && LSPMass==1100.0) isbench=true;
      //if (StopMass==1300.0 && LSPMass==1000.0) isbench=true;
      if (StopMass==1000.0 && LSPMass==50.0) isbench=true;
      if (StopMass==1000.0 && LSPMass==450.0) isbench=true;
      if (StopMass==700.0 && LSPMass==450.0) isbench=true;
      if (StopMass==450.0 && LSPMass==300.0) isbench=true;
      //if (isbench && met>250.0 && passLeptVeto && ntopjets>=ntopcut && nbotjets>=nbcut)
      //if (isbench && met>250.0 && passLeptVeto && ntopjets>=ntopcut && nbotjets==nbcut)
      if (isbench && met>250.0 && passLeptVeto && ntopjets==ntopcut && nbotjets==nbcut)
      {
	double metcutc=250.0;
	double mt2cutc=200.0;
	double HTcutc=300.0;
	double HTtopcutc=0.0;

	//if (StopMass==1900.0 && LSPMass==1100.0) corrWeight=thisweight*20136.0/20202.0;
	//if (StopMass==1300.0 && LSPMass==1000.0) corrWeight=thisweight*0.0460525/0.00163547*20136.0/42472.0;
	if (StopMass==1000.0 && LSPMass==450.0) corrWeight=thisweight*17885.0/15813.0;
	if (StopMass==700.0 && LSPMass==450.0) corrWeight=thisweight/0.00615134*0.0670476*17885.0/24735.0;
	if (StopMass==450.0 && LSPMass==300.0) corrWeight=thisweight/0.00615134*0.948333*17885.0/347312.0;

	for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
	{
	  mt2cutc=minmt2cut+mt2binc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
  	  {
	    if (met>metcutc+metbinc*50.0 && mt2>mt2cutc)
  	    {
	      //if (StopMass==1900.0 && LSPMass==100.0) nSEvents[metbinc][mt2binc]+=corrWeight;
	      //if (StopMass==1900.0 && LSPMass==1100.0) nSEvents2[metbinc][mt2binc]+=corrWeight;
	      //if (StopMass==1300.0 && LSPMass==1000.0) nSEvents3[metbinc][mt2binc]+=corrWeight;
	      if (StopMass==1000.0 && LSPMass==50.0) nSEvents[metbinc][mt2binc]+=corrWeight;
	      if (StopMass==1000.0 && LSPMass==450.0) nSEvents2[metbinc][mt2binc]+=corrWeight;
	      if (StopMass==700.0 && LSPMass==450.0) nSEvents3[metbinc][mt2binc]+=corrWeight;
	      if (StopMass==450.0 && LSPMass==300.0) nSEvents4[metbinc][mt2binc]+=corrWeight;
	      //if (thisweight>10.0) std::cout << "Warning: thisweight = " << thisweight << std::endl;
	      //if (nSEvents[metbinc][mt2binc]>1000.0) std::cout << "Warning: nSEvents[" << metbinc << "][" << mt2binc << "] = " << nSEvents[metbinc][mt2binc] << std::endl;
	    }
  	  }

	  for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	  {
	    if (HT>minHTcut+HTbinc*50.0 && mt2>mt2cutc)
  	    {
	      nSEventsMT2HT[mt2binc][HTbinc]+=corrWeight;
	    }
	  }

  	}

	for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	{
	  HTcutc=minHTcut+HTbinc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
	  {
	    if (met>metcutc+metbinc*50.0 && HT>HTcutc)
	    //if (met<metcutc+metbinc*50.0 && HT<HTcutc)
	    //if (met>metcutc+metbinc*50.0 && HT<HTcutc)
	    //if (met<metcutc+metbinc*50.0 && HT>HTcutc)
	    {
	      //if (StopMass==1900.0 && LSPMass==100.0) nSEventsMETHT[metbinc][HTbinc]+=corrWeight;
	      //if (StopMass==1900.0 && LSPMass==1100.0) nSEvents2METHT[metbinc][HTbinc]+=corrWeight;
	      //if (StopMass==1300.0 && LSPMass==1000.0) nSEvents3METHT[metbinc][HTbinc]+=corrWeight;
	      if (StopMass==1000.0 && LSPMass==50.0) nSEventsMETHT[metbinc][HTbinc]+=corrWeight;
	      if (StopMass==1000.0 && LSPMass==450.0) nSEvents2METHT[metbinc][HTbinc]+=corrWeight;
	      if (StopMass==700.0 && LSPMass==450.0) nSEvents3METHT[metbinc][HTbinc]+=corrWeight;
	      //if (StopMass==450.0 && LSPMass==300.0) nSEvents4METHT[metbinc][HTbinc]+=corrWeight;
	    }
	  }
	}

	for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
	{
	  HTtopcutc=minHTtopcut+HTbinc*50.0;
	  for (int metbinc=0;metbinc<nmetbin;++metbinc)
	  {
	    if (met>metcutc+metbinc*50.0 && HTtop>HTtopcutc)
	    {
	      nSEventsMETHTtop[metbinc][HTbinc]+=corrWeight;
	    }
	  }
	}


      }

    }//end of inner loop
  }//end of Samples loop


  //std::cout << "nBGEvents[0] = " << nBGEvents[0] << std::endl;
  //std::cout << "nSEvents[0] = " << nSEvents[0] << std::endl;

  // debug (cout)
  //  for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
  //  {
//    for (int metbinc=0;metbinc<nmetbin;++metbinc)
//   {
//      std::cout << "B[" << metbinc << "] = " << nBGEvents[metbinc][mt2binc] << std::endl;
//    }
//    std::cout << std::endl;
//
//    for (int metbinc=0;metbinc<nmetbin;++metbinc)
//    {
//      std::cout << "S[" << metbinc << "] = " << nSEvents[metbinc][mt2binc] << std::endl;
//    }
//    std::cout << std::endl;
//
//    for (int metbinc=0;metbinc<nmetbin;++metbinc)
//    {
//      if (nBGEvents[metbinc][mt2binc]>0)
//      {
//	const double sb=nSEvents[metbinc][mt2binc]/nBGEvents[metbinc][mt2binc];
//	std::cout << "S/B[" << metbinc << "] = " << sb << std::endl;
//      }
//      else std::cout << "S/B[" << metbinc << "] = 0" << std::endl;
//    }
//  }

  FSLHistgram myFSLHistgram;
  myFSLHistgram.BookHistgram( (dir_out + "FSLHistgram.root").c_str() );
  for (int mt2binc=0;mt2binc<nmt2bin;++mt2binc)
  {
    for (int metbinc=0;metbinc<nmetbin;++metbinc)
    {
      myFSLHistgram.h2_B->SetBinContent(metbinc+1,mt2binc+1,nBGEvents[metbinc][mt2binc]);
      myFSLHistgram.h2_S->SetBinContent(metbinc+1,mt2binc+1,nSEvents[metbinc][mt2binc]);
      myFSLHistgram.h2_nmuCS->SetBinContent(metbinc+1,mt2binc+1,nmuCSEvents[metbinc][mt2binc]);
      const double sb=nSEvents[metbinc][mt2binc]/std::sqrt(nSEvents[metbinc][mt2binc]+nBGEvents[metbinc][mt2binc]+0.1*0.1*nSEvents[metbinc][mt2binc]*nSEvents[metbinc][mt2binc]+0.3*0.3*nBGEvents[metbinc][mt2binc]*nBGEvents[metbinc][mt2binc]+1.5);
      const double sb2=nSEvents2[metbinc][mt2binc]/std::sqrt(nSEvents2[metbinc][mt2binc]+nBGEvents[metbinc][mt2binc]+0.1*0.1*nSEvents2[metbinc][mt2binc]*nSEvents2[metbinc][mt2binc]+0.3*0.3*nBGEvents[metbinc][mt2binc]*nBGEvents[metbinc][mt2binc]+1.5);
      const double sb3=nSEvents3[metbinc][mt2binc]/std::sqrt(nSEvents3[metbinc][mt2binc]+nBGEvents[metbinc][mt2binc]+0.1*0.1*nSEvents3[metbinc][mt2binc]*nSEvents3[metbinc][mt2binc]+0.3*0.3*nBGEvents[metbinc][mt2binc]*nBGEvents[metbinc][mt2binc]+1.5);
      const double sb4=nSEvents4[metbinc][mt2binc]/std::sqrt(nSEvents4[metbinc][mt2binc]+nBGEvents[metbinc][mt2binc]+0.1*0.1*nSEvents4[metbinc][mt2binc]*nSEvents4[metbinc][mt2binc]+0.3*0.3*nBGEvents[metbinc][mt2binc]*nBGEvents[metbinc][mt2binc]+1.5);
      if (nBGEvents[metbinc][mt2binc]>0)
      {
	//const double sb=nSEvents[metbinc][mt2binc]/nBGEvents[metbinc][mt2binc];
	const double qstat=2.0*(std::sqrt(nSEvents[metbinc][mt2binc]+nBGEvents[metbinc][mt2binc])-std::sqrt(nBGEvents[metbinc][mt2binc]));
	myFSLHistgram.h2_SOverB->SetBinContent(metbinc+1,mt2binc+1,sb);
	myFSLHistgram.h2_SOverB2->SetBinContent(metbinc+1,mt2binc+1,sb2);
	myFSLHistgram.h2_SOverB3->SetBinContent(metbinc+1,mt2binc+1,sb3);
	myFSLHistgram.h2_SOverB4->SetBinContent(metbinc+1,mt2binc+1,sb4);
	myFSLHistgram.h2_Q->SetBinContent(metbinc+1,mt2binc+1,qstat);
      }
      else
      {
	//myFSLHistgram.h2_SOverB->SetBinContent(metbinc+1,mt2binc+1,0.0);
	myFSLHistgram.h2_SOverB->SetBinContent(metbinc+1,mt2binc+1,sb);
	myFSLHistgram.h2_SOverB2->SetBinContent(metbinc+1,mt2binc+1,sb2);
	myFSLHistgram.h2_SOverB3->SetBinContent(metbinc+1,mt2binc+1,sb3);
	myFSLHistgram.h2_SOverB4->SetBinContent(metbinc+1,mt2binc+1,sb4);
	myFSLHistgram.h2_Q->SetBinContent(metbinc+1,mt2binc+1,0.0);
      }
    }

    for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
    {
      myFSLHistgram.h2_BMT2HT->SetBinContent(mt2binc+1,HTbinc+1,nBGEventsMT2HT[mt2binc][HTbinc]);
      myFSLHistgram.h2_SMT2HT->SetBinContent(mt2binc+1,HTbinc+1,nSEventsMT2HT[mt2binc][HTbinc]);
      myFSLHistgram.h2_nmuCSMT2HT->SetBinContent(mt2binc+1,HTbinc+1,nmuCSEventsMT2HT[mt2binc][HTbinc]);
      const double sb=nSEventsMT2HT[mt2binc][HTbinc]/std::sqrt(nSEventsMT2HT[mt2binc][HTbinc]+nBGEventsMT2HT[mt2binc][HTbinc]+0.1*0.1*nSEventsMT2HT[mt2binc][HTbinc]*nSEventsMT2HT[mt2binc][HTbinc]+0.3*0.3*nBGEventsMT2HT[mt2binc][HTbinc]*nBGEventsMT2HT[mt2binc][HTbinc]+1.5);
      if (nBGEventsMT2HT[mt2binc][HTbinc]>0)
      {
	//const double sb=nSEvents[mt2binc][HTbinc]/nBGEvents[mt2binc][HTbinc];
	//const double qstat=2.0*(std::sqrt(nSEvents[mt2binc][HTbinc]+nBGEvents[mt2binc][HTbinc])-std::sqrt(nBGEvents[mt2binc][HTbinc]));
	myFSLHistgram.h2_SOverBMT2HT->SetBinContent(mt2binc+1,HTbinc+1,sb);
	//myFSLHistgram.h2_Q->SetBinContent(mt2binc+1,HTbinc+1,qstat);
      }
      else
      {
	//myFSLHistgram.h2_SOverB->SetBinContent(mt2binc+1,HTbinc+1,0.0);
	myFSLHistgram.h2_SOverBMT2HT->SetBinContent(mt2binc+1,HTbinc+1,sb);
	//myFSLHistgram.h2_Q->SetBinContent(mt2binc+1,HTbinc+1,0.0);
      }
    }

  }

  for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
  {
    for (int metbinc=0;metbinc<nmetbin;++metbinc)
    {
      myFSLHistgram.h2_BMETHT->SetBinContent(metbinc+1,HTbinc+1,nBGEventsMETHT[metbinc][HTbinc]);
      myFSLHistgram.h2_SMETHT->SetBinContent(metbinc+1,HTbinc+1,nSEventsMETHT[metbinc][HTbinc]);
      myFSLHistgram.h2_nmuCSMETHT->SetBinContent(metbinc+1,HTbinc+1,nmuCSEventsMETHT[metbinc][HTbinc]);
      const double sb=nSEventsMETHT[metbinc][HTbinc]/std::sqrt(nSEventsMETHT[metbinc][HTbinc]+nBGEventsMETHT[metbinc][HTbinc]+0.1*0.1*nSEventsMETHT[metbinc][HTbinc]*nSEventsMETHT[metbinc][HTbinc]+0.3*0.3*nBGEventsMETHT[metbinc][HTbinc]*nBGEventsMETHT[metbinc][HTbinc]+1.5);
      const double sb2=nSEvents2METHT[metbinc][HTbinc]/std::sqrt(nSEvents2METHT[metbinc][HTbinc]+nBGEventsMETHT[metbinc][HTbinc]+0.1*0.1*nSEvents2METHT[metbinc][HTbinc]*nSEvents2METHT[metbinc][HTbinc]+0.3*0.3*nBGEventsMETHT[metbinc][HTbinc]*nBGEventsMETHT[metbinc][HTbinc]+1.5);
      const double sb3=nSEvents3METHT[metbinc][HTbinc]/std::sqrt(nSEvents3METHT[metbinc][HTbinc]+nBGEventsMETHT[metbinc][HTbinc]+0.1*0.1*nSEvents3METHT[metbinc][HTbinc]*nSEvents3METHT[metbinc][HTbinc]+0.3*0.3*nBGEventsMETHT[metbinc][HTbinc]*nBGEventsMETHT[metbinc][HTbinc]+1.5);
      if (nBGEventsMETHT[metbinc][HTbinc]>0)
      {
	//const double sb=nSEventsMETHT[metbinc][HTbinc]/nBGEventsMETHT[metbinc][HTbinc];
	const double qstat=2.0*(std::sqrt(nSEventsMETHT[metbinc][HTbinc]+nBGEventsMETHT[metbinc][HTbinc])-std::sqrt(nBGEventsMETHT[metbinc][HTbinc]));
	myFSLHistgram.h2_SOverBMETHT->SetBinContent(metbinc+1,HTbinc+1,sb);
	myFSLHistgram.h2_SOverBMETHT2->SetBinContent(metbinc+1,HTbinc+1,sb2);
	myFSLHistgram.h2_SOverBMETHT3->SetBinContent(metbinc+1,HTbinc+1,sb3);
	myFSLHistgram.h2_QMETHT->SetBinContent(metbinc+1,HTbinc+1,qstat);
      }
      else
      {
	//myFSLHistgram.h2_SOverBMETHT->SetBinContent(metbinc+1,HTbinc+1,0.0);
	myFSLHistgram.h2_SOverBMETHT->SetBinContent(metbinc+1,HTbinc+1,sb);
	myFSLHistgram.h2_SOverBMETHT2->SetBinContent(metbinc+1,HTbinc+1,sb2);
	myFSLHistgram.h2_SOverBMETHT3->SetBinContent(metbinc+1,HTbinc+1,sb3);
	myFSLHistgram.h2_QMETHT->SetBinContent(metbinc+1,HTbinc+1,0.0);
      }
    }
  }

  for (int HTbinc=0;HTbinc<nHTbin;++HTbinc)
  {
    for (int metbinc=0;metbinc<nmetbin;++metbinc)
    {
      myFSLHistgram.h2_BMETHTtop->SetBinContent(metbinc+1,HTbinc+1,nBGEventsMETHTtop[metbinc][HTbinc]);
      myFSLHistgram.h2_SMETHTtop->SetBinContent(metbinc+1,HTbinc+1,nSEventsMETHTtop[metbinc][HTbinc]);
      myFSLHistgram.h2_nmuCSMETHTtop->SetBinContent(metbinc+1,HTbinc+1,nmuCSEventsMETHTtop[metbinc][HTbinc]);
      // S/sqrt(S+B + 0.1*S*0.1*S + 0.3*B*0.3*B + 1.5)
      const double sb=nSEventsMETHTtop[metbinc][HTbinc]/std::sqrt(nSEventsMETHTtop[metbinc][HTbinc]+nBGEventsMETHTtop[metbinc][HTbinc]+0.1*nSEventsMETHTtop[metbinc][HTbinc]*0.1*nSEventsMETHTtop[metbinc][HTbinc]+0.3*0.3*nBGEventsMETHTtop[metbinc][HTbinc]*nBGEventsMETHTtop[metbinc][HTbinc]+1.5);
      if (nBGEventsMETHTtop[metbinc][HTbinc]>0)
      {
	//S/B
	//const double sb=nSEventsMETHTtop[metbinc][HTbinc]/nBGEventsMETHTtop[metbinc][HTbinc];
	const double qstat=2.0*(std::sqrt(nSEventsMETHTtop[metbinc][HTbinc]+nBGEventsMETHTtop[metbinc][HTbinc])-std::sqrt(nBGEventsMETHTtop[metbinc][HTbinc]));
	myFSLHistgram.h2_SOverBMETHTtop->SetBinContent(metbinc+1,HTbinc+1,sb);
	myFSLHistgram.h2_QMETHTtop->SetBinContent(metbinc+1,HTbinc+1,qstat);
      }
      else
      {
	myFSLHistgram.h2_SOverBMETHTtop->SetBinContent(metbinc+1,HTbinc+1,sb);
	//myFSLHistgram.h2_SOverBMETHTtop->SetBinContent(metbinc+1,HTbinc+1,0.0);
	myFSLHistgram.h2_QMETHTtop->SetBinContent(metbinc+1,HTbinc+1,0.0);
      }
    }
  }


  (myFSLHistgram.oFile)->Write();
  (myFSLHistgram.oFile)->Close();

  return ;
}

void LoopSSAllMC( SSSampleWeight& mySSSampleWeight )
{
  //clock to monitor the run time
  size_t t0 = clock();
  std::vector<SSSampleInfo>::iterator iter_SSSampleInfos;
  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file

  SSDataCard mySSDataCard;
  std::cout << "Let's do sensitivity study: " << std::endl;

  //begin to loop over all events
  for(iter_SSSampleInfos = mySSSampleWeight.SSSampleInfos.begin(); iter_SSSampleInfos != mySSSampleWeight.SSSampleInfos.end(); iter_SSSampleInfos++)
  {    
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_SSSampleInfos).chain);

    double thisweight = (*iter_SSSampleInfos).weight;
    std::cout <<"Sample Type: "<< (*iter_SSSampleInfos).Tag << "; Weight: " << thisweight << std::endl;

    while(tr.getNextEvent())
    {
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      //searchbin variables
      int ntopjets = tr.getVar<int>("nTop");
      int nbotjets = tr.getVar<int>("nBot");
      double mt2 = tr.getVar<double>("mt2");
      double met = tr.getVar<double>("met");
      double ht = tr.getVar<double>("ht");
      //closure plots variables
      //int njets30 = tr.getVar<int>("nJets30");
      //int njets50 = tr.getVar<int>("nJets50");
      //double ht = tr.getVar<double>("ht");

      bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
      if(!passLeptVeto) continue;

      int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met );
      //int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met, ht );

      if(searchbin_id<0)
      { 
        std::cout << "Search bin id problem! Please check!!" << std::endl; 
        std::cout << "NTop : " << ntopjets << "; NBot : " << nbotjets << "; MT2 : " << mt2 << "; MET : " << met << "; HT : " << ht << std::endl; 
        continue; 
      }
      //To separate the LL and Hadtau

      if( ((*iter_SSSampleInfos).Tag).find("TTJets") != std::string::npos )
      { 
        bool isLL = tr.getVar<bool>("isLL");
        double TTJetsSF = 0.8;
        if(isLL) 
        {
          mySSDataCard.DC_sb_MC_LL[searchbin_id] += thisweight*TTJetsSF; 
        }
        else 
        { 
          mySSDataCard.DC_sb_MC_HadTau[searchbin_id] += thisweight*TTJetsSF;
          mySSDataCard.DC_sb_MC_HadTau_NMCforsystunc[searchbin_id]++;
          //std::cout << "NTop: " << ntopjets << ",NBoT: " << nbotjets << ",MT2: " << mt2 << ",MET: " << met << std::endl;
        }
      }
      if( ((*iter_SSSampleInfos).Tag).find("ST_tW") != std::string::npos )
      { 
        bool isLL = tr.getVar<bool>("isLL");
        if(isLL)
        { 
          mySSDataCard.DC_sb_MC_LL[searchbin_id] += thisweight;
        }
        else 
        {
          mySSDataCard.DC_sb_MC_HadTau[searchbin_id] += thisweight;
          mySSDataCard.DC_sb_MC_HadTau_NMCforsystunc[searchbin_id]++;
          //std::cout << "NTop: " << ntopjets << ",NBoT: " << nbotjets << ",MT2: " << mt2 << ",MET: " << met << std::endl;
        }
      }
      if( ((*iter_SSSampleInfos).Tag).find("WJetsToLNu_HT") != std::string::npos )
      {
        bool isLL = tr.getVar<bool>("isLL");
        if(isLL)
        {
          mySSDataCard.DC_sb_MC_LL[searchbin_id] += thisweight;
        }
        else 
        {
          mySSDataCard.DC_sb_MC_HadTau[searchbin_id] += thisweight;
          mySSDataCard.DC_sb_MC_HadTau_NMCforsystunc[searchbin_id]++;
          //std::cout << "NTop: " << ntopjets << ",NBoT: " << nbotjets << ",MT2: " << mt2 << ",MET: " << met << std::endl;
        }
      }
      if( ((*iter_SSSampleInfos).Tag).find("ZJetsToNuNu_HT") != std::string::npos )
      {
        double RNorm = 0.783;  
        int nJets30 = tr.getVar<int>("nJets30"); double njetRWF = 1;
        nJets30<8 ? njetRWF = NJetRweightingFactor[nJets30-1] : njetRWF = NJetRweightingFactor[7];    
        mySSDataCard.DC_sb_MC_Zinv[searchbin_id] += thisweight*njetRWF*RNorm; 
        //effective number of event from Joe, we just have 2 type of weight, ZJetsToNuNu_HT-400To600 and ZJetsToNuNu_HT-600ToInf
        mySSDataCard.DC_sb_MC_Zinv_cs_up[searchbin_id]+=thisweight*njetRWF*RNorm; mySSDataCard.DC_sb_MC_Zinv_cs_dn[searchbin_id]+=thisweight*njetRWF*RNorm*thisweight*njetRWF*RNorm;
      }
      if( ((*iter_SSSampleInfos).Tag).find("QCD_HT") != std::string::npos )
      { 
        mySSDataCard.DC_sb_MC_QCD[searchbin_id] += thisweight; 
      }
      //TTZ, becareful to the negative prediction
      if( ((*iter_SSSampleInfos).Tag).find("TTZTo") != std::string::npos )
      { 
        bool isNegativeWeight = false; isNegativeWeight = tr.getVar<bool>("isNegativeWeight");
        bool isGenZLep = false; isGenZLep = tr.getVar<bool>("isGenZLep");
        bool isGenWLep = false; isGenWLep = tr.getVar<bool>("isGenWLep");
        if((!isGenZLep) && (!isGenWLep))
        {
          isNegativeWeight ? mySSDataCard.DC_sb_MC_TTZ[searchbin_id] -= thisweight : mySSDataCard.DC_sb_MC_TTZ[searchbin_id] += thisweight;
          isNegativeWeight ? mySSDataCard.DC_sb_MC_TTZ_cs[searchbin_id]-- : mySSDataCard.DC_sb_MC_TTZ_cs[searchbin_id]++;
        }
      }
      if( ((*iter_SSSampleInfos).Tag).find("TTWJetsTo") != std::string::npos )
      { 
        bool isNegativeWeight = false; isNegativeWeight = tr.getVar<bool>("isNegativeWeight");
        bool isGenZLep = false; isGenZLep = tr.getVar<bool>("isGenZLep");
        bool isGenWLep = false; isGenWLep = tr.getVar<bool>("isGenWLep");
        if((!isGenZLep) && (!isGenWLep))
        {
          isNegativeWeight ? mySSDataCard.DC_sb_MC_TTZ[searchbin_id] -= thisweight : mySSDataCard.DC_sb_MC_TTZ[searchbin_id] += thisweight;
          isNegativeWeight ? mySSDataCard.DC_sb_MC_TTZ_cs[searchbin_id]-- : mySSDataCard.DC_sb_MC_TTZ_cs[searchbin_id]++;
        }
      }
    }//end of inner loop
  }//end of Samples loop

  //mySSDataCard.printDC_AllFiles("_45BinsLUMI2016Moriond");
  //mySSDataCard.printDC_AllFiles("_45BinsLUMI2016ICHEP");
  //mySSDataCard.printDC_AllFiles("_37BinsLUMI2016Moriond");
  //mySSDataCard.printDC_AllFiles("_37BinsLUMI2016ICHEP");
  //mySSDataCard.printDC_AllFiles("_126BinsLUMI2016ICHEP");
  //mySSDataCard.printDC_AllFiles("_69BinsLUMI2016ICHEP");
  //mySSDataCard.printDC_AllFiles("_59BinsLUMI2016ICHEP");
  //mySSDataCard.printDC_AllFiles("_59BinsLUMI2017Moriond");
  mySSDataCard.printDC_AllFiles("_84BinsLUMI2017Moriond");
  return ;
}

void LoopSignalCard( std::string RunMode )
{
  TChain *chain= new TChain("stopTreeMaker/SSTree");
  if(RunMode.find("T1tttt") != std::string::npos)
  { 
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T1tttt_FastSim_scan_stopFlatNtuples.root"); 
  }
  else if(RunMode.find("T2tt") != std::string::npos)
  { 
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T2tt_FastSim_scan_150to250_stopFlatNtuples.root");
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T2tt_FastSim_scan_250to350_stopFlatNtuples.root"); 
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T2tt_FastSim_scan_350to400_stopFlatNtuples.root"); 
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T2tt_FastSim_scan_400to1200_stopFlatNtuples.root"); 
  }
  else if(RunMode.find("T5ttcc") != std::string::npos)
  {
    chain->Add("root://cmseos.fnal.gov//store/group/lpcsusyhad/hua/Skimmed_2015Nov15/Sensitivity_MC_v11_p2b/SSTrimAndSlimmed_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T5ttcc_FastSim_scan_stopFlatNtuples.root");
  }
  else { std::cout << "bad RunMode for signal card!" << std::endl; return ; }

  TFile *pre_trim_file = 0;
  if(RunMode.find("T1tttt") != std::string::npos) pre_trim_file = new TFile("SignalScanBeforeBaseline/SSSignalScan_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T1tttt_FastSim_scan_stopFlatNtuples.root");
  else if(RunMode.find("T2tt") != std::string::npos) pre_trim_file = new TFile("SignalScanBeforeBaseline/SSSignalScan_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T2tt_FastSim_scan_stopFlatNtuples.root");
  else if(RunMode.find("T5ttcc") != std::string::npos) pre_trim_file = new TFile("SignalScanBeforeBaseline/SSSignalScan_Spring16_80X_Nov_2016_Ntp_v11p0_new_IDs_SMS-T5ttcc_FastSim_scan_stopFlatNtuples.root");
  else { std::cout<<"pre_trim file NOT provided for RunMode : "<<RunMode.c_str()<<std::endl; return; }

  //clock to monitor the run time
  size_t t0 = clock();

  std::cout << "Let's generate signal Data Card!" << std::endl;
  //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
  NTupleReader tr(chain);
  double lumi = LUMI, nevents = chain->GetEntries();
  std::vector<SignalDataCard> mySignalDataCardVec;
  std::vector<SignalDataCard>::iterator iter_mySignalDataCardVec;
  std::cout << RunMode << "; NEvent: " << nevents << std::endl;
  //begin loop of events in SSTree
  while(tr.getNextEvent())
  {
    if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

    bool passLeptVeto = tr.getVar<bool>("passLeptVeto");
    if(!passLeptVeto) continue;

    //searchbin variables
    int ntopjets = tr.getVar<int>("nTop");
    int nbotjets = tr.getVar<int>("nBot");
    double mt2 = tr.getVar<double>("mt2");
    double met = tr.getVar<double>("met");
    double ht = tr.getVar<double>("ht");
    //define the search bin id
    int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met );
    //int searchbin_id = mySearchBins.find_Binning_Index( nbotjets , ntopjets , mt2, met, ht );
    if(searchbin_id<0) continue;

    double SusyMotherMass = tr.getVar<double>("SusyMotherMass");
    double SusyLSPMass    = tr.getVar<double>("SusyLSPMass");
    //if((int)SusyMotherMass%25>1.1 || (int)SusyLSPMass%25>1.1){ std::cout << "Werid mass point not in 25 interveral GeV range!" << " SusyMotherMass : " << SusyMotherMass << "; SusyLSPMass : " << SusyLSPMass << std::endl; }

    bool thisMassPoint = false;
    for(iter_mySignalDataCardVec = mySignalDataCardVec.begin(); iter_mySignalDataCardVec != mySignalDataCardVec.end(); iter_mySignalDataCardVec++)
    {
      thisMassPoint = ((*iter_mySignalDataCardVec).MMass == (int) SusyMotherMass) && ((*iter_mySignalDataCardVec).DMass == (int) SusyLSPMass);
      if(thisMassPoint)
      {
        //(*iter_mySignalDataCardVec).DC_sb_MC_Signal[searchbin_id]+=thisweight;
        (*iter_mySignalDataCardVec).DC_sb_MC_Signal_cs[searchbin_id]++;
        break;
      }
    }
    if(!thisMassPoint)
    {
      SignalDataCard oneSignalDataCard;
      oneSignalDataCard.MMass = (int)SusyMotherMass; oneSignalDataCard.DMass = (int)SusyLSPMass;
      //oneSignalDataCard.DC_sb_MC_Signal[searchbin_id]+=thisweight;
      oneSignalDataCard.DC_sb_MC_Signal_cs[searchbin_id]++;
      mySignalDataCardVec.push_back(oneSignalDataCard);
    }
  }//end of event loop in SSTree
  
  //begin loop of mass point
  for(iter_mySignalDataCardVec = mySignalDataCardVec.begin(); iter_mySignalDataCardVec != mySignalDataCardVec.end(); iter_mySignalDataCardVec++)
  {
    int SusyMotherMass = (*iter_mySignalDataCardVec).MMass;
    int SusyLSPMass = (*iter_mySignalDataCardVec).DMass;

    //get xsec from xSec.h
    double xsec = 0, xsec_err = 0;
    if(RunMode.find("T1tttt") != std::string::npos)
    {
      (*iter_mySignalDataCardVec).DC_SignalType = "T1tttt";
      if( xSecMap_glgl.find(SusyMotherMass) != xSecMap_glgl.end() )
      {
        xsec = xSecMap_glgl.find(SusyMotherMass)->second;
        xsec_err = xSecErrMap_glgl.find(SusyMotherMass)->second * xsec;
      }
    }
    else if(RunMode.find("T2tt") != std::string::npos)
    {
      (*iter_mySignalDataCardVec).DC_SignalType = "T2tt";
      if( xSecMap.find(SusyMotherMass) != xSecMap.end() )
      {
        xsec = xSecMap.find(SusyMotherMass)->second;
        xsec_err = xSecErrMap.find(SusyMotherMass)->second * xsec;
      }
    }
    if(!(xsec>0)){ std::cout << "mass point not in the xSec Map??!!" << " SusyMotherMass : " << SusyMotherMass << "; SusyLSPMass : " << SusyLSPMass << std::endl; continue; }

    //get total event from signal scan file
    //char tmpStr[100];
    //sprintf(tmpStr, "totEntries_%d_%d", SusyMotherMass, SusyLSPMass);
    TH2D * h2_thisSig_totEntries = (TH2D*) pre_trim_file->Get("h_totEvt_xSusyMotherMass_ySusyLSPMass");
    int binx = h2_thisSig_totEntries->GetXaxis()->FindBin(SusyMotherMass); int biny = h2_thisSig_totEntries->GetYaxis()->FindBin(SusyLSPMass);
    double thisSig_totEntries = h2_thisSig_totEntries->GetBinContent(binx,biny);
    if(!(thisSig_totEntries>0)){ std::cout << "mass point do not have entry in signal scan, strange!" << " SusyMotherMass : " << SusyMotherMass << "; SusyLSPMass : " << SusyLSPMass << std::endl; continue; }

    //calculate average weight
    double thisweight = xsec*lumi/thisSig_totEntries;
    (*iter_mySignalDataCardVec).DC_all_MC_Signal_avgweight = thisweight;
    
    for(int i=0; i<NSB; i++)
    {
      (*iter_mySignalDataCardVec).DC_sb_MC_Signal_avgweight[i] = thisweight;
      (*iter_mySignalDataCardVec).DC_sb_MC_Signal[i] = (*iter_mySignalDataCardVec).DC_sb_MC_Signal_cs[i] * thisweight;    
    }

    if( h2_thisSig_totEntries ) delete h2_thisSig_totEntries;

    (*iter_mySignalDataCardVec).print_thisSignalDC();
  }
  pre_trim_file->Close();
}

int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    std::cerr <<"Please give at least 3 arguments " << "RunMode " << " " << "runListMCBG " << " " << "runListMCSG" << "" << "runListMCCS" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./QCD RunMode runlist_QCDMC.txt runlist_Data.txt" << std::endl;
    return -1;
  }

  std::string RunMode = argv[1];
  std::string inputFileList_MC_BG = argv[2];
  std::string inputFileList_MC_SG = argv[3];
  std::string inputFileList_CS = argv[4];

  //std::cout << "The valid run modes are: SSCS SSAllMC SignalCardT1tttt SignalCardT2tt" << std::endl;
  std::cout << "The run mode we have right now is: " << RunMode << std::endl;

  if( RunMode == "SSCS" )
  {
    double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
    double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2

    SSSampleWeight mySSSampleWeightCS;
    mySSSampleWeightCS.SSSampleInfo_push_back( "_TTJets_SingleLeptFromT_"   , 831.76*0.5*TTbar_SingleLept_BR, 53057043, LUMI, 1, inputFileList_CS.c_str() );
    mySSSampleWeightCS.SSSampleInfo_push_back( "_TTJets_SingleLeptFromTbar_", 831.76*0.5*TTbar_SingleLept_BR, 60494823, LUMI, 1, inputFileList_CS.c_str() );
    mySSSampleWeightCS.SSSampleInfo_push_back( "TTJets_DiLept"              , 831.76*TTbar_DiLept_BR        , 30498962, LUMI, 1, inputFileList_CS.c_str() );
    
    LoopSSCS( mySSSampleWeightCS );
    return 0;
  }
  else if( RunMode == "DSB" )
  {
    double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
    double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2
    SSSampleWeight mySSSampleWeightAllMC;
    
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_TTJets_SingleLeptFromT_"   , 831.76*0.5*TTbar_SingleLept_BR, 53057043, LUMI, 0.8, inputFileList_MC_BG.c_str() ); // use 0.8
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_TTJets_SingleLeptFromTbar_", 831.76*0.5*TTbar_SingleLept_BR, 60494823, LUMI, 0.8, inputFileList_MC_BG.c_str() ); // use 0.8
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTJets_DiLept"              , 831.76*TTbar_DiLept_BR        , 30682233, LUMI, 0.8, inputFileList_MC_BG.c_str() ); // use 0.8
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "ST_tW_top"               ,   35.6,    998400    , LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "ST_tW_antitop"           ,   35.6,    985000    , LUMI, 1, inputFileList_MC_BG.c_str() );
    
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-200To400"   ,   359.7,      19591498, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-400To600"   ,   48.91,       7432746, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-600To800"   ,   12.05,      18088165, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-800To1200"  ,   5.501,       7854734, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-1200To2500" ,   1.329,       7023857, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-2500ToInf"  , 0.03216,       2507809, LUMI, 1.21, inputFileList_MC_BG.c_str() );

    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-200To400"  ,    77.67,      25035015, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-400To600"  ,    10.73,       9290017, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-600To800"  ,  0.853*3,       5712221, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-800To1200" ,  0.394*3,       1944423, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-1200To2500", 0.0974*3,        513471, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-2500ToInf" ,0.00230*3,        405752, LUMI, 1.23, inputFileList_MC_BG.c_str() );

    //Be careful! TTZ has negative weight issue!!
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTZToLLNuNu"             , 0.2529, 1744167 - 635909, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTZToQQ"                 , 0.5297,  550282 - 199118, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTWJetsToLNu"            , 0.2043,   191474 - 61199, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTWJetsToQQ"             , 0.4062,  631804 - 201494, LUMI, 1, inputFileList_MC_BG.c_str() );

    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT300to500"  , 366800  , 54706298, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT500to700"  , 29370   , 63337753, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT700to1000" , 6524    , 45453945, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT1000to1500", 1064    ,  15316362, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT1500to2000", 121.5   ,  11650581, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT2000toInf" , 25.42   ,  6007777, LUMI, 1, inputFileList_MC_BG.c_str() );

    SSSampleWeight mySSSampleWeightSignal;
    // xsec are in https://github.com/susy2015/SusyAnaTools/blob/master/Tools/xSec.h
    mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T2tt", 0.00615134, 17885, LUMI, 1, inputFileList_MC_SG.c_str() ); // t2tt 1000 50
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T2tt", 0.00615134, 15813, LUMI, 1, inputFileList_MC_SG.c_str() ); // t2tt 1000 450
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T2tt", 0.0670476, 24735, LUMI, 1, inputFileList_MC_SG.c_str() ); // t2tt 700 450
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T2tt", 0.948333, 347312, LUMI, 1, inputFileList_MC_SG.c_str() ); // t2tt 450 300
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T1tttt", 0.00163547, 0, LUMI, 1, inputFileList_MC_SG.c_str() ); // t1tttt 1900 50
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T1tttt", 0.00163547, 20136, LUMI, 1, inputFileList_MC_SG.c_str() ); // t1tttt 1900 100
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T1tttt", 0.00163547, 20202, LUMI, 1, inputFileList_MC_SG.c_str() ); // t1tttt 1900 1100
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T1tttt", 0.0460525, 42472, LUMI, 1, inputFileList_MC_SG.c_str() ); // t1tttt 1300 1000
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T5ttcc", 1.0, 100000, LUMI, 1, inputFileList_MC_SG.c_str() ); // t5ttcc 1700 0
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T5ttcc", 1.0, 100000, LUMI, 1, inputFileList_MC_SG.c_str() ); // t5ttcc 1800 200
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T5ttcc", 1.0, 100000, LUMI, 1, inputFileList_MC_SG.c_str() ); // t5ttcc 1800 1000
    //mySSSampleWeightSignal.SSSampleInfo_push_back( "SMS-T5ttcc", 1.0, 100000, LUMI, 1, inputFileList_MC_SG.c_str() ); // t5ttcc 1300 1100

    LoopDSB( mySSSampleWeightAllMC, mySSSampleWeightSignal);
    return 0;
  }
  else if( RunMode == "SSAllMC" )
  {
    double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
    double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2

    //sample needed in the basic check loop
    SSSampleWeight mySSSampleWeightAllMC;
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_TTJets_SingleLeptFromT_"   , 831.76*0.5*TTbar_SingleLept_BR, 53057043, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_TTJets_SingleLeptFromTbar_", 831.76*0.5*TTbar_SingleLept_BR, 60494823, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTJets_DiLept"              , 831.76*TTbar_DiLept_BR        , 30682233, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "ST_tW_top"               ,   35.6,    998400    , LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "ST_tW_antitop"           ,   35.6,    985000    , LUMI, 1, inputFileList_MC_BG.c_str() );
 
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-200To400"   ,   359.7,      19591498, LUMI, 1.21, inputFileList_MC_BG.c_str() );   
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-400To600"   ,   48.91,       7432746, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-600To800"   ,   12.05,      18088165, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-800To1200"  ,   5.501,       7854734, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-1200To2500" ,   1.329,       7023857, LUMI, 1.21, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_WJetsToLNu_HT-2500ToInf"  , 0.03216,       2507809, LUMI, 1.21, inputFileList_MC_BG.c_str() );

    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-200To400"  ,    77.67,      25035015, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-400To600"  ,    10.73,       9290017, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-600To800"  ,  0.853*3,       5712221, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-800To1200" ,  0.394*3,       1944423, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-1200To2500", 0.0974*3,        513471, LUMI, 1.23, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "_ZJetsToNuNu_HT-2500ToInf" ,0.00230*3,        405752, LUMI, 1.23, inputFileList_MC_BG.c_str() );

    //Be careful! TTZ has negative weight issue!!
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTZToLLNuNu"             , 0.2529, 1744167 - 635909, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTZToQQ"                 , 0.5297,  550282 - 199118, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTWJetsToLNu"            , 0.2043,   191474 - 61199, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "TTWJetsToQQ"             , 0.4062,  631804 - 201494, LUMI, 1, inputFileList_MC_BG.c_str() );

    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT300to500"  , 366800  , 54706298, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT500to700"  , 29370   , 19199088, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT700to1000" , 6524    , 29294258, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT1000to1500", 1064    , 15316362, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT1500to2000", 121.5   ,  7803965, LUMI, 1, inputFileList_MC_BG.c_str() );
    mySSSampleWeightAllMC.SSSampleInfo_push_back( "QCD_HT2000toInf" , 25.42   ,  6007777, LUMI, 1, inputFileList_MC_BG.c_str() );

    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "SMS-T1tttt_mGluino-1200_mLSP-800",  0.0856418,  147194, LUMI, 1, inputFileList_MC_SG.c_str() );
    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "SMS-T1tttt_mGluino-1500_mLSP-100",  0.0141903,  103140, LUMI, 1, inputFileList_MC_SG.c_str() );
    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "SMS-T2tt_mStop-500_mLSP-325"     ,  0.51848  ,  388207, LUMI, 1, inputFileList_MC_SG.c_str() );
    //mySSSampleWeightAllMC.SSSampleInfo_push_back( "SMS-T2tt_mStop-850_mLSP-100"     ,  0.0189612,  240685, LUMI, 1, inputFileList_MC_SG.c_str() );

    LoopSSAllMC( mySSSampleWeightAllMC );
    return 0;
  }
  else if( RunMode == "SignalCardT2tt" || RunMode == "SignalCardT1tttt" )
  {
    LoopSignalCard( RunMode );
  }
  else
  {
    std::cout << "Invalide RunMode!!" << std::endl;
    return 0;
  }

  return 0;
}
