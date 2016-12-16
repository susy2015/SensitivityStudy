#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TList.h"
#include "TF1.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TMatrixDSym.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/searchBins.h"

#include "SBGeometry.h"
#include "SSReWeighting.h"

#include "DC_sb_LL_Header.h"
std::string dir_out = "";

//############finish the definition of class AccRecoEffs######################
//baseline cut function definition
static BaselineVessel *myBaselineVessel;
void mypassBaselineFunc(NTupleReader& tr)
{
  (*myBaselineVessel)(tr);
}

SBGeometry mySBGeometry;

class SSDataCard
{
 public:
  //fake variables for all background, LL, HadTau, Zinv, QCD, TTZ, Rare for rate, stat unc, and syst unc
  //rate and cs: basic input from MCs
  double DC_sb_MC_Data[NSB] = {0};
  double DC_sb_MC_LL[NSB] = {0}, DC_sb_MC_HadTau[NSB] = {0}, DC_sb_MC_Zinv[NSB] = {0}, DC_sb_MC_QCD[NSB] = {0}, DC_sb_MC_TTZ[NSB] = {0}, DC_sb_MC_Rare[NSB] = {0}; 
  double DC_sb_MC_LL_cs[NSB] = {0}, DC_sb_MC_HadTau_cs[NSB] = {0}, DC_sb_MC_Zinv_cs[NSB] = {0}, DC_sb_MC_QCD_cs[NSB] = {0}, DC_sb_MC_TTZ_cs[NSB] = {0}, DC_sb_MC_Rare_cs[NSB] = {0};
  //avgweight and stat unc can be easily faked from rate and cs, but systunc is tricky
  double DC_sb_MC_LL_avgweight[NSB] = {0}, DC_sb_MC_HadTau_avgweight[NSB] = {0}, DC_sb_MC_Zinv_avgweight[NSB] = {0}, DC_sb_MC_QCD_avgweight[NSB] = {0}, DC_sb_MC_TTZ_avgweight[NSB] = {0}, DC_sb_MC_Rare_avgweight[NSB] = {0};
  double DC_sb_MC_LL_statunc[NSB] = {0}, DC_sb_MC_HadTau_statunc[NSB] = {0}, DC_sb_MC_Zinv_statunc[NSB] = {0}, DC_sb_MC_QCD_statunc[NSB] = {0}, DC_sb_MC_TTZ_statunc[NSB] = {0}, DC_sb_MC_Rare_statunc[NSB] = {0};
  double DC_sb_MC_LL_systunc[NSB] = {0}, DC_sb_MC_HadTau_systunc[NSB] = {0}, DC_sb_MC_Zinv_systunc[NSB] = {0}, DC_sb_MC_QCD_systunc[NSB] = {0}, DC_sb_MC_TTZ_systunc[NSB] = {0}, DC_sb_MC_Rare_systunc[NSB] = {0};
 
  //special hadtau syst mistag, and MC N for LL non closure syst unc, and up dn cs for zinv
  double DC_sb_NMC_LL_cs[NSB] = {0};
  double DC_sb_MC_HadTau_systunc_mistag[NSB] = {0}, DC_sb_MC_HadTau_NMCforsystunc[NSB] = {0};
  double DC_sb_MC_Zinv_cs_up[NSB]={0}, DC_sb_MC_Zinv_cs_dn[NSB] = {0};
  void printDC_AllFiles(std::string sbtag);
 private:
  void fake_avg_uncs();
};

void SSDataCard::fake_avg_uncs()
{
  std::cout << "Faking syst uncs in Data card!" << std::endl;
  //SBGeometry mySBGeometry;
  //SearchBins mySearchBins("SB_69_2016");
  SearchBins mySearchBins("SB_59_2016");

  for(int i=0;i<NSB;i++)
  {
    DC_sb_MC_Data[i] = DC_sb_MC_LL[i] + DC_sb_MC_HadTau[i] + DC_sb_MC_Zinv[i] + /*DC_sb_MC_QCD[i] +*/ DC_sb_MC_TTZ[i] /*+ DC_sb_MC_Rare[i]*/; 
    //set LL CS numbers from header file
    DC_sb_MC_LL_cs[i] = head_DC_sb_MC_LL_cs[i]; DC_sb_NMC_LL_cs[i] = head_DC_sb_NMC_LL_cs[i];
    //set zinv numbers from up dn variables 
    DC_sb_MC_Zinv_cs_dn[i]>0 ? DC_sb_MC_Zinv_cs[i] = DC_sb_MC_Zinv_cs_up[i]*DC_sb_MC_Zinv_cs_up[i]/DC_sb_MC_Zinv_cs_dn[i] : DC_sb_MC_Zinv_cs[i] = 0;

    //set avgweight from cs and rate, need for LL Zinv and TTZ, no need for qcd and hadtau
    DC_sb_MC_LL_cs[i]>0 ? DC_sb_MC_LL_avgweight[i] = DC_sb_MC_LL[i]/DC_sb_MC_LL_cs[i] : DC_sb_MC_LL_avgweight[i] = DC_sb_MC_LL_avgweight[i-1];
    DC_sb_MC_Zinv_cs[i]>0 ? DC_sb_MC_Zinv_avgweight[i] = DC_sb_MC_Zinv[i]/DC_sb_MC_Zinv_cs[i] : DC_sb_MC_Zinv_avgweight[i] = DC_sb_MC_Zinv_avgweight[i-1];
    DC_sb_MC_TTZ_cs[i]>0 ? DC_sb_MC_TTZ_avgweight[i] = DC_sb_MC_TTZ[i]/DC_sb_MC_TTZ_cs[i] : DC_sb_MC_TTZ_avgweight[i] = DC_sb_MC_TTZ_avgweight[i-1];
    if(DC_sb_MC_LL_avgweight[i]<=0) DC_sb_MC_LL_avgweight[i] = DC_sb_MC_LL_avgweight[i-1];
    if(DC_sb_MC_Zinv_avgweight[i]<=0) DC_sb_MC_Zinv_avgweight[i] = DC_sb_MC_Zinv_avgweight[i-1];
    if(DC_sb_MC_TTZ_avgweight[i]<=0) DC_sb_MC_TTZ_avgweight[i] = DC_sb_MC_TTZ_avgweight[i-1];    

    //set stat unc. stat unc need for QCD and hadtau, no need for LL,Zinv and TTZ
    DC_sb_MC_HadTau[i]>0 ? DC_sb_MC_HadTau_statunc[i] = 0.5*std::sqrt(DC_sb_MC_HadTau[i])/DC_sb_MC_HadTau[i] : DC_sb_MC_HadTau_statunc[i] = 0;
    if(DC_sb_MC_HadTau_statunc[i]>1) DC_sb_MC_HadTau_statunc[i] = 0.998;
    //seems we still need LL stat unc added in herer
    round(DC_sb_MC_LL_cs[i])>0 ? DC_sb_MC_LL_statunc[i] = 1/std::sqrt(round(DC_sb_MC_LL_cs[i])) : DC_sb_MC_LL_statunc[i] = 0;
    DC_sb_MC_QCD[i]>0 ? DC_sb_MC_QCD_statunc[i] = std::sqrt(DC_sb_MC_QCD[i])/DC_sb_MC_QCD[i] : DC_sb_MC_QCD_statunc[i] = 0;

    //syst unc setting, need for all BGs, need to be fixed!
    //LL just set for closure unc, relation proved by Florent's plot
    DC_sb_NMC_LL_cs[i]>0 ? DC_sb_MC_LL_systunc[i] = 2.69/std::sqrt(DC_sb_NMC_LL_cs[i]) : DC_sb_MC_LL_systunc[i] = DC_sb_MC_LL_systunc[i-1];
    //if(DC_sb_MC_LL_systunc[i]>1) DC_sb_MC_LL_systunc[i]=0.998;
    //Had Tau
    DC_sb_MC_HadTau_NMCforsystunc[i]>0 ? DC_sb_MC_HadTau_systunc[i] = 1.66/std::sqrt(DC_sb_MC_HadTau_NMCforsystunc[i]) : DC_sb_MC_HadTau_systunc[i] = DC_sb_MC_HadTau_systunc[i-1];
    //if(DC_sb_MC_HadTau_systunc[i]>1) DC_sb_MC_HadTau_systunc[i]=0.998;
    SearchBins::searchBinDef outBinDef; mySearchBins.find_BinBoundaries( i, outBinDef );
		outBinDef.bJet_lo_>=3 ? DC_sb_MC_HadTau_systunc_mistag[i] = 0.1 : DC_sb_MC_HadTau_systunc_mistag[i] = 0.05; 
    //SBBoundaries outBinDef; mySBGeometry.SBIDToBinBoundaries( i, outBinDef );
    //outBinDef.nbot_lo>=3 ? DC_sb_MC_HadTau_systunc_mistag[i] = 0.1 : DC_sb_MC_HadTau_systunc_mistag[i] = 0.05; 

    //Zinv
		DC_sb_MC_Zinv_cs[i]>0 ? DC_sb_MC_Zinv_systunc[i] = 10*std::sqrt(DC_sb_MC_Zinv_cs[i])/DC_sb_MC_Zinv_cs[i] : DC_sb_MC_Zinv_systunc[i] = 2.0;
    if(DC_sb_MC_Zinv_systunc[i]>2.0) DC_sb_MC_Zinv_systunc[i]=2.0; 
    //if( DC_sb_MC_Zinv_systunc[i]>1 ) DC_sb_MC_Zinv_systunc[i]=0.998;
    //QCD, 150% of prediction
    DC_sb_MC_QCD_systunc[i] = 1.5;
    DC_sb_MC_TTZ_systunc[i] = 0.3*std::abs(DC_sb_MC_TTZ[i]);
  }
  return ;
}

void SSDataCard::printDC_AllFiles(std::string sbtag)
{
  fake_avg_uncs();
  
  std::ofstream Datafile (("_Data"+sbtag+".txt").c_str());
  if (Datafile.is_open())
  { 
    Datafile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    Datafile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    Datafile << "sample = signal\n";
    Datafile << "channel = "; for(int i=0;i<NSB;i++){ Datafile << "bin" << i+1 << " "; } Datafile << "\n";
    
    Datafile << "\n# Predicted central numbers (need from all backgrounds)\n";
    Datafile << "rate = "; for(int i=0;i<NSB;i++){ Datafile << round(DC_sb_MC_Data[i]) << " "; } Datafile << "\n";
    
    Datafile.close();
  }
  else std::cout << "Unable to open Datafile";

  std::ofstream LLfile (("_LL"+sbtag+".txt").c_str());
  if (LLfile.is_open())
  {
    LLfile << "# The words after \"#\" are comments. No need to remove them.\n";
    LLfile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    LLfile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    LLfile << "sample = lostle     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    LLfile << "channel = "; for(int i=0;i<NSB;i++){ LLfile << "bin" << i+1 << " "; } LLfile << "\n";

    LLfile << "\n# Predicted central numbers (need from all backgrounds)\n";
    LLfile << "rate = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL[i] << " "; } LLfile << "\n";

    LLfile << "\n# Control sample events for lost lepton; Raw MC yields for Zinv and ttZ; No need from had. tau and QCD (will be ignored).\n";
    LLfile << "# cs_event = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL_cs[i] << " "; } LLfile << "\n";
    //round up to integer
    LLfile << "cs_event = "; for(int i=0;i<NSB;i++){ LLfile << round(DC_sb_MC_LL_cs[i]) << " "; } LLfile << "\n";

    LLfile << "\n# Average weight for lost lepton, Zinv and ttZ: avg_weight x cs_event = rate. Will be ignored for had. tau and QCD.\n";
    LLfile << "avg_weight = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL_avgweight[i] << " "; } LLfile << "\n";

    LLfile << "\n# All the uncertainties in absolute numbers for both stat and syst uncertainties. Needed from QCD and hadronic tau; No need from lost lepton, Zinv and ttZ (will be ignored). For symmetric uncertainties, put the up and dn to be the same.\n";
    LLfile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL_statunc[i]/*"No Need!!"*/ << " "; /*break;*/ } LLfile << "\n";
    LLfile << "stat_unc_dn = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL_statunc[i]/*"No Need!!"*/ << " "; /*break;*/ } LLfile << "\n";

    LLfile << "\n# List of all the systematical uncertainties. Do not need combine them. The \"pdf\", \"blackhole\", \"darkmatter\", \"susy\" are keywords to label the different systematic sources (use meaningful names or make comments).\n";
    LLfile << "syst_unc_closure_up = "; for(int i=0;i<NSB;i++){ LLfile << DC_sb_MC_LL_systunc[i] << " "; } LLfile << "\n";
    LLfile << "syst_unc_closure_dn = "; for(int i=0;i<NSB;i++){ if(DC_sb_MC_LL_systunc[i]>0.999) LLfile << "0.998 "; else LLfile << DC_sb_MC_LL_systunc[i] << " "; } LLfile << "\n";
    LLfile << "syst_unc_dimuon_up = " ; for(int i=0;i<NSB;i++){ LLfile << 0.003 << " "; } LLfile << "\n";
    LLfile << "syst_unc_dimuon_dn = " ; for(int i=0;i<NSB;i++){ LLfile << 0.003 << " "; } LLfile << "\n";
    LLfile << "syst_unc_diele_up = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.012 << " "; } LLfile << "\n";
    LLfile << "syst_unc_diele_dn = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.012 << " "; } LLfile << "\n";
    LLfile << "syst_unc_purity_up = " ; for(int i=0;i<NSB;i++){ LLfile << 0.00 << " "; } LLfile << "\n";
    LLfile << "syst_unc_purity_dn = " ; for(int i=0;i<NSB;i++){ LLfile << 0.02 << " "; } LLfile << "\n";
    LLfile << "syst_unc_mt_up = "     ; for(int i=0;i<NSB;i++){ LLfile << 0.01 << " "; } LLfile << "\n";
    LLfile << "syst_unc_mt_dn = "     ; for(int i=0;i<NSB;i++){ LLfile << 0.01 << " "; } LLfile << "\n";
    LLfile << "syst_unc_acc_up = "    ; for(int i=0;i<NSB;i++){ LLfile << 0.10 << " "; } LLfile << "\n";
    LLfile << "syst_unc_acc_dn = "    ; for(int i=0;i<NSB;i++){ LLfile << 0.10 << " "; } LLfile << "\n";
    LLfile << "syst_unc_muiso_up = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_muiso_dn = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_eiso_up = "   ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_eiso_dn = "   ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_mureco_up = " ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_mureco_dn = " ; for(int i=0;i<NSB;i++){ LLfile << 0.05 << " "; } LLfile << "\n";
    LLfile << "syst_unc_ereco_up = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.08 << " "; } LLfile << "\n";
    LLfile << "syst_unc_ereco_dn = "  ; for(int i=0;i<NSB;i++){ LLfile << 0.08 << " "; } LLfile << "\n";
    LLfile << "syst_unc_isotrk_up = " ; for(int i=0;i<NSB;i++){ LLfile << 0.08 << " "; } LLfile << "\n";
    LLfile << "syst_unc_isotrk_dn = " ; for(int i=0;i<NSB;i++){ LLfile << 0.08 << " "; } LLfile << "\n";
    LLfile.close();
  }
  else std::cout << "Unable to open LLfile";

  std::ofstream HadTaufile (("_HadTau"+sbtag+".txt").c_str());
  if (HadTaufile.is_open())
  {
    HadTaufile << "# The words after \"#\" are comments. No need to remove them.\n";
    HadTaufile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    HadTaufile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    HadTaufile << "sample = hadtau     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    HadTaufile << "channel = "; for(int i=0;i<NSB;i++){ HadTaufile << "bin" << i+1 << " "; } HadTaufile << "\n";

    HadTaufile << "\n# Predicted central numbers (need from all backgrounds)\n";
    HadTaufile << "rate = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau[i] << " "; } HadTaufile << "\n";

    HadTaufile << "\n# Control sample events for lost lepton; Raw MC yields for Zinv and ttZ; No need from had. tau and QCD (will be ignored).\n";
    HadTaufile << "cs_event = "; for(int i=0;i<NSB;i++){ HadTaufile << /*DC_sb_MC_HadTau_cs[i]*/"No Need!!" << " "; break; } HadTaufile << "\n";

    HadTaufile << "\n# Average weight for lost lepton, Zinv and ttZ: avg_weight x cs_event = rate. Will be ignored for had. tau and QCD.\n";
    HadTaufile << "avg_weight = "; for(int i=0;i<NSB;i++){ HadTaufile << /*DC_sb_MC_HadTau_avgweight[i]*/"No Need!!" << " "; break; } HadTaufile << "\n";

    HadTaufile << "\n# All the uncertainties in absolute numbers for both stat and syst uncertainties. Needed from QCD and hadronic tau; No need from lost lepton, Zinv and ttZ (will be ignored). For symmetric uncertainties, put the up and dn to be the same.\n";
    HadTaufile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_statunc[i] << " "; } HadTaufile << "\n";
    HadTaufile << "stat_unc_dn = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_statunc[i] << " "; } HadTaufile << "\n";

    HadTaufile << "\n# List of all the systematical uncertainties. Do not need combine them. The \"pdf\", \"blackhole\", \"darkmatter\", \"susy\" are keywords to label the different systematic sources (use meaningful names or make comments).\n";
    HadTaufile << "#Acc\n";
    HadTaufile << "syst_unc_pdf_up = "    ; for(int i=0;i<NSB;i++){ HadTaufile << 0.05 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_pdf_dn = "    ; for(int i=0;i<NSB;i++){ HadTaufile << 0.05 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_MCScale_up = "; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_MCScale_dn = "; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Mt cut (met scale variation)\n";
    HadTaufile << "syst_unc_Mt_up = "     ; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_Mt_dn = "     ; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#tau->mu contamination factor(statistical precision)\n";
    HadTaufile << "syst_unc_taumu_up = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_taumu_dn = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Isolated track veto\n";
    HadTaufile << "syst_unc_isotrk_up = " ; for(int i=0;i<NSB;i++){ HadTaufile << 0.05 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_isotrk_dn = " ; for(int i=0;i<NSB;i++){ HadTaufile << 0.05 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Muon Effciency\n";
    HadTaufile << "syst_unc_mureco_up = " ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_mureco_dn = " ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_muiso_up = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_muiso_dn = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#B mis-tag rate of tau\n";
    HadTaufile << "syst_unc_mistag_up = " ; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_systunc_mistag[i] << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_mistag_dn = " ; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_systunc_mistag[i] << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Tau Template (JES)\n";
    HadTaufile << "syst_unc_temp_up = "   ; for(int i=0;i<NSB;i++){ HadTaufile << 0.02 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_temp_dn = "   ; for(int i=0;i<NSB;i++){ HadTaufile << 0.02 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Closure\n";
    HadTaufile << "syst_unc_closure_up = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_systunc[i] << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_closure_dn = "; for(int i=0;i<NSB;i++){ if(DC_sb_MC_HadTau_systunc[i]>0.999) HadTaufile << "0.998 "; else HadTaufile << DC_sb_MC_HadTau_systunc[i] << " "; } HadTaufile << "\n";
    HadTaufile << "\n#lostlepton overlap\n";
    HadTaufile << "syst_unc_llovr_up = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0.024 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_llovr_dn = "  ; for(int i=0;i<NSB;i++){ HadTaufile << 0.024 << " "; } HadTaufile << "\n";
    HadTaufile << "\n#Trigger Eff\n";
    HadTaufile << "syst_unc_trg_up = "    ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    HadTaufile << "syst_unc_trg_dn = "    ; for(int i=0;i<NSB;i++){ HadTaufile << 0.01 << " "; } HadTaufile << "\n";
    //special piece for HadTau data card
    HadTaufile << "\n#For Hongxuan, Input for estimate hadtau syst unc\n";
    HadTaufile << "# exp_yield              = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau[i] << " "; } HadTaufile << "\n";
    HadTaufile << "# exp_yield_rel_stat_unc = "; for(int i=0;i<NSB;i++){ HadTaufile << 1/std::sqrt(DC_sb_MC_HadTau_NMCforsystunc[i]) << " "; } HadTaufile << "\n";
    HadTaufile << "# 1.12 alpha             = "; for(int i=0;i<NSB;i++){ HadTaufile << 1.12/std::sqrt(DC_sb_MC_HadTau_NMCforsystunc[i]) << " "; } HadTaufile << "\n";
    HadTaufile << "# NMC                    = "; for(int i=0;i<NSB;i++){ HadTaufile << DC_sb_MC_HadTau_NMCforsystunc[i] << " "; } HadTaufile << "\n";

    HadTaufile.close();
  }
  else std::cout << "Unable to open HadTaufile";

  std::ofstream Zinvfile (("_Zinv"+sbtag+".txt").c_str());
  if (Zinvfile.is_open())
  {
    Zinvfile << "# The words after \"#\" are comments. No need to remove them.\n";
    Zinvfile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    Zinvfile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    Zinvfile << "sample = zinv     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    Zinvfile << "channel = "; for(int i=0;i<NSB;i++){ Zinvfile << "bin" << i+1 << " "; } Zinvfile << "\n";
    
    Zinvfile << "\n# Predicted central numbers (need from all backgrounds)\n";
    Zinvfile << "rate = "; for(int i=0;i<NSB;i++){ Zinvfile << DC_sb_MC_Zinv[i] << " "; } Zinvfile << "\n";
    
    Zinvfile << "\n# Control sample events for lost lepton; Raw MC yields for Zinv and ttZ; No need from had. tau and QCD (will be ignored).\n";
    Zinvfile << "#cs_event = "; for(int i=0;i<NSB;i++){ Zinvfile << DC_sb_MC_Zinv_cs[i] << " "; } Zinvfile << "\n";
    Zinvfile << "cs_event = "; for(int i=0;i<NSB;i++){ Zinvfile << round(DC_sb_MC_Zinv_cs[i]) << " "; } Zinvfile << "\n";

    Zinvfile << "\n# Average weight for lost lepton, Zinv and ttZ: avg_weight x cs_event = rate. Will be ignored for had. tau and QCD.\n";
    Zinvfile << "avg_weight = "; for(int i=0;i<NSB;i++){ Zinvfile << DC_sb_MC_Zinv_avgweight[i] << " "; } Zinvfile << "\n";
    
    Zinvfile << "\n# All the uncertainties in absolute numbers for both stat and syst uncertainties. Needed from QCD and hadronic tau; No need from lost lepton, Zinv and ttZ (will be ignored). For symmetric uncertainties, put the up and dn to be the same.\n";
    Zinvfile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ Zinvfile << /*DC_sb_MC_Zinv_statunc[i]*/"No Need!!" << " "; break; } Zinvfile << "\n";
    Zinvfile << "stat_unc_dn = "; for(int i=0;i<NSB;i++){ Zinvfile << /*DC_sb_MC_Zinv_statunc[i]*/"No Need!!" << " "; break; } Zinvfile << "\n";
    
    Zinvfile << "\n# List of all the systematical uncertainties. Do not need combine them. The \"pdf\", \"blackhole\", \"darkmatter\", \"susy\" are keywords to label the different systematic sources (use meaningful names or make comments).\n";
    Zinvfile << "syst_unc_norm_up = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0.3 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_norm_dn = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0.3 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_shape_central_up = "; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_shape_central_dn = "; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_shape_stat_up = "   ; for(int i=0;i<NSB;i++){ Zinvfile << DC_sb_MC_Zinv_systunc[i] << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_shape_stat_dn = "   ; for(int i=0;i<NSB;i++){ if(DC_sb_MC_Zinv_systunc[i]>0.999) Zinvfile << "0.998 "; else Zinvfile << DC_sb_MC_Zinv_systunc[i] << " ";} Zinvfile << "\n";
    Zinvfile << "syst_unc_jeu_up = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_jeu_dn = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_meu_up = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_meu_dn = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_scale_up = "        ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_scale_dn = "        ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_pdf_up = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_pdf_dn = "          ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_trig_up = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_trig_dn = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_btag_up = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_btag_dn = "         ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_bmistag_up = "      ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile << "syst_unc_bmistag_dn = "      ; for(int i=0;i<NSB;i++){ Zinvfile << 0 << " "; } Zinvfile << "\n";
    Zinvfile.close();
  }
  else std::cout << "Unable to open Zinvfile";
  
  std::ofstream QCDfile (("_QCD"+sbtag+".txt").c_str());
  if (QCDfile.is_open())
  {
    QCDfile << "# The words after \"#\" are comments. No need to remove them.\n";
    QCDfile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    QCDfile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    QCDfile << "sample = qcd     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    QCDfile << "channel = "; for(int i=0;i<NSB;i++){ QCDfile << "bin" << i+1 << " "; } QCDfile << "\n";

    QCDfile << "\n# Predicted central numbers (need from all backgrounds)\n";
    QCDfile << "rate = "; for(int i=0;i<NSB;i++){ QCDfile << DC_sb_MC_QCD[i] << " "; } QCDfile << "\n";

    QCDfile << "\n# Control sample events for lost lepton; Raw MC yields for Zinv and ttZ; No need from had. tau and QCD (will be ignored).\n";
    QCDfile << "cs_event = "; for(int i=0;i<NSB;i++){ QCDfile << /*DC_sb_MC_QCD_cs[i]*/"No Need!!" << " "; break; } QCDfile << "\n";

    QCDfile << "\n# Average weight for lost lepton, Zinv and ttZ: avg_weight x cs_event = rate. Will be ignored for had. tau and QCD.\n";
    QCDfile << "avg_weight = "; for(int i=0;i<NSB;i++){ QCDfile << /*DC_sb_MC_QCD_avgweight[i]*/"No Need!!" << " "; break; } QCDfile << "\n";

    QCDfile << "\n# All the uncertainties in absolute numbers for both stat and syst uncertainties. Needed from QCD and hadronic tau; No need from lost lepton, Zinv and ttZ (will be ignored). For symmetric uncertainties, put the up and dn to be the same.\n";
    QCDfile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ QCDfile << DC_sb_MC_QCD_statunc[i] << " "; } QCDfile << "\n";
    QCDfile << "stat_unc_dn = "; for(int i=0;i<NSB;i++){ QCDfile << DC_sb_MC_QCD_statunc[i] << " "; } QCDfile << "\n";

    QCDfile << "\n# List of all the systematical uncertainties. Do not need combine them. The \"pdf\", \"blackhole\", \"darkmatter\", \"susy\" are keywords to label the different systematic sources (use meaningful names or make comments).\n";
    QCDfile << "syst_unc_all_up = "; for(int i=0;i<NSB;i++) { QCDfile << DC_sb_MC_QCD_systunc[i] << " "; } QCDfile << "\n";
    QCDfile << "syst_unc_all_dn = "; for(int i=0;i<NSB;i++) { if(DC_sb_MC_QCD_systunc[i]>0.999) QCDfile << "0.998 "; else QCDfile << DC_sb_MC_QCD_systunc[i] << " "; } QCDfile << "\n";
    QCDfile.close();
  }
  else std::cout << "Unable to open QCDfile";

  std::ofstream TTZfile (("_TTZ"+sbtag+".txt").c_str());
  if (TTZfile.is_open())
  {
    TTZfile << "# The words after \"#\" are comments. No need to remove them.\n";
    TTZfile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    TTZfile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    TTZfile << "sample = ttz     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    TTZfile << "channel = "; for(int i=0;i<NSB;i++){ TTZfile << "bin" << i+1 << " "; } TTZfile << "\n";

    TTZfile << "\n# Predicted central numbers (need from all backgrounds)\n";
    TTZfile << "# rate = "; for(int i=0;i<NSB;i++){ TTZfile << DC_sb_MC_TTZ[i] << " "; } TTZfile << "\n";
    TTZfile << "rate = "; for(int i=0;i<NSB;i++){ if(DC_sb_MC_TTZ[i]>0) TTZfile << DC_sb_MC_TTZ[i] << " "; else TTZfile << "0 "; } TTZfile << "\n";

    TTZfile << "\n# Control sample events for lost lepton; Raw MC yields for Zinv and ttZ; No need from had. tau and QCD (will be ignored).\n";
    TTZfile << "#cs_event = "; for(int i=0;i<NSB;i++){ TTZfile << DC_sb_MC_TTZ_cs[i] << " "; } TTZfile << "\n";
    TTZfile << "cs_event = "; for(int i=0;i<NSB;i++){ if(DC_sb_MC_TTZ_cs[i]>0) TTZfile << round(DC_sb_MC_TTZ_cs[i]) << " "; else TTZfile << "0 ";} TTZfile << "\n";

    TTZfile << "\n# Average weight for lost lepton, Zinv and ttZ: avg_weight x cs_event = rate. Will be ignored for had. tau and QCD.\n";
    TTZfile << "avg_weight = "; for(int i=0;i<NSB;i++){ TTZfile << DC_sb_MC_TTZ_avgweight[i] << " "; } TTZfile << "\n";

    TTZfile << "\n# All the uncertainties in absolute numbers for both stat and syst uncertainties. Needed from QCD and hadronic tau; No need from lost lepton, Zinv and ttZ (will be ignored). For symmetric uncertainties, put the up and down to be the same.\n";
    TTZfile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ TTZfile << /*DC_sb_MC_TTZ_statunc[i]*/"No Need!!" << " "; break; } TTZfile << "\n";
    TTZfile << "stat_unc_down = "; for(int i=0;i<NSB;i++){ TTZfile << /*DC_sb_MC_TTZ_statunc[i]*/"No Need!!" << " "; break; } TTZfile << "\n";

    TTZfile << "\n# List of all the systematical uncertainties. Do not need combine them. The \"pdf\", \"blackhole\", \"darkmatter\", \"susy\" are keywords to label the different systematic sources (use meaningful names or make comments).\n";
    TTZfile << "syst_unc_pdf_up   = "  ; for(int i=0;i<NSB;i++){ TTZfile << 0 << " "; } TTZfile << "\n";
    TTZfile << "syst_unc_pdf_down = "  ; for(int i=0;i<NSB;i++){ TTZfile << 0 << " "; } TTZfile << "\n";
    TTZfile << "syst_unc_scale_up   = "; for(int i=0;i<NSB;i++){ TTZfile << 0 << " "; } TTZfile << "\n";
    TTZfile << "syst_unc_scale_down = "; for(int i=0;i<NSB;i++){ TTZfile << 0 << " "; } TTZfile << "\n";
    TTZfile << "syst_unc_rate_up   = " ; for(int i=0;i<NSB;i++){ TTZfile << DC_sb_MC_TTZ_systunc[i] << " "; } TTZfile << "\n";
    TTZfile << "syst_unc_rate_down = " ; for(int i=0;i<NSB;i++){ if(DC_sb_MC_TTZ_systunc[i]>0.999) TTZfile << "0.998 "; else TTZfile << DC_sb_MC_TTZ_systunc[i] << " ";} TTZfile << "\n";
    TTZfile.close();
  }
  else std::cout << "Unable to open TTZfile";

  /*
  std::ofstream Rarefile (("_Rare"+sbtag+".txt").c_str());
  if (Rarefile.is_open())
  {
    Rarefile << "rate = "; for(int i=0;i<NSB;i++) { Rarefile << DC_sb_MC_Rare[i] << " "; } Rarefile << "\n";
    Rarefile << "syst_unc_all = "; for(int i=0;i<NSB;i++) { Rarefile << DC_sb_MC_Rare_systunc[i] << " "; } Rarefile << "\n";
    Rarefile.close();
  }
  else std::cout << "Unable to open Rarefile";
  */
  return ;
}

class SignalDataCard
{
 public:
  std::string DC_SignalType;
  int MMass;
  int DMass;
  double DC_all_MC_Signal_avgweight;
  double DC_sb_MC_Signal[NSB] = {0}, DC_sb_MC_Signal_cs[NSB] = {0};
  double DC_sb_MC_Signal_avgweight[NSB] = {0};
  double DC_sb_MC_Signal_statunc[NSB] = {0}, DC_sb_MC_Signal_systunc[NSB] = {0};
  void print_thisSignalDC();
 private:
  void fake_avg_uncs();
};

void SignalDataCard::fake_avg_uncs()
{
  std::cout << "Faking syst uncs in Data card! Mother Mass: " << MMass << "; Daughter Mass: " << DMass << std::endl;
  for(int i=0;i<NSB;i++)
  {
    DC_sb_MC_Signal_cs[i]>0 ? DC_sb_MC_Signal_avgweight[i] = DC_sb_MC_Signal[i]/DC_sb_MC_Signal_cs[i] : DC_sb_MC_Signal_avgweight[i] = DC_all_MC_Signal_avgweight;
    DC_sb_MC_Signal_cs[i]>0 ? DC_sb_MC_Signal_statunc[i] = 1/std::sqrt(DC_sb_MC_Signal_cs[i]) : DC_sb_MC_Signal_statunc[i] = 0;
    if(DC_SignalType=="T1tttt"){ MMass-DMass>400 ? DC_sb_MC_Signal_systunc[i]=0.05 : DC_sb_MC_Signal_systunc[i]=0.15; }
    else if(DC_SignalType=="T2tt"){ MMass-DMass>200 ? DC_sb_MC_Signal_systunc[i]=0.05 : DC_sb_MC_Signal_systunc[i]=0.15; }
    else DC_sb_MC_Signal_systunc[i]=0.1;
  }
}
void SignalDataCard::print_thisSignalDC()
{
  fake_avg_uncs();
  std::ofstream Signalfile (("signal_"+std::to_string(MMass)+"_"+std::to_string(DMass)+".txt").c_str());
  if (Signalfile.is_open())
  {
    Signalfile << "# The words after \"#\" are comments. No need to remove them.\n";
    Signalfile << "luminosity = " << LUMI << "     # in pb-1 (FIXED)\n";
    Signalfile << "channels = " << NSB << "     # total number of channels -> following our search bin definition (FIXED)\n";
    Signalfile << "sample = signal     # name of the background: hadtau, lostle, zinv, qcd, ttz\n";
    Signalfile << "channel = "; for(int i=0;i<NSB;i++){ Signalfile << "bin" << i+1 << " "; } Signalfile << "\n";

    Signalfile << "rate = "; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal[i] << " "; } Signalfile << "\n";
    Signalfile << "cs_event = "; for(int i=0;i<NSB;i++){ Signalfile << round(DC_sb_MC_Signal_cs[i]) << " "; } Signalfile << "\n";
    Signalfile << "avg_weight = "; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal_avgweight[i] << " "; } Signalfile << "\n";
    Signalfile << "contam = "; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    
    Signalfile << "stat_unc_up = "; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal_statunc[i]/*"No Need!!"*/ << " "; } Signalfile << "\n";
    Signalfile << "stat_unc_dn = "; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal_statunc[i]/*"No Need!!"*/ << " "; } Signalfile << "\n";

    Signalfile << "lumi_unc_up = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.027 << " "; } Signalfile << "\n";
    Signalfile << "lumi_unc_dn = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.027 << " "; } Signalfile << "\n";
    Signalfile << "bTagSF_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.1 << " "; } Signalfile << "\n";
    Signalfile << "bTagSF_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.1 << " "; } Signalfile << "\n";
    Signalfile << "mistagSF_up = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.02 << " "; } Signalfile << "\n";
    Signalfile << "mistagSF_dn = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.02 << " "; } Signalfile << "\n";
    Signalfile << "pdfUnc_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "pdfUnc_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "scaleUnc_up = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "scaleUnc_dn = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "isrUnc_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal_systunc[i] << " "; } Signalfile << "\n";
    Signalfile << "isrUnc_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << DC_sb_MC_Signal_systunc[i] << " "; } Signalfile << "\n";
    Signalfile << "metMag_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "metMag_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "metPhi_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "metPhi_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "jetJEC_up = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "jetJEC_dn = "                ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "trigUnc_up = "               ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "trigUnc_dn = "               ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "lepVetoUnc_up = "            ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "lepVetoUnc_dn = "            ; for(int i=0;i<NSB;i++){ Signalfile << 0 << " "; } Signalfile << "\n";
    Signalfile << "genTopSF_up = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "genTopSF_dn = "              ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "mistaggenTopSF_up = "        ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "mistaggenTopSF_dn = "        ; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "data_vs_MC_recoTop_unc_up = "; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile << "data_vs_MC_recoTop_unc_dn = "; for(int i=0;i<NSB;i++){ Signalfile << 0.05 << " "; } Signalfile << "\n";
    Signalfile.close();
  }
  else std::cout << "Unable to open Signalfile";
}

class SSCSHistgram
{
 public:
  void BookHistgram(const char *);

  TFile *oFile;
  //NTop, NBot plots
  TH2D *h_ss_ntopnbot_MC_MuCS, *h_ss_ntopnbot_MC_ElCS;
  //MET MT2 plots after top bot
  TH2D *h_ss_metmt2_MC_MuCS[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_metmt2_MC_ElCS[NTOPJETS_BINS][NBOTJETS_BINS];
};

void SSCSHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  h_ss_ntopnbot_MC_MuCS = new TH2D("h_ss_ntopnbot_MC_MuCS","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);
  h_ss_ntopnbot_MC_ElCS = new TH2D("h_ss_ntopnbot_MC_ElCS","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);

  for(int i=0;i<NTOPJETS_BINS;i++)
  {
    for(int j=0;j<NBOTJETS_BINS;j++)
    { 
      std::string ntnbtag = "NT"+std::to_string(i+1)+"NB"+std::to_string(j+1);
      /*
      if(i==NTOPJETS_BINS-1)
      {
        h_ss_metmt2_MC_MuCS[i][j] = new TH2D(("h_ss_metmt2_MC_MuCS"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_ElCS[i][j] = new TH2D(("h_ss_metmt2_MC_ElCS"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      }
      else if(i!=NTOPJETS_BINS-1 && j==NBOTJETS_BINS-1)
      {
        h_ss_metmt2_MC_MuCS[i][j] = new TH2D(("h_ss_metmt2_MC_MuCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_ElCS[i][j] = new TH2D(("h_ss_metmt2_MC_ElCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      }
      else
      {
        h_ss_metmt2_MC_MuCS[i][j] = new TH2D(("h_ss_metmt2_MC_MuCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
        h_ss_metmt2_MC_ElCS[i][j] = new TH2D(("h_ss_metmt2_MC_ElCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
      }
      */
      int ij=i*NBOTJETS_BINS+j;
      int metsize = mySBGeometry.NMETBINS[ij], mt2size = mySBGeometry.NMT2BINS[ij];
      double metbins_edge[mySBGeometry.metbins_edge.at(ij).size()], mt2bins_edge[mySBGeometry.mt2bins_edge.at(ij).size()];
      std::copy ( mySBGeometry.metbins_edge.at(ij).begin(), mySBGeometry.metbins_edge.at(ij).end(), metbins_edge );
      std::copy ( mySBGeometry.mt2bins_edge.at(ij).begin(), mySBGeometry.mt2bins_edge.at(ij).end(), mt2bins_edge );
      
      //if(i==NTOPJETS_BINS-1 || j==NBOTJETS_BINS-1)
      //{
      //  h_ss_metmt2_MC_MuCS[i][j] = new TH2D(("h_ss_metmt2_MC_MuCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //  h_ss_metmt2_MC_ElCS[i][j] = new TH2D(("h_ss_metmt2_MC_ElCS"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //}
      //else
      //{
        h_ss_metmt2_MC_MuCS[i][j] = new TH2D(("h_ss_metmt2_MC_MuCS"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
        h_ss_metmt2_MC_ElCS[i][j] = new TH2D(("h_ss_metmt2_MC_ElCS"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
      //}
    }
  }
  return ;
}

class SSHistgram
{
 public:
  void BookHistgram(const char *);

  TFile *oFile;
  //NTop, NBot plots
  TH2D *h_ss_ntopnbot_MC_AllBG;
  TH2D *h_ss_ntopnbot_MC_T1tttt_mGluino1200_mLSP800, *h_ss_ntopnbot_MC_T1tttt_mGluino1500_mLSP100;
  TH2D *h_ss_ntopnbot_MC_T2tt_mStop500_mLSP325, *h_ss_ntopnbot_MC_T2tt_mStop850_mLSP100; 
  //MET MT2 plots after top bot
  TH2D *h_ss_metmt2_MC_AllBG[NTOPJETS_BINS][NBOTJETS_BINS];
  TH2D *h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];
  TH2D *h_ss_metmt2_MC_T2tt_mStop500_mLSP325[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_metmt2_MC_T2tt_mStop850_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];
};

void SSHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  h_ss_ntopnbot_MC_AllBG = new TH2D("h_ss_ntopnbot_MC_AllBG","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);
  h_ss_ntopnbot_MC_T1tttt_mGluino1200_mLSP800 = new TH2D("h_ss_ntopnbot_MC_T1tttt_mGluino1200_mLSP800","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);
  h_ss_ntopnbot_MC_T1tttt_mGluino1500_mLSP100 = new TH2D("h_ss_ntopnbot_MC_T1tttt_mGluino1500_mLSP100","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);
  h_ss_ntopnbot_MC_T2tt_mStop500_mLSP325 = new TH2D("h_ss_ntopnbot_MC_T2tt_mStop500_mLSP325","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);
  h_ss_ntopnbot_MC_T2tt_mStop850_mLSP100 = new TH2D("h_ss_ntopnbot_MC_T2tt_mStop850_mLSP100","",NTOPJETS_BINS,mySBGeometry.ntopbins_edge,NBOTJETS_BINS,mySBGeometry.nbotbins_edge);

  for(int i=0;i<NTOPJETS_BINS;i++)
  {
    for(int j=0;j<NBOTJETS_BINS;j++)
    {  
      std::string ntnbtag = "NT"+std::to_string(i+1)+"NB"+std::to_string(j+1);
      /*
      if(i==NTOPJETS_BINS-1)
      {
        h_ss_metmt2_MC_AllBG[i][j] = new TH2D(("h_ss_metmt2_MC_AllBG"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",1,metbins_edge[0],metbins_edge[MET_BINS],1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      }
      else if(i!=NTOPJETS_BINS-1 && j==NBOTJETS_BINS-1)
      {
        h_ss_metmt2_MC_AllBG[i][j] = new TH2D(("h_ss_metmt2_MC_AllBG"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
        h_ss_metmt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      }
      else
      {
        h_ss_metmt2_MC_AllBG[i][j] = new TH2D(("h_ss_metmt2_MC_AllBG"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
        h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
        h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
        h_ss_metmt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
        h_ss_metmt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,MT2_BINS,mt2bins_edge);
      }
      */
      int ij=i*NBOTJETS_BINS+j;
      int metsize = mySBGeometry.NMETBINS[ij], mt2size = mySBGeometry.NMT2BINS[ij];
      double metbins_edge[mySBGeometry.metbins_edge.at(ij).size()], mt2bins_edge[mySBGeometry.mt2bins_edge.at(ij).size()];
      std::copy ( mySBGeometry.metbins_edge.at(ij).begin(), mySBGeometry.metbins_edge.at(ij).end(), metbins_edge );
      std::copy ( mySBGeometry.mt2bins_edge.at(ij).begin(), mySBGeometry.mt2bins_edge.at(ij).end(), mt2bins_edge );

      //if(i==NTOPJETS_BINS-1 /*|| j==NBOTJETS_BINS-1*/)
      //{
      //  h_ss_metmt2_MC_AllBG[i][j] = new TH2D(("h_ss_metmt2_MC_AllBG"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //  h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //  h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //  h_ss_metmt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //  h_ss_metmt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",MET_BINS,metbins_edge,1,mt2bins_edge[0],mt2bins_edge[MT2_BINS]);
      //}
      //else
      //{
        h_ss_metmt2_MC_AllBG[i][j] = new TH2D(("h_ss_metmt2_MC_AllBG"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
        h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
        h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
        h_ss_metmt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
        h_ss_metmt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH2D(("h_ss_metmt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",metsize,metbins_edge,mt2size,mt2bins_edge);
      //} 
    }
  }
  return ;
}

#define SSAUXBGBin 5

class FSLHistgram
{
 public:
  void BookHistgram(const char *);
  
  TFile *oFile;
  TH2D *h2_B;
  TH2D *h2_S;
  TH2D *h2_SOverB;
  TH2D *h2_SOverB2;
  TH2D *h2_SOverB3;
  TH2D *h2_SOverB4;
  TH2D *h2_Q;
  TH2D *h2_nmuCS;

  TH2D *h2_BMETHT;
  TH2D *h2_SMETHT;
  TH2D *h2_SOverBMETHT;
  TH2D *h2_SOverBMETHT2;
  TH2D *h2_SOverBMETHT3;
  TH2D *h2_QMETHT;
  TH2D *h2_nmuCSMETHT;

  TH2D *h2_BMETHTtop;
  TH2D *h2_SMETHTtop;
  TH2D *h2_SOverBMETHTtop;
  TH2D *h2_QMETHTtop;
  TH2D *h2_nmuCSMETHTtop;

  TH2D *h2_BMT2HT;
  TH2D *h2_SMT2HT;
  TH2D *h2_SOverBMT2HT;
  TH2D *h2_nmuCSMT2HT;

};
void FSLHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");

  h2_B = new TH2D("h2_B","h2_B",17,250.0,1100.0,16,200.0,1000.0);
  h2_S = new TH2D("h2_S","h2_S",17,250.0,1100.0,16,200.0,1000.0);
  h2_SOverB = new TH2D("h2_SOverB","h2_SoverB",17,250.0,1100.0,16,200.0,1000.0);
  h2_SOverB2 = new TH2D("h2_SOverB2","h2_SoverB2",17,250.0,1100.0,16,200.0,1000.0);
  h2_SOverB3 = new TH2D("h2_SOverB3","h2_SoverB3",17,250.0,1100.0,16,200.0,1000.0);
  h2_SOverB4 = new TH2D("h2_SOverB4","h2_SoverB4",17,250.0,1100.0,16,200.0,1000.0);
  h2_Q = new TH2D("h2_Q","h2_Q",17,250.0,1100.0,16,200.0,1000.0);
  h2_nmuCS = new TH2D("h2_nmuCS","h2_nmuCS",17,250.0,1100.0,16,200.0,1000.0);

  h2_BMETHT = new TH2D("h2_BMETHT","h2_BMETHT",17,250.0,1100.0,34,300.0,2000.0);
  h2_SMETHT = new TH2D("h2_SMETHT","h2_SMETHT",17,250.0,1100.0,34,300.0,2000.0);
  h2_SOverBMETHT = new TH2D("h2_SOverBMETHT","h2_SoverBMETHT",17,250.0,1100.0,34,300.0,2000.0);
  h2_SOverBMETHT2 = new TH2D("h2_SOverBMETHT2","h2_SoverBMETHT2",17,250.0,1100.0,34,300.0,2000.0);
  h2_SOverBMETHT3 = new TH2D("h2_SOverBMETHT3","h2_SoverBMETHT3",17,250.0,1100.0,34,300.0,2000.0);
  h2_QMETHT = new TH2D("h2_QMETHT","h2_QMETHT",17,250.0,1100.0,34,300.0,2000.0);
  h2_nmuCSMETHT = new TH2D("h2_nmuCSMETHT","h2_nmuCSMETHT",17,250.0,1100.0,34,300.0,2000.0);

  h2_BMETHTtop = new TH2D("h2_BMETHTtop","h2_BMETHTtop",17,250.0,1100.0,34,0.0,1700.0);
  h2_SMETHTtop = new TH2D("h2_SMETHTtop","h2_SMETHTtop",17,250.0,1100.0,34,0.0,1700.0);
  h2_SOverBMETHTtop = new TH2D("h2_SOverBMETHTtop","h2_SoverBMETHTtop",17,250.0,1100.0,34,0.0,1700.0);
  h2_QMETHTtop = new TH2D("h2_QMETHTtop","h2_QMETHTtop",17,250.0,1100.0,34,0.0,1700.0);
  h2_nmuCSMETHTtop = new TH2D("h2_nmuCSMETHTtop","h2_nmuCSMETHTtop",17,250.0,1100.0,34,0.0,1700.0);

  h2_BMT2HT = new TH2D("h2_BMT2HT","h2_BMT2HT",12,200.0,800.0,34,300.0,2000.0);
  h2_SMT2HT = new TH2D("h2_SMT2HT","h2_SMT2HT",12,200.0,800.0,34,300.0,2000.0);
  h2_SOverBMT2HT = new TH2D("h2_SOverBMT2HT","h2_SoverBMT2HT",12,200.0,800.0,34,300.0,2000.0);
  //h2_QMT2HT = new TH2D("h2_QMT2HT","h2_QMT2HT",17,250.0,1100.0,34,300.0,2000.0);
  h2_nmuCSMT2HT = new TH2D("h2_nmuCSMT2HT","h2_nmuCSMT2HT",12,200.0,800.0,34,300.0,2000.0);


}

class SSAUX1DHistgram
{
 public:
  void BookHistgram(const char *);
  
  TFile *oFile;
  //MET MT2 plots after top bot
  TH1D *h_ss_aux_met_MC_AllBG[NTOPJETS_BINS][NBOTJETS_BINS][SSAUXBGBin];
  TH1D *h_ss_aux_met_MC_T1tttt_mGluino1200_mLSP800[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_aux_met_MC_T1tttt_mGluino1500_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];
  TH1D *h_ss_aux_met_MC_T2tt_mStop500_mLSP325[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_aux_met_MC_T2tt_mStop850_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];

  TH1D *h_ss_aux_mt2_MC_AllBG[NTOPJETS_BINS][NBOTJETS_BINS][SSAUXBGBin];
  TH1D *h_ss_aux_mt2_MC_T1tttt_mGluino1200_mLSP800[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_aux_mt2_MC_T1tttt_mGluino1500_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];
  TH1D *h_ss_aux_mt2_MC_T2tt_mStop500_mLSP325[NTOPJETS_BINS][NBOTJETS_BINS], *h_ss_aux_mt2_MC_T2tt_mStop850_mLSP100[NTOPJETS_BINS][NBOTJETS_BINS];

};

void SSAUX1DHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");

  for(int i=0;i<NTOPJETS_BINS;i++)
  {
    for(int j=0;j<NBOTJETS_BINS;j++)
    {
      std::string ntnbtag = "NT"+std::to_string(i+1)+"NB"+std::to_string(j+1);
      
      int ij=i*NBOTJETS_BINS+j;
      int metsize = mySBGeometry.NMETBINS[ij], mt2size = mySBGeometry.NMT2BINS[ij];
      double metbins_edge[mySBGeometry.metbins_edge.at(ij).size()], mt2bins_edge[mySBGeometry.mt2bins_edge.at(ij).size()];
      std::copy ( mySBGeometry.metbins_edge.at(ij).begin(), mySBGeometry.metbins_edge.at(ij).end(), metbins_edge );
      std::copy ( mySBGeometry.mt2bins_edge.at(ij).begin(), mySBGeometry.mt2bins_edge.at(ij).end(), mt2bins_edge );
     
      for(int k=0;k<SSAUXBGBin;k++)
      { 
        std::string smalltag;

        if (k == 0) smalltag = "LL";
        else if (k == 1) smalltag = "HadTau";
        else if (k == 2) smalltag = "Zinv";
        else if (k == 3) smalltag = "QCD";
        else smalltag = "TTZ";
        
				h_ss_aux_met_MC_AllBG[i][j][k] = new TH1D(("h_ss_aux_met_MC_AllBG"+ntnbtag+"_"+smalltag).c_str(),"",(metbins_edge[metsize]+500-metbins_edge[0])/50,metbins_edge[0],metbins_edge[metsize]+500);
        h_ss_aux_mt2_MC_AllBG[i][j][k] = new TH1D(("h_ss_aux_mt2_MC_AllBG"+ntnbtag+"_"+smalltag).c_str(),"",(mt2bins_edge[mt2size]+500-mt2bins_edge[0])/50,mt2bins_edge[0],mt2bins_edge[mt2size]+500);

        h_ss_aux_met_MC_AllBG[i][j][k]->SetFillColor(k+2);
        h_ss_aux_mt2_MC_AllBG[i][j][k]->SetFillColor(k+2);
        h_ss_aux_met_MC_AllBG[i][j][k]->SetLineColor(k+2);
        h_ss_aux_mt2_MC_AllBG[i][j][k]->SetLineColor(k+2);
      }

      h_ss_aux_met_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH1D(("h_ss_aux_met_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",(metbins_edge[metsize]+500-metbins_edge[0])/50,metbins_edge[0],metbins_edge[metsize]+500);
      h_ss_aux_met_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH1D(("h_ss_aux_met_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",(metbins_edge[metsize]+500-metbins_edge[0])/50,metbins_edge[0],metbins_edge[metsize]+500);
      h_ss_aux_met_MC_T2tt_mStop500_mLSP325[i][j] = new TH1D(("h_ss_aux_met_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",(metbins_edge[metsize]+500-metbins_edge[0])/50,metbins_edge[0],metbins_edge[metsize]+500);
      h_ss_aux_met_MC_T2tt_mStop850_mLSP100[i][j] = new TH1D(("h_ss_aux_met_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",(metbins_edge[metsize]+500-metbins_edge[0])/50,metbins_edge[0],metbins_edge[metsize]+500);

      h_ss_aux_mt2_MC_T1tttt_mGluino1200_mLSP800[i][j] = new TH1D(("h_ss_aux_mt2_MC_T1tttt_mGluino1200_mLSP800"+ntnbtag).c_str(),"",(mt2bins_edge[mt2size]+500-mt2bins_edge[0])/50,mt2bins_edge[0],mt2bins_edge[mt2size]+500);
      h_ss_aux_mt2_MC_T1tttt_mGluino1500_mLSP100[i][j] = new TH1D(("h_ss_aux_mt2_MC_T1tttt_mGluino1500_mLSP100"+ntnbtag).c_str(),"",(mt2bins_edge[mt2size]+500-mt2bins_edge[0])/50,mt2bins_edge[0],mt2bins_edge[mt2size]+500);
      h_ss_aux_mt2_MC_T2tt_mStop500_mLSP325[i][j] = new TH1D(("h_ss_aux_mt2_MC_T2tt_mStop500_mLSP325"+ntnbtag).c_str(),"",(mt2bins_edge[mt2size]+500-mt2bins_edge[0])/50,mt2bins_edge[0],mt2bins_edge[mt2size]+500);
      h_ss_aux_mt2_MC_T2tt_mStop850_mLSP100[i][j] = new TH1D(("h_ss_aux_mt2_MC_T2tt_mStop850_mLSP100"+ntnbtag).c_str(),"",(mt2bins_edge[mt2size]+500-mt2bins_edge[0])/50,mt2bins_edge[0],mt2bins_edge[mt2size]+500);
    }
  }
}

//##########functions to calculate Delta_R and Delta Phi###############
double DeltaPhi(double phi1, double phi2) 
{
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) 
{
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}
