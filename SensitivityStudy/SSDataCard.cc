#include "SSDataCard.h"

void SSDataCard::fake_avg_uncs()
{
  std::cout << "Faking syst uncs in Data card!" << std::endl;
  //SBGeometry mySBGeometry;
  SearchBins mySearchBins(sb_tag);
  
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
