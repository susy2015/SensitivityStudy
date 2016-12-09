#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string>

#include "SusyAnaTools/Tools/searchBins.h"

#include "ConstantsSnippet.h"

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
