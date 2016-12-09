#include "SignalDataCard.h"

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
