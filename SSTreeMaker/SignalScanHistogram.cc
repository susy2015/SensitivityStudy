#include "SignalScanHistogram.h"

void SignalScanHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
  //closure plots on different variables
  h_totEvt_xSusyMotherMass_ySusyLSPMass = new TH2D("h_totEvt_xSusyMotherMass_ySusyLSPMass","",(SusyMotherMass_max-SusyMotherMass_min)/25,SusyMotherMass_min,SusyMotherMass_max,(SusyLSPMass_max-SusyLSPMass_min)/25,SusyMotherMass_min,SusyLSPMass_max);

  return ;
}
