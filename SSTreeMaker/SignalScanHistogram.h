#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

class SignalScanHistgram
{
 public:
  void BookHistgram(const char *);
  TFile *oFile;
  //closure plots on different variables and search bins
  TH2D *h_totEvt_xSusyMotherMass_ySusyLSPMass;
 private:
  double SusyMotherMass_max = 2500;
  double SusyMotherMass_min = 100;
  double SusyLSPMass_max = 1800;
  double SusyLSPMass_min = 0;  
};

