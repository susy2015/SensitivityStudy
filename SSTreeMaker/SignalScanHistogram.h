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
  double SusyMotherMass_max = 3000;
  double SusyMotherMass_min = 0;
  double SusyLSPMass_max = 2000;
  double SusyLSPMass_min = 0;  
};

