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
#include "SusyAnaTools/Tools/searchBins.h"

#include "SSReWeighting.h"
#include "SSDataCard.h"
#include "SignalDataCard.h"

#include "ConstantsSnippet.h"

std::string dir_out = "";

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
