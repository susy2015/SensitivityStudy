#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <set>

#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include "SSAUX1DPlots.h"

int main(int argc, char* argv[])
{

  if (argc < 1)
  {
    std::cerr <<"Please give at least 1 argument " << "TargetDirName" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./SSAUX1DPlots SSAUX1DPlots20160517" << std::endl;
    return -1;
  }
  std::string DirName = argv[1];

  std::vector<Plotting_Parameter> myPlotting_Paramete = 
  { 
    //T1tttt_1200800
    { "T1tttt_mGluino1200_mLSP800","NT1NB1","met" },
    { "T1tttt_mGluino1200_mLSP800","NT1NB2","met" },
    { "T1tttt_mGluino1200_mLSP800","NT1NB3","met" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB1","met" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB2","met" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB3","met" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB1","met" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB2","met" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB3","met" },
    { "T1tttt_mGluino1200_mLSP800","NT1NB1","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT1NB2","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT1NB3","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB1","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB2","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT2NB3","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB1","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB2","mt2" },
    { "T1tttt_mGluino1200_mLSP800","NT3NB3","mt2" },
    //T1tttt_1500100
    { "T1tttt_mGluino1500_mLSP100","NT1NB1","met" },
    { "T1tttt_mGluino1500_mLSP100","NT1NB2","met" },
    { "T1tttt_mGluino1500_mLSP100","NT1NB3","met" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB1","met" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB2","met" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB3","met" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB1","met" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB2","met" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB3","met" },
    { "T1tttt_mGluino1500_mLSP100","NT1NB1","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT1NB2","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT1NB3","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB1","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB2","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT2NB3","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB1","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB2","mt2" },
    { "T1tttt_mGluino1500_mLSP100","NT3NB3","mt2" },
    //T2tt_500325
    { "T2tt_mStop500_mLSP325","NT1NB1","met" },
    { "T2tt_mStop500_mLSP325","NT1NB2","met" },
    { "T2tt_mStop500_mLSP325","NT1NB3","met" },
    { "T2tt_mStop500_mLSP325","NT2NB1","met" },
    { "T2tt_mStop500_mLSP325","NT2NB2","met" },
    { "T2tt_mStop500_mLSP325","NT2NB3","met" },
    { "T2tt_mStop500_mLSP325","NT3NB1","met" },
    { "T2tt_mStop500_mLSP325","NT3NB2","met" },
    { "T2tt_mStop500_mLSP325","NT3NB3","met" },
    { "T2tt_mStop500_mLSP325","NT1NB1","mt2" },
    { "T2tt_mStop500_mLSP325","NT1NB2","mt2" },
    { "T2tt_mStop500_mLSP325","NT1NB3","mt2" },
    { "T2tt_mStop500_mLSP325","NT2NB1","mt2" },
    { "T2tt_mStop500_mLSP325","NT2NB2","mt2" },
    { "T2tt_mStop500_mLSP325","NT2NB3","mt2" },
    { "T2tt_mStop500_mLSP325","NT3NB1","mt2" },
    { "T2tt_mStop500_mLSP325","NT3NB2","mt2" },
    { "T2tt_mStop500_mLSP325","NT3NB3","mt2" },
    //T2tt_850100
    { "T2tt_mStop850_mLSP100","NT1NB1","met" },
    { "T2tt_mStop850_mLSP100","NT1NB2","met" },
    { "T2tt_mStop850_mLSP100","NT1NB3","met" },
    { "T2tt_mStop850_mLSP100","NT2NB1","met" },
    { "T2tt_mStop850_mLSP100","NT2NB2","met" },
    { "T2tt_mStop850_mLSP100","NT2NB3","met" },
    { "T2tt_mStop850_mLSP100","NT3NB1","met" },
    { "T2tt_mStop850_mLSP100","NT3NB2","met" },
    { "T2tt_mStop850_mLSP100","NT3NB3","met" },
    { "T2tt_mStop850_mLSP100","NT1NB1","mt2" },
    { "T2tt_mStop850_mLSP100","NT1NB2","mt2" },
    { "T2tt_mStop850_mLSP100","NT1NB3","mt2" },
    { "T2tt_mStop850_mLSP100","NT2NB1","mt2" },
    { "T2tt_mStop850_mLSP100","NT2NB2","mt2" },
    { "T2tt_mStop850_mLSP100","NT2NB3","mt2" },
    { "T2tt_mStop850_mLSP100","NT3NB1","mt2" },
    { "T2tt_mStop850_mLSP100","NT3NB2","mt2" },
    { "T2tt_mStop850_mLSP100","NT3NB3","mt2" }
  };


  SSAUX1DPlots mySSAUX1DPlots;
  mySSAUX1DPlots.Initialization(DirName);
  //mySSAUX1DPlots.SSAUX1DPlotsLoop("","NBNT","met");
  std::vector<Plotting_Parameter>::iterator iter_p;

  for( iter_p = myPlotting_Paramete.begin() ; iter_p != myPlotting_Paramete.end() ; iter_p ++)
  {
    mySSAUX1DPlots.SSAUX1DPlotsLoop(
                                    (*iter_p).hist_tag,
                                    (*iter_p).ntnb_tag,
                                    (*iter_p).var_tag
                                   );
  }
  myPlotting_Paramete.clear();
  return 0;
}

