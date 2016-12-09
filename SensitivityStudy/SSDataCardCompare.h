#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"

#include "SusyAnaTools/Tools/searchBins.h"

#include "ConstantsSnippet.h"

SearchBins mySearchBins("SB_45_2015");

class DCComparePlots
{
 public:
  std::string target_DIR;
  void PrintPlotsName();
  void DCComparePlotsLoop(
                          std::string sample_type,
                          std::string var_type
                         );
 private:
  void GetArrayfromDC(
                      std::ifstream &DC,
                      std::string vartype,
                      double (&value)[NSB]
                     );
};

void DCComparePlots::GetArrayfromDC(
                                    std::ifstream &DC,
                                    std::string vartype,
                                    double (&value)[NSB]
                                   )
{
  std::string line;
  while (std::getline(DC, line))
  {
    if( (line.find(vartype) != std::string::npos) && line.at(0)!='#' )
    {
      //std::cout << line << std::endl;
      std::string numstr = line.substr( line.find("=")+1 );
      //std::cout << numstr << std::endl;
      std::stringstream ssin(numstr); int i=0;
      while (ssin.good() && i < NSB)
      {
        std::string arr; ssin >> arr; 
        //std::cout << arr << " ";
        value[i] = std::stod (arr);
        ++i; 
      }
      //std::cout << std::endl;
      break;
    }
  }
  //for(int i=0;i<NSB;i++){ std::cout << value[i] << " "; } std::cout << std::endl;
  return ;
}

void DCComparePlots::DCComparePlotsLoop(
                                        std::string sample_type,
                                        std::string var_type
                                       )
{
  std::string realDCname, fakeDCname;
  if(sample_type=="LL"){ realDCname = "DataCard_real_45Ref/lostle.txt"; fakeDCname = "DataCard_fake/_LL_45BinsLUMI2015Moriond.txt"; }
  else if(sample_type=="HadTau"){ realDCname = "DataCard_real_45Ref/hadtau.txt"; fakeDCname = "DataCard_fake/_HadTau_45BinsLUMI2015Moriond.txt"; }
  else if(sample_type=="Zinv"){ realDCname = "DataCard_real_45Ref/zinv.txt"; fakeDCname = "DataCard_fake/_Zinv_45BinsLUMI2015Moriond.txt"; }
  else if(sample_type=="TTZ"){ realDCname = "DataCard_real_45Ref/ttz.txt"; fakeDCname = "DataCard_fake/_TTZ_45BinsLUMI2015Moriond.txt"; }
  else { std::cout << "Bad sample type!" << std::endl; return; }

  std::ifstream realDC(realDCname.c_str());
  std::ifstream fakeDC(fakeDCname.c_str());
  double realvalue[NSB] = {0}, fakevalue[NSB] = {0};
  GetArrayfromDC(realDC,var_type,realvalue);
  GetArrayfromDC(fakeDC,var_type,fakevalue);
  for(int i=0;i<NSB;i++){ std::cout << realvalue[i] << " "; } std::cout << std::endl;
  for(int i=0;i<NSB;i++){ std::cout << fakevalue[i] << " "; } std::cout << std::endl;

  TH1D * h_fake = new TH1D(("h_ss_DCP_fake_"+sample_type+"_"+var_type).c_str(),"",NSB+1,0,NSB+1);
  TH1D * h_real = new TH1D(("h_ss_DCP_real_"+sample_type+"_"+var_type).c_str(),"",NSB+1,0,NSB+1);

  for(int i=0;i<NSB;i++)
  {
    h_fake->SetBinContent(i,fakevalue[i]);
    h_real->SetBinContent(i,realvalue[i]);
  }

  //Set Style for h_real and h_fake
  h_fake->SetMarkerStyle(21);
  h_fake->SetMarkerColor(kBlue);
  h_fake->SetMarkerSize(0.9);
  h_fake->SetLineColor(h_fake->GetMarkerColor());
  h_fake->SetLineWidth(3);
  h_fake->GetXaxis()->SetRangeUser(0,NSB);
  h_fake->GetXaxis()->SetTitle("Search Bin");
  h_fake->GetYaxis()->SetTitleOffset(0.6);
  h_fake->GetYaxis()->SetTitleFont(42);
  h_fake->GetYaxis()->SetTitleSize(0.065);
  h_fake->GetYaxis()->SetLabelSize(0.04);
  h_fake->GetYaxis()->SetLabelFont(42);
  h_fake->GetYaxis()->SetTitle(var_type.c_str());

  h_real->SetMarkerStyle(20);
  h_real->SetMarkerColor(kRed);
  h_real->SetMarkerSize(0.9);
  h_real->SetLineColor(h_real->GetMarkerColor());
  h_real->SetLineWidth(3);
  h_real->GetXaxis()->SetRangeUser(0,NSB);
  h_real->GetXaxis()->SetTitle("Search Bin");
  h_real->GetYaxis()->SetTitleOffset(0.6);
  h_real->GetYaxis()->SetTitleFont(42);
  h_real->GetYaxis()->SetTitleSize(0.065);
  h_real->GetYaxis()->SetLabelSize(0.04);
  h_real->GetYaxis()->SetLabelFont(42);
  h_real->GetYaxis()->SetTitle(var_type.c_str());

  // Ratio plots
  TH1* h_ratio;
  h_ratio = static_cast<TH1*>(h_real->Clone("Ratio"));
  h_ratio->Divide(h_fake);
	h_ratio->SetMarkerSize(1);
  h_ratio->GetYaxis()->SetTitle("#frac{Real}{Fake}");
  //h_ratio->GetYaxis()->SetRangeUser(0.0,5.1);
  h_ratio->SetTitle("");
  h_ratio->SetStats(0);
  h_ratio->SetLineWidth(1);
  h_ratio->GetYaxis()->SetTitleSize(0.15);
  h_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_ratio->GetYaxis()->SetTitleFont(42);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelFont(42);
  h_ratio->GetXaxis()->SetLabelOffset(0.01);
  h_ratio->GetXaxis()->SetLabelFont(42);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->GetXaxis()->SetTitleSize(0.16);
  h_ratio->GetXaxis()->SetTitleFont(42);
  h_ratio->GetXaxis()->SetTitleOffset(0.6);

  const std::string titre = "Data Card Comparison: " + sample_type + " " + var_type;
  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);

   //Create Legend
	TLegend* leg = new TLegend(0.55,0.75,0.90,0.90);
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  //leg->SetFillStyle();
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("Data Card Comparison");
  leg->AddEntry(h_real,"2015 Moriond Data Card","P");
  leg->AddEntry(h_fake,"Fake Data Card","P");


	//Draw plots on Canvas
  TCanvas *c = new TCanvas("c","",50,50,800,600); 
	gStyle->SetOptStat(0);

  TPad *pad = (TPad*) c->GetPad(0); 
  pad->Clear();
  pad->Divide(2, 1);

  double divRatio = 0.20;
  double labelRatio = (1-divRatio)/divRatio;
  double small = 0;

  pad->cd(1); 
  TPad *pad1 = (TPad*) pad->GetPad(1);
  pad1->SetPad("", "", 0, divRatio, 1.0, 1.0, kWhite);
  pad1->SetBottomMargin(0.005);
  pad1->SetBorderMode(0);
  
	h_real->Draw("P");
	h_fake->SetFillColor(kBlue-4);
	h_fake->SetFillStyle(3001);
	h_fake->Draw("P same");
  double maxdrawsb = 1.0; maxdrawsb = *std::max_element(realvalue,realvalue+NSB);
  mySearchBins.drawSBregionDef(0.0, maxdrawsb);
  title->Draw("same");
  leg->Draw("same");
  c->Update(); 
  
	pad->cd(2);
  TPad *pad2 = (TPad*) pad->GetPad(2);
  pad2->SetPad("ratio", "", 0, 0, 1.0, divRatio, kWhite);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(small);
  pad2->SetBorderMode(0);

  h_ratio->SetMaximum(1.8);
  h_ratio->SetMinimum(0.4);

  TLine *tl_one = new TLine();
  tl_one->SetLineStyle(2);
  tl_one->SetLineColor(1);
  tl_one->SetLineWidth(2);
  
  h_ratio->Draw("P");
  tl_one->DrawLine(0,1.,NSB,1.);

  c->SaveAs( TString("DCPValidation/_") + TString(sample_type) + TString("_") + TString(var_type) + TString(".png") );
  c->SaveAs( TString("DCPValidation/_") + TString(sample_type) + TString("_") + TString(var_type) + TString(".pdf") );
  c->SaveAs( TString("DCPValidation/_") + TString(sample_type) + TString("_") + TString(var_type) + TString(".C") );

  return ;
}

struct DC_Compare_Parameter
{
  std::string sample_type;//LL, HadTau, Zinv and TTZ
  std::string var_type;//variable name in the data we are interested
};
