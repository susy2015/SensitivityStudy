#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdlib.h>

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
#include "THStack.h"

//#include "SusyAnaTools/Tools/searchBins.h"

#include "SSReWeighting.h"

class SSAUX1DPlots
{
 public:
  std::string target_DIR;

  TFile * fin;
  TList * list;

  std::string lumi_str;  

  double scale = 1;

  void Initialization(std::string dir); 
  void PrintPlotsName();
  void SSAUX1DPlotsLoop(
                        TString hist_tag,
                        TString ntnb_tag,
                        TString var_tag
                       );
};

void SSAUX1DPlots::Initialization(std::string dir)
{
  target_DIR = dir;
  system( ("mkdir " + dir).c_str() );

  fin = TFile::Open("SSAUX1DAllMC.root");
  list = fin->GetListOfKeys();
  //convert lumi from double pb-1 to string, fb-1
  std::ostringstream strs;
  strs << (LUMI/1000);
  lumi_str = strs.str();
}

void SSAUX1DPlots::PrintPlotsName()
{
  for(int i  = 0 ; i < list->GetSize() ; i++)
  {
    std::cout<<"Name: "<< list->At(i)->GetName() << "("<< i <<")"<<std::endl;
  }
  return ;
}

void SSAUX1DPlots::SSAUX1DPlotsLoop( 
                                     TString hist_tag,
                                     TString ntnb_tag,
                                     TString var_tag
                                   )
{ 
  TH1D * h;
  THStack * hs_BG = new THStack("hs_BG","");

  TLegend* leg = new TLegend(0.6,0.6,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);

  int NHist = list->GetSize();
  std::string nametag;

	//Get plots from root file
  for(int i  = 0 ; i < NHist ; i++)
  {
    nametag = hist_tag;
    if( TString(list->At(i)->GetName()).Contains( ntnb_tag ) && TString(list->At(i)->GetName()).Contains( var_tag ) )
    {
      //choose one signal point
      if( TString(list->At(i)->GetName()).Contains( hist_tag ) )
      { 
        h = (TH1D*)fin->Get(list->At(i)->GetName())->Clone();
        leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), hist_tag, "l");
      }
      if( TString(list->At(i)->GetName()).Contains( "AllBG" ) )
      {
        hs_BG->Add( (TH1D*)fin->Get(list->At(i)->GetName()) );
        std::string smalltag;
        if( TString(list->At(i)->GetName()).Contains( "LL" ) ) { smalltag = "LL"; leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), smalltag.c_str(), "f"); }
        if( TString(list->At(i)->GetName()).Contains( "HadTau" ) ) { smalltag = "HadTau"; leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), smalltag.c_str(), "f"); }
        if( TString(list->At(i)->GetName()).Contains( "Zinv" ) ) { smalltag = "Zinv"; leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), smalltag.c_str(), "f"); }
        if( TString(list->At(i)->GetName()).Contains( "QCD" ) ) { smalltag = "QCD"; leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), smalltag.c_str(), "f"); }
        if( TString(list->At(i)->GetName()).Contains( "TTZ" ) ) { smalltag = "TTZ"; leg->AddEntry( (TH1D*)fin->Get(list->At(i)->GetName()), smalltag.c_str(), "f"); }
      }
    }
    else continue;
	}

  if( var_tag=="met" )
  { 
    h->GetXaxis()->SetTitle("MET[GeV]");
    //hs_BG->GetXaxis()->SetTitle("MET[GeV]");
  }
  else if( var_tag=="mt2" )
  {
    h->GetXaxis()->SetTitle("MT2[GeV]");
    //hs_BG->GetXaxis()->SetTitle("MT2[GeV]");
  }
  else
  {
    std::cout << "Unknow plot name??!!" << std::endl;
  }

  if     ( var_tag.Contains( "met" ) ) nametag+=" MET";
  else if( var_tag.Contains( "mt2" ) ) nametag+=" MT2";
  else nametag+="";

  if     ( ntnb_tag.Contains( "NT1NB1" ) ) nametag += " : 1 Top, 1 Bot";
  else if( ntnb_tag.Contains( "NT1NB2" ) ) nametag += " : 1 Top, 2 Bot";
  else if( ntnb_tag.Contains( "NT1NB3" ) ) nametag += " : 1 Top, >=3 Bot";
  else if( ntnb_tag.Contains( "NT2NB1" ) ) nametag += " : 2 Top, 1 Bot";
  else if( ntnb_tag.Contains( "NT2NB2" ) ) nametag += " : 2 Top, 2 Bot";
  else if( ntnb_tag.Contains( "NT2NB3" ) ) nametag += " : 2 Top, >=3 Bot";
  else if( ntnb_tag.Contains( "NT3NB1" ) ) nametag += " : >= 3 Top, 1 Bot";
  else if( ntnb_tag.Contains( "NT3NB2" ) ) nametag += " : >= 3 Top, 2 Bot";
  else if( ntnb_tag.Contains( "NT3NB3" ) ) nametag += " : >= 3 Top, >=3 Bot";
  else nametag += "";
  //Create LUMI stamp
  const std::string titre= nametag + " 8 fb-1";

  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);

  h->SetMarkerSize(2);
  //Draw plots on Canvas
  TCanvas *c = new TCanvas("c","",50,50,800,600);
  c->SetLogy();
  //HistStyle::init();
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");
  
  hs_BG->Draw("hist");
  h->Draw("same e0");
  title->Draw("same");
  leg->Draw("same");

  c->SaveAs( target_DIR + TString("/_") + hist_tag + "_" + ntnb_tag + "_" + var_tag + TString(".png") );
  c->SaveAs( target_DIR + TString("/_") + hist_tag + "_" + ntnb_tag + "_" + var_tag + TString(".pdf") );
  c->SaveAs( target_DIR + TString("/_") + hist_tag + "_" + ntnb_tag + "_" + var_tag + TString(".C") );

  return ;
}

struct Plotting_Parameter
{
  TString hist_tag;
  TString ntnb_tag;
  TString var_tag;
};
