
//void Plot_eff()
{
  TFile f_eff("FSLHistgram.root");

  const std::string titre="CMS Simulation";

  bool do_c1=false;
  bool do_c2=false;
  bool do_c3=false;
  bool do_c4=false;
  bool do_c5=false;
 
  do_c1=true;
  do_c2=true;
  do_c3=true;
  //do_c4=true;
  do_c5=true;
 
  gStyle->SetPaintTextFormat("2.2f");

  if (do_c1)
  {

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();
  //c1->SetLogy();

  TH2F *h2_B_ori_c1 = (TH2F*)f_eff.Get("h2_BMETHT");
  h2_B_ori_c1->Draw("colztext");
  h2_B_ori_c1->SetTitle("");
  h2_B_ori_c1->SetXTitle("MET cut [GeV]");
  h2_B_ori_c1->SetYTitle("HT cut [GeV]");
  h2_B_ori_c1->SetStats(0);

  c1->SaveAs( "c1.C" );
  }

  if (do_c2)
  {

  TCanvas *c2 = new TCanvas("c2", "c2",0,51,1920,1004);
  c2->SetFillColor(0);
  c2->cd();
  //c2->SetLogy();

  TH2F *h2_B_ori_c2 = (TH2F*)f_eff.Get("h2_SMETHT");
  h2_B_ori_c2->Draw("colztext");
  h2_B_ori_c2->SetTitle("");
  h2_B_ori_c2->SetXTitle("MET cut [GeV]");
  h2_B_ori_c2->SetYTitle("HT cut [GeV]");
  h2_B_ori_c2->SetStats(0);

  c2->SaveAs( "c2.C" );
  }

  gStyle->SetPaintTextFormat("2.4f");

  if (do_c3)
  {

  TCanvas *c3 = new TCanvas("c3", "c3",0,51,1920,1004);
  c3->SetFillColor(0);
  c3->cd();
  //c3->SetLogy();

  TH2F *h2_B_ori_c3 = (TH2F*)f_eff.Get("h2_SOverBMETHT");
  h2_B_ori_c3->Draw("colztext");
  h2_B_ori_c3->SetTitle("");
  h2_B_ori_c3->SetXTitle("MET cut [GeV]");
  h2_B_ori_c3->SetYTitle("HT cut [GeV]");
  h2_B_ori_c3->SetStats(0);

  c3->SaveAs( "c3.C" );
  }

  if (do_c4)
  {

  TCanvas *c4 = new TCanvas("c4", "c4",0,51,1920,1004);
  c4->SetFillColor(0);
  c4->cd();
  //c4->SetLogy();

  TH2F *h2_B_ori_c4 = (TH2F*)f_eff.Get("h2_QMETHT");
  h2_B_ori_c4->Draw("colztext");
  h2_B_ori_c4->SetTitle("");
  h2_B_ori_c4->SetXTitle("MET cut [GeV]");
  h2_B_ori_c4->SetYTitle("HT cut [GeV]");
  h2_B_ori_c4->SetStats(0);

  c4->SaveAs( "c4.C" );
  }

  gStyle->SetPaintTextFormat("2.2f");

  if (do_c5)
  {

  TCanvas *c5 = new TCanvas("c5", "c5",0,51,1920,1004);
  c5->SetFillColor(0);
  c5->cd();
  //c5->SetLogy();

  TH2F *h2_B_ori_c5 = (TH2F*)f_eff.Get("h2_nmuCSMETHT");
  h2_B_ori_c5->Draw("colztext");
  h2_B_ori_c5->SetTitle("");
  h2_B_ori_c5->SetXTitle("MET cut [GeV]");
  h2_B_ori_c5->SetYTitle("HT cut [GeV]");
  h2_B_ori_c5->SetStats(0);

  c5->SaveAs( "c5.C" );
  }

  //gPad->Modified();
  //gPad->Update();
  //return;

}
