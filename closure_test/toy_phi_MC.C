#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>

double GetPhiInRange(double phi)
{
  double result = phi;
  if (result < 0)
  {
    result = result + 2 * TMath::Pi();
  }
  if (result > 2 * TMath::Pi())
  {
    result = result - 2 * TMath::Pi();
  }
  return result;
}

double Poly2(double *x, double *par)
{
  return par[0] * (1 + 2.0 * par[1] * TMath::Cos(2.0 * x[0]));
}

double Poly3(double *x, double *par)
{
  return par[0] * (1 + 2.0 * par[1] * TMath::Power(TMath::Cos(2.0 * x[0]), 1.0) + par[2] * TMath::Power(TMath::Sin(3.0 * x[0]), 1.0));
}

double v2(double *x, double *par)
{
  return 1 + 2 * 0.1 * TMath::Power(TMath::Cos(2.0 * x[0]), 1.0);
}

double calcrestmom(TLorentzVector mother, TLorentzVector daug)
{
  TVector3 boost_vector = mother.BoostVector();
  daug.Boost(-boost_vector);
  return daug.P();
}

double calcv2(double phi){
  return cos(2*phi);
}

void toy_phi_MC()
{
  TF1 *f2 = new TF1("f2", "1+2.0*0.5*TMath::Cos(2.0*x)", 0, 2 * TMath::Pi());
  TF1 *f3 = new TF1("f3", "1+2.0*0.5*TMath::Cos(2.0*x)", 0, 2 * TMath::Pi());
  TF1 *f4 = new TF1("f4", "1+2.0*0.5*TMath::Cos(2.0*x)", 0, 2 * TMath::Pi());

  float pt1, eta1, pt2, eta2, pt3, eta3;
  float phi1, phi2, phi3;

  TRandom3 rand;

  TH1F *hPhiProton = new TH1F("hPhiProton", ";#phi (rad);counts", 72, 0.0, 2 * TMath::Pi());
  hPhiProton->SetMarkerStyle(20);
  hPhiProton->SetMarkerColor(kBlack);
  TH1F *hPhiDeuteron = new TH1F("hPhiDeuteron", ";#phi (rad);counts", 72, 0.0, 2 * TMath::Pi());
  hPhiDeuteron->SetMarkerStyle(20);
  hPhiDeuteron->SetMarkerColor(kRed);
  TH1F *hPhiHelium = new TH1F("hPhiHelium", ";#phi (rad);counts", 72, 0.0, 2 * TMath::Pi());
  hPhiHelium->SetMarkerStyle(20);
  hPhiHelium->SetMarkerColor(kBlue);

  TH1F *hV2Proton = new TH1F("hV2Proton", ";cos(2(#phi - #Psi_{2}));counts", 200, -1., 1.);
  hV2Proton->SetMarkerStyle(20);
  hV2Proton->SetMarkerColor(kBlack);
  TH1F *hV2Deuteron = new TH1F("hV2Deuteron", ";cos(2(#phi - #Psi_{2}));counts", 200, -1., 1.);
  hV2Deuteron->SetMarkerStyle(20);
  hV2Deuteron->SetMarkerColor(kRed);
  TH1F *hV2Helium = new TH1F("hV2Helium", ";cos(2(#phi - #Psi_{2}));counts", 200, -1., 1.);
  hV2Helium->SetMarkerStyle(20);
  hV2Helium->SetMarkerColor(kBlue);

  TLorentzVector mother, deuteron, helium3;
  TLorentzVector d1, d2, d3;
  TLorentzVector d1dummy, d2dummy, d3dummy;

  for (int i = 0; i < 5000000; i++)
  {
    pt1 = rand.Uniform(0.0, 10.0);
    eta1 = rand.Uniform(-0.8, 0.8);
    pt2 = rand.Uniform(0.0, 10.0);
    eta2 = rand.Uniform(-0.8, 0.8);
    pt3 = rand.Uniform(0.0, 10.0);
    eta3 = rand.Uniform(-0.8, 0.8);

    phi1 = f2->GetRandom(0.0, 2.0 * TMath::Pi());
    phi2 = f3->GetRandom(0.0, 2.0 * TMath::Pi());
    phi3 = f4->GetRandom(0.0, 2.0 * TMath::Pi());

    d1.SetPtEtaPhiM(pt1, eta1, phi1, 0.938);
    d2.SetPtEtaPhiM(pt2, eta2, phi2, 0.938);
    d3.SetPtEtaPhiM(pt3, eta3, phi3, 0.938);

    d1dummy.SetPtEtaPhiM(pt1, eta1, phi1, 0.938);
    d2dummy.SetPtEtaPhiM(pt2, eta2, phi2, 0.938);
    d3dummy.SetPtEtaPhiM(pt3, eta3, phi3, 0.938);

    mother = d1 + d2;
    deuteron.SetPtEtaPhiM(mother.Pt(), mother.Eta(), mother.Phi(), 1.8756);

    helium3 = d1 + d2 + d3;
    TLorentzVector helium3_vec;
    helium3_vec.SetPtEtaPhiM(helium3.Pt(), helium3.Eta(), helium3.Phi(), 2.80839);

    double proton1restmom = calcrestmom(deuteron, d1dummy);
    double proton2restmom = calcrestmom(deuteron, d2dummy);

    double p1 = calcrestmom(helium3_vec, d1dummy);
    double p2 = calcrestmom(helium3_vec, d2dummy);
    double p3 = calcrestmom(helium3_vec, d3dummy);

    hPhiProton->Fill(GetPhiInRange(d1.Phi()));
    hV2Proton->Fill(calcv2(d1.Phi()));

    if (TMath::Abs(proton1restmom - proton2restmom) < 0.3)
    {
      hPhiDeuteron->Fill(GetPhiInRange(deuteron.Phi()));
      hV2Deuteron->Fill(calcv2(deuteron.Phi()));
    }

    if (TMath::Abs(p1 - p2) < 0.3 && TMath::Abs(p1 - p3) < 0.3 && TMath::Abs(p2 - p3) < 0.3)
    {
      hPhiHelium->Fill(GetPhiInRange(helium3.Phi()));
      hV2Helium->Fill(calcv2(helium3.Phi()));
    }
  }

  // phi distributions

  hPhiProton->Sumw2();
  hPhiDeuteron->Sumw2();
  hPhiHelium->Sumw2();

  hPhiProton->Scale(1.0 / hPhiProton->Integral());
  hPhiDeuteron->Scale(1.0 / hPhiDeuteron->Integral());
  hPhiHelium->Scale(1.0 / hPhiHelium->Integral());

  hPhiProton->SetStats(0);
  hPhiProton->GetXaxis()->SetTitle("#phi (rad)");
  hPhiProton->GetYaxis()->SetTitle("normalized counts");

  TCanvas *cPhi = new TCanvas("cPhi", "phi distributions", 800, 600);
  cPhi->DrawFrame(-0.01, -0.01, 6.29, 0.05, ";#phi (rad);counts");
  hPhiProton->Draw("same");
  hPhiDeuteron->Draw("same");
  hPhiHelium->Draw("same");

  TLegend *legPhi = new TLegend(0.25, 0.71, 0.38, 0.86, "", "brNDC");
  legPhi->AddEntry(hPhiProton, "proton", "pl");
  legPhi->AddEntry(hPhiDeuteron, "deuteron", "lp");
  legPhi->AddEntry(hPhiHelium, "helium3", "lp");
  legPhi->Draw();

  // v2 distributions
  hV2Proton->Sumw2();
  hV2Deuteron->Sumw2();
  hV2Helium->Sumw2();

  hV2Proton->Scale(1.0 / hV2Proton->Integral());
  hV2Deuteron->Scale(1.0 / hV2Deuteron->Integral());
  hV2Helium->Scale(1.0 / hV2Helium->Integral());

  hV2Proton->SetStats(0);
  hV2Proton->GetXaxis()->SetTitle("#phi (rad)");
  hV2Proton->GetYaxis()->SetTitle("normalized counts");

  TCanvas *cV2 = new TCanvas("cV2", "v2 distributions", 800, 600);
  cV2->DrawFrame(-1.01, -0.01, 1.001, 0.05, ";cos(2(#phi - #Psi_{2}));counts");
  hV2Proton->Draw("same");
  hV2Deuteron->Draw("same");
  hV2Helium->Draw("same");

  TLegend *legV2 = new TLegend(0.25, 0.71, 0.38, 0.86, "", "brNDC");
  legV2->AddEntry(hV2Proton, "proton", "pl");
  legV2->AddEntry(hV2Deuteron, "deuteron", "lp");
  legV2->AddEntry(hV2Helium, "helium3", "lp");
  legV2->Draw();

  // Calculate mean and error for each histogram
  double meanProton = hV2Proton->GetMean();
  double errorProton = hV2Proton->GetMeanError();
  double meanDeuteron = hV2Deuteron->GetMean();
  double errorDeuteron = hV2Deuteron->GetMeanError();
  double meanHelium = hV2Helium->GetMean();
  double errorHelium = hV2Helium->GetMeanError();

  // Add TPaveText to display the mean and error
  TPaveText *paveV2 = new TPaveText(0.65, 0.65, 0.85, 0.85, "NDC");
  paveV2->SetFillColor(0);
  paveV2->SetTextAlign(12);
  paveV2->SetBorderSize(0);
  paveV2->AddText(Form("Proton: %.3f #pm %.3f", meanProton, errorProton));
  paveV2->AddText(Form("Deuteron: %.3f #pm %.3f", meanDeuteron, errorDeuteron));
  paveV2->AddText(Form("Helium3: %.3f #pm %.3f", meanHelium, errorHelium));
  paveV2->Draw();

  TFile *f = new TFile("toy_phi_MC.root", "recreate");
  hPhiProton->Write();
  hPhiDeuteron->Write();
  hPhiHelium->Write();
  cPhi->Write();
  hV2Proton->Write();
  hV2Deuteron->Write();
  hV2Helium->Write();
  cV2->Write();
}
