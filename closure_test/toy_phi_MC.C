#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLorentzVector.h>


double GetPhiInRange(double phi)
{
  double result = phi;
  while (result < 0) {
    result = result + 2. * TMath::Pi() / 2;
  }
  while (result > 2. * TMath::Pi() / 2) {
    result = result - 2. * TMath::Pi() / 2;
  }
  return result;
}

Double_t Poly2(Double_t *x, Double_t *par)
{
  return  par[0]*(1+2.0*par[1]*TMath::Cos(2.0*x[0]));
}

Double_t Poly3(Double_t *x, Double_t *par)
{
  return  par[0]*(1 +2.0*par[1]*TMath::Power(TMath::Cos(2.0*x[0]),1.0) + par[2]*TMath::Power(TMath::Sin(3.0*x[0]),1.0));
}

Double_t v2(Double_t *x, Double_t *par)
{
  return  1+ 2*0.1*TMath::Power(TMath::Cos(2.0*x[0]),1.0);
}

Double_t calcrestmom(TLorentzVector mother, TLorentzVector daug){
  double beta = mother.Beta();
  double betax = mother.Px()/mother.E();
  double betay = mother.Py()/mother.E();
  double betaz = mother.Pz()/mother.E();
  daug.Boost(-betax, -betay, -betaz);
  return daug.P();
}

void toy_phi_MC()
{
  TF1 *f2=new TF1("f2","1+2.0*0.5*TMath::Cos(2.0*x)",0,TMath::Pi());
  TF1 *f3=new TF1("f3","1+2.0*0.5*TMath::Cos(2.0*x)",0,TMath::Pi());
  TF1 *f4=new TF1("f4","1+2.0*0.5*TMath::Cos(2.0*x)",0,TMath::Pi());

  float pt1, eta1, pt2, eta2, pt3, eta3;
  float phi1, phi2, phi3;

  TRandom3 rand;

  TH1F *h1 = new TH1F("h1","h1", 36, 0.0, TMath::Pi());
  TH1F *h2 = new TH1F("h2","h2", 36, 0.0, TMath::Pi());
  TH1F *h3 = new TH1F("h3","h3", 36, 0.0, TMath::Pi());

  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(4);

  TLorentzVector mother, deuteron, helium3;
  TLorentzVector d1, d2, d3;
  TLorentzVector d1dummy, d2dummy, d3dummy;

  for(int i = 0; i < 5000000; i++)
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

    h1->Fill(GetPhiInRange(d1.Phi()));

    if (TMath::Abs(proton1restmom - proton2restmom) < 0.3) {
      h2->Fill(GetPhiInRange(deuteron.Phi()));
    }

    if (TMath::Abs(p1 - p2) < 0.3 && TMath::Abs(p1 - p3) < 0.3 && TMath::Abs(p2 - p3) < 0.3) {
      h3->Fill(GetPhiInRange(helium3.Phi()));
    }
  }

  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();

  h1->Scale(1.0 / h1->Integral());
  h2->Scale(1.0 / h2->Integral());
  h3->Scale(1.0 / h3->Integral());

  h1->SetStats(0);
  h1->GetXaxis()->SetTitle("#phi (rad)");
  h1->GetYaxis()->SetTitle("normalized counts");
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  
  TLegend *leg = new TLegend(0.404873, 0.694853, 0.705029, 0.894853);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetTextColor(1);
  leg->SetBorderSize(0);
  leg->AddEntry(h1, "proton", "pl");
  leg->AddEntry(h2, "deuteron", "lp");
  leg->AddEntry(h3, "helium3", "lp");
  leg->Draw();

}
