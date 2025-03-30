#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

void closure() {

  gStyle->SetOptStat(0); // Disable statistics box
  // Parameters
  const int nTracksPerEvent = 1000000000; // Number of events
  const double v2 = 0.2;      // v2 parameter
  const int nBins = 100;      // Number of bins in the histogram
  const double phiMin = 0.0;  // Minimum phi value
  const double phiMax = TMath::TwoPi(); // Maximum phi value
  float ptMin = 0.5;
  float ptMax = 0.6;
  float etaMin = -0.8;
  float etaMax = 0.8;

  // Random number generator
  TRandom3 rand(0); // Seed with 0 for random initialization

  // Create a histogram
  TH1D *hPhi = new TH1D("hPhi", "Azimuthal Distribution;#phi;Counts", nBins, phiMin, phiMax);
  TH1D *hPt = new TH1D("hPt", "#it{p}_{T} Distribution;#it{p}_{T} (Gev/#it{c});Counts", 100, ptMin, ptMax);
  TH1D *hEta = new TH1D("hEta", "#eta Distribution;#eta;Counts", 100, etaMin, etaMax);

  struct Particle {
    double pX;
    double pY;
    double pZ;
    long long int id;

    void print() const{
      std::cout << "pX: " << pX << ", pY: " << pY << ", pZ: " << pZ << std::endl;
    }

    void printPtEtaPhi() const {
      double pt = TMath::Sqrt(pX * pX + pY * pY);
      double phi = TMath::ATan2(pY, pX);
      double eta = 0.5 * TMath::Log((TMath::Sqrt(pt * pt + pZ * pZ) + pZ) / (TMath::Sqrt(pt * pt + pZ * pZ) - pZ));
      std::cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << std::endl;
    }

  };

  std::vector<Particle> vParticles;

  // Generate random values for phi, pt, and eta
  for (int i = 0; i < nTracksPerEvent; ++i) {

    long long int id = 0;

    while (true) {
      // Generate a random phi according to the previous distribution
      double phi = rand.Uniform(phiMin, phiMax);
      double r = rand.Uniform(0, 1);
      double fPhi = 1 + 2 * v2 * TMath::Cos(phi);
      if (r < fPhi / (1 + 2 * v2)) {
        // Fill pt and eta histograms
        double pt = rand.Uniform(ptMin, ptMax);
        double eta = rand.Uniform(etaMin, etaMax);
        hPt->Fill(pt);
        hEta->Fill(eta);
        hPhi->Fill(phi);
        vParticles.push_back({pt * TMath::Cos(phi), pt * TMath::Sin(phi), pt * TMath::SinH(eta), id++});
        break;
      }
    }
  }

  auto canCoalesce = [](const Particle& p1, const Particle& p2) {
    // Check if the two protons can coalesce
    double deltaX = p1.pX - p2.pX;
    double deltaY = p1.pY - p2.pY;
    double deltaZ = p1.pZ - p2.pZ;
    return (deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) < 0.150; // Example threshold
  };

  std::vector<long long int> usedIds;
  long long int helium_id = 0;
  std::vector<Particle> vHelium;

  for (const auto& p : vParticles) {
    if (std::find(usedIds.begin(), usedIds.end(), p.id) != usedIds.end()) {
      continue;
    }
    for (const auto& q : vParticles) {
      if (std::find(usedIds.begin(), usedIds.end(), q.id) != usedIds.end()) {
        continue;
      }
      for (const auto& r : vParticles) {
        if (std::find(usedIds.begin(), usedIds.end(), r.id) != usedIds.end()) {
          continue;
        }
        if (&p != &q && &p != &r && &q != &r) {
          // Process the combination of p, q, and r
          // std::cout << "Combination: " << std::endl;
          // p.print();
          // q.print();
          // r.print();
          // std::cout << "-----------------" << std::endl;

          if (canCoalesce(p, q) && canCoalesce(p, r) && canCoalesce(q, r)) {
            // Create a new struct with summed px, py, pz
            Particle newParticle;
            newParticle.pX = p.pX + q.pX + r.pX;
            newParticle.pY = p.pY + q.pY + r.pY;
            newParticle.pZ = p.pZ + q.pZ + r.pZ;
            newParticle.id = helium_id++;

            // Add their ids to a vector of used ids
            usedIds.push_back(p.id);
            usedIds.push_back(q.id);
            usedIds.push_back(r.id);

            vHelium.push_back(newParticle);

            // Print the new Particle and used ids
            std::cout << "Helium: " << std::endl;
            newParticle.printPtEtaPhi();

          }
        }
      }
    }

  }



  // Draw pt and eta histograms
  TCanvas *c2 = new TCanvas("c2", "#it{p}_{T} and #eta Distributions", 800, 600);
  c2->Divide(3, 1);

  c2->cd(1);

  hPhi->Draw();
  hPhi->SetMinimum(0);

  // Add a text box to display the value of v2
  TPaveText *pave = new TPaveText(0.73, 0.22, 0.86, 0.32, "NDC");
  pave->SetFillColor(0);
  pave->SetTextAlign(12);
  pave->SetBorderSize(0);
  pave->AddText(Form("v_{2} = %.2f", v2));
  pave->Draw();

  // Draw a horizontal line at the average of the counts
  double averageCounts = hPhi->Integral() / nBins;
  TLine *avgLine = new TLine(phiMin, averageCounts, phiMax, averageCounts);
  avgLine->SetLineColor(kRed);
  avgLine->SetLineWidth(2);
  avgLine->SetLineStyle(2); // Dashed line
  avgLine->Draw();

  c2->cd(2);
  hPt->Draw();
  hPt->SetMinimum(0);

  c2->cd(3);
  hEta->Draw();
  hEta->SetMinimum(0);

  // Save histograms to the output file
  TFile *outputFile = new TFile("closure_test.root", "RECREATE");
  hPhi->Write();
  hPt->Write();
  hEta->Write();
  c2->SaveAs("pt_eta_distributions.png");



}
