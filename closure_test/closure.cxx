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
#include <TString.h>

void closure(int suffix = 0, float pt_min_arg = 0.5, float pt_max_arg = 0.6)
{

  gStyle->SetOptStat(0); // Disable statistics box
  // Parameters
  const long long nTracks = 100000000;  // Number of events
  const double v2 = 0.2;                // v2 parameter
  const int nBins = 100;                // Number of bins in the histogram
  const double phiMin = 0.0;            // Minimum phi value
  const double phiMax = TMath::TwoPi(); // Maximum phi value
  float ptMin = pt_min_arg;
  float ptMax = pt_max_arg;
  float etaMin = -0.8;
  float etaMax = 0.8;
  double coalescence_parameter = 0.30;

  // Random number generator
  TRandom3 rand(0); // Seed with 0 for random initialization

  // Proton Hisograms
  TH1D *hPhiProton = new TH1D("hPhiProton", "Azimuthal Distribution;#phi;Counts", nBins, phiMin, phiMax);
  TH1D *hPtProton = new TH1D("hPtProton", "#it{p}_{T} Distribution;#it{p}_{T} (Gev/#it{c});Counts", 100, ptMin, ptMax);
  TH1D *hEtaProton = new TH1D("hEtaProton", "#eta Distribution;#eta;Counts", 100, etaMin, etaMax);

  // Helium histograms
  TH1D *hPhiHelium = new TH1D("hPhiHelium", "Azimuthal Distribution;#phi;Counts", nBins, phiMin, phiMax);
  TH1D *hPtHelium = new TH1D("hPtHelium", "#it{p}_{T} Distribution;#it{p}_{T} (Gev/#it{c});Counts", 100, 3 * ptMin, 3 * ptMax);
  TH1D *hEtaHelium = new TH1D("hEtaHelium", "#eta Distribution;#eta;Counts", 100, etaMin, etaMax);

  struct Particle
  {
    double pX;
    double pY;
    double pZ;
    long long int id;

    void print() const
    {
      std::cout << "pX: " << pX << ", pY: " << pY << ", pZ: " << pZ << std::endl;
    }

    void printPtEtaPhi() const
    {
      double pt = TMath::Sqrt(pX * pX + pY * pY);
      double phi = TMath::ATan2(pY, pX);
      double eta = 0.5 * TMath::Log((TMath::Sqrt(pt * pt + pZ * pZ) + pZ) / (TMath::Sqrt(pt * pt + pZ * pZ) - pZ));
      std::cout << "pt: " << pt << ", eta: " << eta << ", phi: " << phi << std::endl;
    }
  };

  auto canCoalesce = [&](const Particle &p1, const Particle &p2)
  {
    // Check if the two protons can coalesce
    double deltaX = p1.pX - p2.pX;
    double deltaY = p1.pY - p2.pY;
    double deltaZ = p1.pZ - p2.pZ;
    return TMath::Sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ) < coalescence_parameter; // Example threshold
  };

  std::vector<Particle> vParticles;
  long long int helium_id = 0;

  // Define the number of subsets for pt, eta, and phi
  const int nPtBins = 5;
  const int nEtaBins = 10;
  const int nPhiBins = 20;

  // Calculate the step sizes for each range
  double ptStep = (ptMax - ptMin) / nPtBins;
  double etaStep = (etaMax - etaMin) / nEtaBins;
  double phiStep = (phiMax - phiMin) / nPhiBins;

  // Loop over subsets
  for (int iPt = 0; iPt < nPtBins; ++iPt)
  {
    double ptStart = ptMin + iPt * ptStep;
    double ptEnd = ptStart + ptStep;

    for (int iEta = 0; iEta < nEtaBins; ++iEta)
    {
      double etaStart = etaMin + iEta * etaStep;
      double etaEnd = etaStart + etaStep;

      for (int iPhi = 0; iPhi < nPhiBins; ++iPhi)
      {
        double phiStart = phiMin + iPhi * phiStep;
        double phiEnd = phiStart + phiStep;

        // printf("Subset: pt [%.2f, %.2f], eta [%.2f, %.2f], phi [%.2f, %.2f]\n", ptStart, ptEnd, etaStart, etaEnd, phiStart, phiEnd);

        // Vector to store particles in the current subset
        std::vector<Particle> subsetParticles;

        // Generate particles within this subset
        for (long i = 0; i < nTracks; ++i)
        {
          long long int id = 0;

          while (true)
          {
            // Generate random phi, pt, and eta within the subset
            double phi = rand.Uniform(phiStart, phiEnd);
            double r = rand.Uniform(0, 1);
            double fPhi = 1 + 2 * v2 * TMath::Cos(phi);
            if (r < fPhi / (1 + 2 * v2))
            {
              double pt = rand.Uniform(ptStart, ptEnd);
              double eta = rand.Uniform(etaStart, etaEnd);

              // Fill histograms
              hPtProton->Fill(pt);
              hEtaProton->Fill(eta);
              hPhiProton->Fill(phi);

              // Add particle to the subset vector
              subsetParticles.push_back({pt * TMath::Cos(phi), pt * TMath::Sin(phi), pt * TMath::SinH(eta), id++});
              break;
            }
          }
        }

        // Perform coalescence within the current subset
        std::vector<long long int> usedIds;
        for (const auto &p : subsetParticles)
        {
          if (std::find(usedIds.begin(), usedIds.end(), p.id) != usedIds.end())
          {
            continue;
          }
          for (const auto &q : subsetParticles)
          {
            if (std::find(usedIds.begin(), usedIds.end(), q.id) != usedIds.end())
            {
              continue;
            }
            for (const auto &r : subsetParticles)
            {
              if (std::find(usedIds.begin(), usedIds.end(), r.id) != usedIds.end())
              {
                continue;
              }
              if (&p != &q && &p != &r && &q != &r)
              {
                if (canCoalesce(p, q) && canCoalesce(p, r) && canCoalesce(q, r))
                {
                  Particle newParticle;
                  newParticle.pX = p.pX + q.pX + r.pX;
                  newParticle.pY = p.pY + q.pY + r.pY;
                  newParticle.pZ = p.pZ + q.pZ + r.pZ;
                  newParticle.id = helium_id++;

                  usedIds.push_back(p.id);
                  usedIds.push_back(q.id);
                  usedIds.push_back(r.id);

                  double pt_helium = TMath::Sqrt(newParticle.pX * newParticle.pX + newParticle.pY * newParticle.pY);
                  double p_helium = TMath::Sqrt(pt_helium * pt_helium + newParticle.pZ * newParticle.pZ);
                  double eta_helium = 0.5 * TMath::Log((p_helium + newParticle.pZ) / (p_helium - newParticle.pZ));
                  double phi_helium = TMath::ATan2(newParticle.pY, newParticle.pX);

                  hPtHelium->Fill(pt_helium);
                  hEtaHelium->Fill(eta_helium);
                  hPhiHelium->Fill(phi_helium);

                  // std::cout << "Helium: " << std::endl;
                  // newParticle.printPtEtaPhi();
                }
              }
            }
          }
        }

        // Clear the subset vector for the next subset
        subsetParticles.clear();
      }
    }
  }

  // Draw pt and eta histograms
  TCanvas *cProton = new TCanvas("cProton", "#it{p}_{T} and #eta Distributions", 800, 600);
  cProton->Divide(3, 1);

  cProton->cd(1);

  hPhiProton->Draw();
  hPhiProton->SetMinimum(0);

  // Add a text box to display the value of v2
  TPaveText *pave = new TPaveText(0.73, 0.22, 0.86, 0.32, "NDC");
  pave->SetFillColor(0);
  pave->SetTextAlign(12);
  pave->SetBorderSize(0);
  pave->AddText(Form("v_{2} = %.2f", v2));
  pave->Draw();

  // Draw a horizontal line at the average of the counts
  double averageCounts = hPhiProton->Integral() / nBins;
  TLine *avgLine = new TLine(phiMin, averageCounts, phiMax, averageCounts);
  avgLine->SetLineColor(kRed);
  avgLine->SetLineWidth(2);
  avgLine->SetLineStyle(2); // Dashed line
  avgLine->Draw();

  cProton->cd(2);
  hPtProton->Draw();
  hPtProton->SetMinimum(0);

  cProton->cd(3);
  hEtaProton->Draw();
  hEtaProton->SetMinimum(0);

  cProton->SaveAs("pt_eta_distributions_proton.pdf");

  TCanvas *cHelium = new TCanvas("cHelium", "#it{p}_{T} and #eta Distributions", 800, 600);
  cHelium->Divide(3, 1);
  cHelium->cd(1);
  hPhiHelium->Draw();
  hPhiHelium->SetMinimum(0);
  cHelium->cd(2);
  hPtHelium->Draw();
  hPtHelium->SetMinimum(0);
  cHelium->cd(3);
  hEtaHelium->Draw();
  hEtaHelium->SetMinimum(0);
  // Draw a horizontal line at the average of the counts
  double averageCountsHelium = hPhiHelium->Integral() / nBins;
  TLine *avgLineHelium = new TLine(phiMin, averageCountsHelium, phiMax, averageCountsHelium);
  avgLineHelium->SetLineColor(kRed);
  avgLineHelium->SetLineWidth(2);
  avgLineHelium->SetLineStyle(2); // Dashed line
  avgLineHelium->Draw();
  // Draw the canvas
  cHelium->SaveAs("pt_eta_distributions_helium.pdf");

  // Save histograms to the output file
  TFile *outputFile = new TFile(Form("closure_test_%d.root", suffix), "RECREATE");
  hPhiProton->Write();
  hPtProton->Write();
  hEtaProton->Write();
  cProton->Write();

  hPhiHelium->Write();
  hPtHelium->Write();
  hEtaHelium->Write();
  cHelium->Write();
  outputFile->Close();
}
