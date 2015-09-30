#define AAZZAnalysis_cxx
#include "AAZZAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TLorentzVector.h>

#include <time.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <stdlib.h>

void AAZZAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L AAZZAnalysis.C
//      Root > AAZZAnalysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   Double_t const mass_Z  = 91.1876;

   Double_t z1_xi=mass_Z-10., z1_xf=mass_Z+10., z2_xi=0., z2_xf=120., zz_xi=0., zz_xf=5000., mpps_xi=300., mpps_xf=3000., pT_i=0., pT_f=500.;
   Double_t z1_bin=80, z2_bin=80, zz_bin=120, mpps_bin =70, pT_bin=100, PPSDetX_bin=30, PPSDetY_bin=30;

   TH1D* hPtEle = new TH1D("hPtEle","hPtEle",250,0,3000);   hPtEle->SetFillColor(4);   hPtEle->GetXaxis()->SetTitle("p_{T}(GeV)");
   TH1D* hmassEle = new TH1D("hmassEle","hmassEle",250,0,2);   hmassEle->SetFillColor(4);   hmassEle->GetXaxis()->SetTitle("mass(GeV)");
   TH1D* hEtaEle = new TH1D("hEtaEle","hEtaEle",100,-4,4);   hEtaEle->SetFillColor(4);   hEtaEle->GetXaxis()->SetTitle("#eta");
   TH1D* hPhiEle = new TH1D("hPhiEle","hPhiEle",100,-4,4);   hPhiEle->SetFillColor(4);   hPhiEle->GetXaxis()->SetTitle("#phi");
   TH1D* hPtMu = new TH1D("hPtMu","hPtMu",250,0,3000);   hPtMu->SetFillColor(4);   hPtMu->GetXaxis()->SetTitle("p_{T}(GeV)");
   TH1D* hmassMu = new TH1D("hmassMu","hmassMu",200,0,0.2);   hmassMu->SetFillColor(4);   hmassMu->GetXaxis()->SetTitle("mass(GeV)");
   TH1D* hEtaMu = new TH1D("hEtaMu","hEtaMu",100,-4,4);   hEtaMu->SetFillColor(4);   hEtaMu->GetXaxis()->SetTitle("#eta");
   TH1D* hPhiMu = new TH1D("hPhiMu","hPhiMu",100,-4,4);   hPhiMu->SetFillColor(4);   hPhiMu->GetXaxis()->SetTitle("#phi");
   TH1D* hVtz = new TH1D("hVtz","hVtz",250,-25,25);   hVtz->SetFillColor(4);   hVtz->GetXaxis()->SetTitle("Vtz (cm)");
   TH1D* hVtxDiffZ = new TH1D("hVtxDiffZ","hVtxDiffZ",500,-0.05,0.05);   hVtxDiffZ->SetFillColor(4);   hVtxDiffZ->GetXaxis()->SetTitle("VtxDiffZ (cm)");
   TH1D* hExtratracks1mm = new TH1D("hExtratracks1mm","hExtratracks1mm",70,0,70); hExtratracks1mm->SetFillColor(4);   hExtratracks1mm->GetXaxis()->SetTitle("Extratracks1mm");
   TH1D* hExtratracks2mm = new TH1D("hExtratracks2mm","hExtratracks2mm",70,0,70); hExtratracks2mm->SetFillColor(4);   hExtratracks2mm->GetXaxis()->SetTitle("Extratracks2mm");
   TH1D* hExtratracks5mm = new TH1D("hExtratracks5mm","hExtratracks5mm",70,0,70); hExtratracks5mm->SetFillColor(4);   hExtratracks5mm->GetXaxis()->SetTitle("Extratracks5mm");
   TH1D* hExtratracks1cm = new TH1D("hExtratracks1cm","hExtratracks1cm",70,0,70); hExtratracks1cm->SetFillColor(4);   hExtratracks1cm->GetXaxis()->SetTitle("Extratracks1cm");
   //-------------------- Histograms for particles reconstruction --------------------------------------------
   TString type = "Eventos/", unit = "GeV", scale, scale_zz, scale_mpps;
   std::stringstream sac;
   
   //---------------------- Z's ---------------------
   sac << (z1_xf - z1_xi)/z1_bin;   scale = sac.str();
   TH1D* h_Z1ee = new TH1D("h_Z1ee","Z1 #rightarrow 2e",z1_bin,z1_xi,z1_xf);   h_Z1ee->SetFillColor(0);   
   h_Z1ee->GetXaxis()->SetTitle("m_{2e} (GeV)");   h_Z1ee->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z1uu = new TH1D("h_Z1uu","Z1 #rightarrow 2#mu",z1_bin,z1_xi,z1_xf);   h_Z1uu->SetFillColor(0);
   h_Z1uu->GetXaxis()->SetTitle("m_{2#mu} (GeV)");   h_Z1uu->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z1 = new TH1D("h_Z1","Z1 #rightarrow ll",z1_bin,z1_xi,z1_xf);   h_Z1->SetFillColor(0);
   h_Z1->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z1->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z1Cand = new TH1D("h_Z1Cand","Z1 #rightarrow ll",z1_bin,z1_xi,z1_xf);   h_Z1Cand->SetFillColor(0);
   h_Z1Cand->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z1Cand->GetYaxis()->SetTitle(type+scale+unit);
    
   sac.str("");   sac << (z2_xf - z2_xi)/z2_bin;   scale = sac.str();
   TH1D* h_Z2ee = new TH1D("h_Z2ee","Z2 #rightarrow 2e",z2_bin,z2_xi,z2_xf);   h_Z2ee->SetFillColor(0);   
   h_Z2ee->GetXaxis()->SetTitle("m_{2e} (GeV)");   h_Z2ee->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z2uu = new TH1D("h_Z2uu","Z2 #rightarrow 2#mu",z2_bin,z2_xi,z2_xf);   h_Z2uu->SetFillColor(0);   
   h_Z2uu->GetXaxis()->SetTitle("m_{2#mu} (GeV)");   h_Z2uu->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z2 = new TH1D("h_Z2","Z2 #rightarrow ll",z2_bin,z2_xi,z2_xf);   h_Z2->SetFillColor(0);
   h_Z2->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z2->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_Z2Cand = new TH1D("h_Z2Cand","Z2 #rightarrow ll",z2_bin,z2_xi,z2_xf);   h_Z2Cand->SetFillColor(0);
   h_Z2Cand->GetXaxis()->SetTitle("m_{ll} (GeV)");   h_Z2Cand->GetYaxis()->SetTitle(type+scale+unit);
    
   TH1D* h_ee = new TH1D("h_ee","Z #rightarrow 2e",z2_bin,z2_xi,z2_xf);   h_ee->SetFillColor(0);   
   h_ee->GetXaxis()->SetTitle("m_{2e} (GeV)");   h_ee->GetYaxis()->SetTitle(type+scale+unit);
   TH1D* h_uu = new TH1D("h_uu","Z #rightarrow 2#mu",z2_bin,z2_xi,z2_xf);   h_uu->SetFillColor(0);
   h_uu->GetXaxis()->SetTitle("m_{2#mu} (GeV)");   h_uu->GetYaxis()->SetTitle(type+scale+unit);
    
   //------------------------------ Z mass scatter -------------------------------------------------------
   TH2D* Zmass_scatter = new TH2D("Zmass_scatter","Z Mass Scatter", z1_bin,z1_xi,z1_xf, z2_bin,z2_xi,z2_xf);
   Zmass_scatter->GetXaxis()->SetTitle("Z on shell mass (GeV)");
   Zmass_scatter->GetYaxis()->SetTitle("Z off shell mass (GeV)");
   Zmass_scatter->SetMarkerStyle(6);
   
   //-----------------------------------------------------------------------------------------------------
    
   //-------------------- ZZ's histograms -------------------
   sac.str("");   sac << (zz_xf - zz_xi)/zz_bin;   scale_zz = sac.str();
   TH1D* h_ZZ_M_4u = new TH1D("h_ZZ_M_4u","4#mu #rightarrow ZZ (M)",zz_bin, zz_xi, zz_xf);
   h_ZZ_M_4u->SetFillColor(kRed);
   h_ZZ_M_4u->GetXaxis()->SetTitle("m_{4#mu} (GeV)");
   h_ZZ_M_4u->GetYaxis()->SetTitle(type+scale_zz+unit);
   
   TH1D* h_ZZ_M_4e = new TH1D("h_ZZ_M_4e","4e #rightarrow ZZ (M)",zz_bin, zz_xi, zz_xf);
   h_ZZ_M_4e->SetFillColor(kGreen);
   h_ZZ_M_4e->GetXaxis()->SetTitle("m_{4e} (GeV)");
   h_ZZ_M_4e->GetYaxis()->SetTitle(type+scale_zz+unit);
   
   TH1D* h_ZZ_M_eu = new TH1D("h_ZZ_M_eu","2e2#mu/2#mu2e #rightarrow ZZ (M)",zz_bin, zz_xi, zz_xf);
   h_ZZ_M_eu->SetFillColor(kYellow);
   h_ZZ_M_eu->GetXaxis()->SetTitle("m_{2e2#mu/2#mu2e} (GeV)");
   h_ZZ_M_eu->GetYaxis()->SetTitle(type+scale_zz+unit);
    
   TH1D* h_ZZ_M = new TH1D("h_ZZ_M","4l #rightarrow ZZ (M)",zz_bin, zz_xi, zz_xf);
   h_ZZ_M->SetFillColor(kBlue);
   h_ZZ_M->GetXaxis()->SetTitle("m_{4l} (GeV)");
   h_ZZ_M->GetYaxis()->SetTitle(type+scale_zz+unit);
   
   TH1D* h_ZZ_M_signal = new TH1D("h_ZZ_M_signal","4l #rightarrow ZZ (M)",zz_bin, zz_xi, zz_xf);
   h_ZZ_M_signal->SetFillColor(kBlue);
   h_ZZ_M_signal->GetXaxis()->SetTitle("m_{4l} (GeV)");
   h_ZZ_M_signal->GetYaxis()->SetTitle(type+scale_zz+unit);
   
   TH1D* h_ZZ_pT_4u = new TH1D("h_ZZ_pT_4u","4#mu #rightarrow ZZ (p_{T})",pT_bin, pT_i, pT_f);
   h_ZZ_pT_4u->SetFillColor(kRed);
   h_ZZ_pT_4u->GetXaxis()->SetTitle("pT_{4#mu} (GeV^{2})");
    
   TH1D* h_ZZ_pT_4e = new TH1D("h_ZZ_pT_4e","4e #rightarrow ZZ (p_{T})",pT_bin, pT_i, pT_f);
   h_ZZ_pT_4e->SetFillColor(kGreen);
   h_ZZ_pT_4e->GetXaxis()->SetTitle("pT_{4e} (GeV^{2})");
    
   TH1D* h_ZZ_pT_eu = new TH1D("h_ZZ_pT_eu","2e2#mu/2#mu2e #rightarrow ZZ (p_{T})",pT_bin, pT_i, pT_f);
   h_ZZ_pT_eu->SetFillColor(kYellow);
   h_ZZ_pT_eu->GetXaxis()->SetTitle("pT_{2e2#mu/2#mu2e} (GeV^{2})");
    
   TH1D* h_ZZ_pT = new TH1D("h_ZZ_pT","4l #rightarrow ZZ (p_{T})",pT_bin, pT_i, pT_f);
   h_ZZ_pT->SetFillColor(kBlue);
   h_ZZ_pT->GetXaxis()->SetTitle("pT_{4l} (GeV^{2})");
   
   TH1D* h_ZZ_pT_signal = new TH1D("h_ZZ_pT_signal","4l #rightarrow ZZ (p_{T})",pT_bin, pT_i, pT_f);
   h_ZZ_pT_signal->SetFillColor(kBlue);
   h_ZZ_pT_signal->GetXaxis()->SetTitle("pT_{4l} (GeV^{2})");
    
   TH2D* h_ZZ_pTvsM = new TH2D("h_ZZ_pTvsM","4l #rightarrow ZZ (p_{T} vs m)",zz_bin, zz_xi, zz_xf, pT_bin, pT_i, pT_f);
   h_ZZ_pTvsM->SetFillColor(kBlue);
   h_ZZ_pTvsM->GetXaxis()->SetTitle("m_{4l} (GeV)");
   h_ZZ_pTvsM->GetYaxis()->SetTitle("p_{T}^{4l} (5GeV^{-1})");
    
   TH1D* h_ZZ_Eta_4u = new TH1D("h_ZZ_Eta_4u","4#mu #rightarrow ZZ (#eta)",100, -10, 10);
   h_ZZ_Eta_4u->SetFillColor(kRed);
   h_ZZ_Eta_4u->GetXaxis()->SetTitle("#eta_{4#mu}");
    
   TH1D* h_ZZ_Eta_4e = new TH1D("h_ZZ_Eta_4e","4e #rightarrow ZZ (#eta)",100, -10, 10);
   h_ZZ_Eta_4e->SetFillColor(kGreen);
   h_ZZ_Eta_4e->GetXaxis()->SetTitle("#eta_{4e}");
    
   TH1D* h_ZZ_Eta_eu = new TH1D("h_ZZ_Eta_eu","2e2#mu/2#mu2e #rightarrow ZZ (#eta)",100, -10, 10);
   h_ZZ_Eta_eu->SetFillColor(kYellow);
   h_ZZ_Eta_eu->GetXaxis()->SetTitle("#eta_{2e2#mu/2#mu2e}");
    
   TH1D* h_ZZ_Eta = new TH1D("h_ZZ_Eta","4l #rightarrow ZZ (#eta)",100, -10, 10);
   h_ZZ_Eta->SetFillColor(kBlue);
   h_ZZ_Eta->GetXaxis()->SetTitle("#eta");
    
   TH1D* h_ZZ_Phi_4u = new TH1D("h_ZZ_Phi_4u","4#mu #rightarrow ZZ (#phi)",100, -4., 4.);
   h_ZZ_Phi_4u->SetFillColor(kRed);
   h_ZZ_Phi_4u->GetXaxis()->SetTitle("#phi_{4#mu}");
    
   TH1D* h_ZZ_Phi_4e = new TH1D("h_ZZ_Phi_4e","4e #rightarrow ZZ (#phi)",100, -4., 4.);
   h_ZZ_Phi_4e->SetFillColor(kGreen);
   h_ZZ_Phi_4e->GetXaxis()->SetTitle("#phi_{4e}");
    
   TH1D* h_ZZ_Phi_eu = new TH1D("h_ZZ_Phi_eu","2e2#mu/2#mu2e #rightarrow ZZ (#phi)",100, -4., 4.);
   h_ZZ_Phi_eu->SetFillColor(kYellow);
   h_ZZ_Phi_eu->GetXaxis()->SetTitle("#phi_{2e2#mu/2#mu2e}");
    
   TH1D* h_ZZ_Phi = new TH1D("h_ZZ_Phi","4l #rightarrow ZZ (#phi)",100, -4., 4.);
   h_ZZ_Phi->SetFillColor(kBlue);
   h_ZZ_Phi->GetXaxis()->SetTitle("#phi");
    
   TH2D* h_ZZ_PhivsEta = new TH2D("h_ZZ_PhivsEta","4l #rightarrow ZZ (#phi vs #eta)",100, -10, 10, 100, -4., 4.);
   h_ZZ_PhivsEta->GetXaxis()->SetTitle("#eta");
   h_ZZ_PhivsEta->GetYaxis()->SetTitle("#phi");
    
   TH1D* h_ZZ_Y_4u = new TH1D("h_ZZ_Y_4u","4#mu #rightarrow ZZ (y)",200, -1.5, 1.5);
   h_ZZ_Y_4u->SetFillColor(kRed);
   h_ZZ_Y_4u->GetXaxis()->SetTitle("y_{4#mu}");
   
   TH1D* h_ZZ_Y_4e = new TH1D("h_ZZ_Y_4e","4e #rightarrow ZZ (y)",200, -1.5, 1.5);
   h_ZZ_Y_4e->SetFillColor(kGreen);
   h_ZZ_Y_4e->GetXaxis()->SetTitle("y_{4e}");
    
   TH1D* h_ZZ_Y_eu = new TH1D("h_ZZ_Y_eu","2e2#mu/2#mu2e #rightarrow ZZ (y)",200, -1.5, 1.5);
   h_ZZ_Y_eu->SetFillColor(kYellow);
   h_ZZ_Y_eu->GetXaxis()->SetTitle("y_{2e2#mu/2#mu2e}");
    
   TH1D* h_ZZ_Y = new TH1D("h_ZZ_Y","4l #rightarrow ZZ (y)",200, -1.5, 1.5);
   h_ZZ_Y->SetFillColor(kBlue);
   h_ZZ_Y->GetXaxis()->SetTitle("y");

   // User variables
   Int_t ZZ4u=0, ZZ4e=0, ZZ2u2e=0, NZZ=0;
   Int_t ZZ4u_sel_mass=0, ZZ4e_sel_mass=0, ZZ2u2e_sel_mass=0, NZZ_sel_mass=0, N_sel_vertex=0;

   Int_t n_events = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      ++n_events;
      //===================================== On shell Z =======================================================================
      
      Int_t e1=-1, e2=-1, e3=-1, e4=-1, 
            u1=-1, u2=-1, u3=-1, u4=-1; 
      Int_t ZZuu=0, ZZee=0, ZZue=0, ZZeu=0;
      //Int_t lep1=-1, lep2=-1;

      TLorentzVector Z_uu, Z_ee, Zuu_off, Zee_off;
      TLorentzVector Z_on, Z_off;
      TLorentzVector q11, q12, q21, q22;
      Z_on.SetPxPyPzE(0.,0.,0.,0.); Z_off.SetPxPyPzE(0.,0.,0.,0.);

      Double_t Z_diff_mass_min = 999.;
      for( Int_t i = 0; i < n_mu; i++ )
      for( Int_t j=i+1; j < n_mu; j++ ) {
	 if( mu_charge[i]*mu_charge[j] < 0 ) {
	    if( mu_pt[i] < 5. || mu_pt[j] < 5. ||
		fabs( mu_eta[i] ) > 2.4 || fabs( mu_eta[j] ) > 2.4 ) continue;

	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(mu_px[i] + mu_px[j],
		              mu_py[i] + mu_py[j],
		              mu_pz[i] + mu_pz[j],
		              mu_energy[i] + mu_energy[j]);

	    // Select Z candidate with mass closest to nominal
	    if( fabs( Z_cand.M() - mass_Z ) < Z_diff_mass_min ) { 
	       Z_diff_mass_min = fabs( Z_cand.M() - mass_Z ); u1 = i; u2 = j;
	       Z_uu = Z_cand;
	    }
	 }
      }
     
      Z_diff_mass_min = 999.;
      for( Int_t i = 0; i < n_e; i++ ) 
      for( Int_t j=i+1; j < n_e; j++ ){
	 if( e_charge[i]*e_charge[j] < 0 ){
	    if( e_pt[i] < 7. || e_pt[j] < 7. ||
		fabs( e_eta[i] ) > 2.5 || fabs( e_eta[j] ) > 2.5 ) continue;
	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(e_px[i] + e_px[j],
   		              e_py[i] + e_py[j],
		              e_pz[i] + e_pz[j],
		              e_energy[i] + e_energy[j]);

	    // Select Z candidate with mass closest to nominal
	    if( fabs( Z_cand.M() - mass_Z ) < Z_diff_mass_min ) { 
	       Z_diff_mass_min = fabs( Z_cand.M() - mass_Z ); e1 = i; e2 = j;
	       Z_ee = Z_cand;
	    }
	 }
      }
    
      //============================================ Define Z1 and searches Z2 ====================================================== 
     
      //******* Case 1: On shell Z, built with muons, more close to nominal Z mass *********************
      if(  ( u1!=-1 && u2!=-1 && ( e1==-1 || e2==-1 ) ) ||
           ( u1!=-1 && u2!=-1 && e1!=-1 && e2!=-1 && fabs( Z_uu.M() - mass_Z ) <= fabs( Z_ee.M() - mass_Z ) ) ) {
         Z_on = Z_uu;
         Int_t lep1=-1, lep2=-1;
         if( mu_charge[u1] < 0 && mu_charge[u2] > 0 ){ lep1 = u1; lep2 = u2; }
         else if( mu_charge[u1] > 0 && mu_charge[u2] < 0 ){ lep1 = u2; lep2 = u1; }
         q11.SetPtEtaPhiM(mu_pt[lep1], mu_eta[lep1], mu_phi[lep1], mu_mass[lep1]);
         q12.SetPtEtaPhiM(mu_pt[lep2], mu_eta[lep2], mu_phi[lep2], mu_mass[lep2]);
         //--------------------- Selection of 1st lepton of off shell Z -----------------------------------
         Double_t max_ptu=0;
         for( Int_t i = 0; i < n_mu; i++ ){
            if( i != u1 && i != u2 ){
               if( mu_pt[i] < 5. || fabs( mu_eta[i] ) > 2.4 ) continue;
               if( mu_pt[i] > max_ptu ){
                  max_ptu = mu_pt[i]; u3 = i;
               }
            }
         }
         Double_t max_pte=0;
         for(Int_t i=0; i< n_e; i++){
            if( e_pt[i] < 7. || fabs( e_eta[i] ) > 2.5 ) continue;
            if( e_pt[i] > max_pte ){
               max_pte = e_pt[i]; e3 = i;
            }
         }
         //-------------------- Selection of 2nd lepton of off shell Z --------------------------------------------------------------------------
	 max_ptu=0;
	 for( Int_t i = 0; i < n_mu; i++ ){
	    if( u3 != -1 && mu_charge[i]*mu_charge[u3] < 0 &&
                            i!=u1 && i!=u2 && i!=u3 ){
	       if( mu_pt[i] < 5. || fabs( mu_eta[i] ) > 2.4 ) continue;
	       if( mu_pt[i] > max_ptu ){
		  max_ptu = mu_pt[i]; u4 = i;
	       }
	    }
	 }
         max_pte=0;
         for(Int_t i=0; i<n_e; i++){
            if( e3 != -1 && e_charge[i]*e_charge[e3] < 0 && i!=e3 ){
               if( e_pt[i] < 7. || fabs( e_eta[i] ) > 2.5 ) continue;
               if( e_pt[i] > max_pte ){
                  max_pte = e_pt[i]; e4 = i;
               }
            }
         }

         //------- Choice of Z2 Case 1 -------------------
         if( ( u3!=-1 && u4!=-1 && ( e3==-1 || e4==-1 ) ) || ( u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && mu_pt[u3]+mu_pt[u4] >= e_pt[e3]+e_pt[e4] ) ){
	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(mu_px[u3] + mu_px[u4],
		              mu_py[u3] + mu_py[u4],
		              mu_pz[u3] + mu_pz[u4],
		              mu_energy[u3] + mu_energy[u4]);
            Z_off = Z_cand;
            ++ZZuu;
            //-------------------------------------------------------------------------------------------------------------------
            Int_t lep1=-1, lep2=-1;
            if( mu_charge[u3] < 0 && mu_charge[u4] > 0 ){ lep1 = u3; lep2 = u4; }
            else if( mu_charge[u3] > 0 && mu_charge[u4] < 0 ){ lep1 = u4; lep2 = u3; }
            q21.SetPtEtaPhiM(mu_pt[lep1], mu_eta[lep1], mu_phi[lep1], mu_mass[lep1]);
            q22.SetPtEtaPhiM(mu_pt[lep2], mu_eta[lep2], mu_phi[lep2], mu_mass[lep2]);
         }
         //-------------------------------------------------------------------------------------------------------------------
         if( ( e3!=-1 && e4!=-1 && ( u3==-1 || u4==-1 ) ) || (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && mu_pt[u3]+mu_pt[u4] < e_pt[e3]+e_pt[e4]) ){
	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(e_px[e3] + e_px[e4],
		              e_py[e3] + e_py[e4],
		              e_pz[e3] + e_pz[e4],
		              e_energy[e3] + e_energy[e4]);
            Z_off = Z_cand;
            ++ZZue;
            //-------------------------------------------------------------------------------------------------------------------------
            Int_t lep1=-1, lep2=-1;
            if( e_charge[e3] < 0 && e_charge[e4] > 0 ){ lep1 = e3; lep2 = e4; }
            else if( e_charge[e3] > 0 && e_charge[e4] < 0 ){ lep1 = e4; lep2 = e3; }
            q21.SetPtEtaPhiM(e_pt[lep1], e_eta[lep1], e_phi[lep1], e_mass[lep1]);
            q22.SetPtEtaPhiM(e_pt[lep2], e_eta[lep2], e_phi[lep2], e_mass[lep2]);
         }
      } //End Case 1 (Z(uu) on-shell)
      //************* Case 2: On shell Z, built with electrons
      else if( ( e1!=-1 && e2!=-1 && ( u1==-1 || u2==-1 ) ) ||
               ( u1!=-1 && u2!=-1 && e1!=-1 && e2!=-1 && fabs( Z_uu.M() - mass_Z ) > fabs( Z_ee.M() - mass_Z) ) ){
         Z_on = Z_ee;
         //------------------------------------------------------------------------------------------------------------------------
         Int_t lep1=-1, lep2=-1;
         if( e_charge[e1] < 0 && e_charge[e2] > 0 ){ lep1 = e1; lep2 = e2; }
         else if( e_charge[e1] > 0 && e_charge[e2] < 0 ){ lep1 = e2; lep2 = e1; }
         q11.SetPtEtaPhiM(e_pt[lep1], e_eta[lep1], e_phi[lep1], e_mass[lep1]);
         q12.SetPtEtaPhiM(e_pt[lep2], e_eta[lep2], e_phi[lep2], e_mass[lep2]);
         //------------------------------------ Selection of 1st lepton of off shell Z ---------------------------------------------------
         Double_t max_ptu=0; 
         for( Int_t i = 0; i < n_mu; i++ ) {
            if( mu_pt[i] < 5. || fabs( mu_eta[i] ) > 2.4 ) continue;
            if( mu_pt[i] > max_ptu ){
               max_ptu = mu_pt[i]; u3 = i;
            }
         }
         Double_t max_pte=0;
         for(Int_t i=0; i< n_e; i++){
            if( i != e1 && i != e2 ){
               if( e_pt[i] < 7. || fabs( e_eta[i] ) > 2.5 ) continue;
               if( e_pt[i] > max_pte ){
                  max_pte = e_pt[i]; e3 = i;
               }
            }
         }
         //------------------------- Selection of 2nd lepton of off shell Z --------------------------------------------------------------------------
         max_ptu=0;
         for( Int_t i = 0; i < n_mu; i++ ) {
            if( u3 != -1 && mu_charge[i]*mu_charge[u3] < 0 && i!=u3 ){
               if( mu_pt[i] < 5. || fabs( mu_eta[i] ) > 2.4 ) continue;
               if( mu_pt[i] > max_ptu ){
                  max_ptu = mu_pt[i]; u4 = i;
               }
            }
         }
         max_pte=0;
         for(Int_t i=0; i<n_e; i++){
            if( e3 != -1 && e_charge[i]*e_charge[e3] < 0 && i!=e1 && i!=e2 && i!=e3 ){
               if( e_pt[i] < 7. || fabs( e_eta[i] ) > 2.5 ) continue;
               if( e_pt[i] > max_pte ) {
                  max_pte = e_pt[i]; e4 = i;
               }
            }
         }
         //------- Choice of off shell Z Case 2 -------------------
         if( ( u3!=-1 && u4!=-1 && ( e3==-1 || e4==-1 ) ) || (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && mu_pt[u3]+mu_pt[u4] >= e_pt[e3]+e_pt[e4]) ){ 
	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(mu_px[u3] + mu_px[u4],
		              mu_py[u3] + mu_py[u4],
		              mu_pz[u3] + mu_pz[u4],
		              mu_energy[u3] + mu_energy[u4]);
            Z_off = Z_cand;
            ++ZZeu;
            //-------------------------------------------------------------------------------------------------------------------
            Int_t lep1=-1, lep2=-1;
            if( mu_charge[u3] < 0 && mu_charge[u4] > 0 ){ lep1 = u3; lep2 = u4; }
            else if( mu_charge[u3] > 0 && mu_charge[u4] < 0 ){ lep1 = u4; lep2 = u3; }
            q21.SetPtEtaPhiM(mu_pt[lep1], mu_eta[lep1], mu_phi[lep1], mu_mass[lep1]);
            q22.SetPtEtaPhiM(mu_pt[lep2], mu_eta[lep2], mu_phi[lep2], mu_mass[lep2]);
         }
         //-------------------------------------------------------------------------------------------------------------------    
         if( ( e3!=-1 && e4!=-1 && ( u3==-1 || u4==-1 ) ) || (u3!=-1 && u4!=-1 && e3!=-1 && e4!=-1 && mu_pt[u3]+mu_pt[u4] < e_pt[e3]+e_pt[e4]) ){
	    TLorentzVector Z_cand; 
	    Z_cand.SetPxPyPzE(e_px[e3] + e_px[e4],
		              e_py[e3] + e_py[e4],
		              e_pz[e3] + e_pz[e4],
		              e_energy[e3] + e_energy[e4]);
            Z_off = Z_cand;
            ++ZZee;
            //-------------------------------------------------------------------------------------------------------------------------
            Int_t lep1=-1, lep2=-1;
            if( e_charge[e3] < 0 && e_charge[e4] > 0 ){ lep1 = e3; lep2 = e4; }
            else if( e_charge[e3] > 0 && e_charge[e4] < 0 ){ lep1 = e4; lep2 = e3; }
            q21.SetPtEtaPhiM(e_pt[lep1], e_eta[lep1], e_phi[lep1], e_mass[lep1]);
            q22.SetPtEtaPhiM(e_pt[lep2], e_eta[lep2], e_phi[lep2], e_mass[lep2]);
         }
      } // End Case 2 (Z(ee) on-shell)
      //==================================================================================================================================
    
      //----------    ZZ       -----------------
      Double_t M_ZZ = -999;
      Double_t pT_ZZ = -999;
      Double_t Eta_ZZ = -999;
      Double_t Phi_ZZ = -999;
      Double_t Y_ZZ = -999;
      Double_t E_ZZ = -999;
      Double_t Pz_ZZ = -999;
      cout << "Event " << n_events << ": " << endl;
      cout << "ZZuu | ZZee | ZZue | ZZeu = " << ZZuu << " | " << ZZee << " | " << ZZue << " | " << ZZeu << endl;
      if( ZZuu != 0 || ZZee != 0 || ZZue != 0 || ZZeu != 0 ){ 

         if(ZZuu != 0)                   ++ZZ4u;
         else if(ZZee != 0)              ++ZZ4e;
         else if(ZZue != 0 || ZZeu != 0) ++ZZ2u2e;
         ++NZZ; 

         if ( Z_off.M() < z2_xi || Z_off.M() > z2_xf ) continue;
         if ( Z_on.M() < z1_xi || Z_on.M() > z1_xf ) continue;

         h_Z1->Fill(Z_on.M()); h_Z1uu->Fill(Z_on.M()); h_uu->Fill(Z_on.M());  
         h_Z2->Fill(Z_off.M()); h_Z2uu->Fill(Z_off.M()); h_uu->Fill(Z_off.M());

         M_ZZ = ((Z_on+Z_off).M());
         pT_ZZ = ((Z_on+Z_off).Pt());
         Eta_ZZ = ((Z_on+Z_off).Eta());
         E_ZZ = ( (Z_on+Z_off).E() );
         Pz_ZZ = ( (Z_on+Z_off).Pz() ); 

         Phi_ZZ = TMath::ATan2( (Z_on+Z_off).Py() , (Z_on+Z_off).Px() );
         if ( Phi_ZZ >= TMath::Pi() ) Phi_ZZ -= 2*TMath::Pi(); // [-pi,pi)
         if ( Phi_ZZ < -TMath::Pi() ) Phi_ZZ += 2*TMath::Pi();
         Y_ZZ = 0.5*log( (E_ZZ + Pz_ZZ)/(E_ZZ - Pz_ZZ) );

	 if(ZZuu != 0){
	    ++ZZ4u_sel_mass;
	    h_ZZ_M_4u->Fill(M_ZZ);
	    h_ZZ_pT_4u->Fill(pT_ZZ);
	    h_ZZ_Eta_4u->Fill(Eta_ZZ);
	    h_ZZ_Phi_4u->Fill(Phi_ZZ);
	    h_ZZ_Y_4u->Fill(Y_ZZ);
	 } else if(ZZee != 0){
	    ++ZZ4e_sel_mass;
	    h_ZZ_M_4e->Fill(M_ZZ);
	    h_ZZ_pT_4e->Fill(pT_ZZ);
	    h_ZZ_Eta_4e->Fill(Eta_ZZ);
	    h_ZZ_Phi_4e->Fill(Phi_ZZ);
	    h_ZZ_Y_4e->Fill(Y_ZZ);
         } else if(ZZue != 0 || ZZeu != 0){
            ++ZZ2u2e_sel_mass;
	    h_ZZ_M_eu->Fill(M_ZZ);
	    h_ZZ_pT_eu->Fill(pT_ZZ);
	    h_ZZ_Eta_eu->Fill(Eta_ZZ);
	    h_ZZ_Phi_eu->Fill(Phi_ZZ);
	    h_ZZ_Y_eu->Fill(Y_ZZ);
         }
         ++NZZ_sel_mass;

         h_ZZ_M->Fill(M_ZZ);
         h_ZZ_pT->Fill(pT_ZZ);
         h_ZZ_Eta->Fill(Eta_ZZ);
         h_ZZ_Phi->Fill(Phi_ZZ);
         h_ZZ_PhivsEta->Fill(Eta_ZZ,Phi_ZZ);
         h_ZZ_Y->Fill(Y_ZZ);
         h_ZZ_pTvsM->Fill( M_ZZ , pT_ZZ );

         Zmass_scatter->Fill(Z_on.M(),Z_off.M());
      } 
      // End of ZZ analysis 
      //=====================================================================    

   }
   cout << "Number of entries analized: " << n_events << endl; 
   cout << "----------------------------------------------------------------------------" << endl;
   cout << "ZZ4u= " << ZZ4u << "  |ZZ4e= " << ZZ4e << "  |ZZ2u2e= " << ZZ2u2e << endl;
   cout << "NZZ= " << NZZ << endl;

   TString extension, name, path;
   extension = ".root";
   name = "GammaGammaZZ";
   path = "root/";
   TFile f(path+name+extension,"RECREATE");
   
   //Double_t scale = norm/n_events;
    
   hPtEle->Write();
   hEtaEle->Write();
   hPhiEle->Write();
   hPtMu->Write();
   hEtaMu->Write();
   hPhiMu->Write();
   hmassEle->Write();
   hmassMu->Write();

   hVtz->Write();
   hVtxDiffZ->Write();
   
   hExtratracks1mm->Write();
   hExtratracks2mm->Write();
   hExtratracks5mm->Write();
   hExtratracks1cm->Write();

   h_ee->Write();
   h_uu->Write();
   h_Z1ee->Write();
   h_Z2ee->Write();
   h_Z1uu->Write();
   h_Z2uu->Write();
   h_Z1->Write();
   h_Z2->Write();
   Zmass_scatter->Write();
   
   h_ZZ_M_4u->Write();
   h_ZZ_M_4e->Write();
   h_ZZ_M_eu->Write();
   h_ZZ_M->Write();
   h_ZZ_M_signal->Write();
   h_ZZ_pT_4u->Write();
   h_ZZ_pT_4e->Write();
   h_ZZ_pT_eu->Write();
   h_ZZ_pT->Write();
   h_ZZ_pT_signal->Write();
   h_ZZ_Eta_4u->Write();
   h_ZZ_Eta_4e->Write();
   h_ZZ_Eta_eu->Write();
   h_ZZ_Eta->Write();
   h_ZZ_Phi_4u->Write();
   h_ZZ_Phi_4e->Write();
   h_ZZ_Phi_eu->Write();
   h_ZZ_Phi->Write();
   h_ZZ_Y_4u->Write();
   h_ZZ_Y_4e->Write();
   h_ZZ_Y_eu->Write();
   h_ZZ_Y->Write();
   h_ZZ_pTvsM->Write();
   h_ZZ_PhivsEta->Write();

   f.Close();

   //----------------- Information saved in txt --------------------
   time_t date;
   struct tm* timeinfo;
   time(&date);
   timeinfo = localtime(&date);
   
   FILE* ZZ_analysis; 
   ZZ_analysis = fopen ("out_ZZ_analysis.txt","a"); 
   fprintf(ZZ_analysis,"%s\n\n","================================================================================================================================");
   fprintf(ZZ_analysis,"%s","ZZ -> 4l analysis\n");
   fprintf(ZZ_analysis,"%s\n","--------------------------------------------------------------------------------------------------------------------------------");
   fprintf(ZZ_analysis,"%s",asctime(timeinfo));
   fprintf(ZZ_analysis,"%s\n","--------------------------------------------------------------------------------------------------------------------------------");
   fprintf(ZZ_analysis,"%s","No smearing\n");
   fprintf(ZZ_analysis,"%s\n","------------------------------------------------------------------------------------------------------------------");
   fprintf(ZZ_analysis,"%s\t\t %d\n","N events:           ",n_events);
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ4u:               ",ZZ4u);
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ4e:               ",ZZ4e);
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ2u2e:             ",ZZ2u2e);
   fprintf(ZZ_analysis,"%s\t\t %d\n","NZZ:                ",NZZ);
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ4u_sel_mass:      ",ZZ4u_sel_mass);
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ4e_sel_mass:      ",ZZ4e_sel_mass); 
   fprintf(ZZ_analysis,"%s\t\t %d\n","ZZ2u2e_sel_mass:    ",ZZ2u2e_sel_mass); 
   fprintf(ZZ_analysis,"%s\t\t %d\n","NZZ_sel_mass:       ",NZZ_sel_mass);  
   fprintf(ZZ_analysis,"%s\n\n","================================================================================================================================");
   fprintf(ZZ_analysis,"%s\n","");
   fclose(ZZ_analysis);
   //--------------------------------------------------------------------------

   cout << "ROOT file created: " << path+name+extension << endl;       
   cout << "____________________________________________________________________________" << endl;
   cout << "============================================================================\n" << endl;
}
