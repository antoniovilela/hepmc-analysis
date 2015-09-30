#ifndef ANALYSIS_SD_DIJETS_HEPMC
#define ANALYSIS_SD_DIJETS_HEPMC

#include "TH1D.h"
#include "TH2D.h"
#include "Math/LorentzVector.h"

namespace math {
   typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVectorD;
   typedef XYZTLorentzVectorD XYZTLorentzVector; 
}

struct Histos {
   static TH1D* hist_NEvents;
   static TH1D* hist_t1;
   static TH1D* hist_pt2_1;
   static TH1D* hist_xigen1;
   static TH1D* hist_beta1;
   static TH1D* hist_log_beta1;
   static TH1D* hist_t2;
   static TH1D* hist_pt2_2;
   static TH1D* hist_xigen2;
   static TH1D* hist_beta2;
   static TH1D* hist_log_beta2;
   static TH2D* hist_Q2_vs_x1;
   static TH2D* hist_Q2_vs_x2;
   static TH2D* hist_PDFQ2_vs_x1;
   static TH2D* hist_PDFQ2_vs_x2;
   static TH2D* hist_Q2_vs_beta1;
   static TH2D* hist_Q2_vs_beta2;
   static TH2D* hist_PDFQ2_vs_beta1;
   static TH2D* hist_PDFQ2_vs_beta2;
   static TH2D* hist_PDFQ2_vs_Q2;
   static TH2D* hist_PT_vs_Q;
   static TH2D* hist_PT_vs_PDFQ;
   static TH2D* hist_x2_vs_x1;
   static TH2D* hist_PDFx1_vs_x1;
   static TH2D* hist_PDFx2_vs_x2;
};

#endif

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/PdfInfo.h"
#include "HepMC/IO_GenEvent.h"

#include "TFile.h"

#include <iostream>
#include <string>
#include <fstream>

TH1D* Histos::hist_NEvents = NULL;
TH1D* Histos::hist_t1 = NULL;
TH1D* Histos::hist_pt2_1 = NULL;
TH1D* Histos::hist_xigen1 = NULL;
TH1D* Histos::hist_beta1 = NULL;
TH1D* Histos::hist_log_beta1 = NULL;
TH1D* Histos::hist_t2 = NULL;
TH1D* Histos::hist_pt2_2 = NULL;
TH1D* Histos::hist_xigen2 = NULL;
TH1D* Histos::hist_beta2 = NULL;
TH1D* Histos::hist_log_beta2 = NULL;
TH2D* Histos::hist_Q2_vs_x1 = NULL;
TH2D* Histos::hist_Q2_vs_x2 = NULL;
TH2D* Histos::hist_PDFQ2_vs_x1 = NULL;
TH2D* Histos::hist_PDFQ2_vs_x2 = NULL;
TH2D* Histos::hist_Q2_vs_beta1 = NULL;
TH2D* Histos::hist_Q2_vs_beta2 = NULL;
TH2D* Histos::hist_PDFQ2_vs_beta1 = NULL;
TH2D* Histos::hist_PDFQ2_vs_beta2 = NULL;
TH2D* Histos::hist_PDFQ2_vs_Q2 = NULL;
TH2D* Histos::hist_PT_vs_Q = NULL;
TH2D* Histos::hist_PT_vs_PDFQ = NULL;
TH2D* Histos::hist_x2_vs_x1 = NULL;
TH2D* Histos::hist_PDFx1_vs_x1 = NULL;
TH2D* Histos::hist_PDFx2_vs_x2 = NULL;

int main() { 

   // Configuration 
   // Output file name
   //std::string outputFilename( "analysis-SD-Dijets.root" );
   std::string outputFilename( "analysis-SD-Dijets-FPMC_SDDijets_13TeV_PTMIN_30.root" );
   // Input HepMC file
   HepMC::IO_GenEvent ascii_in(
      "/afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/FPMC/HepMC/fpmc_hepmc_test/FPMC_SDDijets_13TeV_PTMIN_30.hepmc",
      std::ios::in);
   
   bool debug = true;

   // HepPDT particle table
   /*const char infile[] = "/afs/cern.ch/work/a/antoniov/Workspace/Generators/HepPDT-3.04.01-gcc463/data/particle.tbl"; 

   std::ifstream pdfile( infile );
   if( !pdfile ) { 
      std::cerr << ">>> Cannot open " << infile << std::endl;
      exit(-1);
   }
   */

   TH1::SetDefaultSumw2(true);

   // Initialize histos 
   Histos::hist_NEvents = new TH1D("hist_NEvents","",1,0,1);
   Histos::hist_NEvents->GetXaxis()->SetBinLabel(1,"All");

   Histos::hist_t1     = new TH1D("hist_t1","t proton",100,0.,4.);
   Histos::hist_pt2_1  = new TH1D("hist_pt2_1","pt2 proton",100,0.,4.);
   Histos::hist_xigen1 = new TH1D("hist_xigen1","#xi proton",100,0.,0.1);
   Histos::hist_beta1  = new TH1D("hist_beta1","",100,0.,1.);
   Histos::hist_log_beta1  = new TH1D("hist_log_beta1","",200,-6.,0.);

   Histos::hist_t2     = new TH1D("hist_t2","t proton",100,0.,4.);
   Histos::hist_pt2_2  = new TH1D("hist_pt2_2","pt2 proton",100,0.,4.);
   Histos::hist_xigen2 = new TH1D("hist_xigen2","#xi proton",100,0.,0.1);
   Histos::hist_beta2  = new TH1D("hist_beta2","",100,0.,1.);
   Histos::hist_log_beta2  = new TH1D("hist_log_beta2","",200,-6.,0.);

   Histos::hist_Q2_vs_x1 = new TH2D("hist_Q2_vs_x1","",200,-6.,0.,200,0.,8.);
   Histos::hist_Q2_vs_x2 = new TH2D("hist_Q2_vs_x2","",200,-6.,0.,200,0.,8.);

   Histos::hist_PDFQ2_vs_x1 = new TH2D("hist_PDFQ2_vs_x1","",200,-6.,0.,200,0.,8.);
   Histos::hist_PDFQ2_vs_x2 = new TH2D("hist_PDFQ2_vs_x2","",200,-6.,0.,200,0.,8.);
    
   Histos::hist_Q2_vs_beta1 = new TH2D("hist_Q2_vs_beta1","",200,-6.,0.,200,0.,8.);
   Histos::hist_Q2_vs_beta2 = new TH2D("hist_Q2_vs_beta2","",200,-6.,0.,200,0.,8.);
   Histos::hist_PDFQ2_vs_beta1 = new TH2D("hist_PDFQ2_vs_beta1","",200,-6.,0.,200,0.,8.);
   Histos::hist_PDFQ2_vs_beta2 = new TH2D("hist_PDFQ2_vs_beta2","",200,-6.,0.,200,0.,8.);
 
   Histos::hist_PDFQ2_vs_Q2 = new TH2D("hist_PDFQ2_vs_Q2","",200,0.,8.,200,0.,8.);
   Histos::hist_PT_vs_Q     = new TH2D("hist_PT_vs_Q","",200,0.,400.,200,0.,400.);
   Histos::hist_PT_vs_PDFQ  = new TH2D("hist_PT_vs_PDFQ","",200,0.,400.,200,0.,400.);

   Histos::hist_x2_vs_x1    = new TH2D("hist_x2_vs_x1","",200,-6.,0.,200,-6.,0.);
   Histos::hist_PDFx1_vs_x1 = new TH2D("hist_PDFx1_vs_x1","",200,-6.,0.,200,-6.,0.);
   Histos::hist_PDFx2_vs_x2 = new TH2D("hist_PDFx2_vs_x2","",200,-6.,0.,200,-6.,0.);
  
   //=========================================================
   //
   int event_number = 0;
   HepMC::GenEvent* evt = ascii_in.read_next_event();
   //evt->print();
   //
   char selection_name[20];

   while ( evt ) {
      ++event_number;

      // Access PDF info
      HepMC::PdfInfo const* pdf_info = evt->pdf_info();

      if ( event_number % 50 == 1 ) {

	 std::cout << "Processing Event " << event_number
	           << " HepMC # " << evt->event_number() << std::endl;

	 if(debug){ 
	    int idx = 0;
	    for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it,++idx) {
	       std::cout << idx << "  " << (*it)->status() << "  " << (*it)->pdg_id() << std::endl;
	    }
	 }
      }

      sprintf(selection_name,"All");
      Histos::hist_NEvents->Fill( 0 );

      // Find inital partons
      HepMC::GenParticle* beam1 = 0;
      HepMC::GenParticle* beam2 = 0;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 101 ) { beam1 = *it; break; }
      }
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 102 ) { beam2 = *it; break; }
      }

      HepMC::GenParticle* pom1 = 0;
      HepMC::GenParticle* pom2 = 0;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->pdg_id() == 990 && (*it)->status() == 3 ) {
	    if( pom1 == 0 ) pom1 = *it;
	    else if( pom2 == 0) pom2 = *it;
	 }
	 if( pom1 != 0 && pom2 != 0) break;
      }

      HepMC::GenParticle* parton1 = 0;
      HepMC::GenParticle* parton2 = 0;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 121 ) { parton1 = *it; break; }
      }
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 122 ) { parton2 = *it; break; }
      }

      HepMC::GenParticle* outParton1 = 0;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 123 ) { outParton1 = *it; break; }
      }

      math::XYZTLorentzVector totMomentum(0,0,0,0);
      totMomentum += parton1->momentum();
      totMomentum += parton2->momentum();
      math::XYZTLorentzVector tMomentum(0,0,0,0);
      tMomentum += parton1->momentum();
      tMomentum -= outParton1->momentum();
      math::XYZTLorentzVector uMomentum(0,0,0,0);
      uMomentum += parton2->momentum();
      uMomentum -= outParton1->momentum();
      double s = totMomentum.mag2();
      double t = tMomentum.mag2();
      double u = uMomentum.mag2();
      //double scale = totMomentum.mass();
      double scale = sqrt(2*s*t*u/(s*s + t*t + u*u));
      double beam1_pz = beam1->momentum().pz(); // positive
      double beam2_pz = beam2->momentum().pz(); // negative 
      double x1 = parton1->momentum().pz()/beam1_pz;
      double x2 = parton2->momentum().pz()/beam2_pz;

      if(pdf_info) std::cout << "x1, x2, PDF scale, event scale= " << pdf_info->x1() << ", " << pdf_info->x2() 
	                     << ", " << pdf_info->scalePDF() << ", " << evt->event_scale() << std::endl;
      std::cout << "hard partonic x1, x2, scale, pT= " << x1 << ", " << x2 << ", " << scale << ", " << outParton1->momentum().perp() << std::endl;

      // Look for protons
      HepMC::GenParticle* proton1 = 0;
      HepMC::GenParticle* proton2 = 0;
      double pz1max = 0.;
      double pz2min = 0.;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 double pz = (*it)->momentum().pz();
	 if(((*it)->pdg_id() == 2212)&&((*it)->status() == 1)&&(pz > 0.75*beam1_pz)){
	    if(pz > pz1max){proton1 = *it;pz1max=pz;}
	 } else if(((*it)->pdg_id() == 2212)&&((*it)->status() == 1)&&(pz < 0.75*beam2_pz)){
	    if(pz < pz2min){proton2 = *it;pz2min=pz;}
	 }		
      }

      Histos::hist_Q2_vs_x1->Fill( log10(x1), log10(scale*scale) );
      Histos::hist_Q2_vs_x2->Fill( log10(x2), log10(scale*scale) );
      double scalePDF = -1;
      if(pdf_info){
         scalePDF = pdf_info->scalePDF();
	 Histos::hist_PDFQ2_vs_x1->Fill( log10(x1), log10(scalePDF*scalePDF) );
	 Histos::hist_PDFQ2_vs_x2->Fill( log10(x2), log10(scalePDF*scalePDF) );
	 Histos::hist_PDFQ2_vs_Q2->Fill( log10(scale*scale), log10(scalePDF*scalePDF) );
      }
      Histos::hist_PT_vs_Q->Fill( scale, outParton1->momentum().perp() );
      Histos::hist_PT_vs_PDFQ->Fill( scalePDF, outParton1->momentum().perp() );
      Histos::hist_x2_vs_x1->Fill( log10(x1), log10(x2) );
      double PDFx1 = -1;
      double PDFx2 = -1;
      if(pdf_info) {
         PDFx1 = pdf_info->x1(); 
         PDFx2 = pdf_info->x2(); 
	 Histos::hist_PDFx1_vs_x1->Fill( log10(x1), log10(PDFx1) );
	 Histos::hist_PDFx2_vs_x2->Fill( log10(x2), log10(PDFx2) );
      }

      if(proton1 && pom1){
	 math::XYZTLorentzVector gam1Momentum( pom1->momentum() );
	 double gam1Q2 = gam1Momentum.mag2();
	 math::XYZTLorentzVector t1Momentum(0,0,0,0);
	 t1Momentum += proton1->momentum();
	 t1Momentum -= beam1->momentum();
	 double t1 = t1Momentum.mag2();
	 double xigen1 = 1 - proton1->momentum().pz()/beam1_pz;
	 double beta1 = -1;
	 double PDFbeta1 = -1;
	 if(pdf_info) { PDFbeta1 = pdf_info->x1(); beta1 = PDFbeta1; }
         else         { beta1 = x1/xigen1; } 
	 Histos::hist_t1->Fill(-t1);
	 Histos::hist_pt2_1->Fill( proton1->momentum().perp2() );
	 Histos::hist_xigen1->Fill(xigen1);
	 Histos::hist_beta1->Fill(beta1);
	 Histos::hist_log_beta1->Fill( log10(beta1) );
	 if(pdf_info) Histos::hist_PDFQ2_vs_beta1->Fill( log10(beta1), log10(scalePDF*scalePDF) );
	 Histos::hist_Q2_vs_beta1->Fill( log10(beta1), log10(scale*scale) );
	 std::cout << "Gamma 1 Q2= " << gam1Q2 << std::endl;
	 std::cout << "Proton 1 xi, PDFbeta, beta, t, pt2= " << xigen1 << "  " << PDFbeta1 << "  " << beta1
                                                             << "  " << t1 << "  " << proton1->momentum().perp2() << std::endl;
	 std::cout << "Proton 1 pt, eta, phi= " << proton1->momentum().perp() 
                                                << "  " << proton1->momentum().pseudoRapidity() 
                                                << "  " << proton1->momentum().phi() << std::endl;
      }	

      if(proton2 && pom2){
	 math::XYZTLorentzVector gam2Momentum( pom2->momentum() );
	 double gam2Q2 = gam2Momentum.mag2();
	 math::XYZTLorentzVector t2Momentum(0,0,0,0);
	 t2Momentum += proton2->momentum();
	 t2Momentum -= beam2->momentum();
	 double t2 = t2Momentum.mag2();
	 double xigen2 = 1 - proton2->momentum().pz()/beam2_pz;
	 double beta2 = -1;
	 double PDFbeta2 = -1;
	 if(pdf_info) { PDFbeta2 = pdf_info->x1(); beta2 = PDFbeta2; }
         else         { beta2 = x2/xigen2; } 
	 Histos::hist_t2->Fill(-t2);
	 Histos::hist_pt2_2->Fill( proton2->momentum().perp2() );
	 Histos::hist_xigen2->Fill(xigen2);
	 Histos::hist_beta2->Fill(beta2);
	 Histos::hist_log_beta2->Fill( log10(beta2) );
	 if(pdf_info) Histos::hist_PDFQ2_vs_beta2->Fill( log10(beta2), log10(scalePDF*scalePDF) );
	 Histos::hist_Q2_vs_beta2->Fill( log10(beta2), log10(scale*scale) );
	 std::cout << "Gamma 2 Q2= " << gam2Q2 << std::endl;
	 std::cout << "Proton 2 xi, PDFbeta, beta, t, pt2= " << xigen2 << "  " << PDFbeta2 << "  " << beta2
                                                             << "  " << t2 << "  " << proton2->momentum().perp2() << std::endl;
	 std::cout << "Proton 2 pt, eta, phi= " << proton2->momentum().perp() 
                                                << "  " << proton2->momentum().pseudoRapidity() 
                                                << "  " << proton2->momentum().phi() << std::endl;
      }

      delete evt;
      ascii_in >> evt;
   }

   // Save histograms into file
   TFile file(outputFilename.c_str(),"RECREATE");
   Histos::hist_NEvents->Write();
   Histos::hist_t1->Write();
   Histos::hist_pt2_1->Write();
   Histos::hist_xigen1->Write();
   Histos::hist_beta1->Write();
   Histos::hist_log_beta1->Write();
   Histos::hist_t2->Write();
   Histos::hist_pt2_2->Write();
   Histos::hist_xigen2->Write();
   Histos::hist_beta2->Write();
   Histos::hist_log_beta2->Write();
   Histos::hist_Q2_vs_x1->Write();
   Histos::hist_Q2_vs_x2->Write();
   Histos::hist_PDFQ2_vs_x1->Write();
   Histos::hist_PDFQ2_vs_x2->Write();
   Histos::hist_Q2_vs_beta1->Write();
   Histos::hist_Q2_vs_beta2->Write(); 
   Histos::hist_PDFQ2_vs_beta1->Write();
   Histos::hist_PDFQ2_vs_beta2->Write(); 
   Histos::hist_PDFQ2_vs_Q2->Write(); 
   Histos::hist_PT_vs_Q->Write();
   Histos::hist_PT_vs_PDFQ->Write();
   Histos::hist_x2_vs_x1->Write();
   Histos::hist_PDFx1_vs_x1->Write();
   Histos::hist_PDFx2_vs_x2->Write();
   file.Close();

}
