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
   static TH1D* hist_mp_1;
   static TH1D* hist_t2;
   static TH1D* hist_pt2_2;
   static TH1D* hist_xigen2;
   static TH1D* hist_beta2;
   static TH1D* hist_log_beta2;
   static TH1D* hist_mp_2;
   static TH1D* hist_theta_1;
   static TH1D* hist_theta_x_1;
   static TH1D* hist_theta_y_1;
   static TH1D* hist_theta_2;
   static TH1D* hist_theta_x_2;
   static TH1D* hist_theta_y_2;
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
   static TH2D* hist_log_xi1_vs_log_t1;
   static TH2D* hist_log_xi2_vs_log_t2;
   static TH2D* hist_t1_vs_pt2;
   static TH2D* hist_t2_vs_pt2;
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
TH1D* Histos::hist_mp_1 = NULL;
TH1D* Histos::hist_t2 = NULL;
TH1D* Histos::hist_pt2_2 = NULL;
TH1D* Histos::hist_xigen2 = NULL;
TH1D* Histos::hist_beta2 = NULL;
TH1D* Histos::hist_log_beta2 = NULL;
TH1D* Histos::hist_mp_2 = NULL;
TH1D* Histos::hist_theta_1 = NULL;
TH1D* Histos::hist_theta_x_1 = NULL;
TH1D* Histos::hist_theta_y_1 = NULL;
TH1D* Histos::hist_theta_2 = NULL;
TH1D* Histos::hist_theta_x_2 = NULL;
TH1D* Histos::hist_theta_y_2 = NULL;
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
TH2D* Histos::hist_log_xi1_vs_log_t1 = NULL;
TH2D* Histos::hist_log_xi2_vs_log_t2 = NULL;
TH2D* Histos::hist_t1_vs_pt2 = NULL;
TH2D* Histos::hist_t2_vs_pt2 = NULL;

static double m_p = 0.938272046; // PDG

int main() { 

   // Configuration 
   // Output file name
   //std::string outputFilename( "analysis-SD-Dijets.root" );
   std::string outputFilename( "analysis-SD-Dijets_8TeV_tmp.root" );
   //std::string outputFilename( "analysis-SD-Dijets-FPMC_SDDijets_13TeV_PTMIN_30.root" );
   // Input HepMC file
   HepMC::IO_GenEvent ascii_in(
      //"/tmp/antoniov/FPMC_SDDijets_SD1_13TeV_PTMIN_20-v1_1.hepmc",
      "/tmp/antoniov/FPMC_SDDijets_SD1_8TeV_PTMIN_20-v1_1.hepmc",
       std::ios::in
   );
 
   //int maxEvents = -1;  
   int maxEvents = 50000;
  
   bool debug = true;

   bool select_parton_pt = true;
   double min_pt = 20.;
 
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
   Histos::hist_NEvents = new TH1D("hist_NEvents","",2,0,2);
   Histos::hist_NEvents->GetXaxis()->SetBinLabel(1,"All");
   Histos::hist_NEvents->GetXaxis()->SetBinLabel(2,"PartonPt");

   Histos::hist_t1     = new TH1D("hist_t1","t proton",100,0.,4.);
   Histos::hist_pt2_1  = new TH1D("hist_pt2_1","pt2 proton",100,0.,4.);
   Histos::hist_xigen1 = new TH1D("hist_xigen1","#xi proton",100,0.,0.1);
   Histos::hist_beta1  = new TH1D("hist_beta1","",100,0.,1.);
   Histos::hist_log_beta1  = new TH1D("hist_log_beta1","",200,-6.,0.);
   Histos::hist_mp_1       = new TH1D("hist_mp_1","",200,-2.,2.);

   Histos::hist_t2     = new TH1D("hist_t2","t proton",100,0.,4.);
   Histos::hist_pt2_2  = new TH1D("hist_pt2_2","pt2 proton",100,0.,4.);
   Histos::hist_xigen2 = new TH1D("hist_xigen2","#xi proton",100,0.,0.1);
   Histos::hist_beta2  = new TH1D("hist_beta2","",100,0.,1.);
   Histos::hist_log_beta2  = new TH1D("hist_log_beta2","",200,-6.,0.);
   Histos::hist_mp_2       = new TH1D("hist_mp_2","",200,-2.,2.);

   Histos::hist_theta_1 = new TH1D("hist_theta_1","",200,-150.e-06,150.e-06);
   Histos::hist_theta_x_1 = new TH1D("hist_theta_x_1","",200,-150.e-06,150.e-06);
   Histos::hist_theta_y_1 = new TH1D("hist_theta_y_1","",200,-150.e-06,150.e-06);

   Histos::hist_theta_2 = new TH1D("hist_theta_2","",200,-150.e-06,150.e-06);
   Histos::hist_theta_x_2 = new TH1D("hist_theta_x_2","",200,-150.e-06,150.e-06);
   Histos::hist_theta_y_2 = new TH1D("hist_theta_y_2","",200,-150.e-06,150.e-06);

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

   Histos::hist_log_xi1_vs_log_t1 = new TH2D("hist_log_xi1_vs_log_t1","",200,-4.,0.,200,-4.,0.);
   Histos::hist_log_xi2_vs_log_t2 = new TH2D("hist_log_xi2_vs_log_t2","",200,-4.,0.,200,-4.,0.);
  
   Histos::hist_t1_vs_pt2 = new TH2D("hist_t1_vs_pt2","",200,0.,4.,200,0.,4.);
   Histos::hist_t2_vs_pt2 = new TH2D("hist_t2_vs_pt2","",200,0.,4.,200,0.,4.);

   //=========================================================
   //
   int event_number = 0;
   HepMC::GenEvent* evt = ascii_in.read_next_event();
   //evt->print();
   //
   TString selection_name;

   while ( evt && ( maxEvents > 0 && event_number < maxEvents ) ) {

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

      selection_name = "All";
      Histos::hist_NEvents->Fill( selection_name, 1. );

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
      HepMC::GenParticle* outParton2 = 0;
      for(HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
	 if( (*it)->status() == 124 ) { outParton2 = *it; break; }
      }

      if( debug ) {
         std::cout << "Beam proton 1 px,py,pz,E= " << std::setprecision(6) << beam1->momentum().px() << ", " << beam1->momentum().py() << ", "
                                                   << beam1->momentum().pz() << ", " << beam1->momentum().e() << std::endl;
         std::cout << "Beam proton 2 px,py,pz,E= " << std::setprecision(6) << beam2->momentum().px() << ", " << beam2->momentum().py() << ", "
                                                   << beam2->momentum().pz() << ", " << beam2->momentum().e() << std::endl;
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
      double beam1_px = beam1->momentum().px();
      double beam1_py = beam1->momentum().py();
      double beam1_pz = beam1->momentum().pz(); // positive
      double beam1_e  = beam1->momentum().e();
      double beam2_px = beam2->momentum().px();
      double beam2_py = beam2->momentum().py();
      double beam2_pz = beam2->momentum().pz(); // negative 
      double beam2_e  = beam2->momentum().e();
      //---
      beam1_px = 0.; beam1_py = 0.;
      //beam1_e  = beam1_e; 
      //beam1_e  = 6500.; 
      //beam1_e  = 4000.; 
      beam1_pz = sqrt(beam1_e*beam1_e - m_p*m_p); 

      beam2_px = 0.; beam2_py = 0.;
      //beam2_e  = beam2_e;
      //beam2_e  = 6500.;
      //beam2_e  = 4000.;
      beam2_pz = sqrt(beam2_e*beam2_e - m_p*m_p);

      if( debug ) {
         std::cout << "Corr. beam proton 1 px,py,pz,E= " << std::setprecision(10) << beam1_px << ", " << beam1_py << ", "
                                                         << beam1_pz << ", " << beam1_e << std::endl;
         std::cout << "Corr. beam proton 2 px,py,pz,E= " << std::setprecision(10) << beam2_px << ", " << beam2_py << ", "
                                                         << beam2_pz << ", " << beam2_e << std::endl;
      }
      //---
      double x1 = parton1->momentum().pz()/beam1_pz;
      double x2 = parton2->momentum().pz()/beam2_pz;

      if(pdf_info) std::cout << "x1, x2, PDF scale, event scale= " << pdf_info->x1() << ", " << pdf_info->x2() 
	                     << ", " << pdf_info->scalePDF() << ", " << evt->event_scale() << std::endl;
      std::cout << "hard partonic x1, x2, scale, pT= " << x1 << ", " << x2 << ", " << scale << ", " << outParton1->momentum().perp() << std::endl;

      if( select_parton_pt && ( outParton1->momentum().perp() < min_pt || outParton2->momentum().perp() < min_pt ) ) continue;

      selection_name = "PartonPt";
      Histos::hist_NEvents->Fill( selection_name, 1. );

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
         if( debug ) {
            std::cout << "Proton 1 px,py,pz,E= " << std::setprecision(10) << proton1->momentum().px() << ", " << proton1->momentum().py() << ", "
                                                 << proton1->momentum().pz() << ", " << proton1->momentum().e() << std::endl;
         }
         double proton1_px = proton1->momentum().px();
         double proton1_py = proton1->momentum().py();
         double proton1_pz = proton1->momentum().pz();
         double proton1_e  = proton1->momentum().e();
         proton1_e = sqrt(proton1_px*proton1_px + proton1_py*proton1_py + proton1_pz*proton1_pz + m_p*m_p);
         //proton1_pz = sqrt(proton1_e*proton1_e -  proton1_px*proton1_px - proton1_py*proton1_py - m_p*m_p);
         math::XYZTLorentzVector proton1Momentum(proton1_px,proton1_py,proton1_pz,proton1_e);
         if( debug ) {
            std::cout << "Corr. proton 1 px,py,pz,E= " << std::setprecision(10) << proton1Momentum.px() << ", " << proton1Momentum.py() << ", "
                                                       << proton1Momentum.pz() << ", " << proton1Momentum.energy() << std::endl;
         }

	 math::XYZTLorentzVector gam1Momentum( pom1->momentum() );
	 double gam1Q2 = gam1Momentum.mag2();
	 math::XYZTLorentzVector t1Momentum(0,0,0,0);
	 t1Momentum += proton1Momentum;
	 t1Momentum -= beam1->momentum();
	 double t1 = t1Momentum.mag2();
	 double xigen1 = 1. - proton1Momentum.pz()/beam1_pz;
	 double beta1 = -1;
	 double PDFbeta1 = -1;
	 if(pdf_info) { PDFbeta1 = pdf_info->x1(); beta1 = PDFbeta1; }
         else         { beta1 = x1/xigen1; } 
         // theta_X, theta_y
         double theta = proton1Momentum.theta();
         double phi = proton1Momentum.phi();
         double theta_x = theta*cos( phi );
         double theta_y = theta*sin( phi );
	 Histos::hist_t1->Fill(-t1);
	 Histos::hist_pt2_1->Fill( proton1Momentum.perp2() );
	 Histos::hist_xigen1->Fill(xigen1);
	 Histos::hist_beta1->Fill(beta1);
	 Histos::hist_log_beta1->Fill( log10(beta1) );
         Histos::hist_mp_1->Fill( proton1Momentum.mass() ); 
         Histos::hist_theta_1->Fill( theta );
         Histos::hist_theta_x_1->Fill( theta_x );
         Histos::hist_theta_y_1->Fill( theta_y );
         Histos::hist_log_xi1_vs_log_t1->Fill( log10(-t1), log10(xigen1) );
         Histos::hist_t1_vs_pt2->Fill( proton1Momentum.perp2(), -t1 );
	 if(pdf_info) Histos::hist_PDFQ2_vs_beta1->Fill( log10(beta1), log10(scalePDF*scalePDF) );
	 Histos::hist_Q2_vs_beta1->Fill( log10(beta1), log10(scale*scale) );
	 std::cout << "Gamma 1 Q2= " << gam1Q2 << std::endl;
	 std::cout << "Proton 1 xi, PDFbeta, beta, t, pt2= " << xigen1 << "  " << PDFbeta1 << "  " << beta1
                                                             << "  " << t1 << "  " << proton1Momentum.perp2() << std::endl;
	 std::cout << "Proton 1 pt, eta, phi= " << proton1Momentum.pt() 
                                                << "  " << proton1Momentum.eta() 
                                                << "  " << proton1Momentum.phi() << std::endl;
	 std::cout << "Proton 1 theta, phi, theta_x, theta_y= " << theta << "  " << phi << "  " << theta_x << "  " << theta_y << std::endl;
      }	

      if(proton2 && pom2){
         if( debug ) {
            std::cout << "Proton 2 px,py,pz,E= " << proton2->momentum().px() << ", " << proton2->momentum().py() << ", "
                                                 << proton2->momentum().pz() << ", " << proton2->momentum().e() << std::endl;
         }
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
   Histos::hist_mp_1->Write();
   Histos::hist_t2->Write();
   Histos::hist_pt2_2->Write();
   Histos::hist_xigen2->Write();
   Histos::hist_beta2->Write();
   Histos::hist_log_beta2->Write();
   Histos::hist_mp_2->Write();
   Histos::hist_theta_1->Write();
   Histos::hist_theta_x_1->Write();
   Histos::hist_theta_y_1->Write();
   Histos::hist_theta_2->Write();
   Histos::hist_theta_x_2->Write();
   Histos::hist_theta_y_2->Write();
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
   Histos::hist_log_xi1_vs_log_t1->Write();
   Histos::hist_log_xi2_vs_log_t2->Write();
   Histos::hist_t1_vs_pt2->Write();
   Histos::hist_t2_vs_pt2->Write();
   file.Close();

}
