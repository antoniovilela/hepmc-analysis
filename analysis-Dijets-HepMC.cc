
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IO_GenEvent.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TDataType.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <math.h>
#include <algorithm>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector> 
#include <map> 
//#include <cstdio>

using namespace std;
using namespace fastjet;

const double M_P = 0.9382720813;

class IsFinalState {
public:
  /// returns true if the GenParticle does not decay
  bool operator()( const HepMC::GenParticle* p ) { 
    if ( !p->end_vertex() && p->status()==1 ) return 1;
    return 0;
  }
};

//===================================
// Seleção de eventos
//===================================
class SelectEvent {
public:
  /// check this event for goodness
  bool operator()( const HepMC::GenEvent* evt ) { 
    return true;
  }
};

//========================================
// Análise principal
//=========================================

void print_help(std::vector<std::string> const& required_parameters)
{
   std::stringstream oss;
   oss << "Usage: analysis-Dijets-HepMC.exe ";
   std::vector<std::string>::const_iterator it_par = required_parameters.begin();
   std::vector<std::string>::const_iterator pars_end = required_parameters.end();
   for(; it_par != pars_end; ++it_par){
      oss << "--" <<  *it_par << " <" << *it_par << "> ";
   }
   oss << std::endl;
   std::cout << oss.str(); 
}

int main(int argc, char** argv) { 

   vector<string> required_parameters_;
   required_parameters_.push_back("hepmc");
   required_parameters_.push_back("out");
   
   // Read command line parameters 
   std::vector<string> command_line_parameters_;
   if(argc > 1){
      command_line_parameters_.resize(argc - 1);
      std::copy(argv + 1, argv + argc, command_line_parameters_.begin());
   }

   // Help option
   std::vector<std::string> help_str; help_str.push_back("-h"); help_str.push_back("--help");
   for(std::vector<std::string>::const_iterator it_help = help_str.begin();
                                                it_help != help_str.end(); ++it_help){
      if( std::find(command_line_parameters_.begin(), command_line_parameters_.end(), *it_help) != command_line_parameters_.end() ){
         print_help(required_parameters_);
         return 0;
      }
   }

   // Read required parameters
   std::map<std::string,std::string> required_parameters_map_;
   for(std::vector<std::string>::const_iterator it_par = required_parameters_.begin();
                                                it_par != required_parameters_.end(); ++it_par){
      std::stringstream par_ss; par_ss << "--"; par_ss << *it_par;
      std::vector<std::string>::const_iterator it_par_key = std::find(command_line_parameters_.begin(), command_line_parameters_.end(), par_ss.str());

      if( it_par_key == command_line_parameters_.end() ){
	 std::stringstream oss;
	 oss << "ERROR: The following parameter was not set: " << *it_par << std::endl;
	 throw runtime_error( oss.str() );
      }

      std::vector<std::string>::const_iterator it_par_value = it_par_key + 1;

      if( it_par_value == command_line_parameters_.end() ||
          std::find(required_parameters_.begin(), required_parameters_.end(), *it_par_value) != required_parameters_.end() ){
         std::stringstream oss;
	 oss << "ERROR: Invalid value for parameter: " << *it_par << std::endl;      
	 throw runtime_error( oss.str() );

      }
      required_parameters_map_[*it_par] = *it_par_value;
   }

   bool debug = true;
   int maxEvents = -1;
   double EBeam_ = 6500.;

   std::string const filein_ = required_parameters_map_["hepmc"]; 

   std::string const outputFileName = required_parameters_map_["out"];
 
   HepMC::IO_GenEvent ascii_in(filein_.c_str(),std::ios::in);

   // Declare another IO_GenEvent for writing out the good events
   /*size_t filein_aux_pos_1 = filein_.rfind("/");
   std::string filein_substr = filein_.substr( ( ( filein_aux_pos_1 < filein_.size() ) ? (filein_aux_pos_1 + 1) : 0 ) , (filein_.size() - filein_aux_pos_1) );
   HepMC::IO_GenEvent ascii_out_4mu( ( filein_substr.substr(0,filein_substr.find(".")) + std::string("_selected_4mu.hepmc") ).c_str(),std::ios::out);
   HepMC::IO_GenEvent ascii_out_4e( ( filein_substr.substr(0,filein_substr.find(".")) + std::string("_selected_4e.hepmc") ).c_str(),std::ios::out);
   HepMC::IO_GenEvent ascii_out_2mu2e( ( filein_substr.substr(0,filein_substr.find(".")) + std::string("_selected_2mu2e.hepmc") ).c_str(),std::ios::out);
   */

   // Build HepPDT particle table
   const char infile[] = "/afs/cern.ch/work/a/antoniov/Workspace/HEP-Generators/HepPDT-3.04.01-x86_64-slc6-gcc481/data/particle.tbl";   
   std::ifstream pdfile( infile );
   if( !pdfile ) { 
      std::cerr << ">>> Cannot open " << infile << std::endl;
      exit(-1);
   }
   HepPDT::ParticleDataTable pdt( "Particle Table" );
   {
      // Initialize table builder
      HepPDT::TableBuilder tb(pdt);
      if( !addParticleTable( pdfile, tb, true ) ) { 
	 std::cout << ">> Error reading PDG pdt file " << std::endl; 
      }
   }
   // Loop over Particle Data Table
   std::ostringstream oss;
   oss << std::setw(15) << "Particle Id"
       << std::setw(22) << "Particle Name"
       << std::setw(15) << "Three-charge" << std::endl; 
   for( HepPDT::ParticleDataTable::const_iterator p = pdt.begin(); p != pdt.end(); ++p ) {
      const HepPDT::ParticleID & id = p->first;
      int pdgId = id.pid();
      int q3 = id.threeCharge();
      //double q = id.charge();
      const std::string& name = id.PDTname();
      oss << std::setw(15) << pdgId
 	  << std::setw(22) << name
 	  << std::setw(15) << q3 << std::endl;
   }
   std::cout << oss.str();

   //================================================================
   // FastJet
   //==================================================================
   double Rparam = 0.5;
   fastjet::Strategy                strategy = fastjet::Best;
   fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
   fastjet::JetDefinition            *jetDef = NULL;
   jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);
   std::vector<fastjet::PseudoJet> fjInputs;
   //==================================================================

   int seed = 12345;
   TRandom3 random( seed );  

   // Declare ROOT TTree
   int const NMUMAX = 50;
   int n_mu = 0;
   double mu_pt[NMUMAX];
   double mu_px[NMUMAX];
   double mu_py[NMUMAX];
   double mu_pz[NMUMAX];
   double mu_eta[NMUMAX];
   double mu_phi[NMUMAX];
   double mu_energy[NMUMAX];
   double mu_mass[NMUMAX];
   double mu_charge[NMUMAX];

   int const NEMAX = 50;
   int n_e = 0;
   double e_pt[NEMAX];
   double e_px[NEMAX];
   double e_py[NEMAX];
   double e_pz[NEMAX];
   double e_eta[NEMAX];
   double e_phi[NEMAX];
   double e_energy[NEMAX];
   double e_mass[NEMAX];
   double e_charge[NEMAX];

   int const NCHGMAX =300;  
   int n_chg = 0;
   int chg_id[NCHGMAX];
   double chg_pt[NCHGMAX];
   double chg_px[NCHGMAX];
   double chg_py[NCHGMAX];
   double chg_pz[NCHGMAX];
   double chg_ch[NCHGMAX];
   double chg_eta[NCHGMAX];
   double chg_phi[NCHGMAX];
   double chg_energy[NCHGMAX];
   double chg_mass[NCHGMAX];

   int const PROTONMAX = 100;
   int n_proton =0;
   double proton_pt[PROTONMAX];
   double proton_px[PROTONMAX];
   double proton_py[PROTONMAX];
   double proton_pz[PROTONMAX];
   double proton_energy[PROTONMAX];
   double proton_xi[PROTONMAX];
   double proton_t[PROTONMAX];
   double proton_weight[PROTONMAX];
   int    proton_acc[PROTONMAX];

   int const JETMAX = 500;
   int n_jet =0;
   double jet_px[JETMAX];
   double jet_py[JETMAX];
   double jet_pz[JETMAX];
   double jet_pt[JETMAX];
   double jet_energy[JETMAX];

   TTree* T = new TTree("T","Tree");
   T->Branch("n_mu", &n_mu,"n_mu/I");
   T->Branch("mu_pt", &mu_pt,"mu_pt[n_mu]/D");
   T->Branch("mu_px", &mu_px,"mu_px[n_mu]/D");
   T->Branch("mu_py", &mu_py,"mu_py[n_mu]/D"); 
   T->Branch("mu_pz", &mu_pz,"mu_pz[n_mu]/D");
   T->Branch("mu_eta", &mu_eta,"mu_eta[n_mu]/D");
   T->Branch("mu_phi", &mu_phi,"mu_phi[n_mu]/D");
   T->Branch("mu_energy", &mu_energy,"mu_energy[n_mu]/D");
   T->Branch("mu_mass", &mu_mass,"mu_mass[n_mu]/D");
   T->Branch("mu_charge", &mu_charge,"mu_charge[n_mu]/D");

   T->Branch("n_e", &n_e,"n_e/I");
   T->Branch("e_pt", &e_pt,"e_pt[n_e]/D");
   T->Branch("e_px", &e_px,"e_px[n_e]/D");
   T->Branch("e_py", &e_py,"e_py[n_e]/D"); 
   T->Branch("e_pz", &e_pz,"e_pz[n_e]/D");
   T->Branch("e_eta", &e_eta,"e_eta[n_e]/D");
   T->Branch("e_phi", &e_phi,"e_phi[n_e]/D");
   T->Branch("e_energy", &e_energy,"e_energy[n_e]/D");
   T->Branch("e_mass", &e_mass,"e_mass[n_e]/D");
   T->Branch("e_charge", &e_charge,"e_charge[n_e]/D");

   T->Branch("n_chg", &n_chg,"n_chg/I");
   T->Branch("chg_id", &chg_id,"chg_id[n_chg]/I");
   T->Branch("chg_ch", &chg_ch,"chg_ch[n_chg]/D");
   T->Branch("chg_px", &chg_px,"chg_px[n_chg]/D");
   T->Branch("chg_py", &chg_py,"chg_py[n_chg]/D");
   T->Branch("chg_pz", &chg_pz,"chg_pz[n_chg]/D");
   T->Branch("chg_pt", &chg_pt,"chg_pt[n_chg]/D");
   T->Branch("chg_eta", &chg_eta,"chg_eta[n_chg]/D");
   T->Branch("chg_phi", &chg_phi,"chg_phi[n_chg]/D");
   T->Branch("chg_energy", &chg_energy,"chg_energy[n_chg]/D");
   T->Branch("chg_mass", &chg_mass,"chg_mass[n_chg]/D");

   //Proton
   T->Branch("n_proton", &n_proton,"n_proton/I");
   T->Branch("proton_px", &proton_px,"proton_px[n_proton]/D");
   T->Branch("proton_py", &proton_py,"proton_py[n_proton]/D");
   T->Branch("proton_pz", &proton_pz,"proton_pz[n_proton]/D");
   T->Branch("proton_pt", &proton_pt,"proton_pt[n_proton]/D");
   T->Branch("proton_energy", &proton_energy,"proton_energy[n_proton]/D");
   T->Branch("proton_xi", &proton_xi,"proton_xi[n_proton]/D");
   T->Branch("proton_t", &proton_t,"proton_t[n_proton]/D");
   T->Branch("proton_weight", &proton_weight,"proton_weight[n_proton]/D");
   T->Branch("proton_acc", &proton_acc,"proton_acc[n_proton]/I");

   //Jet
   T->Branch("n_jet", &n_jet,"n_jet/I");
   T->Branch("jet_px", &jet_px,"jet_px[n_jet]/D");
   T->Branch("jet_py", &jet_py,"jet_py[n_jet]/D");
   T->Branch("jet_pz", &jet_pz,"jet_pz[n_jet]/D");
   T->Branch("jet_pt", &jet_pt,"jet_pt[n_jet]/D");
   T->Branch("jet_energy", &jet_energy,"jet_energy[n_jet]/D");

   // Declarar os histogramas no ROOT 
   TH1F* h_nmu = new TH1F("nmu","nmu",10,0.,10.);
   TH1F* h_mu_pt = new TH1F("mu_pt","mu_pt",1000,0.,1000.);
   TH1F* h_mu_px = new TH1F("mu_px","mu_px",1000,0.,1000.);
   TH1F* h_mu_py = new TH1F("mu_py","mu_py",1000,0.,1000.);
   TH1F* h_mu_pz = new TH1F("mu_pz","mu_pz",1000,0.,1000.);
   TH1F* h_mu_eta = new TH1F("mu_eta","mu_eta",1000,-5.0,5.0);
   TH1F* h_mu_phi = new TH1F("mu_phi","mu_phi",1000,-M_PI,M_PI);
   TH1F* h_mu_energy = new TH1F("mu_energy","mu_energy",1000,0.,1000.);
   TH1F* h_mu_mass = new TH1F("mu_mass","mu_mass",100,0.,100.);
   TH1F* h_mu_charge = new TH1F("mu_charge","mu_charge",10,-10.,10.);

   TH1F* h_ne = new TH1F("ne","ne",10,0.,10.);
   TH1F* h_e_pt = new TH1F("e_pt","e_pt",1000,0.,1000.);
   TH1F* h_e_px = new TH1F("e_px","e_px",1000,0.,1000.);
   TH1F* h_e_py = new TH1F("e_py","e_py",1000,0.,1000.);
   TH1F* h_e_pz = new TH1F("e_pz","e_pz",1000,0.,1000.);
   TH1F* h_e_eta = new TH1F("e_eta","e_eta",1000,-5.0,5.0);
   TH1F* h_e_phi = new TH1F("e_phi","e_phi",1000,-M_PI,M_PI);
   TH1F* h_e_energy = new TH1F("e_energy","e_energy",1000,0.,1000.);
   TH1F* h_e_mass = new TH1F("e_mass","e_mass",100,0.,100.);
   TH1F* h_e_charge = new TH1F("e_charge","e_charge",10,-10.,10.);

   //Particulas carregadas
   TH1F* h_nchg = new TH1F("nchg","nchg",100,0.,100.);
   TH1F* h_chg_id = new TH1F("chg_id","chg_id",4000,-2000.,2000.);
   TH1F* h_chg_ch = new TH1F("chg_ch","chg_ch",1000,-1000.,1000.);
   TH1F* h_chg_pt = new TH1F("chg_pt","chg_pt",1000,0.,1000.);
   TH1F* h_chg_px = new TH1F("chg_px","chg_px",1000,0.,1000.);
   TH1F* h_chg_py = new TH1F("chg_py","chg_py",1000,0.,1000.);
   TH1F* h_chg_pz = new TH1F("chg_pz","chg_pz",1000,0.,1000.);
   TH1F* h_chg_eta = new TH1F("chg_eta","chg_eta",100,-5.0,5.0);
   TH1F* h_chg_phi = new TH1F("chg_phi","chg_phi",100,-M_PI,M_PI);
   TH1F* h_chg_energy = new TH1F("chg_energy","chg_energy",1000,0.,1000.);
   TH1F* h_chg_mass = new TH1F("chg_mass","chg_mass",100,0.,140.);

   TH1F* h_nproton = new TH1F("nproton","nproton",100,0.,50.);
   TH1F* h_proton_pt = new TH1F("proton_pt","proton_pt",100,0.,1000.);
   TH1F* h_proton_px = new TH1F("proton_px","proton_px",100,0.,1000.);
   TH1F* h_proton_py = new TH1F("proton_py","proton_py",100,0.,1000.);
   TH1F* h_proton_pz = new TH1F("proton_pz","proton_pz",100,0.,1000.);
   TH1F* h_proton_energy = new TH1F("proton_energy","proton_energy",1000,0.,1000.);
   TH1F* h_proton_xi = new TH1F("proton_xi","proton_xi",100,0.,1.);
   TH1F* h_proton_t = new TH1F("proton_t","proton_t",100,0.,4.);
   TH1F* h_proton_weight = new TH1F("proton_weight","proton_weight",100,0.,1.);
   TH1F* h_proton_acc = new TH1F("proton_acc","proton_acc",2,0,2);

   //Jets
   TH1F* h_njet = new TH1F("njet","njet",50,0.,50.);
   TH1F* h_jetpt = new TH1F("jetpt","jetpt",400,0.,400.);

   // EVENT LOOP
   SelectEvent select_;
   IsFinalState isfinal_;

   int icount=0;
   int num_events_all=0;
   int num_good_events=0;

   HepMC::GenEvent* evt = ascii_in.read_next_event();
   while ( evt ) {

      ++icount;

      if( maxEvents >= 0 && icount > maxEvents) { delete evt; break; }

      if ( icount ) std::cout << "Processing Event Number " << icount
	                      << " its # " << evt->event_number() << std::endl;

      ++num_events_all;

      // Reset Tree variables per event
      n_mu = 0;
      for(int imu = 0; imu < NMUMAX; ++imu) {
	 mu_pt[imu] =  -999.;
	 mu_px[imu] =  -999.;
	 mu_py[imu] =  -999.; 
	 mu_pz[imu] =  -999.;
	 mu_eta[imu] = -999.;
	 mu_phi[imu] = -999.;
	 mu_energy[imu] = -999.;
	 mu_mass[imu] = -999.;
	 mu_charge[imu] = -999.; }

      n_e = 0;
      for(int ie = 0; ie < NEMAX; ++ie) {
	 e_pt[ie] = -999.;
	 e_eta[ie] = -999.;
	 e_px[ie] = -999.;
	 e_py[ie] = -999.;
	 e_pz[ie] = -999.;
	 e_phi[ie] = -999.; 
	 e_energy[ie] = -999.;
	 e_mass[ie] = -999.;
	 e_charge[ie] = -999.; }

      n_chg = 0;
      for(int ichg = 0; ichg < NCHGMAX; ++ichg) {
	 chg_id[ichg] = -999;
	 chg_ch[ichg] = -999.;
	 chg_pt[ichg] = -999.;
	 chg_px[ichg] = -999.;
	 chg_py[ichg] = -999.;
	 chg_pz[ichg] = -999.;
	 chg_eta[ichg] = -999.;
	 chg_phi[ichg] = -999.;
	 chg_energy[ichg] = -999.;
	 chg_mass[ichg] = -999.; }

      n_proton=0; 
      for(int iproton = 0; iproton < PROTONMAX; ++iproton) {
	 proton_pt[iproton] = -999.;
	 proton_px[iproton] = -999.;
	 proton_py[iproton] = -999.;
	 proton_pz[iproton] = -999.;
	 proton_energy[iproton] = -999.; 
	 proton_xi[iproton] = -999.;
	 proton_t[iproton] = -999.;
	 proton_weight[iproton] = -999.;
	 proton_acc[iproton] = -1; }

      n_jet=0;
      for(int ijet = 0; ijet < JETMAX; ++ijet) {
	 jet_px[ijet] = -999.;
	 jet_py[ijet] = -999.;
	 jet_pz[ijet] = -999.;
	 jet_pt[ijet] = -999.;
	 jet_energy[ijet] = -999.; }
      fjInputs.clear();

      // Event selection
      // FIXME    
      if ( !select_(evt) ) {
         // Finalize event
         delete evt;
         evt = ascii_in.read_next_event();   
         continue;
      }

      for( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin() ; p != evt->particles_end() ; ++p ){
	 // mu+/-
	 // pT >= 5 GeV ; |eta| < 2.4
	 if ( isfinal_(*p) && abs( (*p)->pdg_id() ) == 13 && (*p)->momentum().perp() >= 5.0 && fabs( (*p)->momentum().eta() ) <= 2.4 ) {
	    int pdg_id = (*p)->pdg_id();

	    HepPDT::ParticleData * pd_mu;
	    pd_mu = pdt.particle( HepPDT::ParticleID( pdg_id ) );
	    double mucharge = pd_mu->charge(); 
	    mu_charge[n_mu] = mucharge;
	    mu_pt[n_mu] = (*p)->momentum().perp();
	    mu_px[n_mu] = (*p)->momentum().px();
	    mu_py[n_mu] = (*p)->momentum().py();
	    mu_pz[n_mu] = (*p)->momentum().pz();
	    mu_eta[n_mu] = (*p)->momentum().eta(); 
	    mu_phi[n_mu] = (*p)->momentum().phi(); 
	    mu_energy[n_mu] = (*p)->momentum().e();
	    mu_mass[n_mu] = (*p)->momentum().m();

	    h_mu_px->Fill(mu_px[n_mu]);
	    h_mu_py->Fill(mu_py[n_mu]);
	    h_mu_pz->Fill(mu_pz[n_mu]);
	    h_mu_pt->Fill(mu_pt[n_mu]);
	    h_mu_eta->Fill(mu_eta[n_mu] );
	    h_mu_phi->Fill(mu_phi[n_mu] );
	    h_mu_energy->Fill(mu_energy[n_mu] );
	    h_mu_mass->Fill(mu_mass[n_mu] );
	    h_mu_charge->Fill(mu_charge[n_mu]);
	    ++n_mu;
	 }

	 // e+/-
	 // pT >= 5 GeV ; |eta| < 2.5
	 if ( isfinal_(*p) && abs( (*p)->pdg_id() ) == 11 && (*p)->momentum().perp() >= 5.0 && fabs( (*p)->momentum().eta() ) <= 2.5 ) {
	    int pdg_id = (*p)->pdg_id();

	    HepPDT::ParticleData * pd_e;
	    pd_e = pdt.particle( HepPDT::ParticleID( pdg_id ) );
	    double echarge = pd_e->charge(); 
	    e_charge[n_e] = echarge;
	    e_pt[n_e] = (*p)->momentum().perp();
	    e_px[n_e] = (*p)->momentum().px();
	    e_py[n_e] = (*p)->momentum().py();
	    e_pz[n_e] = (*p)->momentum().pz();
	    e_eta[n_e] = (*p)->momentum().eta(); 
	    e_phi[n_e] = (*p)->momentum().phi(); 
	    e_energy[n_e] = (*p)->momentum().e();   
	    e_mass[n_e] = (*p)->momentum().m();                      

	    h_e_pt->Fill( e_pt[n_e]);
	    h_e_px->Fill( e_px[n_e] );
	    h_e_py->Fill( e_py[n_e] );
	    h_e_pz->Fill( e_pz[n_e] );
	    h_e_eta->Fill( e_eta[n_e] );
	    h_e_phi->Fill( e_phi[n_e] );
	    h_e_energy->Fill( e_pt[n_e] );
	    h_e_mass->Fill( e_pt[n_e] );
	    h_e_charge->Fill(e_charge[n_e]);
	    ++n_e;
	 }

	 if ( abs((*p)->pdg_id() == 2212) && (*p)->status() == 1 && fabs( (*p)->momentum().pz() ) >= 0.50*EBeam_ ){
	    proton_px[n_proton] = (*p)->momentum().px();
	    proton_py[n_proton] = (*p)->momentum().py();
	    proton_pz[n_proton] = (*p)->momentum().pz();
	    proton_pt[n_proton] = (*p)->momentum().perp();
	    proton_energy[n_proton] =(*p)->momentum().e();
 
	    h_proton_px->Fill( proton_px[n_proton] );
	    h_proton_py->Fill( proton_py[n_proton] );
	    h_proton_pz->Fill( proton_pz[n_proton] );
	    h_proton_pt->Fill( proton_pt[n_proton] );
	    h_proton_energy->Fill( proton_energy[n_proton] );

            // Simulate whether proton is accepted
            double proton_xi_tmp = 1. - fabs( proton_pz[n_proton] )/EBeam_;
            TLorentzVector vec_proton(proton_px[n_proton],proton_py[n_proton],proton_pz[n_proton],proton_energy[n_proton]);
            double beam_pz = TMath::Sqrt( EBeam_*EBeam_ - M_P*M_P );
            TLorentzVector vec_beam_1(0.,0., beam_pz,EBeam_);              
            TLorentzVector vec_beam_2(0.,0.,-beam_pz,EBeam_);              
            if( proton_pz[n_proton] < 0) vec_proton -= vec_beam_2;
            else                         vec_proton -= vec_beam_1;   
            double proton_t_tmp = vec_proton.Mag2();  
             
            proton_xi[n_proton]     = proton_xi_tmp;
            proton_t[n_proton]      = proton_t_tmp;
            proton_weight[n_proton] = 1.;
            proton_acc[n_proton]    = 1;
              
	    h_proton_xi->Fill( proton_xi[n_proton] );
	    h_proton_t->Fill( fabs( proton_t[n_proton] ) );
	    h_proton_weight->Fill( proton_weight[n_proton] );
	    h_proton_acc->Fill( proton_acc[n_proton] );

	    if( debug ) { std::cout << "Found proton xi,t,weight,accepted= " << proton_xi[n_proton] << " , " << proton_t[n_proton] << " , " << proton_weight[n_proton] << " , " << proton_acc[n_proton] << std::endl; std::cout << "\t"; (*p)->print(); }

	    ++n_proton;
	 }

	 // Looking for charged particles
	 if ( isfinal_(*p) ){
	    int pdg_id = (*p)->pdg_id();
	    HepPDT::ParticleData * pd;
	    pd = pdt.particle( HepPDT::ParticleID( pdg_id ) );
	    double charge = pd->charge(); 
	    if( charge != 0. && fabs( (*p)->momentum().eta() ) <= 2.5 && (*p)->momentum().perp() >= 0.100) {
	       chg_id[n_chg] = pdg_id;
	       chg_ch[n_chg] = charge;
	       chg_pt[n_chg] = (*p)->momentum().perp();
	       chg_px[n_chg] = (*p)->momentum().px();
	       chg_py[n_chg] = (*p)->momentum().py();
	       chg_pz[n_chg] = (*p)->momentum().pz();
	       chg_eta[n_chg] = (*p)->momentum().eta(); 
	       chg_phi[n_chg] = (*p)->momentum().phi(); 
	       chg_energy[n_chg] = (*p)->momentum().e(); 
	       chg_mass[n_chg] = (*p)->momentum().m(); 

	       h_chg_id->Fill( chg_id[n_chg] );
	       h_chg_ch->Fill( chg_ch[n_chg] );
	       h_chg_pt->Fill( chg_pt[n_chg]);
	       h_chg_px->Fill( chg_px[n_chg]);
	       h_chg_py->Fill( chg_py[n_chg]);
	       h_chg_pz->Fill( chg_pz[n_chg]);
	       h_chg_eta->Fill( chg_eta[n_chg] );
	       h_chg_phi->Fill( chg_phi[n_chg] );
	       h_chg_energy->Fill( chg_energy[n_chg]);
	       h_chg_mass->Fill( chg_mass[n_chg]);
	       ++n_chg;
	    }
	 }

         // Jet particles
	 if ( isfinal_(*p) ){
	    int pdg_id = (*p)->pdg_id();

            if( pdg_id != 12 &&
                pdg_id != 14 &&
                pdg_id != 16 &&
                fabs( (*p)->momentum().eta() ) <= 5.0 ){
               fjInputs.push_back( 
                  fastjet::PseudoJet( (*p)->momentum().px(),
                                      (*p)->momentum().py(),
                                      (*p)->momentum().pz(),
                                      (*p)->momentum().e() ) );
            }
         }

      } // Loop on particles

      //=================================================
      double ptmin = 10.0;
      vector <fastjet::PseudoJet> inclusiveJets;
      vector <fastjet::PseudoJet> sortedJets;
      fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
      inclusiveJets = clustSeq.inclusive_jets( ptmin );
      sortedJets  = sorted_by_pt( inclusiveJets );
      for(unsigned int i = 0; i < sortedJets.size(); i++) {
	 jet_px[i] = sortedJets[i].px();
	 jet_py[i] = sortedJets[i].py();
	 jet_pz[i] = sortedJets[i].pz();
	 jet_pt[i] = sortedJets[i].pt();
	 jet_energy[i] = sortedJets[i].e();

	 h_jetpt->Fill(jet_pt[i]);
         ++n_jet;
      }
      //==================================================	  

      h_nmu->Fill(n_mu);
      h_ne->Fill(n_e);
      h_nproton->Fill(n_proton);
      h_nchg->Fill( n_chg );
      h_njet->Fill( n_jet );

      T->Fill();

      // Select HepMC output events
      /* ... */

      ++num_good_events;
      delete evt;
      evt = ascii_in.read_next_event();
   } // Loop on events

   std::cout << num_good_events << " out of " << num_events_all << " processed events passed the cuts. Finished." << std::endl;

   T->Print();
   // Output file
   TFile* output = new TFile( outputFileName.c_str() ,"RECREATE");
   output->cd();

   // Write TTree and histograms to file
   T->Write();
   h_nmu->Write();
   h_mu_pt->Write();
   h_mu_px->Write();
   h_mu_py->Write();
   h_mu_pz->Write();
   h_mu_eta->Write();
   h_mu_phi->Write();
   h_mu_energy->Write();
   h_mu_mass->Write();
   h_mu_charge->Write();

   h_ne->Write();
   h_e_pt->Write();
   h_e_px->Write();
   h_e_py->Write();
   h_e_pz->Write();
   h_e_eta->Write();
   h_e_phi->Write();
   h_e_energy->Write();
   h_e_mass->Write();
   h_e_charge->Write();

   h_nchg->Write();
   h_chg_id->Write();
   h_chg_ch->Write();
   h_chg_pt->Write();
   h_chg_px->Write();
   h_chg_py->Write();
   h_chg_pz->Write();
   h_chg_eta->Write();
   h_chg_phi->Write();
   h_chg_energy->Write();
   h_chg_mass->Write();

   h_nproton->Write();
   h_proton_px->Write();
   h_proton_py->Write();
   h_proton_pz->Write();
   h_proton_pt->Write();
   h_proton_energy->Write();
   h_proton_xi->Write();
   h_proton_t->Write();
   h_proton_weight->Write();
   h_proton_acc->Write();

   h_njet->Write();
   h_jetpt->Write();

   output->Close();
   delete jetDef;

   return 0;

} // End of main

