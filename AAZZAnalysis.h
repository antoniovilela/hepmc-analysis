//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 28 19:39:57 2015 by ROOT version 5.34/28
// from TTree T/Tree
// found on file: AAZZ_13TeV.root
//////////////////////////////////////////////////////////

#ifndef AAZZAnalysis_h
#define AAZZAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AAZZAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           n_mu;
   Double_t        mu_pt[5];   //[n_mu]
   Double_t        mu_px[5];   //[n_mu]
   Double_t        mu_py[5];   //[n_mu]
   Double_t        mu_pz[5];   //[n_mu]
   Double_t        mu_eta[5];   //[n_mu]
   Double_t        mu_phi[5];   //[n_mu]
   Double_t        mu_energy[5];   //[n_mu]
   Double_t        mu_mass[5];   //[n_mu]
   Double_t        mu_charge[5];   //[n_mu]
   Int_t           n_e;
   Double_t        e_pt[6];   //[n_e]
   Double_t        e_px[6];   //[n_e]
   Double_t        e_py[6];   //[n_e]
   Double_t        e_pz[6];   //[n_e]
   Double_t        e_eta[6];   //[n_e]
   Double_t        e_phi[6];   //[n_e]
   Double_t        e_energy[6];   //[n_e]
   Double_t        e_mass[6];   //[n_e]
   Double_t        e_charge[6];   //[n_e]
   Int_t           n_chg;
   Int_t           chg_id[88];   //[n_chg]
   Double_t        chg_ch[88];   //[n_chg]
   Double_t        chg_px[88];   //[n_chg]
   Double_t        chg_py[88];   //[n_chg]
   Double_t        chg_pz[88];   //[n_chg]
   Double_t        chg_pt[88];   //[n_chg]
   Double_t        chg_eta[88];   //[n_chg]
   Double_t        chg_phi[88];   //[n_chg]
   Double_t        chg_energy[88];   //[n_chg]
   Double_t        chg_mass[88];   //[n_chg]
   Int_t           n_proton;
   Double_t        proton_px[9];   //[n_proton]
   Double_t        proton_py[9];   //[n_proton]
   Double_t        proton_pz[9];   //[n_proton]
   Double_t        proton_pt[9];   //[n_proton]
   Double_t        proton_energy[9];   //[n_proton]
   Int_t           n_jet;
   Double_t        jet_px[1];   //[n_jet]
   Double_t        jet_py[1];   //[n_jet]
   Double_t        jet_pz[1];   //[n_jet]
   Double_t        jet_pt[1];   //[n_jet]
   Double_t        jet_energy[1];   //[n_jet]
   Int_t           n_Z;
   Double_t        Z_pt[2];   //[n_Z]
   Double_t        Z_px[2];   //[n_Z]
   Double_t        Z_py[2];   //[n_Z]
   Double_t        Z_pz[2];   //[n_Z]
   Double_t        Z_eta[2];   //[n_Z]
   Double_t        Z_phi[2];   //[n_Z]
   Double_t        Z_energy[2];   //[n_Z]
   Double_t        Z_mass[2];   //[n_Z]
   Int_t           Z_decay_first_pid[2];   //[n_Z]
   Double_t        Z_decay_first_pt[2];   //[n_Z]
   Double_t        Z_decay_first_px[2];   //[n_Z]
   Double_t        Z_decay_first_py[2];   //[n_Z]
   Double_t        Z_decay_first_pz[2];   //[n_Z]
   Double_t        Z_decay_first_eta[2];   //[n_Z]
   Double_t        Z_decay_first_phi[2];   //[n_Z]
   Double_t        Z_decay_first_energy[2];   //[n_Z]
   Double_t        Z_decay_first_mass[2];   //[n_Z]
   Double_t        Z_decay_first_charge[2];   //[n_Z]
   Int_t           Z_decay_second_pid[2];   //[n_Z]
   Double_t        Z_decay_second_pt[2];   //[n_Z]
   Double_t        Z_decay_second_px[2];   //[n_Z]
   Double_t        Z_decay_second_py[2];   //[n_Z]
   Double_t        Z_decay_second_pz[2];   //[n_Z]
   Double_t        Z_decay_second_eta[2];   //[n_Z]
   Double_t        Z_decay_second_phi[2];   //[n_Z]
   Double_t        Z_decay_second_energy[2];   //[n_Z]
   Double_t        Z_decay_second_mass[2];   //[n_Z]
   Double_t        Z_decay_second_charge[2];   //[n_Z]

   // List of branches
   TBranch        *b_n_mu;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_px;   //!
   TBranch        *b_mu_py;   //!
   TBranch        *b_mu_pz;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_energy;   //!
   TBranch        *b_mu_mass;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_n_e;   //!
   TBranch        *b_e_pt;   //!
   TBranch        *b_e_px;   //!
   TBranch        *b_e_py;   //!
   TBranch        *b_e_pz;   //!
   TBranch        *b_e_eta;   //!
   TBranch        *b_e_phi;   //!
   TBranch        *b_e_energy;   //!
   TBranch        *b_e_mass;   //!
   TBranch        *b_e_charge;   //!
   TBranch        *b_n_chg;   //!
   TBranch        *b_chg_id;   //!
   TBranch        *b_chg_ch;   //!
   TBranch        *b_chg_px;   //!
   TBranch        *b_chg_py;   //!
   TBranch        *b_chg_pz;   //!
   TBranch        *b_chg_pt;   //!
   TBranch        *b_chg_eta;   //!
   TBranch        *b_chg_phi;   //!
   TBranch        *b_chg_energy;   //!
   TBranch        *b_chg_mass;   //!
   TBranch        *b_n_proton;   //!
   TBranch        *b_proton_px;   //!
   TBranch        *b_proton_py;   //!
   TBranch        *b_proton_pz;   //!
   TBranch        *b_proton_pt;   //!
   TBranch        *b_proton_energy;   //!
   TBranch        *b_n_jet;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_n_Z;   //!
   TBranch        *b_Z_pt;   //!
   TBranch        *b_Z_px;   //!
   TBranch        *b_Z_py;   //!
   TBranch        *b_Z_pz;   //!
   TBranch        *b_Z_eta;   //!
   TBranch        *b_Z_phi;   //!
   TBranch        *b_Z_energy;   //!
   TBranch        *b_Z_mass;   //!
   TBranch        *b_Z_decay_first_pid;   //!
   TBranch        *b_Z_decay_first_pt;   //!
   TBranch        *b_Z_decay_first_px;   //!
   TBranch        *b_Z_decay_first_py;   //!
   TBranch        *b_Z_decay_first_pz;   //!
   TBranch        *b_Z_decay_first_eta;   //!
   TBranch        *b_Z_decay_first_phi;   //!
   TBranch        *b_Z_decay_first_energy;   //!
   TBranch        *b_Z_decay_first_mass;   //!
   TBranch        *b_Z_decay_first_charge;   //!
   TBranch        *b_Z_decay_second_pid;   //!
   TBranch        *b_Z_decay_second_pt;   //!
   TBranch        *b_Z_decay_second_px;   //!
   TBranch        *b_Z_decay_second_py;   //!
   TBranch        *b_Z_decay_second_pz;   //!
   TBranch        *b_Z_decay_second_eta;   //!
   TBranch        *b_Z_decay_second_phi;   //!
   TBranch        *b_Z_decay_second_energy;   //!
   TBranch        *b_Z_decay_second_mass;   //!
   TBranch        *b_Z_decay_second_charge;   //!

   AAZZAnalysis(TTree *tree=0);
   virtual ~AAZZAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AAZZAnalysis_cxx
AAZZAnalysis::AAZZAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AAZZ_13TeV.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("AAZZ_13TeV.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

AAZZAnalysis::~AAZZAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AAZZAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AAZZAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AAZZAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n_mu", &n_mu, &b_n_mu);
   fChain->SetBranchAddress("mu_pt", mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_px", mu_px, &b_mu_px);
   fChain->SetBranchAddress("mu_py", mu_py, &b_mu_py);
   fChain->SetBranchAddress("mu_pz", mu_pz, &b_mu_pz);
   fChain->SetBranchAddress("mu_eta", mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_energy", mu_energy, &b_mu_energy);
   fChain->SetBranchAddress("mu_mass", mu_mass, &b_mu_mass);
   fChain->SetBranchAddress("mu_charge", mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("n_e", &n_e, &b_n_e);
   fChain->SetBranchAddress("e_pt", e_pt, &b_e_pt);
   fChain->SetBranchAddress("e_px", e_px, &b_e_px);
   fChain->SetBranchAddress("e_py", e_py, &b_e_py);
   fChain->SetBranchAddress("e_pz", e_pz, &b_e_pz);
   fChain->SetBranchAddress("e_eta", e_eta, &b_e_eta);
   fChain->SetBranchAddress("e_phi", e_phi, &b_e_phi);
   fChain->SetBranchAddress("e_energy", e_energy, &b_e_energy);
   fChain->SetBranchAddress("e_mass", e_mass, &b_e_mass);
   fChain->SetBranchAddress("e_charge", e_charge, &b_e_charge);
   fChain->SetBranchAddress("n_chg", &n_chg, &b_n_chg);
   fChain->SetBranchAddress("chg_id", chg_id, &b_chg_id);
   fChain->SetBranchAddress("chg_ch", chg_ch, &b_chg_ch);
   fChain->SetBranchAddress("chg_px", chg_px, &b_chg_px);
   fChain->SetBranchAddress("chg_py", chg_py, &b_chg_py);
   fChain->SetBranchAddress("chg_pz", chg_pz, &b_chg_pz);
   fChain->SetBranchAddress("chg_pt", chg_pt, &b_chg_pt);
   fChain->SetBranchAddress("chg_eta", chg_eta, &b_chg_eta);
   fChain->SetBranchAddress("chg_phi", chg_phi, &b_chg_phi);
   fChain->SetBranchAddress("chg_energy", chg_energy, &b_chg_energy);
   fChain->SetBranchAddress("chg_mass", chg_mass, &b_chg_mass);
   fChain->SetBranchAddress("n_proton", &n_proton, &b_n_proton);
   fChain->SetBranchAddress("proton_px", proton_px, &b_proton_px);
   fChain->SetBranchAddress("proton_py", proton_py, &b_proton_py);
   fChain->SetBranchAddress("proton_pz", proton_pz, &b_proton_pz);
   fChain->SetBranchAddress("proton_pt", proton_pt, &b_proton_pt);
   fChain->SetBranchAddress("proton_energy", proton_energy, &b_proton_energy);
   fChain->SetBranchAddress("n_jet", &n_jet, &b_n_jet);
   fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_energy", jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("n_Z", &n_Z, &b_n_Z);
   fChain->SetBranchAddress("Z_pt", Z_pt, &b_Z_pt);
   fChain->SetBranchAddress("Z_px", Z_px, &b_Z_px);
   fChain->SetBranchAddress("Z_py", Z_py, &b_Z_py);
   fChain->SetBranchAddress("Z_pz", Z_pz, &b_Z_pz);
   fChain->SetBranchAddress("Z_eta", Z_eta, &b_Z_eta);
   fChain->SetBranchAddress("Z_phi", Z_phi, &b_Z_phi);
   fChain->SetBranchAddress("Z_energy", Z_energy, &b_Z_energy);
   fChain->SetBranchAddress("Z_mass", Z_mass, &b_Z_mass);
   fChain->SetBranchAddress("Z_decay_first_pid", Z_decay_first_pid, &b_Z_decay_first_pid);
   fChain->SetBranchAddress("Z_decay_first_pt", Z_decay_first_pt, &b_Z_decay_first_pt);
   fChain->SetBranchAddress("Z_decay_first_px", Z_decay_first_px, &b_Z_decay_first_px);
   fChain->SetBranchAddress("Z_decay_first_py", Z_decay_first_py, &b_Z_decay_first_py);
   fChain->SetBranchAddress("Z_decay_first_pz", Z_decay_first_pz, &b_Z_decay_first_pz);
   fChain->SetBranchAddress("Z_decay_first_eta", Z_decay_first_eta, &b_Z_decay_first_eta);
   fChain->SetBranchAddress("Z_decay_first_phi", Z_decay_first_phi, &b_Z_decay_first_phi);
   fChain->SetBranchAddress("Z_decay_first_energy", Z_decay_first_energy, &b_Z_decay_first_energy);
   fChain->SetBranchAddress("Z_decay_first_mass", Z_decay_first_mass, &b_Z_decay_first_mass);
   fChain->SetBranchAddress("Z_decay_first_charge", Z_decay_first_charge, &b_Z_decay_first_charge);
   fChain->SetBranchAddress("Z_decay_second_pid", Z_decay_second_pid, &b_Z_decay_second_pid);
   fChain->SetBranchAddress("Z_decay_second_pt", Z_decay_second_pt, &b_Z_decay_second_pt);
   fChain->SetBranchAddress("Z_decay_second_px", Z_decay_second_px, &b_Z_decay_second_px);
   fChain->SetBranchAddress("Z_decay_second_py", Z_decay_second_py, &b_Z_decay_second_py);
   fChain->SetBranchAddress("Z_decay_second_pz", Z_decay_second_pz, &b_Z_decay_second_pz);
   fChain->SetBranchAddress("Z_decay_second_eta", Z_decay_second_eta, &b_Z_decay_second_eta);
   fChain->SetBranchAddress("Z_decay_second_phi", Z_decay_second_phi, &b_Z_decay_second_phi);
   fChain->SetBranchAddress("Z_decay_second_energy", Z_decay_second_energy, &b_Z_decay_second_energy);
   fChain->SetBranchAddress("Z_decay_second_mass", Z_decay_second_mass, &b_Z_decay_second_mass);
   fChain->SetBranchAddress("Z_decay_second_charge", Z_decay_second_charge, &b_Z_decay_second_charge);
   Notify();
}

Bool_t AAZZAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AAZZAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AAZZAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AAZZAnalysis_cxx
