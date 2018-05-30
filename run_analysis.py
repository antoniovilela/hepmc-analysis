#!/usr/bin/python

import math
from ROOT import TFile, TH1F, TF1, TLorentzVector

def run_analysis( fileName ):
    myfile = TFile( fileName )

    # retrieve the ntuple of interest
    mychain = myfile.Get('T')
    entries = mychain.GetEntriesFast()

    mychain.Print()

    n_events_all = 0
    n_events_selected_jet = 0
    n_events_selected_proton = 0
   
    for jentry in xrange(entries):
	# get the next tree in the chain and verify
	ientry = mychain.LoadTree(jentry)
	if ientry < 0:
	    break

	# copy next entry into memory and verify
	nb = mychain.GetEntry(jentry)
	if nb<=0:
	    continue

        n_events_all += 1
 
        select_jet = False
	n_jet = int(mychain.n_jet)
        jet1_pt = -1.
        jet2_pt = -1.
        jet1_vec = None
        jet2_vec = None
        if n_jet >= 2:
	    jet1_pt = float( mychain.jet_pt[0] )
	    jet2_pt = float( mychain.jet_pt[1] )

	    jet1_vec = TLorentzVector(
		float( mychain.jet_px[0] ),
		float( mychain.jet_py[0] ),
		float( mychain.jet_pz[0] ),
		float( mychain.jet_energy[0] ),
	    )
	    jet2_vec = TLorentzVector(
		float( mychain.jet_px[1] ),
		float( mychain.jet_py[1] ),
		float( mychain.jet_pz[1] ),
		float( mychain.jet_energy[1] ),
	    )
            select_jet = True
     
        if not select_jet:
            continue

        n_events_selected_jet += 1

	select_proton = False

	n_proton = int(mychain.n_proton)
	mass_from_protons = -1.
	if n_proton == 2:
	    proton1_xi = float(mychain.proton_xi[0] )
	    proton2_xi = float(mychain.proton_xi[1] )
	    mass_from_protons = 13000.*math.sqrt(proton1_xi*proton2_xi)
	    select_proton = True

        if not select_proton:
            continue

        n_events_selected_proton += 1

	twojet_vec = jet1_vec + jet2_vec;    

	print "Event %d" % (jentry + 1)
	print "Number of jets: %d, Jet pT = %f, %f, Two-jet mass: %f" % (n_jet,jet1_pt,jet2_pt,twojet_vec.M())
	print "Mass (protons) = %f" % mass_from_protons

    # Summary
    print "\nSummary:"
    print "Number of events: %d" % n_events_all
    print "Number of events (Jet selection): %d" % n_events_selected_jet
    print "Number of events (Proton selection): %d" % n_events_selected_proton

if __name__ == '__main__':

    fileName = 'analysis-Dijets_DPEDijets_13TeV_PTMIN_50-v1.root'
    run_analysis( fileName )

    print "\nFinished."

