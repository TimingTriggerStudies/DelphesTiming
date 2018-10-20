make -j8 HAS_PYTHIA8=true DelphesPythia8
rm pythia_delphes.root ; ./DelphesPythia8 cards/TimingCards/CMS_PhaseII_MTD_v0.tcl pp2HardQCD.pythia pythia_delphes.root | & tee log
