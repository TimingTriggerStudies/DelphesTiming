#!/bin/tcsh
setenv PYTHIA8 "/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/pythia8/230-omkpbe4/"
setenv PYTHIA8DATA "$PYTHIA8/share/Pythia8/xmldoc"
setenv LD_LIBRARY_PATH "$PYTHIA8/lib:$LD_LIBRARY_PATH"
setenv PYTHIA8_INCLUDE_DIR "$PYTHIA8/include"
setenv PYTHIA8_LIBRARY "$PYTHIA8/lib"
