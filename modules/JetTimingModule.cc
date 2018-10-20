/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class JetTimingModule
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/JetTimingModule.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

JetTimingModule::JetTimingModule() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

JetTimingModule::~JetTimingModule()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void JetTimingModule::Init()
{
  // read parameters

  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  fVertexArray = ImportArray(GetString("VertexArray", "HighMassVertexRecover/vertices"));
  fItVertexArray = fVertexArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void JetTimingModule::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void JetTimingModule::Process()
{
  const Double_t c_light = 2.99792458E8;
  const double mm_To_m = 1.0e-03;
  const double s_To_ns = 1.0e9;

  // Find the primary vertex
  Candidate *vertex; 
  fItVertexArray->Reset();
  double tmpSumPt2 = 0;
  double PV_X = -999;
  double PV_Y = -999;
  double PV_Z = -999;
  while((vertex = static_cast<Candidate*>(fItVertexArray->Next()))) {
    if (vertex->SumPT2 > tmpSumPt2) {
      tmpSumPt2 = vertex->SumPT2;
      PV_X = vertex->Position.X();
      PV_Y = vertex->Position.Y();
      PV_Z = vertex->Position.Z();
    }
  }

  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  double EnergyWeightedPhotonTime = 0;
  double PhotonTimeEnergySum = 0;
  double ChargedParticleTimeSum = 0;
  double ChargedParticleCount = 0;

  // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    Candidate *constituent = 0;
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    TIter itConstituents(candidate->GetCandidates());

    // cout << "Jet : " << candidateMomentum.Pt() << " " << candidateMomentum.Eta() << " : " << candidateMomentum.Phi() << "\n";
    while((constituent = static_cast<Candidate*>(itConstituents.Next()))) {

      if (! ( constituent->PID == 22 || 
	      (abs(constituent->Charge) == 1 && constituent->Momentum.Pt() > 2.0)
	      )) continue;
      
      // cout << "Consituent : " << constituent->PID << " " 
      // 	   << constituent->Momentum.Pt() << " " 
      // 	   << constituent->Momentum.Eta() << " "
      // 	   << constituent->Momentum.Phi() << " " 
      	   // << " | " 
      	   // << constituent->L << " | "
	// ;
      
      double ConstituentTime = (constituent->Position.T()*1.0e-3/c_light) * s_To_ns;
      double StraightLineDistance = mm_To_m * sqrt(pow(PV_X - constituent->Position.X(),2) + pow(PV_Y - constituent->Position.Y(),2) + pow(PV_Z - constituent->Position.Z(),2) );
      double StraightLineTOF =  (StraightLineDistance / c_light)* s_To_ns;
      double TOFCorrectedTime = ConstituentTime - StraightLineTOF;
      if (constituent->PID != 22) TOFCorrectedTime = ConstituentTime - (constituent->L*1.0e-3/c_light)*1.0e9;

      // cout << ConstituentTime  << " ---> "
      // 	   << TOFCorrectedTime << " "
      // 	   << " | "
      // 	   << "\n"
      	   // << constituent->Position.X() << " " << constituent->Position.Y() << " " << constituent->Position.Z() << " " 
      	   // << " --> " << PV_X << " " << PV_Y << " " << PV_Z << " | " 
      	   // << StraightLineDistance << " " << StraightLineTOF << " "
      	   // << " --> " << TOFCorrectedTime << " "
      	   // << "\n"
      //;   

   
      // cout << "constituent initial position : " << constituent->InitialPosition.X() << " "
      // 	   << constituent->InitialPosition.Y() << " "
      // 	   << constituent->InitialPosition.Z() << " "
      // 	   << constituent->InitialPosition.T() << " "
      // 	   << "\n";
      // cout << "constituent final position : " << constituent->Position.X() << " "
      // 	   << constituent->Position.Y() << " "
      // 	   << constituent->Position.Z() << " "
      // 	   << constituent->Position.T() << " "
      // 	   << "\n";

      if (constituent->PID == 22) {
	EnergyWeightedPhotonTime += constituent->Momentum.E() * TOFCorrectedTime;
	PhotonTimeEnergySum += constituent->Momentum.E();
      }

      if (abs(constituent->Charge) == 1 && constituent->Momentum.Pt() > 2.0) {
	ChargedParticleTimeSum += TOFCorrectedTime;
	ChargedParticleCount += 1.0;
      }     
    }
    double PhotonTime = EnergyWeightedPhotonTime / PhotonTimeEnergySum;
    double ChargedParticleTime = ChargedParticleTimeSum / ChargedParticleCount;

    // cout << "PhotonTime = " << PhotonTime << "\n";
    // cout << "ChargedParticleTime = " << ChargedParticleTime << "\n";
    // cout << "\n\n";

    ((Jet*)candidate)->PhotonTime = PhotonTime;
    ((Jet*)candidate)->ChargedParticleTime = ChargedParticleTime;

    // apply an efficency formula   
    fOutputArray->Add(candidate);
   
  }
}

//------------------------------------------------------------------------------
