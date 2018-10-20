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

#ifndef JetTimingModule_h
#define JetTimingModule_h

/** \class JetTimingModule
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author S. Xie - California Institute of Technology
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class JetTimingModule: public DelphesModule
{
public:

  JetTimingModule();
  ~JetTimingModule();

  void Init();
  void Process();
  void Finish();

private:

  Int_t fIntParam;
  Double_t fDoubleParam;
  
  std::deque <Double_t> fArrayParam;
  
  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!
  TIterator *fItVertexArray; //!

  const TObjArray *fInputArray; //!
  const TObjArray *fVertexArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(JetTimingModule, 1)
};

#endif
