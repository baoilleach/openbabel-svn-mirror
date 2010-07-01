/**********************************************************************
confsearch.cpp - Conformer Searching Routines (see also forcefield.cpp)

Copyright (C) 2010 Noel O'Boyle <baoilleach@gmail.com>

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

***********************************************************************/

#include <openbabel/babelconfig.h>

#include <openbabel/forcefield.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

#include <float.h> // For DBL_MAX

namespace OpenBabel
{
  int OBForceField::FastRotorSearch(unsigned int finalSteps, unsigned int intermediateSteps)
  {
    int origLogLevel = _loglvl;

//    // -----------------------------------------------------------START
//    // We are going to use the GTD to find the most central rotors
//    std::vector<int> gtd;
//    _mol.GetGTDVector(gtd);
//
//    OBRotorList rl;
//    OBBitVec fixed = _constraints.GetFixedBitVec();
//    rl.SetFixAtoms(fixed);
//    rl.Setup(_mol);
//
//    OBRotorIterator ri;
//    OBRotor *rotor;
//    OBBond *bond;
//
//    std::vector <OBRotor *> rotors;
//    std::vector <unsigned int> gtds, idxs;
//    typedef std::pair<unsigned int, int> uint_rot;
//    std::vector <uint_rot> gtd_rotor;
//    unsigned int i;
//    for (rotor = rl.BeginRotor(ri), i = 0; rotor; rotor = rl.NextRotor(ri), ++i) {
//      bond = rotor->GetBond();
//      unsigned int totGTD = gtd.at( bond->GetBeginAtomIdx() - 1 ) + 
//                            gtd.at( bond->GetEndAtomIdx()   - 1 );
//      gtd_rotor.push_back(uint_rot(totGTD, i));
//#ifdef _DEBUG
//      std::cout << "Rotor " << i << " with bond " << bond->GetIdx() << " has GTD of " << totGTD << std::endl;
//#endif
//    }
//    std::sort(gtd_rotor.begin(), gtd_rotor.end());
//    // Now we have a list of rotors ordered by smallest GTD first
//    // ----------------------------------------------------------END


    // Remove all conformers (e.g. from previous conformer generators) except for current conformer
    double *initialCoord = new double [_mol.NumAtoms() * 3]; // initial state
    double *store_initial = new double [_mol.NumAtoms() * 3]; // store the initial state
    memcpy((char*)initialCoord,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    memcpy((char*)store_initial,(char*)_mol.GetCoordinates(),sizeof(double)*3*_mol.NumAtoms());
    std::vector<double *> newConfs(1, initialCoord);
    _mol.SetConformers(newConfs);

    _energies.clear(); // Wipe any energies from previous conformer generators

    OBRotorList rl;
    OBBitVec fixed = _constraints.GetFixedBitVec();
    rl.SetFixAtoms(fixed);
    rl.Setup(_mol);

    OBRotorIterator ri;
    OBRotamerList rotamerlist;
    rotamerlist.SetBaseCoordinateSets(_mol);
    rotamerlist.Setup(_mol, rl);

    OBRotor *rotor;

    // Start with all of the rotors in their 0 position
    // (perhaps instead I should set them randomly?)
    std::vector<int> rotorKey(rl.Size() + 1, 0);

    unsigned int i, j, minj;
    double currentE, minE;

    // This function relies on the fact that Rotors are ordered from the most
    // central to the most peripheral (due to CompareRotors in rotor.cpp)
    for (rotor = rl.BeginRotor(ri), i = 1; rotor; rotor = rl.NextRotor(ri), ++i) {
      minE = DBL_MAX;  
      
      for (j = 0; j < rotor->GetResolution().size(); j++) { // For each rotor position
        _mol.SetCoordinates(store_initial);
        rotorKey[i] = j;
        rotamerlist.SetCurrentCoordinates(_mol, rotorKey);
        SetupPointers();
        
        if (intermediateSteps) {
          _loglvl = OBFF_LOGLVL_NONE;
          SteepestDescent(intermediateSteps); // energy minimization for conformer
          _loglvl = origLogLevel;
        }
        currentE = Energy(false);
        if (currentE < minE) {
          minE = currentE;
          minj = j;
        }
      } // Finished testing all positions of this rotor
      rotorKey[i] = minj;
#ifdef _DEBUG
      std::cout << "Energy now is " << minE << std::endl;
#endif
    }

    if (finalSteps) {
      _loglvl = OBFF_LOGLVL_NONE;
      SteepestDescent(intermediateSteps); // energy minimization for conformer
      _loglvl = origLogLevel;
    }

    return true;
  }
} // end of namespace OpenBabel

//! \file confsearch.cpp
//! \brief Conformer searching routines
