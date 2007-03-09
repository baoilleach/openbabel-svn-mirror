/**********************************************************************
obfragment - Generate coordinate database of ring fragments

Copyright (C) 2007 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  OBConversion conv;
  OBFormat *inFormat, *canFormat;
  OBMol mol;
  ifstream ifs;
  vector<OBMol> fragments;
  map<string, int> index; // index of cansmi
  string currentCAN;
  unsigned int size;
  OBAtom *atom;
  OBBond *bond;
  bool nonRingAtoms, nonRingBonds;
  char buffer[BUFF_SIZE];

  canFormat = conv.FindFormat("can");
  conv.SetOutFormat(canFormat);

  if (argc < 2)
    {
      cout << "Usage: obfragment <file>" << endl;
      return(-1);
    }

  for (int i = 1; i < argc; i++) {
    cerr << " Reading file " << argv[i] << endl;

    inFormat = conv.FormatFromExt(argv[i]);
    if(inFormat==NULL || !conv.SetInFormat(inFormat))
      {
        cerr << " Cannot read file format for " << argv[i] << endl;
        continue; // try next file
      }
    
    ifs.open(argv[i]);
    
    if (!ifs)
      {
        cout << "Cannot read input file: " << argv[i] << endl;
        continue;
      }
    
    
    while(ifs.peek() != EOF && ifs.good())
      {
        conv.Read(&mol, &ifs);
        mol.DeleteHydrogens(); // remove these before we do anything else
        
        do {
          nonRingAtoms = false;
          size = mol.NumAtoms();
          for (unsigned int i = 1; i <= size; ++i)
            {
              atom = mol.GetAtom(i);
              if (!atom->IsInRing()) {
                mol.DeleteAtom(atom);
                nonRingAtoms = true;
                break; // don't know how many atoms there are
              } else {
                atom->SetAtomicNum(6);
                atom->SetFormalCharge(0);
              }
            }
        } while (nonRingAtoms);
        
        if (mol.NumAtoms() < 3)
          continue;
        
        if (mol.NumBonds() == 0)
          continue;
        
        do {
          nonRingBonds = false;
          size = mol.NumBonds();
          for (unsigned int i = 0; i < size; ++i)
            {
              bond = mol.GetBond(i);
              if (!bond->IsInRing()) {
                mol.DeleteBond(bond);
                nonRingBonds = true;
                break; // don't know how many bonds there are
              }
            }        
        } while (nonRingBonds);

        fragments = mol.Separate();
        for (unsigned int i = 0; i < fragments.size(); ++i)
          {
            if (fragments[i].NumAtoms() < 3) // too small to care
              continue;

            currentCAN = conv.WriteString(&fragments[i], true);
            if (index.find(currentCAN) != index.end()) // already got this
              continue;

            index[currentCAN] = 1; // don't ever write this ring fragment again

            cout << "INDEX " << currentCAN << '\n'; // endl causes a flush

            FOR_ATOMS_OF_MOL (atom, fragments[i]) {
              snprintf(buffer, BUFF_SIZE, "%9.3f %9.3f %9.3f\n",
                       atom->x(), atom->y(), atom->z());
              cout << buffer;
            }
            

          }
        fragments.clear();
	
      } // while reading molecules (in this file)
    ifs.close();
    ifs.clear();
  } // while reading files
    
  return(0);
}
