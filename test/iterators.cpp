 /**********************************************************************
 iterators.cpp - tests for iterators

 This file is part of the Open Babel project.
 For more information, see <http://openbabel.sourceforge.net/>

 Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
 Some portions Copyright (C) 2001-2005 Geoffrey R. Hutchison

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

#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>

namespace OpenBabel
{
  bool SafeOpen(std::ifstream &fs, const char *filename);
  bool SafeOpen(std::ofstream &fs, const char *filename);
}

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
   string smilestypes_file = "attype.00.smi";
#endif

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: iterators\n";
      cout << "   Tests Open Babel iterators." << endl;
      return 0;
    }
  
  cout << endl << "# Testing iterators...  \n";
  
  std::ifstream mifs;
  if (!SafeOpen(mifs, smilestypes_file.c_str()))
    {
      cout << "Bail out! Cannot read test data " << smilestypes_file << endl;
      return -1; // test failed
    }

  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  OBMol mol;

  unsigned int currentTest = 0;
  unsigned int counter = 0;

  // run through atom and bond iterators
  while(mifs.peek() != EOF && mifs.good())
    {
      mol.Clear();
      conv.Read(&mol);
      if (mol.Empty())
        continue;

      counter = 0;
      FOR_ATOMS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # atom iterator failed: expected " << mol.NumAtoms()
             << " but found " << counter << "\n";
      else
        cout << "ok " << ++currentTest << "\n";

      counter = 0;
      FOR_DFS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # depth-first atom iterator failed: expected " 
             << mol.NumAtoms() << " but found " << counter << "\n";
      else
        cout << "ok " << ++currentTest << "\n";

      counter = 0;
      FOR_BFS_OF_MOL(atom, mol)
        ++counter;
      if (counter != mol.NumAtoms())
        cout << "not ok " << ++currentTest 
             << " # breadth-first atom iterator failed: expected " 
             << mol.NumAtoms() << " but found " << counter << "\n";
      else
        cout << "ok " << ++currentTest << "\n";

      counter = 0;
      FOR_BONDS_OF_MOL(bond, mol)
        ++counter;
      if (counter != mol.NumBonds())
        cout << "not ok " << ++currentTest 
             << " # bond iterator failed: expected " << mol.NumBonds()
             << " but found " << counter << "\n";
      else
        cout << "ok " << ++currentTest << "\n";

    }

  // output the number of tests run
  cout << "1.." << currentTest << endl;

  // Passed Test
  return 0;
}
