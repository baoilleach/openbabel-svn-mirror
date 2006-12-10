/**********************************************************************
invalid-smiles.cpp - Test SMILES pattern parsing (rejecting invalid patterns)

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

Some portions Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include <fstream>

namespace OpenBabel
{
bool SafeOpen(std::ifstream &fs, const char *filename);
bool SafeOpen(std::ofstream &fs, const char *filename);
}

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string invalid_file = testdatadir + "invalid-smiles.txt";
  string random1_file = testdatadir + "random";
  string random2_file = testdatadir + "random2";
  string random3_file = testdatadir + "random3";
#else
  string invalid_file = "files/invalid-smiles.txt";
  string random1_file = "files/random";
  string random2_file = "files/random2";
  string random3_file = "files/random3";
#endif

void GenerateFormalChargeReference();

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: invalid-smiles" << endl;
      cout << "   Tests Open Babel SMILES parsing - rejecting invalid patterns."
           << endl;
      return 0;
    }

  cout << "# Testing invalid SMILES parsing..." << endl;

  // make sure to kill off all error reporting
  obErrorLog.StopLogging();

  std::ifstream mifs;
  if (!SafeOpen(mifs, invalid_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << invalid_file << endl;
      return -1; // test failed
    }

  unsigned int currentTest = 0;
  OBMol mol;
  OBConversion conv(&mifs, &cout);
  if (! conv.SetInAndOutFormats("SMI","SMI"))
    {
      cout << "Bail out! SMILES format is not loaded" << endl;
      return -1;
    }

  while(mifs.good())
    {
      mol.Clear();
      if (conv.Read(&mol))
        cout << "not ok " << ++currentTest << " # invalid SMILES was parsed "
             << " but should have been rejected" << endl;
      else
        cout << "ok " << ++currentTest << endl;
    }

  mifs.close();
  mifs.clear();

  // random file#1
  if (!SafeOpen(mifs, random1_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << random1_file << endl;
      return -1; // test failed
    }

  mol.Clear();
  if (conv.Read(&mol, &mifs))
    cout << "not ok " << ++currentTest << " # random data was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << "\n";
  
  mifs.close();
  mifs.clear();

  // random2
  if (!SafeOpen(mifs, random2_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << random2_file << endl;
      return -1; // test failed
    }

  mol.Clear();
  if (conv.Read(&mol, &mifs))
    cout << "not ok " << ++currentTest << " # random data #2 was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << "\n";

  mifs.close();
  mifs.clear();

  // random3
  if (!SafeOpen(mifs, random3_file.c_str()))
    {
      cout << "Bail out! Cannot read file " << random3_file << endl;
      return -1; // test failed
    }

  mol.Clear();
  if (conv.Read(&mol, &mifs))
    cout << "not ok " << ++currentTest << " # random data #3 was parsed "
         << " but should have been rejected\n";
  else
    cout << "ok " << ++currentTest << "\n";

  // return number of tests run
  cout << "1.." << currentTest << endl;

  // Passed tests
  return 0;
}
