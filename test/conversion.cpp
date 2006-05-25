/**********************************************************************
conversion.cpp - Unit tests for Open Babel OBConversion class

Copyright (C) 2005-2006 Geoffrey R. Hutchison
 
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

#include "babelconfig.h"
#include "mol.h"
#include "obconversion.h"

#include <stdio.h>
#include <iostream>

using namespace std;
using namespace OpenBabel;

int main(int argc,char *argv[])
{
  // turn off slow sync with C-style output (we don't use it anyway).
  std::ios::sync_with_stdio(false);

  if (argc != 1)
    {
      cout << "Usage: conversion" << endl;
      cout << " Unit tests for OBConversion " << endl;
      return(-1);
    }

  cout << "# Unit tests for OBConversion \n";

  // the number of tests for "prove"
  cout << "1..9\n";

  cout << "ok 1\n"; // for loading tests

  OBMol obMol;
  OBConversion obConversion;
  obConversion.SetInAndOutFormats("smi", "mdl");
  cout << "ok 2\n";

  obConversion.ReadString(&obMol, "C1=CC=CS1");
  cout << "ok 3\n";

  if (obMol.NumAtoms() == 5) {
    cout << "ok 4\n";
  } else {
    cout << "not ok 4\n";
  }

  obMol.AddHydrogens();
  if (obMol.NumAtoms() == 9) {
    cout << "ok 5\n";
  } else {
    cout << "not ok 5\n";
  }

  if ( (obConversion.WriteString(&obMol)).length() > 0)
    cout << "ok 6\n";
  else
    cout << "not ok 6\n";

  // PR#1474265
  obConversion.WriteFile(&obMol, "test.mdl");
  ifstream ifs("test.mdl");
  if (ifs.good())
    cout << "ok 7\n";
  else
    cout << "not ok 7\n";

  // PR#143577
  obConversion.SetInFormat("mdl");
  obConversion.ReadFile(&obMol, "test.mdl");
  if ( remove("test.mdl") != -1)
    cout << "ok 8\n";
  else
    cout << "not ok 8\n";
  
  // gzip input
  // gzip output

  // multi-molecule reading
  // PR#1465586
  // aromatics.smi
  // attype.00.smi

  //ReadFile()
  //Read()
  //WriteString()
  // GetOutputIndex()
  // IsLast

  //ReadString()
  //IsFirstInput
  //Read()

  // splitting
  
  // splitting using gzip-input
  // PR#1357705
  
  // size 0 input
  // PR#1250900
  
  // RegisterFormat
  // FindFormat
  // FormatFromExt
  // FormatFromMIME
  // GetNextFormat
  // GetDefaultFormat

  // BatchFileName
  // IncrementedFileName

  // option handling
  // AddOption
  // IsOption
  // RemoveOption
  // IsOption

  // SetOptions
  // IsOption

  // RegisterOptionParam
  // GetOptionParams

  // GetInStream
  // GetOutStream
  // SetInStream
  // SetOutStream

  // nasty tests
  obConversion.ReadString(&obMol, "");
  obConversion.Read(&obMol);
  cout << "ok 9\n";

  return(0);
}


