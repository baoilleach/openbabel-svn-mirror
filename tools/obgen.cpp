/**********************************************************************
obgen.cpp - test program for SMILES 3D coordinate generation
          - using systematic rotor search

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison
 
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
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//! \brief  Generate rough 3D coordinates for SMILES (or other 0D files).
//
int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  string basename, filename = "", option, option2, ff = "";

  if (argc < 2) {
    cout << "Usage: obgen <filename> [options]" << endl;
    cout << endl;
    cout << "options:      description:" << endl;
    cout << endl;
    cout << "  -ff         select a forcefield" << endl;
    cout << endl;
    FOR_EACH(OBForceField, iter) {
      cout << "              " << iter.ID() << " - " << iter.Description() << endl;
    }
    exit(-1);
  } else {
    basename = filename = argv[1];
    size_t extPos = filename.rfind('.');

    if (extPos!= string::npos) {
      basename = filename.substr(0, extPos);
    }

    for (int i = 2; i < argc; i++) {
      option = argv[i];
      if ((option == "-ff") && (argc > (i+1)))
        ff = argv[i+1];
    }
  }

  // Find Input filetype
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(filename.c_str());
  OBFormat *format_out = conv.FindFormat("pdb");
    
  if (!format_in || !format_out || !conv.SetInAndOutFormats(format_in, format_out)) {
    cerr << program_name << ": cannot read input/output format!" << endl;
    exit (-1);
  }

  ifstream ifs;
  ofstream ofs;

  // Read the file
  ifs.open(filename.c_str());
  if (!ifs) {
    cerr << program_name << ": cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol;

  for (c=1;;c++) {
      mol.Clear();
      if (!conv.Read(&mol, &ifs))
        break;
      if (mol.Empty())
        break;

      OBForceField* pFF = OBForceField::FindForceField(ff);
      if (!pFF) {
        cerr << program_name << ": could not find forcefield '" << ff << "'." <<endl;
        exit (-1);
      }
 
      mol.AddHydrogens(false, true); // hydrogens must be added before Setup(mol) is called
      
      pFF->SetLogFile(&cerr);
      pFF->SetLogLevel(OBFF_LOGLVL_LOW);
      
      if (!pFF->Setup(mol)) {
        cerr << program_name << ": could not setup force field." << endl;
        exit (-1);
      }
     
      pFF->GenerateCoordinates();
      pFF->SystematicRotorSearch();

      pFF->UpdateCoordinates(mol);
      //pFF->ValidateGradients();
      pFF->SetLogLevel(OBFF_LOGLVL_HIGH);
      pFF->Energy();
      

      //char FileOut[32];
      //sprintf(FileOut, "%s_obgen.pdb", basename.c_str());
      //ofs.open(FileOut);
      //conv.Write(&mol, &ofs);
      //ofs.close();
      conv.Write(&mol, &cout);
  } // end for loop

  return(1);
}
