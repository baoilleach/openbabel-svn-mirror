/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "obutil.h"
#include "parsmart.h"
#include "typer.h"
#include "rotor.h"
#include "binary.h"
#include "commandline.h"
#include "version.h"
#include "data.h"

#include <stdio.h>

#ifdef __sgi
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
#endif

using namespace std;
using namespace OpenBabel;

void usage();

// There isn't a great way to do this -- we need to save argv[0] for usage()
static char *program_name;

int main(int argc,char *argv[])
{
  io_type inFileType = UNDEFINED, outFileType = UNDEFINED;
  bool gotInType = false, gotOutType = false, removeHydrogens = false;
  bool addHydrogens = false;
  int arg, inFileArg, outFileArg;
  char *ext;
  OBFileFormat fileFormat;

  // Parse commandline
  program_name = argv[0];
  inFileArg = 0;
  outFileArg = 0;
  for (arg = 1; arg <= argc; arg++)
    {
      if (argv[arg])
	{
	  if (argv[arg][0] == '-')
	    {
	      switch (argv[arg][1])
		{
		case 'v':
		  {
		    cout << "Open Babel " << BABEL_VERSION << " -- " 
			 << __DATE__ << " -- " << __TIME__ << endl;
		    exit(0);
		  }
		case 'd':
		  removeHydrogens = true;
		  break;
		case 'h':
		  addHydrogens = true;
		  break;

		case 'i':
		  gotInType = true;
		  ext = argv[arg] + 2;
		  if (extab.CanReadExtension(ext))
		    inFileType = extab.FilenameToType(ext);
		  else
		    {
		      cerr << program_name << ": Cannot read input format!" << endl;
		      usage();
		    }
		  break;

		case 'o':
		  gotOutType = true;
		  ext = argv[arg] + 2;
		  if (extab.CanWriteExtension(ext))
		    outFileType = extab.FilenameToType(ext);
		  else
		    {
		      cerr << program_name << ": Cannot write output format!" << endl;
		      usage();
		    }

		  break;

		default:
		  usage();
		  break;
		}
	    }
	  else if (inFileArg == 0)
	    inFileArg = arg;
	  else
	    outFileArg = arg;
	}
    }
  if (inFileArg == 0 || outFileArg == 0)
    usage();

  if (!gotInType)
    {
      if (extab.CanReadExtension(argv[inFileArg]))
	inFileType = extab.FilenameToType(argv[inFileArg]);
      else
	{
	  cerr << program_name << ": Cannot read input format!" << endl;
	  usage();
	}
    }
  if (!gotOutType)
    {
      if (extab.CanWriteExtension(argv[outFileArg]))
	outFileType = extab.FilenameToType(argv[outFileArg]);
      else
	{
	  cerr << program_name << ": Cannot write output format!" << endl;
	  usage();
	}
    }

  ifstream inFileStream(argv[inFileArg]);
  ofstream outFileStream(argv[outFileArg]);

  if (!inFileStream)
    {
      cerr << program_name << ": Cannot read input file!" << endl;
      exit (-1);
    }
  if (!outFileStream)
    {
      cerr << program_name << ": Cannot write to output file!" << endl;
      exit (-1);
    }

  // Finally, we can do some work!
  OBMol mol(inFileType, outFileType);

  fileFormat.ReadMolecule(inFileStream, mol, argv[inFileArg]);
  if (removeHydrogens)
    mol.DeleteHydrogens();
  if (addHydrogens)
    mol.AddHydrogens(false, false);
  fileFormat.WriteMolecule(outFileStream,mol);

  return(0);
}

void usage()
{
unsigned int i;

  cout << "Open Babel " << BABEL_VERSION << " -- " << __DATE__ << " -- "
       << __TIME__ << endl;
  cout << "Usage is : " << endl << program_name 
       << " [-i<input-type>] <name> [-o<output-type>] <name>" << endl;
  cout << endl << "Currently supported input types" << endl;
  for (i = 0; i < extab.Count(); i++)
    if (extab.IsReadable(i))
      cout << "\t" << extab.GetExtension(i) << " -- " 
		<< extab.GetDescription(i) << " file" << endl;
  cout << endl << "Currently supported output types" << endl;
  for (i = 0; i < extab.Count(); i++)
    if (extab.IsWritable(i))
      cout << "\t" << extab.GetExtension(i) << " -- " 
		<< extab.GetDescription(i) << " file" << endl;
  cout << "Additional options : " << endl;
  cout << " -d Delete Hydrogens " << endl;
  cout << " -h Add Hydrogens " << endl;

  exit(0);
}
