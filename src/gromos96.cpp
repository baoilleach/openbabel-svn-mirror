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

/* contributed by Walter Scott (wscott@igc.phys.chem.ethz.ch)

   (Actually the routine was copied from write_xyz and and write_pdb and
   then modified...)

   This is a small routine to write a GROMOS96 formatted 
   "position coordinate block" (POSITION) or a
   "reduced position coordinate block" (POSITIONRED)
   The former has name information (atom and residue names) while
   the latter only has coordinates.
   This version does not support the writing of binary
   GROMOS files.

   NOTE 1: the actual formats used in writing out the coordinates
   do not matter, as GROMOS96 uses free formatted reads.
   Each line may not be longer than 80 characters.

   (Note, however, in the POSITION block, the first 24 (twenty four)
   character on each line are ignored when the line is read in by GROMOS)
   Comments lines, beginning with hash (#) may occur within a block and are
   used as delimiters for easier reading.

   NOTE 2: Many programs specify the units of the coordinates (e.g. Angstrom).
   GROMOS96 does NOT, as all physical constants, from K_B to EPS are 
   NOT hardwired into the code, but specified by the user.
   This allows some (mostly Americans) to use GROMOS96 in KCal and
   Angstrom and the rest of us to use kJoule and nm.
   It also makes it easy to use reduced units.

   We get around this by supplying a routine, wr_sco_gr96, which
   will scale the coordinates by a factor before writing.
   This routine is then called with the factor set to 1.0 in 
   write_gr96A, or to 0.1 in write_gr96N depending on the users choice.
   Thus, we always assume that we have read coordinates in Angstrom.
*/

#include "mol.h"
#include "babelconfig.h"

using namespace std;

namespace OpenBabel
{

bool WriteGromos96(ostream &ofs,OBMol &mol,double fac)
{ 
  char type_name[10];
  char res_name[10],padded_name[10];
  char buffer[BUFF_SIZE];
  int res_num;

  sprintf(buffer,"#GENERATED BY OPEN BABEL %s",BABEL_VERSION);
  ofs << buffer << endl;

  /* GROMOS wants a TITLE block, so let's write one*/
  sprintf(buffer,"TITLE\n%s\nEND",mol.GetTitle());
  ofs << buffer << endl;
  ofs << "POSITION" << endl;

  OBAtom *atom;
  OBResidue *res;
  vector<OBNodeBase*>::iterator i;

  for(atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      if (atom->HasResidue())
	{
	  res = atom->GetResidue();
	  strcpy(res_name,(char*)res->GetName().c_str());
	  strcpy(type_name,(char*)res->GetAtomID(atom).c_str());
	  res_num = res->GetNum();
	}
      else
	{
	  strcpy(type_name,etab.GetSymbol(atom->GetAtomicNum()));
	  strcpy(res_name,"UNK");
	  sprintf(padded_name,"%2s",type_name);
	  strcpy(type_name,padded_name);
	  res_num = 1;
	}
      
      sprintf(buffer,"%5d %5s %5s %6d %15.5f %15.5f %15.5f",
	      res_num,res_name,type_name,atom->GetIdx(),
	      atom->x()*fac,atom->y()*fac,atom->z()*fac);
      ofs << buffer << endl;

      if (!(atom->GetIdx()%10))
      {
	sprintf(buffer,"# %d",atom->GetIdx());
	ofs << buffer << endl;
      }
    }

  ofs << "END" << endl;

  return(true);
}

/* these are the routines that babel calls */
bool WriteGromos96A(ostream &ofs,OBMol &mol)
{ 
  return(WriteGromos96(ofs,mol,1.0));
}

/* convert A -> nm */
bool WriteGromos96N(ostream &ofs,OBMol &mol)
{ 
  return(WriteGromos96(ofs,mol,0.1));
}


}
