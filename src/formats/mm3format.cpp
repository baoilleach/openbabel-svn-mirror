/**********************************************************************
Copyright (C) 2002-2003 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "obconversion.h"

using namespace std;
namespace OpenBabel {

class MM3Format : public OBFormat
{
public:
	//Register this format type ID
	MM3Format() {OBConversion::RegisterFormat("MM3",this);}

	virtual const char* Description() //required
	{ return
"MM3 format\n \
Currently writes MM2 files\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://??";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTREADABLE;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
			ret=WriteMolecule(pmol,pConv);
		delete pOb; 
		return ret;
	};
};

//Make an instance of the format class
MM3Format theMM3Format;

////////////////////////////////////////////////////////////////
// \todo Update to write genuine MM3 files -- this is currently MM2 code.

bool MM3Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

  unsigned int i;
  char buffer[BUFF_SIZE];
  OBBond *bond;
  vector<OBEdgeBase*>::iterator j;
  
  sprintf(buffer,"%6d %-20s   MM2 parameters",mol.NumAtoms(),mol.GetTitle());
  ofs << buffer << endl;

  ttab.SetFromType("INT");

  OBAtom *atom;
  string str,str1;
  for(i = 1;i <= mol.NumAtoms(); i++)
  {
    atom = mol.GetAtom(i);
    str = atom->GetType();
    ttab.SetToType("MM3"); ttab.Translate(str1,str);
    sprintf(buffer,"%6d %2s  %12.6f%12.6f%12.6f %5d",
	    i,
	    etab.GetSymbol(atom->GetAtomicNum()),
	    atom->GetX(),
	    atom->GetY(),
	    atom->GetZ(),
	    atoi((char*)str1.c_str()));
    ofs << buffer;
    
    for (bond = atom->BeginBond(j); bond; bond = atom->NextBond(j))
      {
	sprintf(buffer,"%6d", (bond->GetNbrAtom(atom))->GetIdx());
 	ofs << buffer;
      }

    ofs << endl;
  }

  return(true);
}

} //namespace OpenBabel
//   int i;
//   int type_name;
//   int connections = 0;
//   int attachments = 0;
//   char temp_type[5];

//   char ID[60]; /* filename */
//   int METHOD; /* 0 no cojugated pi system, 1 if conjugated pi system */
//   int N; /* #of atoms */
//   int IPRINT; /* Controls amount of printout */
//   int NSTR; /* Restricted motion data  */
//   int INIT; /* Minimize energy  */
//   int NCONST; /* Read in new constants ? */
//   double TMAX; /* Max time */

//   int NCON; /* Number of connected atoms */
//   int NATTACH; /*Number of attached atoms */
//   double DEL; /* Termianation of geometry optimization */
//   int NSYMM;/* Number of symmetry matrices */
//   int NX; /* Number of coordiante calcualtions or replacement cards */
//   int NROT; /* Reorient */
//   int LABEL; /* Change names or atomic weights */
//   int NDC; /* Dipole and charge interaction energy */
//   int NCALC; /* Crystal conversions */
//   int HFORM; /* Heat of formation */
//   int MVDW; /* Approximate van der Waals */
//   int NDRIVE; /* Dihedral driver */
  
  
//   for (i = 0;i < Bonds; i++)
//   {
//     if ((Valence(Start(i)) == 1) || (Valence(End(i)) == 1))
//       attachments ++;
//   }
  
//   connections = Bonds - attachments;  
//   strcpy(ID,OutfileName);

// /*------ CARD 1 -------*/
//   METHOD = 1;
//   N = Atoms;
//   IPRINT = 3;
//   NSTR = 0;
//   INIT = 0;
//   NCONST = 0;
//   TMAX = 999.0;
// /*------ CARD 2 -------*/
//   DEL = 0.00008;
//   NCON = connections;
//   NATTACH = attachments;
//   NSYMM = 0;
//   NX = 0;
//   NROT = 0;
//   LABEL = 0;
//   NDC = 0;
//   NCALC = 0;
//   HFORM = 0;
//   MVDW = 1;
//   NDRIVE = 0;
  
//   fprintf(file1,"%-60s%d%4d %d  %d %d  %d%-5.0f\n",
// 	  ID,
// 	  METHOD,
// 	  N,
// 	  IPRINT,
// 	  NSTR,
// 	  INIT,
// 	  NCONST,
// 	  TMAX);

//   fprintf(file1,"%1d%4d%5s%4.5f%8s%5d%5d%5d%5d%5d%5d%5d%5d%10d%5d\n",
// 	  0,
// 	  NCON,
// 	  "",
// 	  DEL,
// 	  "",
// 	  NATTACH,
// 	  NSYMM,
// 	  NX,
// 	  NROT,
// 	  LABEL,
// 	  NDC,
// 	  NCALC,
// 	  HFORM,
// 	  MVDW,
// 	  NDRIVE);

//   for(i = 0;i < Bonds; i++)
//   {
//     if ((Valence(Start(i)) > 1) && (Valence(End(i)) > 1))
//       fprintf(file1,"%5d%5d\n",
// 	    Start(i),
// 	    End(i));
//   }
  
//   attachments = 0;
  
//   for(i = 0;i < Bonds; i++)
//   {
//     if ((Valence(Start(i)) == 1) || (Valence(End(i)) == 1))
//     {
//       attachments ++;
//       fprintf(file1,"%5d%5d",
// 	    Start(i),
// 	    End(i));

//       if (((attachments % 8) == 0))
// 	fprintf(file1,"\n");
//     }
//   }
  
//   if (((attachments % 8) != 0))
//     fprintf(file1,"\n");
  
//   for (i = 1;i <= Atoms; i++)
//   {
//     get_output_type(i,"MM2",Type(i),temp_type,all_caps);
//     type_name = atoi(temp_type);
//     type_name = update_mm2_types(mol,i,type_name);
//     fprintf(file1,"  %8.5f  %8.5f  %8.5f%5d(%3d)\n",
// 	    X(i),
// 	    Y(i),
// 	    Z(i),
// 	    type_name,
// 	    i);
//   }
