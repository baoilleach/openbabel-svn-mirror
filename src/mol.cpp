/**********************************************************************
mol.cpp - Handle molecules.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2003 by Michael Banck
 
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
#include "babelconfig.h"

#include "mol.h"
#include "rotamer.h"
#include "phmodel.h"
#include "bondtyper.h"
#include "math/matrix3x3.h"
#include "obiter.h"

#include <sstream>

using namespace std;

namespace OpenBabel
{

  extern bool SwabInt;
  extern OBPhModel	phmodel;
  extern OBAromaticTyper  aromtyper;
  extern OBAtomTyper      atomtyper;
  extern OBBondTyper      bondtyper;


  /** \class OBMol
      \brief Molecule Class
 
      The most important class in Open Babel is OBMol, or the molecule class.
      The OBMol class is designed to store all the basic information
      associated with a molecule, to make manipulations on the connection
      table of a molecule facile, and to provide member functions which
      automatically perceive information about a molecule. A guided tour
      of the OBMol class is a good place to start.
 
      An OBMol class can be declared:
      \code
      OBMol mol;
      \endcode

      For example:
      \code
      #include <iostream.h>
      #include "mol.h"
      #include "obconversion.h"
      int main(int argc,char **argv)
      {
      OBConversion conv(&cin,&cout);
      if(conv.SetInAndOutFormats("SDF","MOL2"))
      {	
      OBMol mol;
      if(conv.Read(&mol))
      ...manipulate molecule 
		
      conv->Write(&mol);
      }
      return(1);
      }
      \endcode
 
      will read in a molecule in SD file format from stdin 
      (or the C++ equivalent cin) and write a MOL2 format file out
      to standard out. Additionally, The input and output formats can
      be altered using the OBConversion class

      Once a molecule has been read into an OBMol (or created via other methods)
      the atoms and bonds
      can be accessed by the following methods:
      \code
      OBAtom *atom;
      atom = mol.GetAtom(5); //random access of an atom
      \endcode
      or
      \code
      OBBond *bond;
      bond = mol.GetBond(14); //random access of a bond
      \endcode
      or
      \code
      FOR_ATOMS_IN_MOL(atom, mol) // iterator access (see OBMolAtomIter)
      \endcode
      or
      \code
      FOR_BONDS_IN_MOL(bond, mol) // iterator access (see OBMolBondIter)
      \endcode
      It is important to note that atom arrays currently begin at 1 and bond arrays
      begin at 0. Requesting atom 0 (\code
      OBAtom *atom = mol.GetAtom(0); \endcode
      will result in an error, but
      \code
      OBBond *bond = mol.GetBond(0);
      \endcode
      is perfectly valid.
      Note that this is expected to change in the near future to simplify coding
      and improve efficiency.
 
      The ambiguity of numbering issues and off-by-one errors led to the use
      of iterators in Open Babel. An iterator is essentially just a pointer, but
      when used in conjunction with Standard Template Library (STL) vectors
      it provides an unambiguous way to loop over arrays. OBMols store their
      atom and bond information in STL vectors. Since vectors are template
      based, a vector of any user defined type can be declared. OBMols declare
      vector<OBNodeBase*> and vector<OBEdgeBase*> to store atom and bond information.
      Iterators are then a natural way to loop over the vectors of atoms and bonds.
 
      A variety of predefined iterators have been created to simplify
      common looping requests (e.g., looping over all atoms in a molecule,
      bonds to a given atom, etc.)

      \code
      #include "obiter.h"
      ...
      #define FOR_ATOMS_OF_MOL(a,m)     for( OBMolAtomIter     a(m); a; a++ )
      #define FOR_BONDS_OF_MOL(b,m)     for( OBMolBondIter     b(m); b; b++ )
      #define FOR_NBORS_OF_ATOM(a,p)    for( OBAtomAtomIter    a(p); a; a++ )
      #define FOR_BONDS_OF_ATOM(b,p)    for( OBAtomBondIter    b(p); b; b++ )
      #define FOR_RESIDUES_OF_MOL(r,m)  for( OBResidueIter     r(m); r; r++ )
      #define FOR_ATOMS_OF_RESIDUE(a,r) for( OBResidueAtomIter a(r); a; a++ )
      ...
      \endcode

      These convenience functions can be used like so:
      \code
      #include "obiter.h"
      #include "mol.h"

      OBMol mol;
      double exactMass = 0.0f;
      FOR_ATOMS_OF_MOL(a, mol)
      {
      exactMass +=  a->GetExactMass();
      }
      \endcode

      Note that with these convenience macros, the iterator "a" (or
      whichever name you pick) is declared for you -- you do not need to
      do it beforehand.
  */

  //
  // OBMol member functions
  //
  void  OBMol::SetTitle(const char *title)
  { 
    _title = title;
    Trim(_title);
  }

  void  OBMol::SetTitle(std::string &title)
  { 
    _title = title;
    Trim(_title);
  }

  bool SortVVInt(const vector<int> &a,const vector<int> &b)
  {
    return(a.size() > b.size());
  }

  bool SortAtomZ(const pair<OBAtom*,double> &a, const pair<OBAtom*,double> &b)
  {
    return (a.second < b.second);
  }

  double OBMol::GetTorsion(int a,int b,int c,int d)
  {
    return(CalcTorsionAngle(((OBAtom*)_vatom[a-1])->GetVector(),
                            ((OBAtom*)_vatom[b-1])->GetVector(),
                            ((OBAtom*)_vatom[c-1])->GetVector(),
                            ((OBAtom*)_vatom[d-1])->GetVector()));
  }

  void OBMol::SetTorsion(OBAtom *a,OBAtom *b,OBAtom *c, OBAtom *d, double ang)
  {
    vector<int> tor;
    vector<int> atoms;

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::SetTorsion", obAuditMsg);

    tor.push_back(a->GetCIdx());
    tor.push_back(b->GetCIdx());
    tor.push_back(c->GetCIdx());
    tor.push_back(d->GetCIdx());

    FindChildren(atoms, b->GetIdx(), c->GetIdx());
    int j;
    for (j = 0 ; (unsigned)j < atoms.size() ; j++ )
      atoms[j] = (atoms[j] - 1) * 3;

    double v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
    double c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
    double c1mag,c2mag,radang,costheta,m[9];
    double x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

    //calculate the torsion angle

    v1x = _c[tor[0]]   - _c[tor[1]];
    v2x = _c[tor[1]]   - _c[tor[2]];
    v1y = _c[tor[0]+1] - _c[tor[1]+1];
    v2y = _c[tor[1]+1] - _c[tor[2]+1];
    v1z = _c[tor[0]+2] - _c[tor[1]+2];
    v2z = _c[tor[1]+2] - _c[tor[2]+2];
    v3x = _c[tor[2]]   - _c[tor[3]];
    v3y = _c[tor[2]+1] - _c[tor[3]+1];
    v3z = _c[tor[2]+2] - _c[tor[3]+2];

    c1x = v1y*v2z - v1z*v2y;
    c2x = v2y*v3z - v2z*v3y;
    c1y = -v1x*v2z + v1z*v2x;
    c2y = -v2x*v3z + v2z*v3x;
    c1z = v1x*v2y - v1y*v2x;
    c2z = v2x*v3y - v2y*v3x;
    c3x = c1y*c2z - c1z*c2y;
    c3y = -c1x*c2z + c1z*c2x;
    c3z = c1x*c2y - c1y*c2x;

    c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
    c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
    if (c1mag*c2mag < 0.01)
      costheta = 1.0; //avoid div by zero error
    else
      costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

    if (costheta < -0.999999)
      costheta = -0.999999;
    if (costheta >  0.999999)
      costheta =  0.999999;

    if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0)
      radang = -acos(costheta);
    else
      radang = acos(costheta);

    //
    // now we have the torsion angle (radang) - set up the rot matrix
    //

    //find the difference between current and requested
    rotang = ang - radang;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;
    //normalize the rotation vector
    mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
    x = v2x/mag;
    y = v2y/mag;
    z = v2z/mag;

    //set up the rotation matrix
    m[0]= t*x*x + cs;
    m[1] = t*x*y + sn*z;
    m[2] = t*x*z - sn*y;
    m[3] = t*x*y - sn*z;
    m[4] = t*y*y + cs;
    m[5] = t*y*z + sn*x;
    m[6] = t*x*z + sn*y;
    m[7] = t*y*z - sn*x;
    m[8] = t*z*z + cs;

    //
    //now the matrix is set - time to rotate the atoms
    //
    tx = _c[tor[1]];
    ty = _c[tor[1]+1];
    tz = _c[tor[1]+2];
    vector<int>::iterator i;
    for (i = atoms.begin(),j=*i;i != atoms.end();i++,j=*i)
      {
        _c[j] -= tx;
        _c[j+1] -= ty;
        _c[j+2]-= tz;
        x = _c[j]*m[0] + _c[j+1]*m[1] + _c[j+2]*m[2];
        y = _c[j]*m[3] + _c[j+1]*m[4] + _c[j+2]*m[5];
        z = _c[j]*m[6] + _c[j+1]*m[7] + _c[j+2]*m[8];
        _c[j] = x;
        _c[j+1] = y;
        _c[j+2] = z;
        _c[j] += tx;
        _c[j+1] += ty;
        _c[j+2] += tz;
      }
  }


  double OBMol::GetTorsion(OBAtom *a,OBAtom *b,OBAtom *c,OBAtom *d)
  {
    return(CalcTorsionAngle(a->GetVector(),
                            b->GetVector(),
                            c->GetVector(),
                            d->GetVector()));
  }

  void OBMol::ContigFragList(std::vector<std::vector<int> >&cfl)
  {
    int j;
    OBAtom *atom;
    OBBond *bond;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator k;
    OBBitVec used,curr,next,frag;
    vector<int> tmp;

    used.Resize(NumAtoms()+1);
    curr.Resize(NumAtoms()+1);
    next.Resize(NumAtoms()+1);
    frag.Resize(NumAtoms()+1);

    while ((unsigned)used.CountBits() < NumAtoms())
      {
        curr.Clear();
        frag.Clear();
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  if (!used.BitIsOn(atom->GetIdx()))
            {
	      curr.SetBitOn(atom->GetIdx());
	      break;
            }

        frag |= curr;
        while (!curr.IsEmpty())
	  {
            next.Clear();
            for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
	      {
                atom = GetAtom(j);
                for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
		  if (!used.BitIsOn(bond->GetNbrAtomIdx(atom)))
		    next.SetBitOn(bond->GetNbrAtomIdx(atom));
	      }

            used |= curr;
            used |= next;
            frag |= next;
            curr = next;
	  }

        tmp.clear();
        frag.ToVecInt(tmp);
        cfl.push_back(tmp);
      }

    sort(cfl.begin(),cfl.end(),SortVVInt);
  }

  /*!
  **\brief Fills the OBGeneric OBTorsionData with torsions from the mol
  */
  void OBMol::FindTorsions()
  {
    //if already has data return
    if(HasData(OBGenericDataType::TorsionData))
      return;

    //get new data and attach it to molecule
    OBTorsionData *torsions = new OBTorsionData;
    SetData(torsions);

    OBTorsion torsion;
    vector<OBEdgeBase*>::iterator bi1,bi2,bi3;
    OBBond* bond;
    OBAtom *a,*b,*c,*d;

    //loop through all bonds generating torsions
    for(bond = BeginBond(bi1);bond;bond = NextBond(bi1))
      {
        b = bond->GetBeginAtom();
        c = bond->GetEndAtom();
        if(b->IsHydrogen() || c->IsHydrogen())
	  continue;

        for(a = b->BeginNbrAtom(bi2);a;a = b->NextNbrAtom(bi2))
	  {
            if(a == c)
	      continue;

            for(d = c->BeginNbrAtom(bi3);d;d = c->NextNbrAtom(bi3))
	      {
                if(d == b)
		  continue;
                torsion.AddTorsion(a,b,c,d);
	      }
	  }
        //add torsion to torsionData
        if(torsion.GetSize())
	  torsions->SetData(torsion);
        torsion.Clear();
      }

    return;
  }

  void OBMol::FindLargestFragment(OBBitVec &lf)
  {
    int j;
    OBAtom *atom;
    OBBond *bond;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator k;
    OBBitVec used,curr,next,frag;

    lf.Clear();
    while ((unsigned)used.CountBits() < NumAtoms())
      {
        curr.Clear();
        frag.Clear();
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  if (!used.BitIsOn(atom->GetIdx()))
            {
	      curr.SetBitOn(atom->GetIdx());
	      break;
            }

        frag |= curr;
        while (!curr.IsEmpty())
	  {
            next.Clear();
            for (j = curr.NextBit(-1);j != curr.EndBit();j = curr.NextBit(j))
	      {
                atom = GetAtom(j);
                for (bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
		  if (!used.BitIsOn(bond->GetNbrAtomIdx(atom)))
		    next.SetBitOn(bond->GetNbrAtomIdx(atom));
	      }

            used |= curr;
            used |= next;
            frag |= next;
            curr = next;
	  }

        if (lf.Empty() || lf.CountBits() < frag.CountBits())
	  lf = frag;
      }
  }

  //! locates all atoms for which there exists a path to 'end'
  //! without going through 'bgn'
  //! children must not include 'end'
  void OBMol::FindChildren(vector<OBAtom*> &children,OBAtom *bgn,OBAtom *end)
  {
    OBBitVec used,curr,next;

    used |= bgn->GetIdx();
    used |= end->GetIdx();
    curr |= end->GetIdx();
    children.clear();

    int i;
    OBAtom *atom,*nbr;
    vector<OBEdgeBase*>::iterator j;

    for (;;)
      {
        next.Clear();
        for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
	  {
            atom = GetAtom(i);
            for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
	      if (!used[nbr->GetIdx()])
                {
		  children.push_back(nbr);
		  next |= nbr->GetIdx();
		  used |= nbr->GetIdx();
                }
	  }
        if (next.Empty())
	  break;
        curr = next;
      }
  }

  //! locates all atoms for which there exists a path to 'second'
  //! without going through 'first'
  //! children must not include 'second'
  void OBMol::FindChildren(vector<int> &children,int first,int second)
  {
    int i;
    OBBitVec used,curr,next;

    used.SetBitOn(first);
    used.SetBitOn(second);
    curr.SetBitOn(second);

    OBAtom *atom;
    OBBond *bond;
    vector<OBEdgeBase*>::iterator j;

    while (!curr.IsEmpty())
      {
        next.Clear();
        for (i = curr.NextBit(-1);i != curr.EndBit();i = curr.NextBit(i))
	  {
            atom = GetAtom(i);
            for (j = atom->BeginBonds(),bond=(OBBond *)*j;
		 j != atom->EndBonds();j++,bond=(OBBond *)*j)
	      if (!used.BitIsOn(bond->GetNbrAtomIdx(atom)))
		next.SetBitOn(bond->GetNbrAtomIdx(atom));
	  }

        used |= next;
        curr = next;
      }

    used.SetBitOff(first);
    used.SetBitOff(second);
    used.ToVecInt(children);
  }

  /*!
  **\brief Calculates the graph theoretical distance of each atom.
  ** Vector is indexed from zero
  */
  bool OBMol::GetGTDVector(vector<int> &gtd)
    //calculates the graph theoretical distance for every atom
    //and puts it into gtd
  {
    gtd.clear();
    gtd.resize(NumAtoms());

    int gtdcount,natom;
    OBBitVec used,curr,next;
    OBAtom *atom,*atom1;
    OBBond *bond;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;

    next.Clear();

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        gtdcount = 0;
        used.Clear();
        curr.Clear();
        used.SetBitOn(atom->GetIdx());
        curr.SetBitOn(atom->GetIdx());

        while (!curr.IsEmpty())
	  {
            next.Clear();
            for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom))
	      {
                atom1 = GetAtom(natom);
                for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j))
		  if (!used.BitIsOn(bond->GetNbrAtomIdx(atom1)) && !curr.BitIsOn(bond->GetNbrAtomIdx(atom1)))
		    if (!(bond->GetNbrAtom(atom1))->IsHydrogen())
		      next.SetBitOn(bond->GetNbrAtomIdx(atom1));
	      }

            used |= next;
            curr = next;
            gtdcount++;
	  }
        gtd[atom->GetIdx()-1] = gtdcount;
      }
    return(true);
  }

  /*!
  **\brief Calculates a set of graph invariant indexes using
  ** the graph theoretical distance, number of connected heavy atoms,
  ** aromatic boolean, ring boolean, atomic number, and 
  ** summation of bond orders connected to the atom.
  ** Vector is indexed from zero
  */
  void OBMol::GetGIVector(vector<unsigned int> &vid)
  {
    vid.clear();
    vid.resize(NumAtoms()+1);

    vector<int> v;
    GetGTDVector(v);

    int i;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;
    for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),i++)
      {
        vid[i] =  (unsigned int)v[i];
        vid[i] += (unsigned int)(atom->GetHvyValence()*100);
        vid[i] += (unsigned int)(((atom->IsAromatic()) ? 1 : 0)*1000);
        vid[i] += (unsigned int)(((atom->IsInRing()) ? 1 : 0)*10000);
        vid[i] += (unsigned int)(atom->GetAtomicNum()*100000);
        vid[i] += (unsigned int)(atom->GetImplicitValence()*10000000);
      }
  }

  static bool OBComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  static bool OBComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.first->GetIdx() < b.first->GetIdx());
  }

  //! counts the number of unique symmetry classes in a list
  static void ClassCount(vector<pair<OBAtom*,unsigned int> > &vp,unsigned int &count)
  {
    count = 0;
    vector<pair<OBAtom*,unsigned int> >::iterator k;
    sort(vp.begin(),vp.end(),OBComparePairSecond);
#if 0 // original version

    unsigned int id=0; // [ejk] appease gcc's bogus "might be undef'd" warning
    for (k = vp.begin();k != vp.end();k++)
      {
        if (k == vp.begin())
	  {
            id = k->second;
            k->second = count = 0;
	  }
        else
	  if (k->second != id)
            {
	      id = k->second;
	      k->second = ++count;
            }
	  else
	    k->second = count;
      }
    count++;
#else // get rid of warning, moves test out of loop, returns 0 for empty input

    k = vp.begin();
    if (k != vp.end())
      {
        unsigned int id = k->second;
        k->second = 0;
        ++k;
        for (;k != vp.end(); ++k)
	  {
            if (k->second != id)
	      {
                id = k->second;
                k->second = ++count;
	      }
            else
	      k->second = count;
	  }
        ++count;
      }
    else
      {
        // [ejk] thinks count=0 might be OK for an empty list, but orig code did
        //++count;
      }
#endif
  }

  //! creates a new vector of symmetry classes base on an existing vector
  //! helper routine to GetGIDVector
  static void	CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,vector<pair<OBAtom*,unsigned int> > &vp2)
  {
    int m,id;
    OBAtom *nbr;
    vector<OBEdgeBase*>::iterator j;
    vector<unsigned int>::iterator k;
    vector<pair<OBAtom*,unsigned int> >::iterator i;
    sort(vp1.begin(),vp1.end(),OBComparePairFirst);
    vp2.clear();
    for (i = vp1.begin();i != vp1.end();i++)
      {
        vector<unsigned int> vtmp;
        for (nbr = i->first->BeginNbrAtom(j);nbr;nbr = i->first->NextNbrAtom(j))
	  vtmp.push_back(vp1[nbr->GetIdx()-1].second);
        sort(vtmp.begin(),vtmp.end(),OBCompareUnsigned);
        for (id=i->second,m=100,k = vtmp.begin();k != vtmp.end();k++,m*=100)
	  id += *k * m;

        vp2.push_back(pair<OBAtom*,unsigned int> (i->first,id));
      }
  }

  /*!
  **\brief Calculates a set of symmetry identifiers for a molecule.
  ** Atoms with the same symmetry ID are symmetrically equivalent.
  ** Vector is indexed from zero
  */
  void OBMol::GetGIDVector(vector<unsigned int> &vgid)
  {
    vector<unsigned int> vgi;
    GetGIVector(vgi);  //get vector of graph invariants

    int i;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;
    vector<pair<OBAtom*,unsigned int> > vp1,vp2;
    for (i=0,atom = BeginAtom(j);atom;atom = NextAtom(j),i++)
      vp1.push_back(pair<OBAtom*,unsigned int> (atom,vgi[i]));

    unsigned int nclass1,nclass2; //number of classes
    ClassCount(vp1,nclass1);

    if (nclass1 < NumAtoms())
      {
        for (i = 0;i < 100;i++) //sanity check - shouldn't ever hit this number
	  {
            CreateNewClassVector(vp1,vp2);
            ClassCount(vp2,nclass2);
            vp1 = vp2;
            if (nclass1 == nclass2)
	      break;
            nclass1 = nclass2;
	  }
      }

    vgid.clear();
    sort(vp1.begin(),vp1.end(),OBComparePairFirst);
    vector<pair<OBAtom*,unsigned int> >::iterator k;
    for (k = vp1.begin();k != vp1.end();k++)
      vgid.push_back(k->second);
  }

  unsigned int OBMol::NumHvyAtoms()
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator(i);
    unsigned int count = 0;

    for(atom = this->BeginAtom(i);atom;atom = this->NextAtom(i))
      {
        if(!atom->IsHydrogen())
	  count++;
      }

    return(count);
  }

  unsigned int OBMol::NumRotors()
  {
    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;

    unsigned int count = 0;
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->IsRotor())
	count++;

    return(count);
  }

  //! Returns a pointer to the atom after a safety check
  //! 0 < idx <= NumAtoms
  OBAtom *OBMol::GetAtom(int idx)
  {
    if ((unsigned)idx < 1 || (unsigned)idx > NumAtoms())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Requested Atom Out of Range", obDebug);
        return((OBAtom*)NULL);
      }

    return((OBAtom*)_vatom[idx-1]);
  }

  OBAtom *OBMol::GetFirstAtom()
  {
    return((_vatom.empty()) ? (OBAtom*)NULL : (OBAtom*)_vatom[0]);
  }

  //! Returns a pointer to the bond after a safety check
  //! 0 <= idx < NumBonds
  OBBond *OBMol::GetBond(int idx)
  {
    if (idx < 0 || (unsigned)idx >= NumBonds())
      {
	obErrorLog.ThrowError(__FUNCTION__, "Requested Bond Out of Range", obDebug);
        return((OBBond*)NULL);
      }

    return((OBBond*)_vbond[idx]);
  }

  OBBond *OBMol::GetBond(int bgn, int end)
  {
    return(GetBond(GetAtom(bgn),GetAtom(end)));
  }

  OBBond *OBMol::GetBond(OBAtom *bgn,OBAtom *end)
  {
    OBAtom *nbr;
    vector<OBEdgeBase*>::iterator i;

    for (nbr = bgn->BeginNbrAtom(i);nbr;nbr = bgn->NextNbrAtom(i))
      if (nbr == end)
	return((OBBond *)*i);

    return(NULL); //just to keep the SGI compiler happy
  }

  OBResidue *OBMol::GetResidue(int idx)
  {
    if (idx < 0 || (unsigned)idx >= NumResidues())
      {
	obErrorLog.ThrowError(__FUNCTION__, "Requested Residue Out of Range", obDebug);
        return((OBResidue*)NULL);
      }

    return (_residue[idx]);
  }

  std::vector<OBInternalCoord*> OBMol::GetInternalCoord()
  {
    if (_internals.empty())
      {
        _internals.push_back((OBInternalCoord*)NULL);
        for(unsigned int i = 1; i <= NumAtoms(); i++)
	  {
            _internals.push_back(new OBInternalCoord);
	  }
        CartesianToInternal(_internals, *this);
      }
    return _internals;
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#findSmallestSetOfSmallestRings">blue-obelisk:findSmallestSetOfSmallestRings</a>.
  vector<OBRing*> &OBMol::GetSSSR()
  {
    if (!HasSSSRPerceived())
      FindSSSR();

    if (!HasData(OBGenericDataType::RingData))
      SetData(new OBRingData);

    OBRingData *rd = (OBRingData *) GetData(OBGenericDataType::RingData);
    return(rd->GetData());
  }

  double OBMol::GetMolWt()
  {
    double molwt=0.0;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    bool UseImplicitH = NumHvyAtoms() && (NumBonds()!=0 || NumAtoms()==1);
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
	if(UseImplicitH)
	  {
	    if (! atom->IsHydrogen())
	      molwt += etab.GetMass(1) * atom->ImplicitHydrogenCount();
	  }
	molwt += atom->GetAtomicMass();
      }
    return(molwt);
  }

  double OBMol::GetExactMass()
  {
    double mass=0.0;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    bool UseImplicitH = NumHvyAtoms() && (NumBonds()!=0 || NumAtoms()==1);
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
	if(UseImplicitH)
	  {
	    if (!atom->IsHydrogen())
	      mass += isotab.GetExactMass(1,1) * atom->ImplicitHydrogenCount();
	  }
	mass += atom->GetExactMass();
      }
    return(mass);
  }

  //! Stochoimetric formula in spaced format e.g. C 4 H 6 O 1
  //! No pair data is stored. Normally use without parameters: GetSpacedFormula()
  string OBMol::GetSpacedFormula(int ones, const char* sp)
  {
    //Default ones=0, sp=' '. 
    //Using ones=1 and sp='' will give unspaced formula (and no pair data entry)
    // These are the atomic numbers of the elements in alphabetical order.
    const int NumElements = 112;
    const int alphabetical[NumElements] = {
      89, 47, 13, 95, 18, 33, 85, 79, 5, 56, 4, 107, 83, 97, 35, 6, 20, 48,
      58, 98, 17, 96, 27, 24, 55, 29, 111, 105, 66, 68, 99, 63, 9, 26, 100, 87, 31,
      64, 32, 1, 2, 72, 80, 67, 108, 53, 49, 77, 19, 36, 57, 3, 103, 71, 101,
      12, 25, 42, 109, 7, 11, 41, 60, 10, 28, 102, 93, 8, 76, 15, 91, 82, 46, 
      61, 84, 59, 78, 94, 88, 37, 75, 104, 45, 86, 44, 16, 51, 21, 34, 106, 14, 
      62, 50, 38, 112, 73, 65, 43, 52, 90, 22, 81, 69, 92, 110, 23, 74, 54, 39, 70, 
      30, 40 };

    int atomicCount[NumElements];
    stringstream formula;

    for (int i = 0; i < NumElements; i++)
      atomicCount[i] = 0;

    bool UseImplicitH = (NumBonds()!=0 || NumAtoms()==1);
    bool HasHvyAtoms = NumHvyAtoms()>0;
    FOR_ATOMS_OF_MOL(a, *this)
      {
	int anum = a->GetAtomicNum();
	if(anum==0)
	  continue;
	bool IsHiso = anum == 1 && a->GetIsotope()>=2;
	if(UseImplicitH)
	  {
	    if (anum == 1 && !IsHiso && HasHvyAtoms)
	      continue;	// skip explicit hydrogens except D,T
	    if(anum==1)
	      {
		if (IsHiso && HasHvyAtoms)
		  --atomicCount[0]; //one of the implicit hydrogens is now explicit
	      }
	    else
	      atomicCount[0] += a->ImplicitHydrogenCount() + a->ExplicitHydrogenCount();
	  }
	if (IsHiso)
	  anum = NumElements + a->GetIsotope() - 3; //pseudo AtNo for D, T
	atomicCount[anum - 1]++;
      }
    
    if (atomicCount[5] != 0) // Carbon (i.e. 6 - 1 = 5)
      {
	if (atomicCount[5] > ones)
	  formula << "C" << sp << atomicCount[5] << sp;
	else if (atomicCount[5] == 1)
	  formula << "C";
	
	atomicCount[5] = 0; // So we don't output C twice
	
	// only output H if there's also carbon -- otherwise do it alphabetical
	if (atomicCount[0] != 0) // Hydrogen (i.e., 1 - 1 = 0)
	  {
	    if (atomicCount[0] > ones)
	      formula << "H" << sp << atomicCount[0] << sp;
	    else if (atomicCount[0] == 1)
	      formula << "H";
	    
	    atomicCount[0] = 0;
	  }
      }
    
    for (int j = 0; j < NumElements; j++)
      {
	char DT[4] = {'D',0,'T',0};
	char* symb;
	int alph = alphabetical[j]-1;
	if (atomicCount[ alph ])
	  {
	    if(alph==NumElements-1)
	      symb = DT + 2;//T
	    else if (alph==NumElements-2)
	      symb = DT; //D
	    else
	      symb = etab.GetSymbol(alphabetical[j]);
	    
	    formula << symb << sp; 
	    if(atomicCount[alph] > ones)
	      formula << sp << atomicCount[alph] << sp;
	  }
      }

    return (formula.str());
  }

  //! Stochoimetric formula (e.g., C4H6O).
  //!   This is either set by OBMol::SetFormula() or generated on-the-fly
  //!   using the "Hill order" -- i.e., C first if present, then H if present
  //!   all other elements in alphabetical order.
  string OBMol::GetFormula()
  {
    string attr = "Formula";
    OBPairData *dp = (OBPairData *) GetData(attr);
  
    if (dp != NULL) // we already set the formula
      return dp->GetValue();

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::SetFormula -- Hill order formula",
			  obAuditMsg);

    string sformula = GetSpacedFormula(1,"");

    dp = new OBPairData;
    dp->SetAttribute(attr);
    dp->SetValue( sformula );
    SetData(dp);

    return sformula;
  }

  void OBMol::SetFormula(string molFormula)
  {
    string attr = "Formula";
    OBPairData *dp = (OBPairData *) GetData(attr);
  
    if (dp == NULL)
      {
	dp = new OBPairData;
	dp->SetAttribute(attr);
      }
    dp->SetValue(molFormula);

    SetData(dp);
  }

  void OBMol::SetTotalCharge(int charge)
  {
    SetFlag(OB_TCHARGE_MOL);
    _totalCharge = charge;
  }

  //! Returns the total molecular charge -- if it has not previously been set
  //!  it is calculated from the atomic formal charge information.
  //!  (This may or may not be correct!)
  //!  If you set atomic charges with OBAtom::SetFormalCharge()
  //!   you really should set the molecular charge with OBMol::SetTotalCharge()
  int OBMol::GetTotalCharge()
  {
    if(HasFlag(OB_TCHARGE_MOL))
      return(_totalCharge);
    else // calculate from atomic formal charges (seems the best default)
      {
	obErrorLog.ThrowError(__FUNCTION__,
			      "Ran OpenBabel::GetTotalCharge -- calculated from formal charges", 
			      obAuditMsg);

        OBAtom *atom;
        vector<OBNodeBase*>::iterator i;
        int chg = 0;

        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  chg += atom->GetFormalCharge();
        return (chg);
      }
  }

  void   OBMol::SetTotalSpinMultiplicity(unsigned int spin)
  {
    SetFlag(OB_TSPIN_MOL);
    _totalSpin = spin;
  }

  //! Returns the total spin multiplicity -- if it has not previously been set
  //!  it is calculated from the atomic spin multiplicity information
  //!  assuming the high-spin case (i.e. it simply sums the atomic spins,
  //!  making no attempt to pair spins).
  //!  However, if you set atomic spins with OBAtom::SetSpinMultiplicity()
  //!   you really should set the molecular spin with
  //!   OBMol::SetTotalSpinMultiplicity()
  unsigned int OBMol::GetTotalSpinMultiplicity()
  {
    if (HasFlag(OB_TSPIN_MOL))
      return(_totalSpin);
    else // calculate from atomic spin information (assuming high-spin case)
      {
	obErrorLog.ThrowError(__FUNCTION__,
			      "Ran OpenBabel::GetTotalSpinMultiplicity -- calculating from atomic spins and high spin",
			      obAuditMsg);

        OBAtom *atom;
        vector<OBNodeBase*>::iterator i;
        unsigned int spin = 1;

        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  {
            if (atom->GetSpinMultiplicity() > 1)
	      spin += atom->GetSpinMultiplicity() - 1;
	  }
        return (spin);
      }
  }

  OBMol &OBMol::operator=(const OBMol &source)
    //only atom and bond info is copied from src to dest
    //Conformers are now copied also, MM 2/7/01
    //Residue information are copied, MM 4-27-01
    //All OBGenericData incl OBRotameterList is copied, CM 2006
  {
    OBMol &src = (OBMol &)source;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;
    OBAtom *atom;
    OBBond *bond;

    Clear();
    BeginModify();

    _vatom.reserve(src.NumAtoms());
    _vbond.reserve(src.NumBonds());

    for (atom = src.BeginAtom(i);atom;atom = src.NextAtom(i))
      AddAtom(*atom);
    for (bond = src.BeginBond(j);bond;bond = src.NextBond(j))
      AddBond(*bond);

    this->_title  = src.GetTitle();
    this->_energy = src.GetEnergy();
    this->_dimension = src.GetDimension();
		
    EndModify();

    //Copy Residue information
    unsigned int NumRes = src.NumResidues();
    if (NumRes)
      {
        unsigned int k;
        OBResidue *src_res=NULL;
        OBResidue *res=NULL;
        OBAtom *src_atom=NULL;
        OBAtom *atom=NULL;
        vector<OBAtom*>::iterator ii;
        for (k=0 ; k<NumRes ; k++)
	  {
            res = NewResidue();
            src_res = src.GetResidue(k);
            res->SetName(src_res->GetName());
            res->SetNum(src_res->GetNum());
            res->SetChain(src_res->GetChain());
            res->SetChainNum(src_res->GetChainNum());
            for (src_atom=src_res->BeginAtom(ii) ; src_atom ; src_atom=src_res->NextAtom(ii))
	      {
                atom = GetAtom(src_atom->GetIdx());
                res->AddAtom(atom);
                res->SetAtomID(atom,src_res->GetAtomID(src_atom));
                res->SetHetAtom(atom,src_res->IsHetAtom(src_atom));
                res->SetSerialNum(atom,src_res->GetSerialNum(src_atom));
	      }
	  }
      }

    //Copy conformer information
    if (src.NumConformers() > 1)
      {
        int k,l;
        vector<double*> conf;
        double* xyz = NULL;
        for (k=0 ; k<src.NumConformers() ; k++)
	  {
            xyz = new double [3*src.NumAtoms()];
            for (l=0 ; l<(int) (3*src.NumAtoms()) ; l++)
	      xyz[l] = src.GetConformer(k)[l];
            conf.push_back(xyz);
	  }
        SetConformers(conf);
      }

    /* Now done with other OBGenericData
    //Copy rotamer list
    OBRotamerList *rml = (OBRotamerList *)src.GetData(OBGenericDataType::RotamerList);
    if (rml && rml->NumAtoms() == src.NumAtoms())
    {
    //Destroy old rotamer list if necessary
    if ((OBRotamerList *)GetData(OBGenericDataType::RotamerList))
    {
    DeleteData(OBGenericDataType::RotamerList);
    }

    //Set base coordinates
    OBRotamerList *cp_rml = new OBRotamerList;
    unsigned int k,l;
    vector<double*> bc;
    double *c=NULL;
    double *cc=NULL;
    for (k=0 ; k<rml->NumBaseCoordinateSets() ; k++)
    {
    c = new double [3*rml->NumAtoms()];
    cc = rml->GetBaseCoordinateSet(k);
    for (l=0 ; l<3*rml->NumAtoms() ; l++)
    c[l] = cc[l];
    bc.push_back(c);
    }
    if (rml->NumBaseCoordinateSets())
    cp_rml->SetBaseCoordinateSets(bc,rml->NumAtoms());

    //Set reference array
    unsigned char *ref = new unsigned char [rml->NumRotors()*4];
    if (ref)
    {
    rml->GetReferenceArray(ref);
    cp_rml->Setup((*this),ref,rml->NumRotors());
    delete [] ref;
    }

    //Set Rotamers
    unsigned char *rotamers = new unsigned char [(rml->NumRotors()+1)*rml->NumRotamers()];
    if (rotamers)
    {
    vector<unsigned char*>::iterator kk;
    unsigned int idx=0;
    for (kk = rml->BeginRotamer();kk != rml->EndRotamer();kk++)
    {
    memcpy(&rotamers[idx],(const unsigned char*)*kk,sizeof(unsigned char)*(rml->NumRotors()+1));
    idx += sizeof(unsigned char)*(rml->NumRotors()+1);
    }
    cp_rml->AddRotamers(rotamers,rml->NumRotamers());
    delete [] rotamers;
    }
    SetData(cp_rml);
    }
    */

    //Copy all the OBGenericData, providing the new molecule, this,
    //for those classes like OBRotameterList which contain Atom pointers
    //OBGenericData classes can choose not to be cloned by returning NULL
    vector<OBGenericData*>::iterator itr;
    for(itr=src.BeginData();itr!=src.EndData();++itr)
      {
	OBGenericData* pCopiedData = (*itr)->Clone(this);
	SetData(pCopiedData);
      }
   
    return(*this);
  }

  OBMol &OBMol::operator+=(const OBMol &source)
  {
    OBMol &src = (OBMol &)source;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;
    OBAtom *atom;
    OBBond *bond;

    BeginModify();

    int prevatms = NumAtoms();

    _title += "_" + string(src.GetTitle());

    for (atom = src.BeginAtom(i) ; atom ; atom = src.NextAtom(i))
      AddAtom(*atom);
    for (bond = src.BeginBond(j) ; bond ; bond = src.NextBond(j))
      AddBond(bond->GetBeginAtomIdx() + prevatms, bond->GetEndAtomIdx() + prevatms, bond->GetBO(), bond->GetFlags());
    if(src.GetDimension()<_dimension)
      _dimension = src.GetDimension();
    EndModify();

    return(*this);
  }

  bool OBMol::Clear()
  {
    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Clear Molecule", obAuditMsg);

    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;
    for (i = _vatom.begin();i != _vatom.end();i++)
      {
        DestroyAtom(*i);
        *i = NULL;
      }
    for (j = _vbond.begin();j != _vbond.end();j++)
      {
        DestroyBond(*j);
        *j = NULL;
      }

    _natoms = _nbonds = 0;

    //Delete residues
    unsigned int ii;
    for (ii=0 ; ii<_residue.size() ; ii++)
      {
        delete _residue[ii];
      }
    _residue.clear();

    //clear out the multiconformer data
    vector<double*>::iterator k;
    for (k = _vconf.begin();k != _vconf.end();k++)
      delete [] *k;
    _vconf.clear();

    if (!_vdata.empty()) //clean up generic data
      {
        vector<OBGenericData*>::iterator m;
        for (m = _vdata.begin();m != _vdata.end();m++)
	  delete *m;
        _vdata.clear();
      }

    _c = (double*) NULL;
    _flags = 0;
    _mod = 0;

    return(true);
  }

  void OBMol::BeginModify()
  {
    //suck coordinates from _c into _v for each atom
    if (!_mod && !Empty())
      {
        OBAtom *atom;
        vector<OBNodeBase*>::iterator i;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  {
            atom->SetVector();
            atom->ClearCoordPtr();
	  }

        vector<double*>::iterator j;
        for (j = _vconf.begin();j != _vconf.end();j++)
	  delete [] *j;

        _c = NULL;
        _vconf.clear();

        //Destroy rotamer list if necessary
        if ((OBRotamerList *)GetData(OBGenericDataType::RotamerList))
	  {
            delete (OBRotamerList *)GetData(OBGenericDataType::RotamerList);
            DeleteData(OBGenericDataType::RotamerList);
	  }
      }

    _mod++;
  }

  void OBMol::EndModify(bool nukePerceivedData)
  {
    if (_mod == 0)
      {
        obErrorLog.ThrowError(__FUNCTION__, "_mod is negative - EndModify() called too many times", obDebug);
        return;
      }

    _mod--;

    if (_mod)
      return;

    if (nukePerceivedData)
      _flags = 0;
    _c = NULL;

    /*
      leave generic data alone for now - just nuke it on clear()
      if (HasData("Comment")) delete [] (char*)GetData("Comment");
      _vdata.clear();
    */

    if (Empty())
      return;

    //if atoms present convert coords into array
    double *c = new double [NumAtoms()*3];
    _c = c;

    int idx;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;
    for (idx=0,atom = BeginAtom(j);atom;atom = NextAtom(j),idx++)
      {
        atom->SetIdx(idx+1);
        (atom->GetVector()).Get(&_c[idx*3]);
        atom->SetCoordPtr(&_c);
      }
    _vconf.push_back(c);

    //kekulize structure
    SetAromaticPerceived();
    Kekulize();
    //kekulize();
    UnsetAromaticPerceived();
		
    //    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
    //      atom->UnsetAromatic();

    //    OBBond *bond;
    //    vector<OBEdgeBase*>::iterator k;
    //    for (bond = BeginBond(k);bond;bond = NextBond(k))
    //      bond->UnsetAromatic();

    UnsetImplicitValencePerceived();
  }

  OBAtom *OBMol::CreateAtom(void)
  {
    return new OBAtom;
  }

  OBBond *OBMol::CreateBond(void)
  {
    return new OBBond;
  }

  void OBMol::DestroyAtom(OBNodeBase *atom)
  {
    if (atom)
      {
        delete atom;
        atom = NULL;
      }
  }

  void OBMol::DestroyBond(OBEdgeBase *bond)
  {
    if (bond)
      {
        delete bond;
        bond = NULL;
      }
  }

  //! \brief Get a new atom to add to a molecule
  //!
  //! Also checks bond_queue for any bonds that should be made to the new atom
  OBAtom *OBMol::NewAtom()
  {
    BeginModify();

    OBAtom *obatom = CreateAtom();
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100

    if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBNodeBase*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++)
	  *j = (OBNodeBase*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();i++)
	  if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
	      vb = (OBVirtualBond*)*i;
	      if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms)
		continue;
	      if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
		  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
		  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
		  verase.push_back(*i);
                }
            }

        if (!verase.empty())
	  DeleteData(verase);
      }

    EndModify();

    return(obatom);
  }

  OBResidue *OBMol::NewResidue()
  {
    OBResidue *obresidue = new OBResidue;
    obresidue->SetIdx(_residue.size());
    _residue.push_back(obresidue);
    return(obresidue);
  }

  //! \brief Add an atom to a molecule
  //!
  //! Also checks bond_queue for any bonds that should be made to the new atom
  bool OBMol::AddAtom(OBAtom &atom)
  {
    BeginModify();

    OBAtom *obatom = CreateAtom();
    *obatom = atom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100

    if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBNodeBase*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++)
	  *j = (OBNodeBase*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = (OBNodeBase*)obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();i++)
	  if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
	      vb = (OBVirtualBond*)*i;
	      if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms)
		continue;
	      if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
		  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
		  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
		  verase.push_back(*i);
                }
            }

        if (!verase.empty())
	  DeleteData(verase);
      }

    EndModify();

    return(true);
  }

  bool OBMol::InsertAtom(OBAtom &atom)
  {
    BeginModify();

    OBAtom *obatom = CreateAtom();
    *obatom = atom;
    obatom->SetIdx(_natoms+1);
    obatom->SetParent(this);


#define OBAtomIncrement 100

    if (_vatom.empty() || _natoms+1 >= (signed)_vatom.size())
      {
        _vatom.resize(_natoms+OBAtomIncrement);
        vector<OBNodeBase*>::iterator j;
        for (j = _vatom.begin(),j+=(_natoms+1);j != _vatom.end();j++)
	  *j = (OBNodeBase*)NULL;
      }
#undef OBAtomIncrement

    _vatom[_natoms] = (OBNodeBase*)obatom;
    _natoms++;

    if (HasData(OBGenericDataType::VirtualBondData))
      {
        /*add bonds that have been queued*/
        OBVirtualBond *vb;
        vector<OBGenericData*> verase;
        vector<OBGenericData*>::iterator i;
        for (i = BeginData();i != EndData();i++)
	  if ((*i)->GetDataType() == OBGenericDataType::VirtualBondData)
            {
	      vb = (OBVirtualBond*)*i;
	      if (vb->GetBgn() > _natoms || vb->GetEnd() > _natoms)
		continue;
	      if (obatom->GetIdx() == static_cast<unsigned int>(vb->GetBgn())
		  || obatom->GetIdx() == static_cast<unsigned int>(vb->GetEnd()))
                {
		  AddBond(vb->GetBgn(),vb->GetEnd(),vb->GetOrder());
		  verase.push_back(*i);
                }
            }

        if (!verase.empty())
	  DeleteData(verase);
      }

    EndModify();

    return(true);
  }

  bool OBMol::AddResidue(OBResidue &residue)
  {
    BeginModify();

    OBResidue *obresidue = new OBResidue;
    *obresidue = residue;

    obresidue->SetIdx(_residue.size());

    _residue.push_back(obresidue);

    EndModify();

    return(true);
  }

  bool OBMol::StripSalts()
  {
    vector<vector<int> > cfl;
    vector<vector<int> >::iterator i,max;

    ContigFragList(cfl);
    if (cfl.empty() || cfl.size() == 1)
      return(false);

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::StripSalts", obAuditMsg);

    max = cfl.begin();
    for (i = cfl.begin();i != cfl.end();i++)
      if ((*max).size() < (*i).size())
	max = i;

    vector<int>::iterator j;
    vector<OBNodeBase*> delatoms;
    for (i = cfl.begin();i != cfl.end();i++)
      if (i != max)
	for (j = (*i).begin();j != (*i).end();j++)
	  delatoms.push_back(GetAtom(*j));

    if (!delatoms.empty())
      {
        int tmpflags = _flags & (~(OB_SSSR_MOL));
        BeginModify();
        vector<OBNodeBase*>::iterator k;
        for (k = delatoms.begin();k != delatoms.end();k++)
	  DeleteAtom((OBAtom*)*k);
        EndModify();
        _flags = tmpflags;
      }

    return(true);
  }

  bool OBMol::DeleteNonPolarHydrogens()
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;
    vector<OBNodeBase*> delatoms;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::DeleteHydrogens -- nonpolar",
			  obAuditMsg);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->IsNonPolarHydrogen())
	delatoms.push_back(atom);

    if (delatoms.empty())
      return(true);

    /*
      
     
    int idx1,idx2;
    vector<double*>::iterator j;
    for (idx1=0,idx2=0,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
    if (!atom->IsHydrogen())
    {
    for (j = _vconf.begin();j != _vconf.end();j++)
    memcpy((char*)&((*j)[idx2*3]),(char*)&((*j)[idx1*3]),sizeof(double)*3);
    idx2++;
    }
    */

    IncrementMod();

    for (i = delatoms.begin();i != delatoms.end();i++)
      DeleteAtom((OBAtom *)*i);

    DecrementMod();

    return(true);
  }

  bool OBMol::DeleteHydrogens()
  {
    OBAtom *atom,*nbr;
    vector<OBNodeBase*>::iterator i;
    vector<OBNodeBase*> delatoms,va;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::DeleteHydrogens", obAuditMsg);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->IsHydrogen())
	delatoms.push_back(atom);

    if (delatoms.empty())
      return(true);

    /* decide whether these flags need to be reset
       _flags &= (~(OB_ATOMTYPES_MOL));
       _flags &= (~(OB_HYBRID_MOL));
       _flags &= (~(OB_PCHARGE_MOL)); 
       _flags &= (~(OB_IMPVAL_MOL));
    */

    //find bonds to delete
    vector<OBEdgeBase*> vdb;
    vector<OBEdgeBase*>::iterator j;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (!atom->IsHydrogen())
	for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
	  if (nbr->IsHydrogen())
	    vdb.push_back(*j);

    IncrementMod();
    for (j = vdb.begin();j != vdb.end();j++)
      DeleteBond((OBBond *)*j); //delete bonds
    DecrementMod();

    int idx1,idx2;
    vector<double*>::iterator k;
    for (idx1=0,idx2=0,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
      if (!atom->IsHydrogen())
        {
	  //Update conformer coordinates
	  for (k = _vconf.begin();k != _vconf.end();k++)
	    memcpy((char*)&((*k)[idx2*3]),(char*)&((*k)[idx1*3]),sizeof(double)*3);

	  idx2++;
	  va.push_back(atom);
        }

    for (i = delatoms.begin();i != delatoms.end();i++)
      {
        DestroyAtom(*i);
        _natoms--;
      }

    _vatom.clear();
    for (i = va.begin();i != va.end();i++)
      _vatom.push_back((OBNodeBase*)*i);

    //_atom = va;
    //_atom.resize(_atom.size()+1);
    //_atom[_atom.size()-1] = NULL;
    _natoms = va.size();

    //reset all the indices to the atoms
    for (idx1=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx1++)
      atom->SetIdx(idx1);

    return(true);
  }

  bool OBMol::DeleteHydrogens(OBAtom *atom)
    //deletes all hydrogens attached to the atom passed to the function
  {
    OBAtom *nbr;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator k;
    vector<OBNodeBase*> delatoms;

    for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
      if (nbr->IsHydrogen())
	delatoms.push_back(nbr);

    if (delatoms.empty())
      return(true);

    IncrementMod();
    for (i = delatoms.begin();i != delatoms.end();i++)
      DeleteHydrogen((OBAtom*)*i);
    DecrementMod();

    return(true);
  }


  bool OBMol::DeleteHydrogen(OBAtom *atom)
    //deletes the hydrogen atom passed to the function
  {
    //find bonds to delete
    OBAtom *nbr;
    vector<OBEdgeBase*> vdb;
    vector<OBEdgeBase*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    IncrementMod();
    for (j = vdb.begin();j != vdb.end();j++)
      DeleteBond((OBBond*)*j); //delete bonds
    DecrementMod();

    int idx;
    if (atom->GetIdx() != NumAtoms())
      {
        idx = atom->GetCIdx();
        int size = NumAtoms()-atom->GetIdx();
        vector<double*>::iterator k;
        for (k = _vconf.begin();k != _vconf.end();k++)
	  memmove((char*)&(*k)[idx],(char*)&(*k)[idx+3],sizeof(double)*3*size);

      }

    _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
    DestroyAtom(atom);
    _natoms--;

    //reset all the indices to the atoms
    vector<OBNodeBase*>::iterator i;
    for (idx=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx++)
      atom->SetIdx(idx);

    return(true);
  }

  bool OBMol::AddHydrogens(bool polaronly,bool correctForPH)
  {
    if (!IsCorrectedForPH() && correctForPH)
      CorrectForPH();

    if (HasHydrogensAdded())
      return(true);
    SetHydrogensAdded();

    if(NumBonds()==0 && NumAtoms()!=1)
      {
	obErrorLog.ThrowError(__FUNCTION__,
			      "Did not run OpenBabel::AddHydrogens on molecule with no bonds", obAuditMsg);
	return true;
      }
    if (!polaronly)
      obErrorLog.ThrowError(__FUNCTION__,
			    "Ran OpenBabel::AddHydrogens", obAuditMsg);
    else
      obErrorLog.ThrowError(__FUNCTION__,
			    "Ran OpenBabel::AddHydrogens -- polar only", obAuditMsg);

    //count up number of hydrogens to add
    OBAtom *atom,*h;
    int hcount,count=0;
    vector<pair<OBAtom*,int> > vhadd;
    vector<OBNodeBase*>::iterator i;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
	if (polaronly && !(atom->IsNitrogen() || atom->IsOxygen() ||
			   atom->IsSulfur() || atom->IsPhosphorus()))
	  continue;

	hcount = atom->GetImplicitValence() - atom->GetValence();

	//Jan 05 Implicit valency now left alone; use spin multiplicity for implicit Hs
	int mult = atom->GetSpinMultiplicity();
	if(mult==2) //radical
	  hcount-=1;
	else if(mult==1 || mult==3) //carbene
	  hcount-=2;
	else if(mult>=4) // as in CH, C etc
	  hcount -= mult-1;

	if (hcount < 0)
	  hcount = 0;
	if (hcount)
	  {
	    vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
	    count += hcount;
	  }
      }

    if (count == 0)
      return(true);
    bool hasCoords = HasNonZeroCoords();

    //realloc memory in coordinate arrays for new hydrogens
    double *tmpf;
    vector<double*>::iterator j;
    for (j = _vconf.begin();j != _vconf.end();j++)
      {
        tmpf = new double [(NumAtoms()+count)*3];
        memset(tmpf,'\0',sizeof(double)*(NumAtoms()+count)*3);
        if (hasCoords)
	  memcpy(tmpf,(*j),sizeof(double)*NumAtoms()*3);
        delete []*j;
        *j = tmpf;
      }

    IncrementMod();

    int m,n;
    vector3 v;
    vector<pair<OBAtom*,int> >::iterator k;
    double hbrad = etab.CorrectedBondRad(1,0);


    for (k = vhadd.begin();k != vhadd.end();k++)
      {
        atom = k->first;
        double bondlen = hbrad+etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
        for (m = 0;m < k->second;m++)
	  {
            for (n = 0;n < NumConformers();n++)
	      {
                SetConformer(n);
                if (hasCoords)
		  {
                    atom->GetNewBondVector(v,bondlen);
                    _c[(NumAtoms())*3]   = v.x();
                    _c[(NumAtoms())*3+1] = v.y();
                    _c[(NumAtoms())*3+2] = v.z();
		  }
                else
		  memset((char*)&_c[NumAtoms()*3],'\0',sizeof(double)*3);
	      }
            h = NewAtom();
            h->SetType("H");
            h->SetAtomicNum(1);

            // copy parent atom residue to added hydrogen     REG 6/30/02

            if (atom->HasResidue())
	      {

                string aname;

                aname = "H";

                // Add the new H atom to the appropriate residue list
                atom->GetResidue()->AddAtom(h);

                // Give the new atom a pointer back to the residue
                h->SetResidue(atom->GetResidue());

                atom->GetResidue()->SetAtomID(h,aname);

	      }

            AddBond(atom->GetIdx(),h->GetIdx(),1);
            h->SetCoordPtr(&_c);
	  }
      }

    DecrementMod();
    SetConformer(0);

    //reset atom type and partial charge flags
    _flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

    return(true);
  }

  bool OBMol::AddPolarHydrogens()
  {
    return(AddHydrogens(true));
  }

  bool OBMol::AddHydrogens(OBAtom *atom)
  {
    OBAtom *h;

    //count up number of hydrogens to add
    int hcount,count=0;
    vector<pair<OBAtom*,int> > vhadd;

    hcount = atom->GetImplicitValence() - atom->GetValence();

    //Jan 05 Implicit valency now left alone; use spin multiplicity for implicit Hs
    int mult = atom->GetSpinMultiplicity();
    if(mult==2) //radical
      hcount-=1;
    else if(mult==1 || mult==3) //carbene
      hcount-=2;

    if (hcount < 0)
      hcount = 0;
    if (hcount)
      {
        vhadd.push_back(pair<OBAtom*,int>(atom,hcount));
        count += hcount;
      }

    if (count == 0)
      return(true);

    //realloc memory in coordinate arrays for new hydroges
    double *tmpf;
    vector<double*>::iterator j;
    for (j = _vconf.begin();j != _vconf.end();j++)
      {
        tmpf = new double [(NumAtoms()+count)*3+10];
        memcpy(tmpf,(*j),sizeof(double)*NumAtoms()*3);
        delete []*j;
        *j = tmpf;
      }

    IncrementMod();

    int m,n;
    vector3 v;
    vector<pair<OBAtom*,int> >::iterator k;
    double hbrad = etab.CorrectedBondRad(1,0);

    for (k = vhadd.begin();k != vhadd.end();k++)
      {
        atom = k->first;
        double bondlen = hbrad+etab.CorrectedBondRad(atom->GetAtomicNum(),atom->GetHyb());
        for (m = 0;m < k->second;m++)
	  {
            for (n = 0;n < NumConformers();n++)
	      {
                SetConformer(n);
                atom->GetNewBondVector(v,bondlen);
                _c[(NumAtoms())*3]   = v.x();
                _c[(NumAtoms())*3+1] = v.y();
                _c[(NumAtoms())*3+2] = v.z();
	      }
            h = NewAtom();
            h->SetType("H");
            h->SetAtomicNum(1);
            AddBond(atom->GetIdx(),h->GetIdx(),1);
            h->SetCoordPtr(&_c);
	  }
      }

    DecrementMod();
    SetConformer(0);

    //reset atom type and partial charge flags
    //_flags &= (~(OB_PCHARGE_MOL|OB_ATOMTYPES_MOL));

    return(true);
  }

  bool OBMol::CorrectForPH()
  {
    if (IsCorrectedForPH())
      return(true);
    phmodel.CorrectForPH(*this);

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::CorrectForPH", obAuditMsg);

    return(true);
  }

  //! \brief set spin multiplicity for H-deficient atoms
  bool OBMol::AssignSpinMultiplicity()
  {
    if (HasSpinMultiplicityAssigned())
      return(true);

    SetSpinMultiplicityAssigned();

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::AssignSpinMultiplicity",
			  obAuditMsg);

    OBAtom *atom;
    int diff;
    vector<OBNodeBase*>::iterator k;
    //if there are any explicit Hs on an atom, except if their _isotope!=0,
    //then they consitute all the Hs.
    //Any discrepancy with the expected atom valency is because it is a radical of some sort.
    //So SMILES CC[2H] is interpreted as CH3CH2D; CC[H] is methyl carbene.
    for (atom = BeginAtom(k);atom;atom = NextAtom(k))
      {    
	if ((!atom->IsHydrogen() && atom->ExplicitHydrogenCount(true)!=0)//exclude D,T
	    || atom->HasNoHForced()) 
	  {
	    diff=atom->GetImplicitValence() - (atom->GetHvyValence() + atom->ExplicitHydrogenCount());
	    if (diff)
	      atom->SetSpinMultiplicity(diff+1);//radicals =2; all carbenes =3
	  }
      }

    vector<OBNodeBase*>::iterator i;
    unsigned int spin = 1;

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (atom->GetSpinMultiplicity() > 1)
	  spin += atom->GetSpinMultiplicity() - 1;
      }
    _totalSpin = spin;

    return (true);
  }


  // Not used anywhere internally -- likely predates OBBase code
  // static void ResetVisit(OBMol &mol,vector<int> &visit,int depth)
  // {
  //     OBBond *bond;
  //     vector<OBEdgeBase*>::iterator i;

  //     for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
  //         if (bond->IsAromatic() && visit[bond->GetIdx()] >= depth)
  //             visit[bond->GetIdx()] = 0;
  // }

  static int ValenceSum(OBAtom *atom)
  {
    int count = atom->GetImplicitValence();

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
      if (bond->IsKDouble())
	count++;

    return(count);
  }

  static bool KekulePropagate(OBAtom *atom,vector<int> &visit,vector<int> &ival,int depth)
  {
    int count = 0;
    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
      if (!visit[bond->GetIdx()])
	count++;

    if (!count)
      return(ValenceSum(atom) == ival[atom->GetIdx()]);

    bool result = true;
    OBAtom *nbr;

    if (ValenceSum(atom) >= ival[atom->GetIdx()])
      {
        for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
	  if (nbr->IsAromatic() && !visit[(*i)->GetIdx()])
            {
	      visit[(*i)->GetIdx()] = depth;
	      ((OBBond*)*i)->SetKSingle();
	      result = KekulePropagate(nbr,visit,ival,depth);
	      if (result)
		break;
	      //            if (!result) break;
            }
      }
    else if (count == 1)
      for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
	if (nbr->IsAromatic() && !visit[(*i)->GetIdx()])
	  {
	    visit[(*i)->GetIdx()] = depth;
	    ((OBBond*)*i)->SetKDouble();
	    result = KekulePropagate(nbr,visit,ival,depth);
	    //break;
	    if (result)
	      break;
	  }
    return(result);
  }

  int GetCurrentValence(OBAtom *atom)
  {
    int count = atom->GetImplicitValence();

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = atom->BeginBond(i);bond;bond = atom->NextBond(i))
      {
        if (bond->IsKDouble())
	  count++;
        else if (bond->IsKTriple())
	  count += 2;
        //      else if (bond->IsSingle()) count++;
        //      else if (bond->IsDouble()) count += 2;
        //      else if (bond->IsTriple()) count += 3;
      }
    return(count);
  }

  bool ExpandKekule(OBMol &mol, vector<OBNodeBase*> &va,
		    vector<OBNodeBase*>::iterator i,
		    vector<int> &maxv,bool secondpass)
  {
    if (i == va.end())
      {
        //check to see that the ideal valence has been achieved for all atoms
        vector<OBNodeBase*>::iterator j;
        for (j = va.begin();j != va.end();j++)
	  {
            //let erroneously aromatic carboxylates pass
            if (((OBAtom*)*j)->IsOxygen() && ((OBAtom*)*j)->GetValence() == 1)
	      continue;
            if (GetCurrentValence((OBAtom*)*j) != maxv[(*j)->GetIdx()])
	      {
                //	      cout << " ExpandKekule atom: " << ((OBAtom*)*j)->GetIdx()
                //		   << " valence is " << (GetCurrentValence((OBAtom*)*j))
                //		   << " should be " << maxv[(*j)->GetIdx()] << endl;
                return(false);
	      }
	  }
        return(true);
      }

    //jump to next atom in list if current atom doesn't have any attached
    //aromatic bonds
    OBBond *bond;
    OBAtom *atom = (OBAtom*)*i;
    vector<OBEdgeBase*>::iterator j;
    bool done = true;
    for (bond = atom->BeginBond(j);bond;bond = atom->NextBond(j))
      if (bond->GetBO() == 5)
        {
	  done = false;
	  break;
        }
    if (done)
      return(ExpandKekule(mol,va,i+1,maxv,secondpass));

    //store list of attached aromatic atoms
    OBAtom *nbr;
    vector<OBEdgeBase*> vb;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      if ((*j)->GetBO() == 5)
        {
	  vb.push_back(*j);
	  ((OBBond *)*j)->SetBO(1);
	  ((OBBond *)*j)->SetKSingle();
        }

    //try setting a double bond
    if (GetCurrentValence(atom) < maxv[atom->GetIdx()])
      {
        for (j = vb.begin();j != vb.end();j++)
	  {
            nbr = ((OBBond *)*j)->GetNbrAtom(atom);
            if (GetCurrentValence(nbr) <= maxv[nbr->GetIdx()])
	      {
                ((OBBond*)*j)->SetKDouble();
                ((OBBond*)*j)->SetBO(2);
                if (ExpandKekule(mol,va,i+1,maxv,secondpass))
		  return(true);
                ((OBBond*)*j)->SetKSingle();
                ((OBBond*)*j)->SetBO(1);
	      }
	  }

        if (secondpass && atom->IsNitrogen() && atom->GetFormalCharge() == 0 &&
	    atom->GetImplicitValence() == 2)
	  {
            atom->IncrementImplicitValence();
            if (ExpandKekule(mol,va,i+1,maxv,secondpass))
	      return(true);
            atom->DecrementImplicitValence();
	  }
      }
    else  //full valence - no double bond to set
      {
        if (ExpandKekule(mol,va,i+1,maxv,secondpass))
	  return(true);

        bool trycharge = false;
        if (secondpass && atom->GetFormalCharge() == 0)
	  {
            if (atom->IsNitrogen() && atom->GetHvyValence() == 3)
	      trycharge = true;
            if (atom->IsOxygen() && atom->GetHvyValence() == 2)
	      trycharge = true;
            if (atom->IsSulfur() && atom->GetHvyValence() == 2)
	      trycharge = true;
	  }

        if (trycharge) //attempt to charge up O,N,S to make a valid kekule form
	  {
            maxv[atom->GetIdx()]++;
            atom->SetFormalCharge(1);
            for (j = vb.begin();j != vb.end();j++)
	      {
                nbr = ((OBBond*)*j)->GetNbrAtom(atom);
                if (GetCurrentValence(nbr) <= maxv[nbr->GetIdx()])
		  {
                    ((OBBond*)*j)->SetKDouble();
                    ((OBBond*)*j)->SetBO(2);
                    if (ExpandKekule(mol,va,i+1,maxv,secondpass))
		      return(true);
                    ((OBBond*)*j)->SetKSingle();
                    ((OBBond*)*j)->SetBO(1);
		  }
	      }
            maxv[atom->GetIdx()]--;
            atom->SetFormalCharge(0);
	  }

        if (secondpass && atom->IsNitrogen() && atom->GetFormalCharge() == 0 &&
	    atom->GetImplicitValence() == 2) //try protonating the nitrogen
	  {
            atom->IncrementImplicitValence();
            if (ExpandKekule(mol,va,i+1,maxv,secondpass))
	      return(true);
            atom->DecrementImplicitValence();
	  }
      }

    //failed to find a valid solution - reset attached bonds
    for (j = vb.begin();j != vb.end();j++)
      {
        ((OBBond*)*j)->SetKSingle();
        ((OBBond*)*j)->SetBO(5);
      }

    return(false);
  }

  void CorrectBadResonanceForm(OBMol &mol)
  {
    string s;
    OBBond *b1,*b2,*b3;
    OBAtom *a1,*a2,*a3,*a4;
    vector<vector<int> > mlist;
    vector<vector<int> >::iterator i;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::CorrectBadResonanceForm", obAuditMsg);

    OBSmartsPattern acid;
    acid.Init("[oD1]c[oD1]");

    //carboxylic acid
    if (acid.Match(mol))
      {
        mlist = acid.GetUMapList();
        for (i = mlist.begin();i != mlist.end();i++)
	  {
            a1 = mol.GetAtom((*i)[0]);
            a2 = mol.GetAtom((*i)[1]);
            a3 = mol.GetAtom((*i)[2]);
            b1 = a2->GetBond(a1);
            b2 = a2->GetBond(a3);
            if (!b1 || !b2)
	      continue;
            b1->SetKDouble();
            b2->SetKSingle();
	  }
      }

    //phosphonic acid
    OBSmartsPattern phosphate;
    phosphate.Init("[p]([oD1])([oD1])([oD1])[#6,#8]");
    if (phosphate.Match(mol))
      {
        mlist = phosphate.GetUMapList();
        for (i = mlist.begin();i != mlist.end();i++)
	  {
            a1 = mol.GetAtom((*i)[0]);
            a2 = mol.GetAtom((*i)[1]);
            a3 = mol.GetAtom((*i)[2]);
            a4 = mol.GetAtom((*i)[3]);
            b1 = a1->GetBond(a2);
            b2 = a1->GetBond(a3);
            b3 = a1->GetBond(a4);

            if (!b1 || !b2 || !b3)
	      continue;
            b1->SetKDouble();
            b2->SetKSingle();
            b3->SetKSingle();
	  }
      }

    //amidene and guanidine
    OBSmartsPattern amidene;
    amidene.Init("[nD1]c([nD1])*");
    if (amidene.Match(mol))
      {
        mlist = amidene.GetUMapList();
        for (i = mlist.begin();i != mlist.end();i++)
	  {
            a1 = mol.GetAtom((*i)[0]);
            a2 = mol.GetAtom((*i)[1]);
            a3 = mol.GetAtom((*i)[2]);
            b1 = a2->GetBond(a1);
            b2 = a2->GetBond(a3);
            if (!b1 || !b2)
	      continue;
            b1->SetKDouble();
            b2->SetKSingle();
	  }
      }
  }

  bool OBMol::PerceiveKekuleBonds()
  {
    if (HasKekulePerceived())
      return(true);
    SetKekulePerceived();

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;

    //initialize kekule bonds
    bool done = true;
    bool badResonanceForm = false;
    vector<bool> varo;
    varo.resize(NumAtoms()+1,false);
    for (bond = BeginBond(i);bond;bond = NextBond(i))
      switch (bond->GetBO())
        {
        case 2:
	  bond->SetKDouble();
	  break;
        case 3:
	  bond->SetKTriple();
	  break;
        case 5:

	  bond->SetKSingle();
	  if (bond->IsInRing())
            {
	      varo[bond->GetBeginAtomIdx()] = true;
	      varo[bond->GetEndAtomIdx()]   = true;
	      done = false;
            }
	  else
	    badResonanceForm = true;

	  break;

        default:
	  bond->SetKSingle();
	  break;
        }

    if (badResonanceForm)
      CorrectBadResonanceForm(*this);

    if (done)
      return(true);

    //set the maximum valence for each aromatic atom
    OBAtom *atom,*nbr;
    vector<OBNodeBase*>::iterator j,k;
    vector<int> maxv;
    maxv.resize(NumAtoms()+1);

    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      if (varo[atom->GetIdx()])
        {
	  switch (atom->GetAtomicNum())
            {
            case 6:
	      maxv[atom->GetIdx()] = 4;
	      break;
            case 8:
            case 16:
            case 34:
            case 52:
	      maxv[atom->GetIdx()] = 2;
	      break;
            case 7:
            case 15:
            case 33:
	      maxv[atom->GetIdx()] = 3;
	      break;
            }
	  //correct valence for formal charges
	  if (atom->IsCarbon())
	    maxv[atom->GetIdx()] -= abs(atom->GetFormalCharge());
	  else
	    maxv[atom->GetIdx()] += atom->GetFormalCharge();

	  if (atom->IsNitrogen() || atom->IsSulfur())
	    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
	      if (nbr->IsOxygen() && (*i)->GetBO() == 2)
		maxv[atom->GetIdx()] += 2;
        }

    bool result = true;
    vector<bool> used;
    used.resize(NumAtoms()+1);
    vector<OBNodeBase*> va,curr,next;
    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      if (varo[atom->GetIdx()] && !used[atom->GetIdx()])
        {
	  va.clear();
	  va.push_back(atom);
	  curr.clear();
	  curr.push_back(atom);
	  used[atom->GetIdx()] = true;

	  for (;!curr.empty();)
            {
	      next.clear();
	      for (k = curr.begin();k != curr.end();k++)
		for (nbr = ((OBAtom*)*k)->BeginNbrAtom(i);nbr;nbr = ((OBAtom*)*k)->NextNbrAtom(i))
		  if (varo[nbr->GetIdx()] && !used[nbr->GetIdx()])
		    {
		      used[nbr->GetIdx()] = true;
		      next.push_back(nbr);
		      va.push_back(nbr);
		    }
	      curr = next;
            }

	  //try it first without protonating aromatic nitrogens
	  if (!ExpandKekule(*this,va,va.begin(),maxv,false) &&
	      !ExpandKekule(*this,va,va.begin(),maxv,true))
            {
	      result = false;
	      //	      cerr << " Died on atom " << atom->GetIdx() << endl;
            }
        }

    if (!result)
      {
	//        cerr << "Kekulization Error = " << GetTitle() << endl;
        //exit(0);
      }

    return(result);
  }

  bool OBMol::Kekulize()
  {
    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    // Not quite sure why this is here -GRH 2003
    //  if (NumAtoms() > 255) return(false);

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Kekulize", obAuditMsg);

    for (bond = BeginBond(i);bond;bond = NextBond(i))
      if (bond->IsKSingle())
	bond->SetBO(1);
      else if (bond->IsKDouble())
	bond->SetBO(2);
      else if (bond->IsKTriple())
	bond->SetBO(3);

    return(true);
  }

  bool OBMol::DeleteAtom(OBAtom *atom)
  {
    if (atom->IsHydrogen())
      return(DeleteHydrogen(atom));

    BeginModify();
    //don't need to do anything with coordinates b/c
    //BeginModify() blows away coordinates

    //find bonds to delete
    OBAtom *nbr;
    vector<OBEdgeBase*> vdb;
    vector<OBEdgeBase*>::iterator j;
    for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
      vdb.push_back(*j);

    for (j = vdb.begin();j != vdb.end();j++)
      DeleteBond((OBBond *)*j); //delete bonds

    _vatom.erase(_vatom.begin()+(atom->GetIdx()-1));
    DestroyAtom(atom);
    _natoms--;

    //reset all the indices to the atoms
    int idx;
    vector<OBNodeBase*>::iterator i;
    for (idx=1,atom = BeginAtom(i);atom;atom = NextAtom(i),idx++)
      atom->SetIdx(idx);

    EndModify();

    return(true);
  }

  bool OBMol::DeleteResidue(OBResidue *residue)
  {
    unsigned short idx = residue->GetIdx();
    for ( unsigned short i = idx ; i < _residue.size() ; i++ )
      _residue[i]->SetIdx(i-1);

    _residue.erase(_residue.begin() + idx);

    if (residue)
      {
        delete residue;
        residue = NULL;
      }

    return(true);
  }

  bool OBMol::AddBond(int first,int second,int order,int stereo,int insertpos)
  {
    if (first == second)
      return(false);

    BeginModify();

    if ((unsigned)first <= NumAtoms() && (unsigned)second <= NumAtoms()
	&& !GetBond(first, second))
      //atoms exist and bond doesn't
      {
        OBBond *bond = CreateBond();
        if (!bond)
	  {
            EndModify();
            return(false);
	  }

        OBAtom *bgn,*end;
        bgn = GetAtom(first);
        end = GetAtom(second);
        if (!bgn || !end)
	  {
            obErrorLog.ThrowError(__FUNCTION__, "Unable to add bond - invalid atom index", obDebug);
            return(false);
	  }
        bond->Set(_nbonds,bgn,end,order,stereo);
        bond->SetParent(this);

        //set aromatic flags if it has the appropriate order
        if (order == 5)
	  {
            bond->SetAromatic();
            bgn->SetAromatic();
            end->SetAromatic();
	  }

#define OBBondIncrement 100
        if (_vbond.empty() || _nbonds+1 >= (signed)_vbond.size())
	  {
            _vbond.resize(_nbonds+OBBondIncrement);
            vector<OBEdgeBase*>::iterator i;
            for (i = _vbond.begin(),i+=(_nbonds+1);i != _vbond.end();i++)
	      *i = (OBEdgeBase*)NULL;
	  }
#undef  OBBondIncrement

        _vbond[_nbonds] = (OBEdgeBase*)bond;
        _nbonds++;

        if (insertpos == -1)
	  {
            bgn->AddBond(bond);
            end->AddBond(bond);
	  }
        else
	  {
            if (insertpos >= static_cast<int>(bgn->GetValence())
		) bgn->AddBond(bond);
            else //need to insert the bond for the connectivity order to be preserved
	      {    //otherwise stereochemistry gets screwed up
                vector<OBEdgeBase*>::iterator bi;
                bgn->BeginNbrAtom(bi);
                bi += insertpos;
                bgn->InsertBond(bi,bond);
	      }
            end->AddBond(bond);
	  }
      }
    else //at least one atom doesn't exist yet - add to bond_q
      SetData(new OBVirtualBond(first,second,order,stereo));

    EndModify();
    return(true);
  }

  bool OBMol::AddBond(OBBond &bond)
  {
    return(AddBond(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond.GetBO(),
                   bond.GetFlags()));
  }

  bool OBMol::DeleteBond(OBBond *bond)
  {
    BeginModify();

    (bond->GetBeginAtom())->DeleteBond(bond);
    (bond->GetEndAtom())->DeleteBond(bond);
    _vbond.erase(_vbond.begin() + bond->GetIdx());

    DestroyBond(bond);

    vector<OBEdgeBase*>::iterator i;
    int j;
    for (bond = BeginBond(i),j=0;bond;bond = NextBond(i),j++)
      bond->SetIdx(j);

    _nbonds--;
    EndModify();
    return(true);
  }

  void OBMol::Align(OBAtom *a1,OBAtom *a2,vector3 &p1,vector3 &p2)
  {
    vector<int> children;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::Align", obAuditMsg);

    //find which atoms to rotate
    FindChildren(children,a1->GetIdx(),a2->GetIdx());
    children.push_back(a2->GetIdx());

    //find the rotation vector and angle
    vector3 v1,v2,v3;
    v1 = p2 - p1;
    v2 = a2->GetVector() - a1->GetVector();
    v3 = cross(v1,v2);
    double angle = vectorAngle(v1,v2);

    //find the rotation matrix
    matrix3x3 m;
    m.RotAboutAxisByAngle(v3,angle);

    //rotate atoms
    vector3 v;
    OBAtom *atom;
    vector<int>::iterator i;
    for (i = children.begin();i != children.end();i++)
      {
        atom = GetAtom(*i);
        v = atom->GetVector();
        v -= a1->GetVector();
        v *= m;   //rotate the point
        v += p1;  //translate the vector
        atom->SetVector(v);
      }
    //set a1 = p1
    a1->SetVector(p1);
  }

  void OBMol::ToInertialFrame()
  {
    double m[9];
    for (int i = 0;i < NumConformers();i++)
      ToInertialFrame(i,m);
  }

  void OBMol::ToInertialFrame(int conf,double *rmat)
  {
    unsigned int i;
    double x,y,z;
    double mi;
    double mass = 0.0;
    double center[3],m[3][3];

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::ToInertialFrame", obAuditMsg);

    for (i = 0;i < 3;i++)
      memset(&m[i],'\0',sizeof(double)*3);
    memset(center,'\0',sizeof(double)*3);

    SetConformer(conf);
    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;
    //find center of mass
    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      {
        mi = atom->GetAtomicMass();
        center[0] += mi*atom->x();
        center[1] += mi*atom->y();
        center[2] += mi*atom->z();
        mass += mi;
      }

    center[0] /= mass;
    center[1] /= mass;
    center[2] /= mass;

    //calculate inertial tensor
    for (atom = BeginAtom(j);atom;atom = NextAtom(j))
      {
        x = atom->x()-center[0];
        y = atom->y()-center[1];
        z = atom->z()-center[2];
        mi = atom->GetAtomicMass();

        m[0][0] += mi*(y*y+z*z);
        m[0][1] -= mi*x*y;
        m[0][2] -= mi*x*z;
        m[1][0] -= mi*x*y;
        m[1][1] += mi*(x*x+z*z);
        m[1][2] -= mi*y*z;
        m[2][0] -= mi*x*z;
        m[2][1] -= mi*y*z;
        m[2][2] += mi*(x*x+y*y);
      }

    /* find rotation matrix for moment of inertia */
    ob_make_rmat(m,rmat);

    /* rotate all coordinates */
    double *c = GetConformer(conf);
    for(i=0; i < NumAtoms();i++)
      {
        x = c[i*3]-center[0];
        y = c[i*3+1]-center[1];
        z = c[i*3+2]-center[2];
        c[i*3]   = x*rmat[0] + y*rmat[1] + z*rmat[2];
        c[i*3+1] = x*rmat[3] + y*rmat[4] + z*rmat[5];
        c[i*3+2] = x*rmat[6] + y*rmat[7] + z*rmat[8];
      }
  }

  /*NF
    istream& operator>> (istream &ifs, OBMol &mol)
    {
    bool retcode = OBFileFormat::ReadMolecule(ifs, mol);
 
    if (!retcode)
    {
    if (mol.GetMod())
    mol.EndModify();
    mol.Clear();
    }
 
    return(ifs);
    }
 
    ostream& operator<< (ostream &ofs, OBMol &mol)
    {
    OBFileFormat::WriteMolecule(ofs, mol);
    return(ofs);
    }
  */

  OBMol::OBMol()
  {
    _natoms = _nbonds = 0;
    _mod = 0;
    _energy = 0.0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _vbond.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
  }

  OBMol::OBMol(const OBMol &mol) :
    OBGraphBase()
  {
    _natoms = _nbonds = 0;
    _mod = 0;
    _totalCharge = 0;
    _dimension = 3;
    _vatom.clear();
    _vbond.clear();
    _vdata.clear();
    _title = "";
    _c = (double*)NULL;
    _flags = 0;
    _vconf.clear();
    _autoPartialCharge = true;
    _autoFormalCharge = true;
    //NF  _compressed = false;
    *this = mol;
  }

  OBMol::~OBMol()
  {
    OBAtom    *atom;
    OBBond    *bond;
    OBResidue *residue;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;
    vector<OBResidue*>::iterator r;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      DestroyAtom(atom);
    for (bond = BeginBond(j);bond;bond = NextBond(j))
      DestroyBond(bond);
    for (residue = BeginResidue(r);residue;residue = NextResidue(r))
      delete residue;

    //clear out the multiconformer data
    vector<double*>::iterator k;
    for (k = _vconf.begin();k != _vconf.end();k++)
      delete [] *k;
    _vconf.clear();
  }

  bool OBMol::HasNonZeroCoords()
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (atom->GetVector() != VZero)
	return(true);

    return(false);
  }

  bool OBMol::Has2D()
  {
    bool hasX,hasY;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    hasX = hasY = false;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (!hasX && !IsNearZero(atom->x()))
	  hasX = true;
        if (!hasY && !IsNearZero(atom->y()))
	  hasY = true;

        if (hasX && hasY)
	  return(true);
      }
    return(false);
  }

  bool OBMol::Has3D()
  {
    bool hasX,hasY,hasZ;
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    hasX = hasY = hasZ = false;
    if (this->_c == NULL)
      return(false);
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (!hasX && !IsNearZero(atom->x()))
	  hasX = true;
        if (!hasY && !IsNearZero(atom->y()))
	  hasY = true;
        if (!hasZ && !IsNearZero(atom->z()))
	  hasZ = true;

        if (hasX && hasY && hasZ)
	  return(true);
      }
    return(false);
  }

  bool OBMol::IsChiral()
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if ((atom->IsCarbon() || atom->IsNitrogen()) && atom->GetHvyValence() > 2 && atom->IsChiral())
	return(true);

    return(false);
  }

  //! Renumber the atoms in this molecule according to the order in the supplied
  //! vector. This will return without action if the supplied vector is empty or
  //! does not have the same number of atoms as the molecule.
  void OBMol::RenumberAtoms(vector<OBNodeBase*> &v)
  {
    if (Empty())
      return;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::RenumberAtoms", obAuditMsg);

    OBAtom *atom;
    vector<OBNodeBase*> va;
    vector<OBNodeBase*>::iterator i;

    va = v;

    //make sure all atoms are represented in the vector
    if (va.empty() || va.size() != NumAtoms())
      return;

    OBBitVec bv;
    for (i = va.begin();i != va.end();i++)
      bv |= (*i)->GetIdx();
    
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      if (!bv[atom->GetIdx()])
	va.push_back(atom);

    int j,k;
    double *c;
    double *ctmp = new double [NumAtoms()*3];

    for (j = 0;j < NumConformers();j++)
      {
        c = GetConformer(j);
        for (k=0,i = va.begin();i != va.end();i++,k++)
	  memcpy((char*)&ctmp[k*3],(char*)&c[((OBAtom*)*i)->GetCIdx()],sizeof(double)*3);
        memcpy((char*)c,(char*)ctmp,sizeof(double)*3*NumAtoms());
      }

    for (k=1,i = va.begin();i != va.end();i++,k++)
      (*i)->SetIdx(k);

    delete [] ctmp;

    _vatom.clear();
    for (i = va.begin();i != va.end();i++)
      _vatom.push_back(*i);
  }

  bool WriteTitles(ostream &ofs, OBMol &mol)
  {
    ofs << mol.GetTitle() << endl;
    return true;
  }

  /*! This method adds single bonds between all atoms
    closer than their combined atomic covalent radii,
    then "cleans up" making sure bonded atoms are not
    closer than 0.4A and the atom does not exceed its valence.
    It implements blue-obelisk:rebondFrom3DCoordinates.
  
  */
  void OBMol::ConnectTheDots(void)
  {
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::ConnectTheDots", obAuditMsg);

    int j,k,max;
    bool unset = false;
    OBAtom *atom,*nbr;
    vector<OBNodeBase*>::iterator i;
    vector<pair<OBAtom*,double> > zsortedAtoms;
    vector<double> rad;
    vector<int> zsorted;

    double *c = new double [NumAtoms()*3];
    rad.resize(_natoms);

    for (j = 0, atom = BeginAtom(i) ; atom ; atom = NextAtom(i), j++)
      {
        (atom->GetVector()).Get(&c[j*3]);
        pair<OBAtom*,double> entry(atom, atom->GetVector().z());
        zsortedAtoms.push_back(entry);
      }
    sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);

    max = zsortedAtoms.size();

    for ( j = 0 ; j < max ; j++ )
      {
        atom   = zsortedAtoms[j].first;
        rad[j] = etab.GetCovalentRad(atom->GetAtomicNum());
        zsorted.push_back(atom->GetIdx()-1);
      }

    int idx1, idx2;
    double d2,cutoff,zd;
    for (j = 0 ; j < max ; j++)
      {
        idx1 = zsorted[j];
        for (k = j + 1 ; k < max ; k++ )
	  {
            idx2 = zsorted[k];

            // bonded if closer than elemental Rcov + tolerance
            cutoff = SQUARE(rad[j] + rad[k] + 0.45);

            zd  = SQUARE(c[idx1*3+2] - c[idx2*3+2]);
            if (zd > 25.0 )
	      break; // bigger than max cutoff

            d2  = SQUARE(c[idx1*3]   - c[idx2*3]);
            d2 += SQUARE(c[idx1*3+1] - c[idx2*3+1]);
            d2 += zd;

            if (d2 > cutoff)
	      continue;
            if (d2 < 0.40)
	      continue;

            atom = GetAtom(idx1+1);
            nbr  = GetAtom(idx2+1);

            if (atom->IsConnected(nbr))
	      continue;
            if (atom->IsHydrogen() && nbr->IsHydrogen())
	      continue;

            AddBond(idx1+1,idx2+1,1);
	  }
      }

    // If between BeginModify and EndModify, coord pointers are NULL
    // setup molecule to handle current coordinates

    if (_c == NULL)
      {
        _c = c;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  atom->SetCoordPtr(&_c);
        _vconf.push_back(c);
        unset = true;
      }

    // Cleanup -- delete long bonds that exceed max valence
    OBBond *maxbond, *bond;
    double maxlength;
    vector<OBEdgeBase*>::iterator l;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        while (atom->BOSum() > static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum()))
	       || atom->SmallestBondAngle() < 45.0)
	  {
            maxbond = atom->BeginBond(l);
            maxlength = maxbond->GetLength();
            for (bond = atom->BeginBond(l);bond;bond = atom->NextBond(l))
	      {
                if (bond->GetLength() > maxlength)
		  {
                    maxbond = bond;
                    maxlength = bond->GetLength();
		  }
	      }
            DeleteBond(maxbond);
	  }
      }

    if (unset)
      {
        _c = NULL;
        for (atom = BeginAtom(i);atom;atom = NextAtom(i))
	  atom->ClearCoordPtr();
        _vconf.resize(_vconf.size()-1);
      }

    delete [] c;
  }

  /*! This method uses bond angles and geometries from current
    connectivity to guess atom types and then filling empty valences
    with multiple bonds. It currently has a pass to detect some
    frequent functional groups. It still needs a pass to detect aromatic
    rings to "clean up." */
  void OBMol::PerceiveBondOrders()
  {
    if (Empty())
      return;
    if (_dimension != 3) return; // not useful on non-3D structures

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::PerceiveBondOrders", obAuditMsg);

    OBAtom *atom, *b, *c;
    vector3 v1, v2;
    double angle;//, dist1, dist2;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;//,k;

    //  BeginModify();

    // Pass 1: Assign estimated hybridization based on avg. angles
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        angle = atom->AverageBondAngle();

        if (angle > 155.0)
	  atom->SetHyb(1);
        else if ( angle <= 155.0 && angle > 115)
	  atom->SetHyb(2);
      } // pass 1

    // Make sure upcoming calls to GetHyb() don't kill these temporary values
    SetHybridizationPerceived();

    // Pass 2: look for 5-member rings with torsions <= 7.5 degrees
    //         and 6-member rings with torsions <= 12 degrees
    //         (set all atoms with at least two bonds to sp2)

    vector<OBRing*> rlist;
    vector<OBRing*>::iterator ringit;
    vector<int> path;
    double torsions = 0.0;

    if (!HasSSSRPerceived())
      FindSSSR();
    rlist = GetSSSR();
    for (ringit = rlist.begin(); ringit != rlist.end(); ringit++)
      {
        if ((*ringit)->PathSize() == 5)
	  {
            path = (*ringit)->_path;
            torsions =
	      ( fabs(GetTorsion(path[0], path[1], path[2], path[3])) +
		fabs(GetTorsion(path[1], path[2], path[3], path[4])) +
		fabs(GetTorsion(path[2], path[3], path[4], path[0])) +
		fabs(GetTorsion(path[3], path[4], path[0], path[1])) +
		fabs(GetTorsion(path[4], path[0], path[1], path[2])) ) / 5.0;
            if (torsions <= 7.5)
	      {
                for (unsigned int ringAtom = 0; ringAtom != path.size(); ringAtom++)
		  {
                    b = GetAtom(path[ringAtom]);
		    // if an aromatic ring atom has valence 3, it is already set
		    // to sp2 because the average angles should be 120 anyway
		    // so only look for valence 2
                    if (b->GetValence() == 2)
		      b->SetHyb(2);
		  }
	      }
	  }
        else if ((*ringit)->PathSize() == 6)
	  {
            path = (*ringit)->_path;
            torsions =
	      ( fabs(GetTorsion(path[0], path[1], path[2], path[3])) +
		fabs(GetTorsion(path[1], path[2], path[3], path[4])) +
		fabs(GetTorsion(path[2], path[3], path[4], path[5])) +
		fabs(GetTorsion(path[3], path[4], path[5], path[0])) +
		fabs(GetTorsion(path[4], path[5], path[0], path[1])) +
		fabs(GetTorsion(path[5], path[0], path[1], path[2])) ) / 6.0;
            if (torsions <= 12.0)
	      {
                for (unsigned int ringAtom = 0; ringAtom != path.size(); ringAtom++)
		  {
                    b = GetAtom(path[ringAtom]);
                    if (b->GetValence() == 2 || b->GetValence() == 3)
		      b->SetHyb(2);
		  }
	      }
	  }
      }

    // Pass 3: "Antialiasing" If an atom marked as sp hybrid isn't
    //          bonded to another or an sp2 hybrid isn't bonded
    //          to another (or terminal atoms in both cases)
    //          mark them to a lower hybridization for now
    bool openNbr;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        if (atom->GetHyb() == 2 || atom->GetHyb() == 1)
	  {
            openNbr = false;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	      {
                if (b->GetHyb() < 3 || b->GetValence() == 1)
		  {
                    openNbr = true;
                    break;
		  }
	      }
            if (!openNbr && atom->GetHyb() == 2)
	      atom->SetHyb(3);
            else if (!openNbr && atom->GetHyb() == 1)
	      atom->SetHyb(2);
	  }
      } // pass 3

    // Pass 4: Check for known functional group patterns and assign bonds
    //         to the canonical form
    //      Currently we have explicit code to do this, but a "bond typer"
    //      is in progress to make it simpler to test and debug.
    bondtyper.AssignFunctionalGroupBonds(*this);

    // Pass 5: Check for aromatic rings and assign bonds as appropriate
    // This is just a quick and dirty approximation that marks everything
    //  as potentially aromatic

    // This doesn't work perfectly, but it's pretty decent.
    //  Need to have a list of SMARTS patterns for common rings
    //  which would "break ties" on complicated multi-ring systems
    // (Most of the current problems lie in the interface with the
    //   Kekulize code anyway, not in marking everything as potentially aromatic)

    bool typed; // has this ring been typed?
    unsigned int loop, loopSize;
    for (ringit = rlist.begin(); ringit != rlist.end(); ringit++)
      {
        typed = false;
        loopSize = (*ringit)->PathSize();
        if (loopSize == 5 || loopSize == 6)
	  {
            path = (*ringit)->_path;
            for(loop = 0; loop < loopSize; loop++)
	      {
                atom = GetAtom(path[loop]);
                if(atom->HasBondOfOrder(2) || atom->HasBondOfOrder(3)
		   || atom->GetHyb() != 2)
		  {
                    typed = true;
                    break;
		  }
	      }

            if (!typed)
	      for(loop = 0; loop < loopSize; loop++)
                {
		  //		cout << " set aromatic " << path[loop] << endl;
		  (GetBond(path[loop], path[(loop+1) % loopSize]))->SetBO(5);
		  (GetBond(path[loop], path[(loop+1) % loopSize]))->UnsetKekule();
                }
	  }
      }
    _flags &= (~(OB_KEKULE_MOL));
    Kekulize();

    // Pass 6: Assign remaining bond types, ordered by atom electronegativity
    vector<pair<OBAtom*,double> > sortedAtoms;
    vector<double> rad;
    vector<int> sorted;
    int iter, max;
    double maxElNeg, shortestBond, currentElNeg;

    for (atom = BeginAtom(i) ; atom ; atom = NextAtom(i))
      {
	// if atoms have the same electronegativity, make sure those with shorter bonds
	// are handled first (helps with assignment of conjugated single/double bonds)
	shortestBond = 1.0e5f;
	for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	  {
	    if (b->GetAtomicNum()!=1) shortestBond =
					min(shortestBond,(atom->GetBond(b))->GetLength());
	  }
	pair<OBAtom*,double> entry(atom,
				   etab.GetElectroNeg(atom->GetAtomicNum())*1e6+shortestBond);

	sortedAtoms.push_back(entry);
      }
    sort(sortedAtoms.begin(), sortedAtoms.end(), SortAtomZ);

    max = sortedAtoms.size();

    for (iter = 0 ; iter < max ; iter++ )
      {
        atom = sortedAtoms[iter].first;
        if ( (atom->GetHyb() == 1 || atom->GetValence() == 1)
	     && atom->BOSum() + 2  <= static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum()))
	     )
	  {
            // loop through the neighbors looking for a hybrid or terminal atom
            // (and pick the one with highest electronegativity first)
            // *or* pick a neighbor that's a terminal atom

            if (atom->HasNonSingleBond() ||
		(atom->GetAtomicNum() == 7 && atom->BOSum() + 2 > 3))
	      continue;

            maxElNeg = 0.0;
            shortestBond = 5000.0;
            c = NULL;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	      {
                currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
                if ( (b->GetHyb() == 1 || b->GetValence() == 1)
		     && b->BOSum() + 2 <= static_cast<unsigned int>(etab.GetMaxBonds(b->GetAtomicNum()))
		     && (currentElNeg > maxElNeg ||
			 (IsNear(currentElNeg,maxElNeg)
			  && (atom->GetBond(b))->GetLength() < shortestBond)) )
		  {
                    if (b->HasNonSingleBond() ||
			(b->GetAtomicNum() == 7 && b->BOSum() + 2 > 3))
		      continue;

                    shortestBond = (atom->GetBond(b))->GetLength();
                    maxElNeg = etab.GetElectroNeg(b->GetAtomicNum());
                    c = b; // save this atom for later use
		  }
	      }
            if (c)
	      (atom->GetBond(c))->SetBO(3);
	  }
        else if ( (atom->GetHyb() == 2 || atom->GetValence() == 1)
                  && atom->BOSum() + 1 <= static_cast<unsigned int>(etab.GetMaxBonds(atom->GetAtomicNum())) )
	  {
            // as above
            if (atom->HasNonSingleBond() ||
		(atom->GetAtomicNum() == 7 && atom->BOSum() + 1 > 3))
	      continue;

            maxElNeg = 0.0;
            shortestBond = 5000.0;
            c = NULL;
            for (b = atom->BeginNbrAtom(j); b; b = atom->NextNbrAtom(j))
	      {
                currentElNeg = etab.GetElectroNeg(b->GetAtomicNum());
                if ( (b->GetHyb() == 2 || b->GetValence() == 1)
		     && b->BOSum() + 1 <= static_cast<unsigned int>(etab.GetMaxBonds(b->GetAtomicNum()))
		     && (GetBond(atom, b))->IsDoubleBondGeometry()
		     && (currentElNeg > maxElNeg ||
			 (IsNear(currentElNeg,maxElNeg)
			  // If only the bond length counts, prefer double bonds in the ring
			  && (((atom->GetBond(b))->GetLength() < shortestBond) 
			      && (!atom->IsInRing() || !c || !c->IsInRing() || b->IsInRing()))
			  || (atom->IsInRing() && c && !c->IsInRing() && b->IsInRing()))))
		  {
                    if (b->HasNonSingleBond() ||
			(b->GetAtomicNum() == 7 && b->BOSum() + 1 > 3))
		      continue;

                    shortestBond = (atom->GetBond(b))->GetLength();
                    maxElNeg = etab.GetElectroNeg(b->GetAtomicNum());
                    c = b; // save this atom for later use
		  }
	      }
            if (c)
	      (atom->GetBond(c))->SetBO(2);
	  }
      } // pass 6

    // Now let the atom typer go to work again
    _flags &= (~(OB_HYBRID_MOL));
    _flags &= (~(OB_KEKULE_MOL));
    _flags &= (~(OB_AROMATIC_MOL));
    _flags &= (~(OB_ATOMTYPES_MOL));
    _flags &= (~(OB_IMPVAL_MOL));
    //  EndModify(true); // "nuke" perceived data
  }

  void OBMol::Center()
  {
    int j,size;
    double *c,x,y,z,fsize;

    size = NumAtoms();
    fsize = -1.0/(double)NumAtoms();

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Center", obAuditMsg);

    vector<double*>::iterator i;
    for (i = _vconf.begin();i != _vconf.end();i++)
      {
        c = *i;
        x = y = z = 0.0;
        for (j = 0;j < size;j++)
	  {
            x += c[j*3];
            y += c[j*3+1];
            z += c[j*3+2];
	  }
        x *= fsize;
        y *= fsize;
        z *= fsize;

        for (j = 0;j < size;j++)
	  {
            c[j*3]+=x;
            c[j*3+1]+=y;
            c[j*3+2]+=z;
	  }
      }

  }

  vector3 OBMol::Center(int nconf)
  {
    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Center", obAuditMsg);

    SetConformer(nconf);

    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    double x=0.0,y=0.0,z=0.0;
    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        x += atom->x();
        y += atom->y();
        z += atom->z();
      }

    x /= (double)NumAtoms();
    y /= (double)NumAtoms();
    z /= (double)NumAtoms();

    vector3 vtmp;
    vector3 v(x,y,z);

    for (atom = BeginAtom(i);atom;atom = NextAtom(i))
      {
        vtmp = atom->GetVector() - v;
        atom->SetVector(vtmp);
      }

    return(v);
  }


  /*! this method adds the vector v to all atom positions in all conformers */
  void OBMol::Translate(const vector3 &v)
  {
    for (int i = 0;i < NumConformers();i++)
      Translate(v,i);
  }

  /*! this method adds the vector v to all atom positions in the
    conformer nconf. If nconf == OB_CURRENT_CONFORMER, then the atom
    positions in the current conformer are translated. */
  void OBMol::Translate(const vector3 &v,int nconf)
  {
    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Translate", obAuditMsg);

    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);

    x = v.x();
    y = v.y();
    z = v.z();
    size = NumAtoms();
    for (i = 0;i < size;i++)
      {
        c[i*3  ] += x;
        c[i*3+1] += y;
        c[i*3+2] += z;
      }
  }

  void OBMol::Rotate(const double u[3][3])
  {
    int i,j,k;
    double m[9];
    for (k=0,i = 0;i < 3;i++)
      for (j = 0;j < 3;j++)
	m[k++] = u[i][j];

    for (i = 0;i < NumConformers();i++)
      Rotate(m,i);
  }

  void OBMol::Rotate(const double m[9])
  {
    for (int i = 0;i < NumConformers();i++)
      Rotate(m,i);
  }

  void OBMol::Rotate(const double m[9],int nconf)
  {
    int i,size;
    double x,y,z;
    double *c = (nconf == OB_CURRENT_CONFORMER)? _c : GetConformer(nconf);

    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::Rotate", obAuditMsg);

    size = NumAtoms();
    for (i = 0;i < size;i++)
      {
        x = c[i*3  ];
        y = c[i*3+1];
        z = c[i*3+2];
        c[i*3  ] = m[0]*x + m[1]*y + m[2]*z;
        c[i*3+1] = m[3]*x + m[4]*y + m[5]*z;
        c[i*3+2] = m[6]*x + m[7]*y + m[8]*z;
      }
  }


  void OBMol::SetConformers(vector<double*> &v)
  {
    vector<double*>::iterator i;
    for (i = _vconf.begin();i != _vconf.end();i++)
      delete [] *i;

    _vconf = v;
    _c = (_vconf.empty()) ? NULL : _vconf[0];

  }

  void OBMol::CopyConformer(double *c,int idx)
  {
    //    obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());
    memcpy((char*)_vconf[idx],(char*)c,sizeof(double)*3*NumAtoms());
  }

  // void OBMol::CopyConformer(double *c,int idx)
  // {
  //   obAssert(!_vconf.empty() && (unsigned)idx < _vconf.size());

  //   unsigned int i;
  //   for (i = 0;i < NumAtoms();i++)
  //     {
  //       _vconf[idx][i*3  ] = (double)c[i*3  ];
  //       _vconf[idx][i*3+1] = (double)c[i*3+1];
  //       _vconf[idx][i*3+2] = (double)c[i*3+2];
  //     }
  // }

  void OBMol::DeleteConformer(int idx)
  {
    if (idx < 0 || idx >= (signed)_vconf.size())
      return;

    delete [] _vconf[idx];
    _vconf.erase((_vconf.begin()+idx));
  }

  ///Converts for instance [N+]([O-])=O to N(=O)=O
  bool OBMol::ConvertDativeBonds()
  {
    obErrorLog.ThrowError(__FUNCTION__,
			  "Ran OpenBabel::ConvertDativeBonds", obAuditMsg);

    //Look for + and - charges on adjacent atoms
    OBAtom* patom;
    vector<OBNodeBase*>::iterator i;
    for (patom = BeginAtom(i);patom;patom = NextAtom(i))
      {
	vector<OBEdgeBase*>::iterator itr;
	OBBond *pbond;
	for (pbond = patom->BeginBond(itr);patom->GetFormalCharge() && pbond;pbond = patom->NextBond(itr))
	  {
	    OBAtom* pNbratom = pbond->GetNbrAtom(patom);
	    int chg1 = patom->GetFormalCharge();
	    int chg2 = pNbratom->GetFormalCharge();
	    if((chg1>0 && chg2<0)|| (chg1<0 && chg2>0))
	      {
		//dative bond. Reduce charges and increase bond order
		if(chg1>0)
		  --chg1;
		else 
		  ++chg1;
		patom->SetFormalCharge(chg1);
		if(chg2>0)
		  --chg2;
		else 
		  ++chg2;
		pNbratom->SetFormalCharge(chg2);
		pbond->SetBO(pbond->GetBO()+1);
	      }
	  }		
      }
    return true;
  }

  OBAtom *OBMol::BeginAtom(vector<OBNodeBase*>::iterator &i)
  {
    i = _vatom.begin();
    return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
  }

  OBAtom *OBMol::NextAtom(vector<OBNodeBase*>::iterator &i)
  {
    i++;
    return((i == _vatom.end()) ? (OBAtom*)NULL : (OBAtom*)*i);
  }

  OBBond *OBMol::BeginBond(vector<OBEdgeBase*>::iterator &i)
  {
    i = _vbond.begin();
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

  OBBond *OBMol::NextBond(vector<OBEdgeBase*>::iterator &i)
  {
    i++;
    return((i == _vbond.end()) ? (OBBond*)NULL : (OBBond*)*i);
  }

} // end namespace OpenBabel

//! \file mol.cpp
//! \brief Handle molecules. Implementation of OBMol.
