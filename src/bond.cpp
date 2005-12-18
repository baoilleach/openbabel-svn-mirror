/**********************************************************************
bond.cpp - Handle OBBond class
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#include "mol.h"
#include "typer.h"

using namespace std;

namespace OpenBabel
{

extern OBAromaticTyper  aromtyper;

/** \class OBBond
    \brief Bond class
 
    The OBBond class is straightforward in its data access and
    modification methods. OBBonds store pointers to the atoms on each end
    of the bond. In storing pointers to atoms instead of integer indices,
    the necessity of having to reorder bonds when atoms are shuffled,
    added, or delete is obviated.
*/

// *******************************
// *** OBBond member functions ***
// *******************************

OBBond::OBBond()
{
    _idx=0;
    _order=0;
    _flags=0;
    _bgn=NULL;
    _end=NULL;
    _vdata.clear();
}

OBBond::~OBBond()
{
/*  if (!_vdata.empty())
    {
        vector<OBGenericData*>::iterator m;
        for (m = _vdata.begin();m != _vdata.end();m++)
            delete *m;
        _vdata.clear();
    }
*/
}

void OBBond::Set(int idx,OBAtom *begin,OBAtom *end,int order,int flags)
{
    SetIdx(idx);
    SetBegin(begin);
    SetEnd(end);
    SetBO(order);
    SetFlag(flags);
}

void OBBond::SetBO(int order)
{
    _order = (char)order;
    if (order == 5)
    {
        SetAromatic();
        if (_bgn)
            _bgn->SetAromatic();
        if (_end)
            _end->SetAromatic();
    }
    else
      {
	if (order == 1)
	  SetKSingle();
	else if (order == 2)
	  SetKDouble();
	else if (order == 3)
	  SetKTriple();

        UnsetAromatic();
      }
}

void OBBond::SetLength(OBAtom *fixed, double length)
{
    unsigned int i;
    OBMol *mol = (OBMol*)fixed->GetParent();
    vector3 v1,v2,v3,v4,v5;
    vector<int> children;

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::SetBondLength", obAuditMsg);

    int a = fixed->GetIdx();
    int b = GetNbrAtom(fixed)->GetIdx();

    mol->FindChildren(children,a,b);
    children.push_back(b);

    v1 = GetNbrAtom(fixed)->GetVector();
    v2 = fixed->GetVector();
    v3 = v1 - v2;
    v3.normalize();
    v3 *= length;
    v3 += v2;
    v4 = v3 - v1;

    for ( i = 0 ; i < children.size() ; i++ )
    {
        v1 = mol->GetAtom(children[i])->GetVector();
        v1 += v4;
        mol->GetAtom(children[i])->SetVector(v1);
        /*
        idx = (children[i]-1) * 3;
        c[idx]   += x; 
        c[idx+1] += y;
        c[idx+2] += z;
        */
    }
}

bool OBBond::IsRotor()
{
    return(_bgn->GetHvyValence() > 1 && _end->GetHvyValence() > 1 &&
           _order == 1 && !IsInRing() && _bgn->GetHyb() != 1 &&
           _end->GetHyb() != 1);
}

bool OBBond::IsAmide()
{
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
    {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
    }

    if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
    {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
    }

    if (!a1 || !a2)
        return(false);
    if (GetBO() != 1)
        return(false);

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
        if (bond->IsCarbonyl())
            return(true);

    return(false);
}

bool OBBond::IsPrimaryAmide()
{
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 7)
    {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
    }

    if (_bgn->GetAtomicNum() == 7 && _end->GetAtomicNum() == 6)
    {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
    }

    if (!a1 || !a2)
        return(false);
    if (GetBO() != 1)
        return(false);

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
        if (bond->IsCarbonyl())
            if (a2->GetHvyValence() == 2)
                return(true);

    return(false);
}

//! \todo Implement this properly -- currently always returns false
bool OBBond::IsSecondaryAmide()
{
    return(false);
}

bool OBBond::IsEster()
{
    OBAtom *a1,*a2;
    a1 = a2 = NULL;

    if (_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8)
    {
        a1 = (OBAtom*)_bgn;
        a2 = (OBAtom*)_end;
    }

    if (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6)
    {
        a1 = (OBAtom*)_end;
        a2 = (OBAtom*)_bgn;
    }

    if (!a1 || !a2)
        return(false);
    if (GetBO() != 1)
        return(false);

    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (bond = a1->BeginBond(i);bond;bond = a1->NextBond(i))
        if (bond->IsCarbonyl())
            return(true);

    return(false);
}

bool OBBond::IsCarbonyl()
{
    if (GetBO() != 2)
        return(false);

    if ((_bgn->GetAtomicNum() == 6 && _end->GetAtomicNum() == 8) ||
            (_bgn->GetAtomicNum() == 8 && _end->GetAtomicNum() == 6))
        return(true);

    return(false);
}

bool OBBond::IsSingle()
{
    if (HasFlag(OB_AROMATIC_BOND))
        return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
    {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
    }

    if ((this->GetBondOrder()==1) && !(HasFlag(OB_AROMATIC_BOND)))
        return(true);

    return(false);
}

bool OBBond::IsDouble()
{
    if	(HasFlag(OB_AROMATIC_BOND))
        return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
    {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
    }

    if ((this->GetBondOrder()==2) && !(HasFlag(OB_AROMATIC_BOND)))
        return(true);

    return(false);
}

bool OBBond::IsTriple()
{
    if	(HasFlag(OB_AROMATIC_BOND))
        return(false);

    if (!((OBMol*)GetParent())->HasAromaticPerceived())
    {
        aromtyper.AssignAromaticFlags(*((OBMol*)GetParent()));
    }

    if ((this->GetBondOrder()==3) && !(HasFlag(OB_AROMATIC_BOND)))
        return(true);

    return(false);
}

bool OBBond::IsAromatic() const
{
    if (((OBBond*)this)->HasFlag(OB_AROMATIC_BOND))
        return(true);

    OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();
    if (!mol->HasAromaticPerceived())
    {
        aromtyper.AssignAromaticFlags(*mol);
        if (((OBBond*)this)->HasFlag(OB_AROMATIC_BOND))
            return(true);
    }

    return(false);
}

/*! This method checks if the geometry around this bond looks unsaturated
  by measuring the torsion angles formed by all connected atoms X-start=end-Y
  and checking that they are close to 0 or 180 degrees */
bool OBBond::IsDoubleBondGeometry()
{
  double torsion;
  OBAtom *nbrStart,*nbrEnd;
  vector<OBEdgeBase*>::iterator i,j;
  // We concentrate on sp2 atoms with valence up to 3 and ignore the rest (like sp1 or S,P)
  // As this is called from PerceiveBondOrders, GetHyb() may still be undefined.
  if (_bgn->GetHyb()==1 || _bgn->GetValence()>3||
      _end->GetHyb()==1 || _end->GetValence()>3)
    return(true);

  for (nbrStart = static_cast<OBAtom*>(_bgn)->BeginNbrAtom(i); nbrStart;
       nbrStart = static_cast<OBAtom*>(_bgn)->NextNbrAtom(i))
    {
      if (nbrStart != _end)
	{ 
	  for (nbrEnd = static_cast<OBAtom*>(_end)->BeginNbrAtom(j);
	       nbrEnd; nbrEnd = static_cast<OBAtom*>(_end)->NextNbrAtom(j))
	    {
	      if (nbrEnd != _bgn)
		{
		  torsion=fabs(CalcTorsionAngle(nbrStart->GetVector(),
						static_cast<OBAtom*>(_bgn)->GetVector(),
						static_cast<OBAtom*>(_end)->GetVector(),
						nbrEnd->GetVector()));

		  // >12&&<168 not enough
		  if (torsion > 15.0f  && torsion < 165.0f)
		    {
		      // Geometry does not match a double bond
		      return(false);
		    }

		}
	    }  // end loop for neighbors of end
	}
    } // end loop for neighbors of start
  return(true);
}

void OBBond::SetKSingle()
{
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KSINGLE_BOND;
}

void OBBond::SetKDouble()
{
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KDOUBLE_BOND;
}

void OBBond::SetKTriple()
{
    _flags &= (~(OB_KSINGLE_BOND|OB_KDOUBLE_BOND|OB_KTRIPLE_BOND));
    _flags |= OB_KTRIPLE_BOND;
}

bool OBBond::IsKSingle()
{
    if (_flags & OB_KSINGLE_BOND)
        return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
        ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KSINGLE_BOND) != 0) ? true : false;
}

bool OBBond::IsKDouble()
{
    if (_flags & OB_KDOUBLE_BOND)
        return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
        ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KDOUBLE_BOND) != 0) ? true : false;
}

bool OBBond::IsKTriple()
{
    if (_flags & OB_KTRIPLE_BOND)
        return(true);
    if (!((OBMol*)GetParent())->HasKekulePerceived())
        ((OBMol*)GetParent())->NewPerceiveKekuleBonds();

    return((_flags & OB_KTRIPLE_BOND) != 0) ? true : false;
}

bool OBBond::IsInRing() const
{
    if (((OBBond*)this)->HasFlag(OB_RING_BOND))
        return(true);

    OBMol *mol = (OBMol*)((OBBond*)this)->GetParent();
    if (!mol->HasRingAtomsAndBondsPerceived())
    {
        mol->FindRingAtomsAndBonds();
        if (((OBBond*)this)->HasFlag(OB_RING_BOND))
            return(true);
    }

    return(false);
}

bool OBBond::IsClosure()
{
    OBMol *mol = (OBMol*)GetParent();
    if (!mol)
        return(false);
    if (mol->HasClosureBondsPerceived())
        return(HasFlag(OB_CLOSURE_BOND));

    mol->SetClosureBondsPerceived();

    obErrorLog.ThrowError(__FUNCTION__,
                          "Ran OpenBabel::PerceiveClosureBonds", obAuditMsg);

    OBBond *bond;
    OBAtom *atom,*nbr;
    OBBitVec uatoms,ubonds;
    vector<OBNodeBase*> curr,next;
    vector<OBNodeBase*>::iterator i;
    vector<OBEdgeBase*>::iterator j;

    uatoms.Resize(mol->NumAtoms()+1);
    ubonds.Resize(mol->NumAtoms()+1);

    for (;uatoms.CountBits() < (signed)mol->NumAtoms();)
    {
        if (curr.empty())
            for (atom = mol->BeginAtom(i);atom;atom = mol->NextAtom(i))
                if (!uatoms[atom->GetIdx()])
                {
                    uatoms |= atom->GetIdx();
                    curr.push_back(atom);
                    break;
                }

        for (;!curr.empty();)
        {
            for (i = curr.begin();i != curr.end();i++)
                for (nbr = ((OBAtom*)*i)->BeginNbrAtom(j);nbr;nbr = ((OBAtom*)*i)->NextNbrAtom(j))
                    if (!uatoms[nbr->GetIdx()])
                    {
                        uatoms |= nbr->GetIdx();
                        ubonds |= (*j)->GetIdx();
                        next.push_back(nbr);
                    }

            curr = next;
            next.clear();
        }
    }

    for (bond = mol->BeginBond(j);bond;bond = mol->NextBond(j))
        if (!ubonds[bond->GetIdx()])
            bond->SetClosure();

    return(HasFlag(OB_CLOSURE_BOND));
}

double OBBond::GetEquibLength()
{
    double length;
    OBAtom *begin, *end;
    // CorrectedBondRad will always return a # now
    //  if (!CorrectedBondRad(GetBeginAtom(),rad1)) return(0.0);
    //  if (!CorrectedBondRad(GetEndAtom(),rad2))   return(0.0);

    begin = GetBeginAtom();
    end = GetEndAtom();
    length = etab.CorrectedBondRad(begin->GetAtomicNum(), begin->GetHyb())
             + etab.CorrectedBondRad(end->GetAtomicNum(), end->GetHyb());

    if (IsAromatic())
        length *= 0.93;
    else if (GetBO() == 2)
        length *= 0.91;
    else if (GetBO() == 3)
        length *= 0.87;
    return(length);
}

double OBBond::GetLength()
{
    double	d2;
    OBAtom *begin, *end;
    begin = GetBeginAtom();
    end = GetEndAtom();

    d2 = SQUARE(begin->GetX() - end->GetX());
    d2 += SQUARE(begin->GetY() - end->GetY());
    d2 += SQUARE(begin->GetZ() - end->GetZ());

    return(sqrt(d2));
}

/*Now in OBBase
// OBGenericData methods
bool OBBond::HasData(string &s)
//returns true if the generic attribute/value pair exists
{
    if (_vdata.empty())
        return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(true);

    return(false);
}

bool OBBond::HasData(const char *s)
//returns true if the generic attribute/value pair exists
{
    if (_vdata.empty())
        return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(true);

    return(false);
}

bool OBBond::HasData(unsigned int dt)
//returns true if the generic attribute/value pair exists
{
    if (_vdata.empty())
        return(false);

    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt)
            return(true);

    return(false);
}

OBGenericData *OBBond::GetData(string &s)
//returns the value given an attribute
{
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(*i);

    return(NULL);
}

OBGenericData *OBBond::GetData(const char *s)
//returns the value given an attribute
{
    vector<OBGenericData*>::iterator i;

    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetAttribute() == s)
            return(*i);

    return(NULL);
}

OBGenericData *OBBond::GetData(unsigned int dt)
{
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt)
            return(*i);
    return(NULL);
}

void OBBond::DeleteData(unsigned int dt)
{
    vector<OBGenericData*> vdata;
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();i++)
        if ((*i)->GetDataType() == dt)
            delete *i;
        else
            vdata.push_back(*i);
    _vdata = vdata;
}

void OBBond::DeleteData(vector<OBGenericData*> &vg)
{
    vector<OBGenericData*> vdata;
    vector<OBGenericData*>::iterator i,j;

    bool del;
    for (i = _vdata.begin();i != _vdata.end();i++)
    {
        del = false;
        for (j = vg.begin();j != vg.end();j++)
            if (*i == *j)
            {
                del = true;
                break;
            }
        if (del)
            delete *i;
        else
            vdata.push_back(*i);
    }
    _vdata = vdata;
}

void OBBond::DeleteData(OBGenericData *gd)
{
    vector<OBGenericData*>::iterator i;
    for (i = _vdata.begin();i != _vdata.end();i++)
        if (*i == gd)
        {
            delete *i;
            _vdata.erase(i);
        }

}
*/

} // end namespace OpenBabel

//! \file bond.cpp
//! \brief Handle OBBond class
