/**********************************************************************
molchrg.h - Assign Gasteiger partial charges.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

#ifndef OB_MOLCHRG_H
#define OB_MOLCHRG_H

namespace OpenBabel {

class GasteigerState;

// class introduction in molchrg.cpp
class OBGastChrg
{
  std::vector <GasteigerState*> _gsv;
  void InitialPartialCharges(OBMol &);
  bool GasteigerSigmaChi(OBAtom *,double &,double &,double &);
 public:
  OBGastChrg(){}
  ~OBGastChrg();
  bool AssignPartialCharges(OBMol &);
  void GSVResize(int);
};

class GasteigerState 
{
  //helper class for OBGastChrg
 public:
  GasteigerState();
  ~GasteigerState() {}
  void SetValues(double _a,double _b,double _c,double _q) 
    {a = _a;b = _b;c = _c;denom=a+b+c;q = _q;}
  double a, b, c;
  double denom;
  double chi;
  double q;
};

}

#define OB_GASTEIGER_DENOM  20.02
#define OB_GASTEIGER_DAMP   0.5
#define OB_GASTEIGER_ITERS  6

#endif // OB_MOLCHRG_H
