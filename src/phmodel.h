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

#ifndef _PHMODEL_H
#define _PHMODEL_H

#include "parsmart.h"
#include "data.h"

namespace OpenBabel {

class OBChemTsfm
{
  std::vector<int>                       _vadel;
  std::vector<std::pair<int,int> >            _vele;
  std::vector<std::pair<int,int> >            _vchrg;
  std::vector<std::pair<int,int> >            _vbdel;
  std::vector<std::pair<std::pair<int,int>,int> >  _vbond;
  OBSmartsPattern _bgn,_end;
public:
  OBChemTsfm() {}
  ~OBChemTsfm() {}
  bool Init(std::string&,std::string&);
  bool Apply(OBMol&);
};


class OBPhModel : public OBGlobalDataBase
{
  std::vector<std::vector<int> >                           _mlist;
  std::vector<OBChemTsfm*>                            _vtsfm;
  std::vector<std::pair<OBSmartsPattern*,std::vector<float> > > _vschrg;
 public:
  OBPhModel();
  ~OBPhModel();

  void ParseLine(const char*);
  void AssignSeedPartialCharge(OBMol&);
  void CorrectForPH(OBMol&);
};



} //namespace OpenBabel

#endif //_PHMODEL_H
