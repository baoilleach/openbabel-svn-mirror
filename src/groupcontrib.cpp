/**********************************************************************
groupcontrib.cpp - Handle logP, PSA, MR, and other group-based predictions

Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu
              2007      by Chris Morley

Original version: JOELib2, http://joelib.sf.net

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <vector>
#include <utility>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/groupcontrib.h>
#include <openbabel/locale.h>

using namespace std;

namespace OpenBabel
{

  const char* OBGroupContrib::Description()
  {
   //Adds name of datafile containing SMARTS strings to the description
    static string txt;
    txt =  _descr;
    txt += "\n Datafile: ";
    txt += _filename;
    txt += "\nOBGroupContrib is definable";
    return txt.c_str();
  }

  bool OBGroupContrib::ParseFile()
  {
    OBSmartsPattern *sp;

    // open data file
    ifstream ifs;

    if (OpenDatafile(ifs, _filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, " Could not find contribution data file.", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    vector<string> vs;
    bool heavy = false;
    string ln;
    while(getline(ifs,ln)){
      if(ln[0]=='#') continue;
      if(ln.find(";heavy")!=string::npos)
        heavy=true;
      if(ln[0]==';') continue;
      tokenize(vs, ln);

      if (vs.size() < 2)
        continue;

      sp = new OBSmartsPattern;//causes non-serious memory leak.
      // Could be cured by copying OBSmartsPattern rather than a pointer in vectors
      if (sp->Init(vs[0]))
      {
        if (heavy)
          _contribsHeavy.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
        else
          _contribsHydrogen.push_back(pair<OBSmartsPattern*, double> (sp, atof(vs[1].c_str())));
      }
      else
      {
        delete sp;
        sp = NULL;
        obErrorLog.ThrowError(__FUNCTION__, " Could not parse SMARTS from contribution data file", obInfo);

        // return the locale to the original one
        obLocale.RestoreLocale();

        return false;
      }
    }

    // return the locale to the original one
    obLocale.RestoreLocale();
    return true;
  }


  double OBGroupContrib::Predict(OBBase* pOb, string* param)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(!pmol)
      return 0.0;

    //Need to add hydrogens, so do this to a copy to leave original unchanged
    OBMol mol(*pmol);
    mol.AddHydrogens(false, false);

    //Read in data, unless it has already been done.
    if(_contribsHeavy.empty() && _contribsHydrogen.empty())
      ParseFile();

    vector<vector<int> > _mlist; // match list for atom typing
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*, double> >::iterator i;

    vector<double> atomValues(mol.NumAtoms(), 0.0);

    OBMol tmpmol;
    tmpmol = mol;

    tmpmol.ConvertDativeBonds();

    // atom contributions
    //cout << "atom contributions:" << endl;
    for (i = _contribsHeavy.begin();i != _contribsHeavy.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  atomValues[(*j)[0] - 1] = i->second;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }

    vector<double> hydrogenValues(tmpmol.NumAtoms(), 0.0);
    //hydrogenValues.resize(tmpmol.NumAtoms());

    // hydrogen contributions
    //cout << "hydrogen contributions:" << endl;
    for (i = _contribsHydrogen.begin();i != _contribsHydrogen.end();++i) {
      if (i->first->Match(tmpmol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
	  int Hcount = tmpmol.GetAtom((*j)[0])->GetValence() - tmpmol.GetAtom((*j)[0])->GetHvyValence();
	  hydrogenValues[(*j)[0] - 1] = i->second * Hcount;
	  //cout << (*j)[0] << " = " << i->first->GetSMARTS() << " : " << i->second << endl;
        }
      }
    }

    // total atomic and hydrogen contribution
    double total = 0.0;

    for (int index = 0; index < tmpmol.NumAtoms(); index++) {
      if (tmpmol.GetAtom(index+1)->IsHydrogen())
        continue;

      total += atomValues[index];
      total += hydrogenValues[index];
    }

    /*
    FOR_ATOMS_OF_MOL (a, tmpmol)
      cout << "hydrogens on atom " << a->GetIdx() << ": " << a->GetValence() - a->GetHvyValence() << endl;
    for (int index = 0; index < tmpmol.NumAtoms(); index++)
      cout << "atom " << index << ": " << atomValues[index] << endl;
    for (int index = 0; index < tmpmol.NumAtoms(); index++)
      cout << "hydrogen " << index << ": " << hydrogenValues[index] << endl;
    */

    return total;
  }

//! \file groupcontrib.cpp
//! \brief Handle logP, PSA and other group-based prediction algorithms.

  }//namespace
