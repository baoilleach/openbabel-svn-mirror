/**********************************************************************
groupcontrib.h - Handle group contribution algorithms.
 
Copyright (C) 2007      by Tim Vandermeersch
              2001-2007 by Stephen Jelfs
              2001-2007 by Joerg Kurt Wegner, me@cheminformatics.eu

Original version: JOELib2, http://joelib.sf.net
 
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

#ifndef OB_GROUPCONTRIB_H
#define OB_GROUPCONTRIB_H

#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/descriptor.h>

namespace OpenBabel
{

  /** \class OBGroupContrib groupcontrib.h <openbabel/groupcontrib.h>
      \brief Handle group contribution algorithms.
 
      This is the base class for calculations that use the JOELib2 contribution 
      algorithm. 
    */
class OBAPI OBGroupContrib : public OBDescriptor
{
public:
  //! constructor. Each instance provides an ID and a datafile.
  OBGroupContrib(const char* ID, const char* filename, const char* descr)
    : OBDescriptor(ID, false), _filename(filename), _descr(descr){}

  /*! Predict the logP, MR, TPSA (each instance of OBGroupContrib 
   *  uses different parameters loaded from its own datafile) for 
   *  molecule mol using the group contributions algorithm from JOELib2.
   */
  virtual const char* Description(){ return _descr;}; 
  virtual double Predict(OBBase* pOb); 

 private:
  bool ParseFile();

  const char* _filename;
  const char* _descr;
  std::vector<std::pair<OBSmartsPattern*, double> > _contribsHeavy; //! heavy atom contributions
  std::vector<std::pair<OBSmartsPattern*, double> > _contribsHydrogen; //!  hydrogen contributions
};

/* The classes OBLogp, OBPSA and OBMR have been replaced by instances of
OBGroupContrib with different IDs.
So instead of: 
      OBLogp logP;
      cout << "logP  " << logP.Predict(mol) << endl;
use:
      OBDescriptor* pDesc = OBDescriptor::FindType("logP")
      if(pDesc)
        cout << "logP  " << pDesc->Predict(&mol) << endl;
*/
} // end namespace OpenBabel

#endif // OB_GROUPCONTRIB_H

//! \file groupcontrib.h
//! \brief Handle group contribution algorithms.
