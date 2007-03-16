/**********************************************************************
groupcontrib.h - Handle group contribution algorithms.
 
Copyright (C) 2007 by Tim Vandermeersch
 
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

namespace OpenBabel
{
  // class introduction in groupcontrib.cpp
  class OBAPI OBGroupContrib
  {
    public:
      //! constructor
      OBGroupContrib();
      //! destructor
      ~OBGroupContrib();
      bool ParseFile(const char *filename);
      /*! Predict the logP, MR, TPSA (see derived classes) for molecule mol using the group contributions
       *  algorithm from JOELib2.
       *
       *  \param mol OBMol object for which to predict the logP, MR, TPSA
       *  \return predicted logP
       */
      double GroupContributions(OBMol &mol);
    private:
      std::vector<std::pair<OBSmartsPattern*, double> > _contribsHeavy; //! heavy atom contributions
      std::vector<std::pair<OBSmartsPattern*, double> > _contribsHydrogen; //!  hydrogen contributions
  };

  // class introduction in logp.cpp
  class OBAPI OBLogP: public OBGroupContrib
  {
    public:
      //! constructor
      OBLogP();
      //! destructor
      ~OBLogP();
      /*! Predict the logP for molecule mol using the group contributions
       *  algorithm from JOELib2.
       *
       *  \param mol OBMol object for which to predict the logP
       *  \return predicted logP
       */
      double Predict(OBMol &mol);
  };
  
  // class introduction in psa.cpp
  class OBAPI OBPSA: public OBGroupContrib
  {
    public:
      //! constructor
      OBPSA();
      //! destructor
      ~OBPSA();
      /*! Predict the TPSA (Topological Polar Surface Area) for molecule mol 
       *  using the group contributions algorithm from JOELib2.
       *
       *  \param mol OBMol object for which to predict the TPSA
       *  \return predicted TPSA
       */
      double Predict(OBMol &mol);
  };

  // class introduction in psa.cpp
  class OBAPI OBMR: public OBGroupContrib
  {
    public:
      //! constructor
      OBMR();
      //! destructor
      ~OBMR();
      /*! Predict the MR (Molecular Reractivity) for molecule mol 
       *  using the group contributions algorithm from JOELib2.
       *
       *  \param mol OBMol object for which to predict the TPSA
       *  \return predicted TPSA
       */
      double Predict(OBMol &mol);
  };


} // end namespace OpenBabel

#endif // OB_GROUPCONTRIB_H

//! \file groupcontrib.h
//! \brief Handle group contribution algorithms.
