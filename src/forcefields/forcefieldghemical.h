/**********************************************************************
forcefieldghemical.h - Ghemical force field.
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#include <vector>
#include <string>
#include <map>

#include <openbabel/forcefield.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>

namespace OpenBabel
{
  // Class OBForceFieldMM2
  // class introduction in forcefield.cpp
  class OBAPI OBForceFieldGhemical: public OBForceField
  {
    protected:
      //! \return Parses the parameter file
      bool ParseParamFile();
      //! \return Sets atomtypes to TRPS in _mol
      bool SetTRPSTypes();
      
      double bondunit, bond_cubic, bond_quartic;
      double angleunit, angle_sextic;
      double stretchbendunit;
      double torsionunit;
      double outplanebendunit;
      double a_expterm, b_expterm, c_expterm;
      double dielectric;
      std::vector<OBFFParameter> _ffbondparams; // a = atom 1 of bond
                                                // b = atom 2 of bond
						// dpar1 = lenght
						// dpar2 = force
      std::vector<OBFFParameter> _ffangleparams; // a = atom 1 of angle abc
                                                 // b = atom 2 of angle abc
						 // c = atom 3 of angle abc
						 // dpar1 = angle
						 // dpar2 = force
      std::vector<OBFFParameter> _fftorsionparams; // a = atom 1 of torsion
                                                   // b = atom 2 of torsion
                                                   // c = atom 3 of torsion
                                                   // d = atom 4 of torsion
						   // dpar1 = v1
						   // dpar2 = v2
						   // dpar3 = v3
      std::vector<OBFFParameter> _ffvdwparams;
      std::vector<vector3> forces; // used to hold forces on each atom



    public:
      //! Setup
      bool Setup(OBMol &mol);
      //! Constructor
      OBForceFieldGhemical(std::string ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      std::string Description()
	{ return "Ghemical force field.";};
      
      std::string GetUnit() { return std::string("kJ/mol"); }


      //! Destructor
      virtual ~OBForceFieldGhemical();
      //! Assignment
      OBForceFieldGhemical &operator = (OBForceFieldGhemical &);
      //! Returns total energy
      double Energy();
     //! Returns the bond stretching energy
      double E_Bond();
      //! Returns the angle bending energy
      double E_Angle();
      //! Returns the torsional energy
      double E_Torsion();
      double E_OOP();
      //! Returns energy due to Van der Waals interactions
      double E_VDW();

  }; // class OBForceFieldGhemical

}// namespace OpenBabel

//! \file forcefieldGhemical.h
//! \brief Ghemical force field

