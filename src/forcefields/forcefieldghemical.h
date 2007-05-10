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
  class OBFFBondCalculationGhemical : public OBFFCalculation
  {
    public:
      double kb, r0, rab, delta;
      int bt; // bondtype
      
      void Compute(bool gradients = true);
  };
  
  class OBFFAngleCalculationGhemical : public OBFFCalculation
  {
    public:
      double ka, theta0, theta, delta;
      
      void Compute(bool gradients = true);
  };
  
  class OBFFTorsionCalculationGhemical : public OBFFCalculation
  {
    public:
      double V, s, n, tor;
      double k1, k2, k3;
      int tt; //torsiontype
      
      void Compute(bool gradients = true);
  };

  class OBFFVDWCalculationGhemical : public OBFFCalculation
  {
    public:
      double ka, Ra, kb, Rb, kab, rab;
      bool is14, samering;

      void Compute(bool gradients = true);
  };

  class OBFFElectrostaticCalculationGhemical : public OBFFCalculation
  {
    public:
      double qq, rab;
      
      void Compute(bool gradients = true);
  };

  // Class OBForceFieldGhemical
  // class introduction in forcefieldghemical.cpp
  class OBAPI OBForceFieldGhemical: public OBForceField
  {
    protected:
      //!  Parses the parameter file
      bool ParseParamFile();
      //!  Sets atomtypes to Ghemical types in _mol
      bool SetGhemicalTypes();
      //!  Sets partial charges to Ghemical charges in _mol
      bool SetGhemicalCharges();
      //! fill OBFFXXXCalculation vectors
      bool SetupCalculations();
      //! Same as OBForceField::GetParameter, but takes (bond/angle/torsion) type in account.
      OBFFParameter* GetParameterGhemical(int type, const char* a, const char* b, const char* c, const char* d, std::vector<OBFFParameter> &parameter);
      //! Returns the negative gradient (force) on atom a
      vector3 GetGradient(OBAtom *a, int terms = OBFF_ENERGY);
      
      // OBFFParameter vectors to contain the parameters
      std::vector<OBFFParameter> _ffbondparams; 
      std::vector<OBFFParameter> _ffangleparams; 
      std::vector<OBFFParameter> _fftorsionparams; 
      std::vector<OBFFParameter> _ffvdwparams;
      std::vector<OBFFParameter> _ffchargeparams;

      // OBFFXXXCalculationYYY vectors to contain the calculations
      std::vector<OBFFBondCalculationGhemical>          _bondcalculations;
      std::vector<OBFFAngleCalculationGhemical>         _anglecalculations;
      std::vector<OBFFTorsionCalculationGhemical>       _torsioncalculations;
      std::vector<OBFFVDWCalculationGhemical>           _vdwcalculations;
      std::vector<OBFFElectrostaticCalculationGhemical> _electrostaticcalculations;

    public:
      //! Constructor
      OBForceFieldGhemical(const char* ID, bool IsDefault=true) : OBForceField(ID, IsDefault)
      {
        ParseParamFile();
      }
      
      //! Destructor
      virtual ~OBForceFieldGhemical();
      
      //! Assignment
      OBForceFieldGhemical &operator = (OBForceFieldGhemical &);
      
      //! Get the description for this force field
      const char* Description() 
      { 
        return "Ghemical force field.";
      }
      
      //! Get the unit in wich the energy is expressed
      std::string GetUnit() 
      { 
        return std::string("kJ/mol"); 
      }

      //! Setup
      bool Setup(OBMol &mol);
      
      //! Returns total energy
      double Energy(bool gradients = true);
     //! Returns the bond stretching energy
      double E_Bond(bool gradients = true);
      //! Returns the angle bending energy
      double E_Angle(bool gradients = true);
      //! Returns the torsional energy
      double E_Torsion(bool gradients = true);
      //! Returns energy due to Van der Waals interactions
      double E_VDW(bool gradients = true);
      //! Returns energy due to electrostatic interactions
      double E_Electrostatic(bool gradients = true);
      
      //! Compare and print the numerical and analytical gradients
      bool ValidateGradients();

  }; // class OBForceFieldGhemical

}// namespace OpenBabel

//! \file forcefieldghemical.h
//! \brief Ghemical force field

