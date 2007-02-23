/**********************************************************************
forcefieldghemical.cpp - Ghemical force field.
 
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
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include "forcefieldghemical.h"

using namespace std;

namespace OpenBabel
{
  double OBFFBondCalculationGhemical::GetEnergy()
  {
    double delta2;
    vector3 vab;
     
    vab = a->GetVector() - b->GetVector();
    rab = vab.length();
    delta = rab - r0;
    delta2 = delta * delta;

    energy = kb * delta2;

    return energy;
  }
   
  vector3 OBFFBondCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, gradient;
    double dE;

    if ((atom != a) && (atom != b))
      return  VZero;
     
    da = a->GetVector();
    db = b->GetVector();
    rab = OBForceField::VectorLengthDerivative(da, db);
    delta = rab - r0;
 
    dE = 2.0f * kb * delta;

    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    }
  }

  double OBForceFieldGhemical::E_Bond()
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    double energy;
    
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES  BOND    BOND       IDEAL       FORCE\n");
      OBFFLog(" I    J     TYPE   LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }
 
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); i++) {

      energy += i->GetEnergy();

      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s    %d   %8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        OBFFLog(logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL BOND STRETCHING ENERGY = %8.3f %s\n",  energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
    return energy;
  }
  
  double OBFFAngleCalculationGhemical::GetEnergy()
  {
    double delta2;

    theta = a->GetAngle(b->GetIdx(), c->GetIdx());

    if (theta0 > 170) {
      delta = 1.0f + cos(theta * DEG_TO_RAD);

      energy = ka * delta * RAD_TO_DEG * RAD_TO_DEG;
    } else {
      delta = theta - theta0;
      delta2 = delta * delta;
    
      energy = ka * delta2;
    }
     
    return energy;
  }
  
  vector3 OBFFAngleCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, dc, gradient;
    double dE;

    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()))
      return VZero;
    
    da = a->GetVector();
    db = b->GetVector();
    dc = c->GetVector();
    theta = OBForceField::VectorAngleDerivative(da, db, dc);
    delta = theta - theta0;
 
    if (theta0 > 170) {
      if (atom == a) {
        gradient = ka * da;
        return gradient;
      } else if (atom == b) {
        gradient = ka * db;
        return gradient;
      } else {
        gradient = ka * dc;
        return gradient;
      }
    } else
      dE = 2.0f * ka * delta;

    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else if (atom == b) {
      gradient = dE * db; // - dE/drab * drab/db = - dE/drab * drab/da - dE/drab * drab/dc 
      return gradient;
    } else {
      gradient = dE * dc; // - dE/drab * drab/dc
      return gradient;
    }
  }
 
  double OBForceFieldGhemical::E_Angle()
  {
    vector<OBFFAngleCalculationGhemical>::iterator i;
    double energy;
    
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "A N G L E   B E N D I N G" << endl << endl;
      *logos << "ATOM TYPES       VALENCE     IDEAL      FORCE" << endl;
      *logos << " I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY" << endl;
      *logos << "-----------------------------------------------------------------------------" << endl;
    }
    
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); i++) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s %s  %8.3f   %8.3f     %8.3f   %8.3f   %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).theta, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        OBFFLog(logbuf);
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM {
      sprintf(logbuf, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }
    return energy;
  }
  
  double OBFFTorsionCalculationGhemical::GetEnergy()
  {
    double cosine, cosine2, cosine3;
    double phi1, phi2, phi3;
    vector3 vab, vbc, vcd, abbc, bccd;

    vab = a->GetVector() - b->GetVector();
    vbc = b->GetVector() - c->GetVector();
    vcd = c->GetVector() - d->GetVector();
    abbc = cross(vab, vbc);
    bccd = cross(vbc, vcd);
    tor = RAD_TO_DEG * acos(dot(abbc, bccd) / (abbc.length() * bccd.length()));
    if (dot(abbc, bccd) > 0.0f)
      tor = -tor;
    
    cosine = cos(DEG_TO_RAD * tor);
    cosine2 = cos(2.0f * DEG_TO_RAD * tor);
    cosine3 = cos(3.0f * DEG_TO_RAD * tor);

    phi1 = 1.0f + cosine;
    phi2 = 1.0f - cosine2;
    phi3 = 1.0f + cosine3;

    energy = k1 * phi1 + k2 * phi2 + k3 * phi3;

    return energy;
  }
  
  vector3 OBFFTorsionCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, dc, dd, gradient;
    double dE, sine, sine2, sine3;

    if ((atom->GetIdx() != a->GetIdx()) && (atom->GetIdx() != b->GetIdx()) && (atom->GetIdx() != c->GetIdx()) &&  (atom->GetIdx() != d->GetIdx()))
      return  VZero;
    
    da = a->GetVector();
    db = b->GetVector();
    dc = c->GetVector();
    dd = d->GetVector();
    tor = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
 
    sine = sin(DEG_TO_RAD * tor);
    sine2 = sin(2.0f * DEG_TO_RAD * tor);
    sine3 = sin(3.0f * DEG_TO_RAD * tor);
    dE = -k1 * sine + k2 * 2.0f * sine2 - k3 * 3.0f * sine3;
    
    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else if (atom == b) {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    } else if (atom == c) {
      gradient = dE * dc; // - dE/drab * drab/dc
      return gradient;
    } else {
      gradient = dE * dd; // - dE/drab * drab/dd
      return gradient;
    }
  }

  double OBForceFieldGhemical::E_Torsion() 
  {
    vector<OBFFTorsionCalculationGhemical>::iterator i;
    double energy;
 
    energy = 0.0f;

    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "T O R S I O N A L" << endl << endl;
      *logos << "----ATOM TYPES-----    FORCE              TORSION" << endl;
      *logos << " I    J    K    L     CONSTANT     s       ANGLE    n    ENERGY" << endl;
      *logos << "----------------------------------------------------------------" << endl;
    }
    
    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); i++) {

      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s %s %s    %6.3f    %5.0f   %8.3f   %1.0f   %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).d->GetType(), (*i).V, (*i).s, (*i).tor, (*i).n, (*i).energy);
        *logos << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos <<  "     TOTAL TORSIONAL ENERGY = " << energy << GetUnit() << endl;
    return energy;
  }

  double OBFFVDWCalculationGhemical::GetEnergy()
  {
    double term6, term12;

    rab = a->GetDistance(b);
    
    term12 = rab / ka;
    term6 = rab / kb;

    term12 = term12 * term12 * term12; // ^3
    term12 = term12 * term12; // ^6
    term12 = term12 * term12; // ^12
    term6 = term6 * term6 * term6; // ^3
    term6 = term6 * term6; // ^6
 
    energy = (1.0f / term12) - (1.0f / term6);

    if (is14 && samering)
      energy *= 2.0f;

    return energy;
  }
  
  vector3 OBFFVDWCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, gradient;
    double dE, term6, term7, term12, term13;

    if ((atom != a) && (atom != b))
      return  VZero;
    
    da = a->GetVector();
    db = b->GetVector();
    rab = OBForceField::VectorLengthDerivative(da, db);
    
    term12 = rab / ka;
    term6 = rab / kb;
    term13 = term12;
    term7 = term6;

    term12 = term12 * term12 * term12; // ^3
    term12 = term12 * term12; // ^6
    term12 = term12 * term12; // ^12
    term13 = term13 * term12; // ^13
    term6 = term6 * term6 * term6; // ^3
    term6 = term6 * term6; // ^6
    term7 = term7 * term6; // ^7

    dE = - (12.0f / ka) * (1.0f / term13) + (6.0f / kb) * (1.0f / term7);
    
    if (is14 && samering)
      dE *= 2.0f;

  
    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    }
  }


  double OBForceFieldGhemical::E_VDW()
  {
    vector<OBFFVDWCalculationGhemical>::iterator i;
    double energy;
 
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << endl << "V A N   D E R   W A A L S" << endl << endl;
      *logos << "ATOM TYPES          " << endl;
      *logos << " I    J        Rij       kij       ENERGY" << endl;
      *logos << "-----------------------------------------" << endl;
      //          XX   XX     -000.000  -000.000  -000.000  -000.000
    }
    
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); i++) {
      
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).kab, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos  << "     TOTAL VAN DER WAALS ENERGY = " << energy << GetUnit() << endl;
    return energy;
  }

  double OBFFElectrostaticCalculationGhemical::GetEnergy()
  {
    vector3 vab;

    vab = a->GetVector() - b->GetVector();
    rab = vab.length();

    energy = qq / rab;

    return energy;
  }
  
  vector3 OBFFElectrostaticCalculationGhemical::GetGradient(OBAtom *atom) 
  {
    vector3 da, db, gradient;
    double dE, rab2;

    if ((atom != a) && (atom != b))
      return  VZero;
     
    da = a->GetVector();
    db = b->GetVector();
    rab = OBForceField::VectorLengthDerivative(da, db);
    rab2 = rab * rab;
    
    dE = -qq / rab2;
    
    if (atom == a) {
      gradient = dE * da; // - dE/drab * drab/da
      return gradient;
    } else  {
      gradient = dE * db; // - dE/drab * drab/db
      return gradient;
    }
  }
  
  double OBForceFieldGhemical::E_Electrostatic()
  {
    vector<OBFFElectrostaticCalculationGhemical>::iterator i;
    double energy;
 
    energy = 0.0f;
    
    IF_OBFF_LOGLVL_HIGH {
      *logos << std::endl << "E L E C T R O S T A T I C   I N T E R A C T I O N S" << std::endl << std::endl;
      *logos << "ATOM TYPES          " << std::endl;
      *logos << " I    J           Rij   332.17*QiQj  ENERGY" << std::endl;
      *logos << "-------------------------------------------" << std::endl;
      //            XX   XX     -000.000  -000.000  -000.000  
    }

    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); i++) {
      
      energy += i->GetEnergy();
      
      IF_OBFF_LOGLVL_HIGH {
        sprintf(logbuf, "%s %s   %8.3f  %8.3f  %8.3f", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).qq, (*i).energy);
        *logos  << logbuf << endl;
      }
    }

    IF_OBFF_LOGLVL_MEDIUM
      *logos << "     TOTAL ELECTROSTATIC ENERGY = " << energy << GetUnit() << endl;
    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldGhemical theForceFieldGhemical("Ghemical", true);
  //***********************************************

  OBForceFieldGhemical::~OBForceFieldGhemical()
  {
  }

  OBForceFieldGhemical &OBForceFieldGhemical::operator=(OBForceFieldGhemical &src)
  {
    _mol = src._mol;

    _ffbondparams    = src._ffbondparams;
    _ffangleparams   = src._ffangleparams;
    _fftorsionparams = src._fftorsionparams;
    _ffvdwparams     = src._ffvdwparams;

    _bondcalculations          = src._bondcalculations;
    _anglecalculations         = src._anglecalculations;
    _torsioncalculations       = src._torsioncalculations;
    _vdwcalculations           = src._vdwcalculations;
    _electrostaticcalculations = src._electrostaticcalculations;

    return *this;
  }

  bool OBForceFieldGhemical::Setup(OBMol &mol)
  {
    _mol = mol;
    SetGhemicalTypes();
    
    if (!SetupCalculations())
      return false;
    
    return true;
  }
  
  bool OBForceFieldGhemical::SetupCalculations()
  {
    OBFFParameter *parameter;
    OBAtom *a, *b, *c, *d;
    bool found;
    
    // 
    // Bond Calculations
    //
    OBFFBondCalculationGhemical bondcalc;
    int bondtype;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = bond->GetBondOrder(); 
      if (bond->IsAromatic())
        bondtype = 5;

      bondcalc.a = a;
      bondcalc.b = b;
      bondcalc.bt = bondtype;

      parameter = GetParameterGhemical(bondtype, a->GetType(), b->GetType(), NULL, NULL,  _ffbondparams);
      if (parameter == NULL) {
        parameter = GetParameterGhemical(bondtype, "FFFF", a->GetType(), NULL, NULL, _ffbondparams);

        if (parameter == NULL) {
          parameter = GetParameterGhemical(bondtype, "FFFF", b->GetType(), NULL, NULL, _ffbondparams);

          if (parameter == NULL) {
            //obErrorLog.ThrowError(__FUNCTION__, "Could not find all bond parameters ", obError);
            //return false;
            bondcalc.kb = KCAL_TO_KJ * 500.0f;
            bondcalc.r0 = 1.100f;

            _bondcalculations.push_back(bondcalc);
            continue;
          }
        }
      }
      bondcalc.kb = KCAL_TO_KJ * parameter->dpar2;
      bondcalc.r0 = parameter->dpar1;

      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    OBFFAngleCalculationGhemical anglecalc;
 
    _anglecalculations.clear();
    
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameter = GetParameter(a->GetType(), b->GetType(), c->GetType(), _ffangleparams);
      if (parameter == NULL) {
        parameter = GetParameter("FFFF", b->GetType(), c->GetType(), _ffangleparams);
        if (parameter == NULL) {
          parameter = GetParameter(a->GetType(), b->GetType(), "FFFF", _ffangleparams);
          if (parameter == NULL) {
            parameter = GetParameter("FFFF", b->GetType(), "FFFF", _ffangleparams);
            if (parameter == NULL) {
              anglecalc.ka = KCAL_TO_KJ * 0.020f;
              anglecalc.theta0 = 120.0f;
            
              _anglecalculations.push_back(anglecalc);
            
	      IF_OBFF_LOGLVL_LOW
                *logos  << "COULD NOT FIND PARAMETERS FOR ANGLE " << a->GetType() << "-" << b->GetType() << "-" << c->GetType() << ", USING DEFAULT PARAMETERS" << endl;

              continue;
            }
	  }
        }
      }
      anglecalc.ka = KCAL_TO_KJ * parameter->dpar2;
      anglecalc.theta0 = parameter->dpar1;
      
      _anglecalculations.push_back(anglecalc);
    }
    
    //
    // Torsion Calculations
    //
    OBFFTorsionCalculationGhemical torsioncalc;
    int torsiontype;
    int s;

    _torsioncalculations.clear();
 
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);
      OBBond *bc = _mol.GetBond(b, c);
      torsiontype = bc->GetBondOrder(); 
      if (bc->IsAromatic())
        torsiontype = 5;
      
      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
      torsioncalc.tt = torsiontype;

      parameter = GetParameterGhemical(torsiontype, a->GetType(), b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
      if (parameter == NULL) {
        parameter = GetParameterGhemical(torsiontype, "FFFF", b->GetType(), c->GetType(), d->GetType(), _fftorsionparams);
        if (parameter == NULL) {
          parameter = GetParameterGhemical(torsiontype, a->GetType(), b->GetType(), c->GetType(), "FFFF", _fftorsionparams);
          if (parameter == NULL) {
            parameter = GetParameterGhemical(torsiontype, "FFFF", b->GetType(), c->GetType(), "FFFF", _fftorsionparams);
            if (parameter == NULL) {
              torsioncalc.V = 0.0f;
              torsioncalc.s = 1.0f;
              torsioncalc.n = 1.0f;
              
	      torsioncalc.k1 = 0.0f;
	      torsioncalc.k2 = 0.0f;
	      torsioncalc.k3 = 0.0f;
              _torsioncalculations.push_back(torsioncalc);
              continue;
              //obErrorLog.ThrowError(__FUNCTION__, "Could not find all torsion parameters", obError);
              //return false;
            }
          }
        }
      }
      torsioncalc.V = KCAL_TO_KJ * parameter->dpar1;
      torsioncalc.s = parameter->dpar2;
      torsioncalc.n = parameter->dpar3;

      s = (int) (torsioncalc.s * torsioncalc.n);
      switch(s) {
        case +3:
          torsioncalc.k1 = 0.0f;
          torsioncalc.k2 = 0.0f;
	  torsioncalc.k3 = torsioncalc.V;
	  break;
        case +2:
          torsioncalc.k1 = 0.0f;
	  torsioncalc.k2 = -torsioncalc.V;
          torsioncalc.k3 = 0.0f;
	  break;
        case +1:
	  torsioncalc.k1 = torsioncalc.V;
          torsioncalc.k2 = 0.0f;
          torsioncalc.k3 = 0.0f;
	  break;
        case -1:
	  torsioncalc.k1 = -torsioncalc.V;
          torsioncalc.k2 = 0.0f;
          torsioncalc.k3 = 0.0f;
	  break;
        case -2:
          torsioncalc.k1 = 0.0f;
	  torsioncalc.k2 = torsioncalc.V;
          torsioncalc.k3 = 0.0f;
	  break;
        case -3:
          torsioncalc.k1 = 0.0f;
          torsioncalc.k2 = 0.0f;
	  torsioncalc.k3 = -torsioncalc.V;
	  break;
      }

      _torsioncalculations.push_back(torsioncalc);     
    }
    
    // 
    // VDW Calculations
    //
    OBFFVDWCalculationGhemical vdwcalc;
    OBFFParameter *parameter_a, *parameter_b;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      parameter_a = GetParameter(a->GetType(), _ffvdwparams);
      if (parameter_a == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Ra = 1.5f;
        vdwcalc.ka = 0.042f;
      } else {
        vdwcalc.Ra = parameter_a->dpar1;
        vdwcalc.ka = parameter_a->dpar2;
      }

      parameter_b = GetParameter(b->GetType(), _ffvdwparams);
      if (parameter_b == NULL) { // no vdw parameter -> use hydrogen
        vdwcalc.Rb = 1.5f;
        vdwcalc.kb = 0.042;;
      } else {
        vdwcalc.Rb = parameter_b->dpar1;
        vdwcalc.kb = parameter_b->dpar2;
      }

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;
     
      //this calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      vdwcalc.kab = KCAL_TO_KJ * sqrt(vdwcalc.ka * vdwcalc.kb);
      
      // 1-4 scaling
      vdwcalc.is14 = false;
      FOR_NBORS_OF_ATOM (nbr, a)
        FOR_NBORS_OF_ATOM (nbr2, &*nbr)
          FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
            if (b == &*nbr3) {
              vdwcalc.is14 = true;
              vdwcalc.kab *= 0.5f;
            }
      
      // not sure why this is needed, but validation showed it works...
      if (a->IsInRingSize(6) && b->IsInRingSize(6) && IsInSameRing(a, b))
        vdwcalc.samering = true;
      else if ((a->IsInRingSize(5) || a->IsInRingSize(4)) && (b->IsInRingSize(5) || b->IsInRingSize(4)))
        vdwcalc.samering = true;
      else
        vdwcalc.samering = false;


      vdwcalc.ka = (vdwcalc.Ra + vdwcalc.Rb) * pow(1.0f * vdwcalc.kab , 1.0f / 12.0f);
      vdwcalc.kb = (vdwcalc.Ra + vdwcalc.Rb) * pow(2.0f * vdwcalc.kab , 1.0f / 6.0f);
      
      _vdwcalculations.push_back(vdwcalc);
    }
    
    // 
    // Electrostatic Calculations
    //
    OBFFElectrostaticCalculationGhemical elecalc;

    _electrostaticcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      elecalc.qq = KCAL_TO_KJ * 332.17f * a->GetPartialCharge() * b->GetPartialCharge();
      
      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;
        
	// 1-4 scaling
        FOR_NBORS_OF_ATOM (nbr, a)
          FOR_NBORS_OF_ATOM (nbr2, &*nbr)
            FOR_NBORS_OF_ATOM (nbr3, &*nbr2)
              if (b == &*nbr3)
                elecalc.qq *= 0.5f;

        _electrostaticcalculations.push_back(elecalc);
      }
    }

    return true;
  }

  bool OBForceFieldGhemical::ParseParamFile()
  {
    vector<string> vs;
    char buffer[80];
    
    OBFFParameter parameter;

    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
      {
        ifs2.open(buffer2.c_str());
        ifsP = &ifs2;
      }

    while (ifsP->getline(buffer, 80)) {
      tokenize(vs, buffer);

      if (EQn(buffer, "bond", 4)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter.dpar1 = atof(vs[4].c_str()); // length
        parameter.dpar2 = atof(vs[5].c_str()); // force cte
        if (EQn(vs[3].c_str(), "S", 1))
          parameter.ipar5 = 1;
        if (EQn(vs[3].c_str(), "D", 1))
          parameter.ipar5 = 2;
        if (EQn(vs[3].c_str(), "T", 1))
          parameter.ipar5 = 3;
        if (EQn(vs[3].c_str(), "C", 1))
          parameter.ipar5 = 5;
        _ffbondparams.push_back(parameter);
      }
      if (EQn(buffer, "angle", 5)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter.dpar1 = atof(vs[5].c_str()); // angle
        parameter.dpar2 = atof(vs[6].c_str()); // force cte
        _ffangleparams.push_back(parameter);
      }
      if (EQn(buffer, "torsion", 7)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        parameter._c = vs[3];
        parameter._d = vs[4];
        parameter.dpar1 = atof(vs[6].c_str()); // force cte
        parameter.dpar3 = atof(vs[8].c_str()); // n
        if (EQn(vs[7].c_str(), "+", 1))
          parameter.dpar2 = +1; // s
        if (EQn(vs[7].c_str(), "-", 1))
          parameter.dpar2 = -1; // s
        if (EQn(vs[5].c_str(), "?S?", 3))
          parameter.ipar5 = 1;
        if (EQn(vs[5].c_str(), "?D?", 3))
          parameter.ipar5 = 2;
        if (EQn(vs[5].c_str(), "?T?", 3))
          parameter.ipar5 = 3;
        if (EQn(vs[5].c_str(), "?C?", 3))
          parameter.ipar5 = 5;
        _fftorsionparams.push_back(parameter);
      }
      if (EQn(buffer, "vdw", 3)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter.dpar1 = atof(vs[2].c_str()); // r
        parameter.dpar2 = atof(vs[3].c_str()); // force cte
        _ffvdwparams.push_back(parameter);
      }
      if (EQn(buffer, "charge", 6)) {
        parameter.clear();
        parameter._a = vs[1];
        parameter._b = vs[2];
        if (EQn(vs[3].c_str(), "S", 1))
          parameter.ipar5 = 1;
        if (EQn(vs[3].c_str(), "D", 1))
          parameter.ipar5 = 2;
        parameter.dpar1 = atof(vs[4].c_str()); // charge
        _ffchargeparams.push_back(parameter);
      }
    }
	
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
 
    return 0;
  }
  
  bool OBForceFieldGhemical::SetGhemicalTypes()
  {
    std::vector<std::vector<int> > _mlist; //!< match list for atom typing
    std::vector<std::pair<OBSmartsPattern*,std::string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[80];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/ghemical.prm
    string buffer2, subbuffer;
    ifstream ifs1, ifs2, *ifsP;
    buffer2 = BABEL_DATADIR;
    buffer2 += FILE_SEP_CHAR;
    subbuffer = buffer2;
    subbuffer += BABEL_VERSION;
    subbuffer += FILE_SEP_CHAR;
    subbuffer += "ghemical.prm";
    buffer2 += "ghemical.prm";

    ifs1.open(subbuffer.c_str());
    ifsP= &ifs1;
    if (!(*ifsP))
      {
        ifs2.open(buffer2.c_str());
        ifsP = &ifs2;
      }
   
    while (ifsP->getline(buffer, 80)) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1]))
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse EXTTYP line in atom type table from atomtyp.txt", obInfo);
          return false;
        }
        
        for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
          if (i->first->Match(_mol)) {
            _mlist = i->first->GetMapList();
            for (j = _mlist.begin();j != _mlist.end();++j) {
              _mol.GetAtom((*j)[0])->SetType(i->second);
            }
          }
        }
      }
    }

    SetGhemicalCharges();
 
    IF_OBFF_LOGLVL_LOW {
      *logos  << endl << "A T O M   T Y P E S" << endl << endl;
      *logos << "IDX\tTYPE" << endl;
      
      FOR_ATOMS_OF_MOL (a, _mol)
        *logos << a->GetIdx() << "\t" << a->GetType() << endl;

      *logos  << endl << "C H A R G E S" << endl << endl;
      *logos << "IDX\tCHARGE" << endl;
      
      FOR_ATOMS_OF_MOL (a, _mol)
        *logos << a->GetIdx() << "\t" << a->GetPartialCharge() << endl;
 
    }
    
    // DEBUG (validation)
    //FOR_ATOMS_OF_MOL (a, _mol)
    //  if (atoi(a->GetType()) != 0)
    //    cout << "ATOMTYPE " << atoi(a->GetType()) << endl;
    //  else
    //    cout << "ATOMTYPE " << a->GetType() << endl;
 
    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();

    return true;
  }
  
  bool OBForceFieldGhemical::SetGhemicalCharges()
  {
    OBAtom *a, *b;
    int bondtype;

    _mol.SetAutomaticPartialCharge(false);
    _mol.SetPartialChargesPerceived();

    // set all partial charges to 0.0
    FOR_ATOMS_OF_MOL (atom, _mol)
      atom->SetPartialCharge(0.0f);

    FOR_BONDS_OF_MOL (bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();	
      bondtype = bond->GetBondOrder(); 

      string _a(a->GetType());
      string _b(b->GetType());

      for (unsigned int idx=0; idx < _ffchargeparams.size(); idx++) {
        if (((_a == _ffchargeparams[idx]._a) && (_b == _ffchargeparams[idx]._b)) && (bondtype == _ffchargeparams[idx].ipar5)) {
          a->SetPartialCharge(a->GetPartialCharge() - _ffchargeparams[idx].dpar1);
	  b->SetPartialCharge(b->GetPartialCharge() + _ffchargeparams[idx].dpar1);
	} else if (((_a == _ffchargeparams[idx]._b) && (_b == _ffchargeparams[idx]._a)) && (bondtype == _ffchargeparams[idx].ipar5)) {
          a->SetPartialCharge(a->GetPartialCharge() + _ffchargeparams[idx].dpar1);
	  b->SetPartialCharge(b->GetPartialCharge() - _ffchargeparams[idx].dpar1);
	}
      }
    }

    return true;
  }

  double OBForceFieldGhemical::Energy()
  {
    double energy;
   
    energy = E_Bond();
    energy += E_Angle();
    energy += E_Torsion();
    energy += E_VDW();
    energy += E_Electrostatic();

    IF_OBFF_LOGLVL_HIGH {
      sprintf(logbuf, "\nTOTAL ENERGY = %8.3f %s", energy, GetUnit().c_str());
      OBFFLog(logbuf);
    }

    return energy;
  }
  
  OBFFParameter* OBForceFieldGhemical::GetParameterGhemical(int type, const char* a, const char* b, const char* c, const char* d, 
                                                            vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    if (a == NULL)
      return NULL;

    if (b == NULL) {
      string _a(a);
      for (unsigned int idx=0; idx < parameter.size(); idx++) 
        if ((_a == parameter[idx]._a) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      return NULL;
    }
    if (c == NULL) {
      string _a(a);
      string _b(b);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b)) && (type == parameter[idx].ipar5) || 
            ((_a == parameter[idx]._b) && (_b == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      }
      return NULL;
    }
    if (d == NULL) {
      string _a(a);
      string _b(b);
      string _c(c);
      for (unsigned int idx=0; idx < parameter.size(); idx++) {
        if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c)) && (type == parameter[idx].ipar5)|| 
            ((_a == parameter[idx]._c) && (_b == parameter[idx]._b) && (_c == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
          par = &parameter[idx];
          return par;
        }
      }
      return NULL;
    }
    string _a(a);
    string _b(b);
    string _c(c);
    string _d(d);

    for (unsigned int idx=0; idx < parameter.size(); idx++) {
      if (((_a == parameter[idx]._a) && (_b == parameter[idx]._b) && (_c == parameter[idx]._c) && 
           (_d == parameter[idx]._d)) && (type == parameter[idx].ipar5) || 
          ((_a == parameter[idx]._d) && (_b == parameter[idx]._c) && (_c == parameter[idx]._b) && 
           (_d == parameter[idx]._a)) && (type == parameter[idx].ipar5)) {
        par = &parameter[idx];
        return par;
      }
    }

    return NULL;
  }
  
  vector3 OBForceFieldGhemical::GetGradient(OBAtom *a, int terms)
  {
    vector<OBFFBondCalculationGhemical>::iterator i;
    vector<OBFFAngleCalculationGhemical>::iterator i2;
    vector<OBFFTorsionCalculationGhemical>::iterator i3;
    vector<OBFFVDWCalculationGhemical>::iterator i4;
    vector<OBFFElectrostaticCalculationGhemical>::iterator i5;

    vector3 grad(0.0f, 0.0f, 0.0f);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EBOND))
      for (i = _bondcalculations.begin(); i != _bondcalculations.end(); i++)
        if (((*i).a->GetIdx() == a->GetIdx()) || ((*i).b->GetIdx() == a->GetIdx()))
          grad += i->GetGradient(&*a);
    

    if ((terms & OBFF_ENERGY) || (terms & OBFF_EANGLE))
      for (i2 = _anglecalculations.begin(); i2 != _anglecalculations.end(); i2++)
        if (((*i2).a->GetIdx() == a->GetIdx()) || ((*i2).b->GetIdx() == a->GetIdx()) || ((*i2).c->GetIdx() == a->GetIdx()))
          grad += i2->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_ETORSION))
      for (i3 = _torsioncalculations.begin(); i3 != _torsioncalculations.end(); i3++)
        if (((*i3).a->GetIdx() == a->GetIdx()) || ((*i3).b->GetIdx() == a->GetIdx()) || ((*i3).c->GetIdx() == a->GetIdx()) || ((*i3).d->GetIdx() == a->GetIdx()))
          grad += i3->GetGradient(&*a);
      
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EVDW))
      for (i4 = _vdwcalculations.begin(); i4 != _vdwcalculations.end(); i4++)
        if (((*i4).a->GetIdx() == a->GetIdx()) || ((*i4).b->GetIdx() == a->GetIdx()))
          grad += i4->GetGradient(&*a);
    
    if ((terms & OBFF_ENERGY) || (terms & OBFF_EELECTROSTATIC))
      for (i5 = _electrostaticcalculations.begin(); i5 != _electrostaticcalculations.end(); i5++)
        if (((*i5).a->GetIdx() == a->GetIdx()) || ((*i5).b->GetIdx() == a->GetIdx()))
          grad += i5->GetGradient(&*a);

    return grad;
  }

  bool OBForceFieldGhemical::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    
    OBFFLog("\nV A L I D A T E   G R A D I E N T S\n\n");
    OBFFLog("ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERRROR (%)   \n");
    OBFFLog("----------------------------------------------------------------------------------------\n");
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"
   
    FOR_ATOMS_OF_MOL (a, _mol) {

      // OBFF_ENERGY
      numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      anagrad = GetGradient(&*a, OBFF_ENERGY);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      anagrad = GetGradient(&*a, OBFF_EBOND);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      *logos << logbuf << endl;
      OBFFLog(logbuf);
      
      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      anagrad = GetGradient(&*a, OBFF_EANGLE);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      anagrad = GetGradient(&*a, OBFF_ETORSION);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      anagrad = GetGradient(&*a, OBFF_EVDW);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      anagrad = GetGradient(&*a, OBFF_EELECTROSTATIC);
      err = ValidateGradientError(numgrad, anagrad);

      sprintf(logbuf, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(logbuf);
    }
  }

} // end namespace OpenBabel

//! \file forcefieldghemical.cpp
//! \brief Ghemical force field
