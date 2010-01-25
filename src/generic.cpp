/**********************************************************************
generic.cpp - Handle OBGenericData classes.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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

#include <string>

#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/math/matrix3x3.h>

// needed for msvc to have at least one reference to AtomClass, AliasData in openbabel library
#include <openbabel/atomclass.h>
#include <openbabel/alias.h>

using namespace std;

namespace OpenBabel
{

  /** \class OBGenericData generic.h <openbabel/generic.h>

      OBGenericData is an abstract base class which defines an interface for
      storage, retrieval, and indexing of arbitrary generic data.
      Subclasses of OBGenericData can be used to store custom data
      on a per-atom, per-bond, per-molecule, or per-residue basis.
      Open Babel currently supports a small subset of chemical functionality
      as OBGenericData types, which will expand over time to support additional
      interconversion (e.g., spectroscopy, dynamics, surfaces...)

      For more information on currently supported types, please see
      the developer wiki:
      http://openbabel.sourceforge.net/wiki/Generic_Data

      For your own custom data, either define a custom subclass using 
      an id from the OBGenericDataType::CustomData0 to 
      OBGenericDataType::CustomData15 slots,
      or store your data as a string and use OBPairData for key/value access.
      The latter is <strong>highly</strong> recommended for various text 
      descriptors
      e.g., in QSAR, atom or bond labels, or other textual data.

      <strong>New in Open Babel, version 2.1</strong>
      is the template-based OBPairTemplate,
      which can be used to store arbitrary data types. There are predefined
      types OBPairInteger and OBPairFloatingPoint for storing integers and
      floating-point values without converting to a string representation.

      Also <strong>new</strong> is the "source" or "origin" of a data
      entry, enumerated by DataOrigin. This can be accessed by
      SetOrigin() and GetOrigin(), as well as via "filtering" methods
      in OBBase, allowing you to separate data read in from a file,
      added by a user, or assigned by Open Babel internally.

      While the library and import routines will set DataOrigin correctly, 
      you should try to annotate data added by your code. Typically this would
      either be userAdded or external. The former refers to something the
      user requested as an annotation, while the latter refers to annotations
      your code adds automatically.

      Example code using OBGenericData:

      @code
      if (mol.HasData(OBGenericDataType::UnitCell))
      {
         uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
         sprintf(buffer,
            "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f",
            uc->GetA(), uc->GetB(), uc->GetC(),
            uc->GetAlpha() , uc->GetBeta(), uc->GetGamma());
         ofs << buffer << endl;
      }

      ...

      vector<OBGenericData*>::iterator k;
      vector<OBGenericData*> vdata = mol.GetData();
      for (k = vdata.begin();k != vdata.end();++k)
         if ((*k)->GetDataType() == OBGenericDataType::PairData)
         {
            ofs << ">  <" << (*k)->GetAttribute() << ">" << endl;
            ofs << ((OBPairData*)(*k))->GetValue() << endl << endl;
         }
      @endcode

      Similar code also works for OBGenericData stored in an OBAtom or 
      OBBond or OBResidue. These examples show use of DataOrigin outside
      of the Open Babel library.

      @code
      string atomLabel; // e.g., from the user adding annotation to an atom
      if (!atom.HasData("UserLabel")) // stored textual data as an OBPairData
      {
         OBPairData *label = new OBPairData;
         label->SetAttribute("UserLabel");
         label->SetValue(atomLabel);
         label->SetOrigin(userInput); // set by user, not by Open Babel

         atom.SetData(label);
      }

      ...

      if (bond.HasData("DisplayType")) // e.g. in a visualization tool
      {
         OBPairData *display = dynamic_cast<OBPairData *> bond.GetData("DisplayType");
         if (display->GetValue() == "wireframe")
         {
            ... // display a wireframe view
         }
      }
      @endcode

      When designing a class derived from OBGenericData you must add a 
      Clone() function. For classes used with OBMol this is used when 
      an OBMol object is copied. If your class member variables contain
      pointers to atoms or bonds then it will be necessary to ensure
      that these are updated in Clone() to refer to the new molecule. Without
      these and similar pointers it is more likely that the very simple 
      clone function
      @code
      virtual OBGenericData* Clone(OBBase* parent) const
         {return new MyNewClass(*this);}
      @endcode
      and the compiler generated copy constructor would be sufficient. 

      It is recommended that, if possible, OBGenericData classes do not
      store atom and bond pointers. Using atom and bond indices instead
      would allow the simple version of Clone() above. See 
      OBRotameterData::Clone for an example of a more complicated version.
      For classes which are not intended to support copying, Clone() can 
      return NULL 
      @code
      virtual OBGenericData* Clone(OBBase* parent) const 
         {return NULL;}
      @endcode
      Clone() is a pure virtual function so that you need to decide what
      kind of function you need and include it explicitly.
  **/

  //
  //member functions for OBGenericData class
  //

  OBGenericData::OBGenericData(const std::string attr, const unsigned int type,
                               const DataOrigin  source):
    _attr(attr), _type(type), _source(source)
  { }

  /* Use default copy constructor and assignment operators
     OBGenericData::OBGenericData(const OBGenericData &src)
     {
     _type = src.GetDataType();
     _attr = src.GetAttribute();
     }


     OBGenericData& OBGenericData::operator = (const OBGenericData &src)
     {
     if(this == &src)
     return(*this);

     _type = src._type;
     _attr = src._attr;

     return(*this);
     }
  */

  //
  //member functions for OBCommentData class
  //

  OBCommentData::OBCommentData():
    OBGenericData("Comment", OBGenericDataType::CommentData)
  { }

  OBCommentData::OBCommentData(const OBCommentData &src) :
    OBGenericData(src), _data(src._data)
  {  }

  //
  //member functions for OBExternalBond class
  //
  OBExternalBond::OBExternalBond(OBAtom *atom,OBBond *bond,int idx):
    _idx(idx), _atom(atom), _bond(bond)
  {  }

  OBExternalBond::OBExternalBond(const OBExternalBond &src):
    _idx(src._idx), _atom(src._atom), _bond(src._bond)
  { }

  //
  //member functions for OBExternalBondData class
  //

  OBExternalBondData::OBExternalBondData():
    OBGenericData("ExternalBondData", OBGenericDataType::ExternalBondData,
                  perceived)
  { }

  void OBExternalBondData::SetData(OBAtom *atom,OBBond *bond,int idx)
  {
    OBExternalBond xb(atom,bond,idx);
    _vexbnd.push_back(xb);
  }

  //
  //member functions for OBPairData class
  //

  OBPairData::OBPairData() :
    OBGenericData("PairData", OBGenericDataType::PairData)
  { }

  //
  //member functions for OBVirtualBond class
  //

  OBVirtualBond::OBVirtualBond():
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    _bgn(0), _end(0), _ord(0), _stereo(0)
  {  }

  OBVirtualBond::OBVirtualBond(int bgn,int end,int ord,int stereo):
    OBGenericData("VirtualBondData", OBGenericDataType::VirtualBondData, perceived),
    _bgn(bgn), _end(end), _ord(ord), _stereo(stereo)
  {  }

  //
  // member functions for OBUnitCell class
  //
  OBUnitCell::OBUnitCell():
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    _a(0.0), _b(0.0), _c(0.0), _alpha(0.0), _beta(0.0), _gamma(0.0),
    _spaceGroup( NULL ), _lattice(Undefined)
  {  
    // We should default to P1 space group unless we know differently
    SetSpaceGroup(1);
  }

  OBUnitCell::OBUnitCell(const OBUnitCell &src) :
    OBGenericData("UnitCell", OBGenericDataType::UnitCell),
    _a(src._a), _b(src._b), _c(src._c), 
    _alpha(src._alpha), _beta(src._beta), _gamma(src._gamma),
    _offset(src._offset),
    _v1(src._v1), _v2(src._v2), _v3(src._v3),
    _spaceGroupName(src._spaceGroupName),
    _spaceGroup(src._spaceGroup),
    _lattice(src._lattice)
  {  }

  OBUnitCell & OBUnitCell::operator=(const OBUnitCell &src)
  {
    if(this == &src)
      return(*this);

    _a = src._a;
    _b = src._b;
    _c = src._c;
    _alpha = src._alpha;
    _beta = src._beta;
    _gamma = src._gamma;
    _offset = src._offset;

    _v1 = src._v1;
    _v2 = src._v2;
    _v3 = src._v3;

    _spaceGroup = src._spaceGroup;
    _spaceGroupName = src._spaceGroupName;
    _lattice = src._lattice;

    return(*this);
  }

	/*!
  ** The angles and lengths of the unitcell will be calculated from the
  ** vectors @p v1, @p v2 and @p v3. Those vectors will as well be
  ** stored internally.
  **Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoNotionalCoordinates">blue-obelisk:convertCartesianIntoNotionalCoordinates</a>
  **\brief Sets the vectors, angles and lengths of the unitcell
  **\param v1 The x-vector
  **\param v2 The y-vector
  **\param v3 The z-vector
  **\see OBUnitCell::GetCellVectors
  */
  void OBUnitCell::SetData(const vector3 v1, const vector3 v2, const vector3 v3)
  {
    _v1 = v1;
    _v2 = v2;
    _v3 = v3;

    _a = _v1.length();
    _b = _v2.length();
    _c = _v3.length();

    // For PR#1961604 -- somewhat contrived example
    if (IsNearZero(_a) && !IsNearZero(_c)) {
      _v1 = _v3; // we'll reset _v3 below
      _a = _c;
      _c = 0.0;
    }

    // Sanity checks for 1D or 2D translation
    if (IsNearZero(_b)) { // 1D
      _v2.Set(v1.y(), -v1.x(), v1.z()); // rotate base vector by 90 degrees
      _b = 999.999;
      _v2 = _b * _v2.normalize(); // set to a large displacement
    }

    if (IsNearZero(_c)) { // 2D or 1D
      _v3 = cross(_v1, _v2);
      _c = 999.999;
      _v3 = _c * _v3.normalize(); // set to a large displacement
    }
    
    _alpha = vectorAngle(_v2, _v3);
    _beta =  vectorAngle(_v1, _v3);
    _gamma = vectorAngle(_v1, _v2);
  }

  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertNotionalIntoCartesianCoordinates">blue-obelisk:convertNotionalIntoCartesianCoordinates</a>
  vector<vector3> OBUnitCell::GetCellVectors()
  {
    vector<vector3> v;
    v.reserve(3);

    // no unit cell vectors
    if (IsNegligible(_v1.length(), 1.0, 1.0e-9) &&
        IsNegligible(_v2.length(), 1.0, 1.0e-9) &&
        IsNegligible(_v3.length(), 1.0, 1.0e-9))
      {
        vector3 temp;
        matrix3x3 m = GetOrthoMatrix();

        temp = vector3(1.0, 0.0, 0.0);
        v.push_back(m * temp);
        temp = vector3(0.0, 1.0, 0.0);
        v.push_back(m * temp);
        temp = vector3(0.0, 0.0, 1.0);
        v.push_back(m * temp);
      }
    else
      {
        v.push_back(_v1);
        // we set these above in case we had a 1D or 2D translation vector system
        if (fabs(_b - 999.999) > 1.0e-1)
          v.push_back(_v2);
        if (fabs(_c - 999.999) > 1.0e-1)
          v.push_back(_v3);
      }

    return v;
  }

  matrix3x3 OBUnitCell::GetCellMatrix()
  {
    matrix3x3 m;

    if (IsNegligible(_v1.length(), 1.0, 1.0e-9) &&
        IsNegligible(_v2.length(), 1.0, 1.0e-9) &&
        IsNegligible(_v3.length(), 1.0, 1.0e-9))
      {
        m = GetOrthoMatrix();
      }
    else
      {
        vector3 v1, v2, v3;
        v1 = _v1;
        v2 = _v2;
        v3 = _v3;
        m = matrix3x3(v1,v2,v3);
      }
    return m;
  }

  // Convert from fractional to Cartesian
  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#calculateOrthogonalisationMatrix">blue-obelisk:calculateOrthogonalisationMatrix</a>
  matrix3x3 OBUnitCell::GetOrthoMatrix()
  {
    matrix3x3 m;
    if (IsNearZero(_c) || IsNearZero(_b)) {
      // 1D or 2D unit cell
    }

    // already here, let's not duplicate the work
    m.FillOrth(_alpha, _beta, _gamma, _a, _b, _c);

    return m;
  }

  // Based on code in PyMMLib: http://pymmlib.sf.net/
  //! Matrix to convert from Cartesian to fractional
  //! Implements <a href="http://qsar.sourceforge.net/dicts/blue-obelisk/index.xhtml#convertCartesianIntoFractionalCoordinates">blue-obelisk:convertCartesianIntoFractionalCoordinates</a> 
  matrix3x3 OBUnitCell::GetFractionalMatrix()
  {
    matrix3x3 m;
    double sinAlpha, sinBeta, sinGamma;
    double cosAlpha, cosBeta, cosGamma;
    double v;

    sinAlpha = sin(_alpha * DEG_TO_RAD);
    sinBeta = sin(_beta * DEG_TO_RAD);
    sinGamma = sin(_gamma * DEG_TO_RAD);
    cosAlpha = cos(_alpha * DEG_TO_RAD);
    cosBeta = cos(_beta * DEG_TO_RAD);
    cosGamma = cos(_gamma * DEG_TO_RAD);

    v = sqrt(1 - SQUARE(cosAlpha) - SQUARE(cosBeta) - SQUARE(cosGamma) +
             2 * cosAlpha*cosBeta*cosGamma);

    m.Set(0,0,  1.0 / _a);
    m.Set(0,1,  -cosGamma / (_a * sinGamma) );
    m.Set(0,2,  (cosGamma * cosAlpha - cosBeta) / (_a * v * sinGamma) );
    m.Set(1,0,  0.0);
    m.Set(1,1,  1.0 / (_b * sinGamma) );
    m.Set(1,2,  (cosGamma * cosBeta - cosAlpha) / (_b * v * sinGamma) );
    m.Set(2,0,  0.0);
    m.Set(2,1,  0.0);
    m.Set(2,2,  sinGamma / (_c * v) );

    return m;
  }

  OBUnitCell::LatticeType OBUnitCell::GetLatticeType( int spacegroup )
  {
	  //	1-2 	Triclinic
	  //	3-15 	Monoclinic
	  //	16-74	Orthorhombic
	  //	75-142 	Tetragonal
	  //	143-167 Rhombohedral
	  //	168-194 Hexagonal
	  //	195-230 Cubic

      if ( spacegroup == 0  && _spaceGroup)
          spacegroup = _spaceGroup->GetId();
	  
	  if ( spacegroup <= 0 )
		  return OBUnitCell::Undefined;

	  else if ( spacegroup == 1 ||
              spacegroup == 2 )
		  return OBUnitCell::Triclinic;
	  
	  else if ( spacegroup >= 3 &&
              spacegroup <= 15 )
		  return OBUnitCell::Monoclinic;
	  
	  else if ( spacegroup >= 16 &&
              spacegroup <= 74 )
		  return OBUnitCell::Orthorhombic;
	  
	  else if ( spacegroup >= 75 &&
              spacegroup <= 142 )
		  return OBUnitCell::Tetragonal;
	  
	  else if ( spacegroup >= 143 &&
              spacegroup <= 167 )
		  return OBUnitCell::Rhombohedral;
	  
	  else if ( spacegroup >= 168 &&
              spacegroup <= 194 )
		  return OBUnitCell::Hexagonal;
	  
	  else if ( spacegroup >= 195 &&
              spacegroup <= 230 )
		  return OBUnitCell::Cubic;

	  //just to be extra sure
	  else // ( spacegroup > 230 )
		  return OBUnitCell::Undefined;
  }
  
  OBUnitCell::LatticeType OBUnitCell::GetLatticeType()
  {
    if (_lattice != Undefined)
      return _lattice;
    else if (_spaceGroup != NULL)
      return GetLatticeType(_spaceGroup->GetId());

    unsigned int rightAngles = 0;
    if (IsApprox(_alpha, 90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(_beta,  90.0, 1.0e-3)) rightAngles++;
    if (IsApprox(_gamma, 90.0, 1.0e-3)) rightAngles++;

    switch (rightAngles)
      {
      case 3:
        if (IsApprox(_a, _b, 1.0e-4) && IsApprox(_b, _c, 1.0e-4))
          _lattice = Cubic;
        else if (IsApprox(_a, _b, 1.0e-4) || IsApprox(_b, _c, 1.0e-4))
          _lattice = Tetragonal;
        else
          _lattice = Orthorhombic;
        break;
      case 2:
        if ( (IsApprox(_alpha, 120.0, 1.0e-3) 
              || IsApprox(_beta, 120.0, 1.0e-3) 
              || IsApprox(_gamma, 120.0f, 1.0e-3))
             && (IsApprox(_a, _b, 1.0e-4) || IsApprox(_b, _c, 1.0e-4)) )
          _lattice = Hexagonal;
        else
          _lattice = Monoclinic;
        break;
      default:
        if (IsApprox(_a, _b, 1.0e-4) && IsApprox(_b, _c, 1.0e-4))
          _lattice = Rhombohedral;
        else
          _lattice = Triclinic;
      }

    return _lattice;
  }

  int OBUnitCell::GetSpaceGroupNumber( std::string name)
  {
    static const char * const spacegroups[] = { 
      "P1", "P-1", "P2", "P2(1)", "C2", "Pm", "Pc", "Cm", "Cc", "P2/m", 
      "P2(1)/m", "C2/m", "P2/c", "P2(1)/c", "C2/c", "P222", "P222(1)", 
      "P2(1)2(1)2", "P2(1)2(1)2(1)", "C222(1)", "C222", "F222", "I222", 
      "I2(1)2(1)2(1)", "Pmm2", "Pmc2(1)", "Pcc2", "Pma2", "Pca2(1)", "Pnc2", 
      "Pmn2(1)", "Pba2", "Pna2(1)", "Pnn2", "Cmm2", "Cmc2(1)", "Ccc2", "Amm2", 
      "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2", "Ima2", "Pmmm", 
      "Pnnn", "Pccm", "Pban", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", 
      "Pbcm", "Pnnm", "Pmmn", "Pbcn", "Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm", 
      "Cccm", "Cmma", "Ccca", "Fmmm", "Fddd", "Immm", "Ibam", "Ibca", "Imma", 
      "P4", "P4(1)", "P4(2)", "P4(3)", "I4", "I4(1)", "P-4", "I-4", "P4/m", 
      "P4(2)/m", "P4/n", "P4(2)/n", "I4/m", "I4(1)/a", "P422", "P42(1)2", 
      "P4(1)22", "P4(1)2(1)2", "P4(2)22", "P4(2)2(1)2", "P4(3)22", "P4(3)2(1)2",
      "I422", "I4(1)22", "P4mm", "P4bm", "P4(2)cm", "P4(2)nm", "P4cc", "P4nc",
      "P4(2)mc", "P4(2)bc", "I4mm", "I4cm", "I4(1)md", "I4(1)cd", "P-42m", 
      "P-42c", "P-42(1)m", "P-42(1)c", "P-4m2", "P-4c2", "P-4b2", "P-4n2", 
      "I-4m2", "I-4c2", "I-42m", "I-42d", "P4/mmm", "P4/mcc", "P4/nbm",
      "P4/nnc", "P4/mbm", "P4/mnc", "P4/nmm", "P4/ncc", "P4(2)/mmc", 
      "P4(2)/mcm", "P4(2)/nbc", "P4(2)/nnm", "P4(2)/mbc", "P4(2)/mnm", 
      "P4(2)/nmc", "P4(2)/ncm", "I4/mmm", "I4/mcm", "I4(1)/amd", "I4(1)/acd",
      "P3", "P3(1)", "P3(2)", "R3", "P-3", "R-3", "P312", "P321", "P3(1)12",
      "P3(1)21", "P3(2)12", "P3(2)21", "R32", "P3m1", "P31m", "P3c1", "P31c",
      "R3m", "R3c", "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m", "R-3c", "P6",
      "P6(1)", "P6(5)", "P6(2)", "P6(4)", "P6(3)", "P-6", "P6/m", "P6(3)/m",
      "P622", "P6(1)22", "P6(5)22", "P6(2)22", "P6(4)22", "P6(3)22", "P6mm",
      "P6cc", "P6(3)cm", "P6(3)mc", "P-6m2", "P-6c2", "P-62m", "P-62c",
      "P6/mmm", "P6/mcc", "P6(3)/mcm", "P6(3)/mmc", "P23", "F23", "I23",
      "P2(1)3", "I2(1)3", "Pm-3", "Pn-3", "Fm-3", "Fd-3", "Im-3", "Pa-3",
      "Ia-3", "P432", "P4(2)32", "F432", "F4(1)32", "I432", "P4(3)32",
      "P4(1)32", "I4(1)32", "P-43m", "F4-3m", "I-43m", "P-43n", "F-43c",
      "I-43d", "Pm-3m", "Pn-3n", "Pm-3n", "Pn-3m", "Fm-3m", "Fm-3c", 
      "Fd-3m", "Fd-3c", "Im-3m", "Ia-3d"
    };

    if (name.length () == 0)
	  {
        if (_spaceGroup != NULL)
          return _spaceGroup->GetId();
        else
          name = _spaceGroupName;
      }
    static const int numStrings = sizeof( spacegroups ) / sizeof( spacegroups[0] );
    for ( int i = 0; i < numStrings; ++i ) {
      if (name == spacegroups[i] ) {
        return i+1;
      }
    }
    return 0; //presumably never reached
  }

  // Helper function -- transform fractional coordinates to ensure they lie in the unit cell
  vector3 transformedFractionalCoordinate(vector3 originalCoordinate)
  {
    // ensure the fractional coordinate is entirely within the unit cell
    vector3 returnValue(originalCoordinate);

    // So if we have -2.08, we take -2.08 - (-2) = -0.08 .... almost what we want
    returnValue.SetX(originalCoordinate.x() - int(originalCoordinate.x()) );
    returnValue.SetY(originalCoordinate.y() - int(originalCoordinate.y()) );
    returnValue.SetZ(originalCoordinate.z() - int(originalCoordinate.z()) );

    if (returnValue.x() < 0.0)
      returnValue.SetX(returnValue.x() + 1.0);
    if (returnValue.y() < 0.0)
      returnValue.SetY(returnValue.y() + 1.0);
    if (returnValue.z() < 0.0)
      returnValue.SetZ(returnValue.z() + 1.0);

    return returnValue;
  }

  void OBUnitCell::FillUnitCell(OBMol *mol)
  {
    const SpaceGroup *sg = GetSpaceGroup(); // the actual space group and transformations for this unit cell
    
    // For each atom, we loop through: convert the coords back to inverse space, apply the transformations and create new atoms
    vector3 uniqueV, newV, updatedCoordinate;
    list<vector3> transformedVectors; // list of symmetry-defined copies of the atom
    list<vector3>::iterator transformIterator, duplicateIterator;
    OBAtom *newAtom;
    list<OBAtom*> atoms; // keep the current list of unique atoms -- don't double-create
    list<vector3> coordinates; // all coordinates to prevent duplicates
    bool foundDuplicate;
    FOR_ATOMS_OF_MOL(atom, *mol)
      atoms.push_back(&(*atom));

    list<OBAtom*>::iterator i;
    for (i = atoms.begin(); i != atoms.end(); ++i) {
      uniqueV = (*i)->GetVector();
      uniqueV *= GetFractionalMatrix();
      uniqueV = transformedFractionalCoordinate(uniqueV);
      coordinates.push_back(uniqueV);
  
      transformedVectors = sg->Transform(uniqueV);
      for (transformIterator = transformedVectors.begin();
           transformIterator != transformedVectors.end(); ++transformIterator) {
        // coordinates are in reciprocal space -- check if it's in the unit cell
        // if not, transform it in place
        updatedCoordinate = transformedFractionalCoordinate(*transformIterator);
        foundDuplicate = false;

        // Check if the transformed coordinate is a duplicate of an atom
        for (duplicateIterator = coordinates.begin();
             duplicateIterator != coordinates.end(); ++duplicateIterator) {
          if (duplicateIterator->distSq(updatedCoordinate) < 1.0e-4) {
            foundDuplicate = true;
            break;
          }
        }
        if (foundDuplicate)
          continue;
        
        coordinates.push_back(updatedCoordinate); // make sure to check the new atom for dupes
        newAtom = mol->NewAtom();
        newAtom->Duplicate(*i);
        newAtom->SetVector(GetOrthoMatrix() * updatedCoordinate);
      } // end loop of transformed atoms
      (*i)->SetVector(GetOrthoMatrix() * uniqueV); // move the atom back into the unit cell
    } // end loop of atoms

    SetSpaceGroup(1); // We've now applied the symmetry, so we should act like a P1 unit cell
  }
  
  double OBUnitCell::GetCellVolume()
  {
    double result = 0.0;
    
    switch ( GetLatticeType() )
      {
      case Triclinic:
        result = _a * _b * _c 
          * sqrt(1
                 - SQUARE(cos( _alpha * DEG_TO_RAD ))
                 - SQUARE(cos( _beta * DEG_TO_RAD ))
                 - SQUARE(cos( _gamma * DEG_TO_RAD ))
                 + 2 * cos( _alpha * DEG_TO_RAD ) * cos( _beta * DEG_TO_RAD ) * cos( _gamma * DEG_TO_RAD )
                 );
        break;
      case Monoclinic:
        result = _a * _b * _c * sin( _beta * DEG_TO_RAD );
        break;
      case Orthorhombic:
        result = _a * _b * _c;
        break;
      case Tetragonal:
        result = _a * _a * _c;
        break;
      case Rhombohedral:
        result = _a * _a * _a
          * sqrt(1
                 - SQUARE(cos( _alpha * DEG_TO_RAD ))
                 - SQUARE(cos( _beta * DEG_TO_RAD ))
                 - SQUARE(cos( _gamma * DEG_TO_RAD ))
                 + 2 * cos( _alpha * DEG_TO_RAD ) * cos( _beta * DEG_TO_RAD ) * cos( _gamma * DEG_TO_RAD )
                 );
        break;
      case Hexagonal:
        result = pow( 3.0, 0.333333333 ) * _a * _a * _c / 2;
        break;
      case Cubic:
        result = _a * _a * _a;
        break;
      default:
        result = 0.0;
      }
    
    return result;
  }
  
  //
  // member functions for OBSymmetryData class
  //
  OBSymmetryData::OBSymmetryData(): 
    OBGenericData("Symmetry", OBGenericDataType::SymmetryData)
  { }

  OBSymmetryData::OBSymmetryData(const OBSymmetryData &src) :
    OBGenericData(src._attr, src._type, src._source), 
    _pointGroup(src._pointGroup), _spaceGroup(src._spaceGroup)
  {  }

  OBSymmetryData & OBSymmetryData::operator=(const OBSymmetryData &src)
  {
    if(this == &src)
      return(*this);

    _pointGroup = src._pointGroup;
    _spaceGroup = src._spaceGroup;
    _source = src._source;

    return(*this);
  }

  OBConformerData::OBConformerData() :
    OBGenericData("Conformers", OBGenericDataType::ConformerData)
  {  }

  OBConformerData::OBConformerData(const OBConformerData &src) :
    OBGenericData("Conformers", OBGenericDataType::ConformerData),
    _vDimension(src._vDimension),
    _vEnergies(src._vEnergies), _vForces(src._vForces),
    _vVelocity(src._vVelocity), _vDisplace(src._vDisplace),
    _vData(src._vData)
  {  }

  OBConformerData & OBConformerData::operator=(const OBConformerData &src)
  {
    if(this == &src)
      return(*this);
    
    _source = src._source;

    _vDimension = src._vDimension;
    _vEnergies = src._vEnergies;
    _vForces = src._vForces;
    _vVelocity = src._vVelocity;
    _vDisplace = src._vDisplace;
    _vData = src._vData;

    return(*this);
  }

  //
  //member functions for OBRingData class
  //

  OBRingData::OBRingData() :
    OBGenericData("RingData", OBGenericDataType::RingData)
  {
    _vr.clear();
  }

  /*!
  **\brief OBRingData copy constructor
  **\param src reference to original OBRingData object (rhs)
  */
  OBRingData::OBRingData(const OBRingData &src)
    :	OBGenericData(src),	//chain to base class
      _vr(src._vr)				//chain to member classes
  {
    //no other memeber data
    //memory management

    vector<OBRing*>::iterator ring;

    for(ring = _vr.begin();ring != _vr.end();++ring)
      {
        OBRing *newring = new OBRing;
        (*newring) = (**ring);	//copy data to new object
        (*ring)    = newring;	//repoint new pointer to new copy of data
      }
  }

  OBRingData::~OBRingData()
  {
    vector<OBRing*>::iterator ring;
    for (ring = _vr.begin();ring != _vr.end();++ring)
      {
        delete *ring;
      }
  }

  /*!
  **\brief OBRingData assignment operator
  **\param src reference to original OBRingData object (rhs)
  **\return reference to changed OBRingData object (lhs)
  */
  OBRingData& OBRingData::operator =(const OBRingData &src)
  {
    //on identity, return
    if(this == &src)
      return(*this);

    //chain to base class
    OBGenericData::operator =(src);

    //member data

    //memory management
    vector<OBRing*>::iterator ring;
    for(ring = _vr.begin();ring != _vr.end();++ring)
      {
        delete &*ring;	//deallocate old rings to prevent memory leak
      }

    _vr.clear();
    _vr = src._vr;	//copy vector properties

    for(ring = _vr.begin();ring != _vr.end();++ring)
      {
        if(*ring == 0)
          continue;

        //allocate and copy ring data
        OBRing *newring = new OBRing;
        (*newring) = (**ring);
        (*ring) = newring;	//redirect pointer
      }
    return(*this);
  }

  OBRing *OBRingData::BeginRing(std::vector<OBRing*>::iterator &i)
  {
    i = _vr.begin();
    return((i == _vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  OBRing *OBRingData::NextRing(std::vector<OBRing*>::iterator &i)
  {
    ++i;
    return((i == _vr.end()) ? (OBRing*)NULL : (OBRing*)*i);
  }

  //
  //member functions for OBAngle class - stores all angles
  //

  /*!
  **\brief Angle default constructor
  */
  OBAngle::OBAngle():
    _vertex(NULL), _termini(NULL, NULL), _radians(0.0)
  {  }

  /*!
  **\brief Angle constructor
  */
  OBAngle::OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b):
    _vertex(vertex), _termini(a, b)
  {
    SortByIndex();
  }

  /*!
  **\brief OBAngle copy constructor
  */
  OBAngle::OBAngle(const OBAngle &src):
    _vertex(src._vertex), _termini(src._termini), _radians(src._radians)
  {  }

  /*!
  **\brief OBAngle assignment operator
  */
  OBAngle& OBAngle::operator = (const OBAngle &src)
  {
    if (this == &src)
      return(*this);

    _vertex         = src._vertex;
    _termini.first  = src._termini.first;
    _termini.second = src._termini.second;
    _radians        = src._radians;

    return(*this);
  }

  /*!
  **\brief Return OBAngle to its original state
  */
  void OBAngle::Clear()
  {
    _vertex         = 0;
    _termini.first  = 0;
    _termini.second = 0;
    _radians        = 0.0;
    return;
  }

  /*!
  **\brief Sets the 3 atoms in the angle
  ** Parameters are pointers to each OBAtom
  */
  void OBAngle::SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b)
  {
    _vertex         = vertex;
    _termini.first  = a;
    _termini.second = b;
    SortByIndex();
    return;
  }

  /*!
  **\brief Sets the 3 atoms in the angle
  **\param atoms a triple of OBAtom pointers, the first must be the vertex
  */
  void OBAngle::SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms)
  {
    _vertex         = atoms.first;
    _termini.first  = atoms.second;
    _termini.second = atoms.third;
    SortByIndex();
    return;
  }

  /*!
  **\brief Retrieves the 3 atom pointer for the angle (vertex first)
  **\return triple of OBAtom pointers 
  */
  triple<OBAtom*,OBAtom*,OBAtom*> OBAngle::GetAtoms()
  {
    triple<OBAtom*,OBAtom*,OBAtom*> atoms;
    atoms.first  = _vertex;
    atoms.second = _termini.first;
    atoms.third  = _termini.second;
    return(atoms);
  }

  /*!
  **\brief sorts atoms in angle by order of indices
  */
  void OBAngle::SortByIndex()
  {
    OBAtom *tmp;

    if(_termini.first->GetIdx() > _termini.second->GetIdx())
      {
        tmp             = _termini.first;
        _termini.first  = _termini.second;
        _termini.second = tmp;
      }
  }

  /*!
  **\brief OBAngle equality operator, is same angle, NOT same value
  **\return boolean equality
  */
  bool OBAngle::operator ==(const OBAngle &other)
  {
    return ((_vertex         == other._vertex)        &&
            (_termini.first  == other._termini.first) &&
            (_termini.second == other._termini.second));
  }

  //
  //member functions for OBAngleData class - stores OBAngle set
  //

  /*!
  **\brief OBAngleData constructor
  */
  OBAngleData::OBAngleData()
    :	OBGenericData("AngleData", OBGenericDataType::AngleData)
  {  }

  /*!
  **\brief OBAngleData copy constructor
  */
  OBAngleData::OBAngleData(const OBAngleData &src)
    :	OBGenericData(src), _angles(src._angles)
  {  }

  /*!
  **\brief OBAngleData assignment operator
  */
  OBAngleData& OBAngleData::operator =(const OBAngleData &src)
  {
    if (this == &src)
      return(*this);

    _source = src._source;
    _angles = src._angles;

    return(*this);
  }

  /*!
  **\brief sets OBAngleData to its original state
  */
  void OBAngleData::Clear()
  {
    _angles.clear();
    return;
  }

  /*!
  **\brief Adds a new angle to OBAngleData
  */
  void OBAngleData::SetData(OBAngle &angle)
  {
    _angles.push_back(angle);
    return;
  }

  /*!
  **\brief Fills an array with the indices of the atoms in the angle (vertex first)
  **\param angles pointer to the pointer to an array of angles atom indices
  **\return True if successful
  */
  bool OBAngleData::FillAngleArray(std::vector<std::vector<unsigned int> > &angles)
  {
    if(_angles.empty())
      return(false);

    vector<OBAngle>::iterator angle;
    
    angles.clear();
    angles.resize(_angles.size());

    unsigned int ct = 0;

    for( angle=_angles.begin(); angle!=_angles.end(); angle++,ct++)
      {
        angles[ct].resize(3);
        angles[ct][0] = angle->_vertex->GetIdx() - 1;
        angles[ct][1] = angle->_termini.first->GetIdx() - 1;
        angles[ct][2] = angle->_termini.second->GetIdx() - 1;
      }
 
    return(true);
  }
  
  /*!
  **\brief Fills an array with the indices of the atoms in the angle (vertex first)
  **\param angles pointer to the pointer to an array of angles atom indices
  **\param size the current number of rows in the array
  **\return int The number of angles
  */
  unsigned int OBAngleData::FillAngleArray(int **angles, unsigned int &size)
  {
    if(_angles.size() > size)
      {
        delete [] *angles;
        *angles = new int[_angles.size()*3];
        size    = (unsigned int)_angles.size();
      }

    vector<OBAngle>::iterator angle;
    int angleIdx = 0;
    for( angle=_angles.begin(); angle!=_angles.end(); ++angle)
      {
        *angles[angleIdx++] = angle->_vertex->GetIdx();
        *angles[angleIdx++] = angle->_termini.first->GetIdx();
        *angles[angleIdx++] = angle->_termini.second->GetIdx();
      }
    return (unsigned int)_angles.size();
  }
  
  //
  //member functions for OBAngleData class - stores OBAngle set
  //

  /*!
  **\brief OBTorsion constructor
  */
  OBTorsion::OBTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d)
  {
    triple<OBAtom*,OBAtom*,double> ad(a,d,0.0);
    _ads.push_back(ad);

    _bc.first  = b;
    _bc.second = c;
  }

  /*!
  **\brief OBTorsion copy constructor
  */
  OBTorsion::OBTorsion(const OBTorsion &src)
    :	_bc(src._bc), _ads(src._ads)
  {}

  /*!
  **\brief Returns all the 4 atom sets in OBTorsion
  */
  vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > OBTorsion::GetTorsions()
  {
    quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> abcd;

    abcd.second = _bc.first;
    abcd.third  = _bc.second;

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > torsions;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;

    for(ad = _ads.begin();ad != _ads.end();++ad)
      {
        abcd.first = ad->first;
        abcd.fourth = ad->second;
        torsions.push_back(abcd);
      }

    return(torsions);
  }

  /*!
  **\brief OBTorsion assignment operator
  */
  OBTorsion& OBTorsion::operator =(const OBTorsion &src)
  {
    if (this == &src)
      return(*this);

    _bc  = src._bc;
    _ads = src._ads;

    return(*this);
  }

  /*!
  **\brief Returns the OBTorsion to its original state
  */
  void OBTorsion::Clear()
  {
    _bc.first  = 0;
    _bc.second = 0;
    _ads.erase(_ads.begin(),_ads.end());
  }

  /*!
  **\brief Sets the angle of a torsion in OBTorsion
  **\param radians the value to assign to the torsion
  **\param index the index into the torsion of the OBTorsion
  **\return boolean success
  */
  bool OBTorsion::SetAngle(double radians,unsigned int index)
  {
    if(index >= _ads.size())
      return(false);

    _ads[index].third = radians;

    return(true);
  }

  /*!
  **\brief Obtains the angle of a torsion in OBTorsion
  **\param radians the value of the angle is set here
  **\param index the index into the torsion of the OBTorsion
  **\return boolean success
  */
  bool OBTorsion::GetAngle(double &radians, unsigned int index)
  {
    if(index >= _ads.size())
      return false;
    radians = _ads[index].third;
    return true;
  }

  unsigned int OBTorsion::GetBondIdx()
  {
    return(_bc.first->GetBond(_bc.second)->GetIdx());
  }

  /*!
  **\brief determines if torsion has only protons on either the a or d end
  **\return boolean 
  */
  bool OBTorsion::IsProtonRotor()
  {
    bool Aprotor = true;
    bool Dprotor = true;
    vector<triple<OBAtom*,OBAtom*,double> >::iterator ad;
    for(ad = _ads.begin();ad != _ads.end() && (Aprotor || Dprotor);++ad)
      {
        if(!ad->first->IsHydrogen())
          Aprotor = false;
        if(!ad->second->IsHydrogen())
          Dprotor = false;
      }
    return (Aprotor || Dprotor);
  }

  /*!
  **\brief adds a new torsion to the OBTorsion object
  */
  bool OBTorsion::AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d)
  {
    if(!Empty() && (b != _bc.first || c != _bc.second))
      return(false);

    if(Empty())
      {
        _bc.first  = b;
        _bc.second = c;
      }

    triple<OBAtom*,OBAtom*,double> ad(a,d,0.0);
    _ads.push_back(ad);

    return(true);
  }

  /*!
  **\brief adds a new torsion to the OBTorsion object
  */
  bool OBTorsion::AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms)
  {
    if(!Empty() && (atoms.second != _bc.first || atoms.third != _bc.second))
      return(false);

    if(Empty())
      {
        _bc.first  = atoms.second;
        _bc.second = atoms.third;
      }

    triple<OBAtom*,OBAtom*,double> ad(atoms.first,atoms.fourth,0.0);
    _ads.push_back(ad);

    return(true);
  }

  //\!brief OBTorsionData ctor
  OBTorsionData::OBTorsionData()
    : OBGenericData("TorsionData", OBGenericDataType::TorsionData)
  {  }

  //
  //member functions for OBTorsionData class - stores OBTorsion set
  //
  OBTorsionData::OBTorsionData(const OBTorsionData &src)
    :	OBGenericData(src), _torsions(src._torsions)
  {  }

  OBTorsionData& OBTorsionData::operator =(const OBTorsionData &src)
  {
    if (this == &src)
      return(*this);

    OBGenericData::operator =(src);

    _source = src._source;
    _torsions = src._torsions;

    return(*this);
  }

  void OBTorsionData::Clear()
  {
    _torsions.clear();
  }

  void OBTorsionData::SetData(OBTorsion &torsion)
  {
    _torsions.push_back(torsion);
  }

  /*!
  **\brief Fills a vector with the indices of the atoms in torsions (ordered abcd)
  **\param torsions reference to the vector of abcd atom sets
  **\return boolean success
  */
  bool OBTorsionData::FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions)
  {
    if(_torsions.empty())
      return(false);

    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > tmpquads,quads;
    vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> >::iterator thisQuad;
    vector<OBTorsion>::iterator torsion;

    //generate set of all 4 atom abcd's from torsion structure
    for (torsion = _torsions.begin();torsion != _torsions.end();++torsion)
      {
        tmpquads = torsion->GetTorsions();
        for(thisQuad = tmpquads.begin();thisQuad != tmpquads.end();++thisQuad)
          quads.push_back(*thisQuad);
      }

    //fill array of torsion atoms

    torsions.clear();
    torsions.resize(quads.size());

    unsigned int ct = 0;

    for (thisQuad = quads.begin();thisQuad != quads.end();++thisQuad,++ct)
      {
        torsions[ct].resize(4);
        torsions[ct][0] = thisQuad->first->GetIdx()-1;
        torsions[ct][1] = thisQuad->second->GetIdx()-1;
        torsions[ct][2] = thisQuad->third->GetIdx()-1;
        torsions[ct][3] = thisQuad->fourth->GetIdx()-1;
      }

    return(true);
  }

  //
  // Member functions for OBChiralDarta
  //
  bool OBChiralData::SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t)
  {
    if (atom4refs.size() != 4)
      {
        obErrorLog.ThrowError(__FUNCTION__, "Incorrect number of atoms atom4refs, should be 4", obDebug);
        return(false);
      }
    switch(t){
    case input: _atom4refs = atom4refs;break;
    case output:_atom4refo = atom4refs;break;
    case calcvolume:_atom4refc = atom4refs;break;
    default: 
      obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
      return(false);
    }
    return (true);
  }

  int OBChiralData::AddAtomRef(unsigned int atomref, atomreftype t)
  {
    switch(t){
    case input: _atom4refs.push_back(atomref);break;
    case output: _atom4refo.push_back(atomref);break;
    case calcvolume:_atom4refc.push_back(atomref);break;
    default:
      obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
      return(false);
    }
              
    return (_atom4refs.size());
  }

  unsigned int OBChiralData::GetAtomRef(int a, atomreftype t)
  {
    switch(t){
    case input: return(_atom4refs[a]);break;
    case output: return(_atom4refo[a]);break;
    case calcvolume: return(_atom4refc[a]);break;
    default:
      obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
      return(false);
    }  
  }
  std::vector<unsigned int> OBChiralData::GetAtom4Refs(atomreftype t) const
  {
    switch (t){
    case output:
      return(_atom4refo);
      break;
    case input:
      return(_atom4refs);
      break;
    case calcvolume:
      return(_atom4refc);
      break;
    default:
      obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
      return(_atom4refo);
    }
  }

  unsigned int OBChiralData::GetSize(atomreftype t) const
  {
    switch (t)
      {
      case output:
        return(unsigned int)_atom4refo.size();
        break;
      case input:
        return(unsigned int)_atom4refs.size();
        break;
      case calcvolume:
        return(unsigned int)_atom4refc.size();
      default:
        obErrorLog.ThrowError(__FUNCTION__, "AtomRefType called is invalid", obDebug);
        return(0);
      }
  }

  // Chiral data is a perceived data type. We might read in some chiral info
  // but this class derives and converts from whatever is read
  OBChiralData::OBChiralData()
    : OBGenericData("ChiralData", OBGenericDataType::ChiralData, perceived)
  {  }

  OBChiralData::OBChiralData(const OBChiralData &src)
    : OBGenericData(src)
  {
    _atom4refs = src._atom4refs;
    _atom4refo = src._atom4refo;
    _atom4refc = src._atom4refc;
    parity     = src.parity;
  }

  OBChiralData & OBChiralData::operator=(const OBChiralData &src)
  {
    if(this == &src)
      return(*this);

    _source = src._source;

    _atom4refs = src._atom4refs;
    _atom4refo = src._atom4refo;
    _atom4refc = src._atom4refc;
    parity=src.parity;
    return(*this);
  }

  void OBChiralData::Clear()
  {
    _atom4refs.clear();
    parity=0;
    _atom4refo.clear();
    _atom4refc.clear();
  }

//
//member functions for OBDOSData class
//

/*!
**\brief Assign the data
**\param fermi The Fermi energy in eV
**\param vEnergies Energy levels in eV
**\param vDensities Density of states in (number of states) / (unit cell)
**\param vIntegration Integrated DOS
*/
void OBDOSData::SetData(double fermi,
                        const std::vector<double> & vEnergies,
                        const std::vector<double> & vDensities,
                        const std::vector<double> & vIntegration)
{
  this->_fermi = fermi;
  this->_vEnergies = vEnergies;
  this->_vIntegration = vIntegration;
  this->_vDensities = vDensities;
}

//
//member functions for OBElectronicTransitionData class
//

/*!
**\brief Assign the basic excitation data
**\param vWavelengths Wavelengths in nm
**\param vForces Oscillator strengths
*/
void OBElectronicTransitionData::SetData(const std::vector<double> & vWavelengths,
                                  const std::vector<double> & vForces)
{
  this->_vWavelengths = vWavelengths;
  this->_vForces = vForces;
}

/*!
**\brief Assign the electronic dipole strengths
**\param vEDipole Electronic dipole moment strength
*/
void OBElectronicTransitionData::SetEDipole(const std::vector<double> & vEDipole)
{
  this->_vEDipole = vEDipole;
}

/*!
**\brief Assign the rotatory strengths (velocity)
**\param vRotatorStrengthsVelocity Vector containing the rotatory strengths
*/
void OBElectronicTransitionData::SetRotatoryStrengthsVelocity(const std::vector<double> & vRotatoryStrengthsVelocity)
{
  this->_vRotatoryStrengthsVelocity = vRotatoryStrengthsVelocity;
}

/*!
**\brief Assign the rotatory strengths (length)
**\param vRotatorStrengthsLength Vector containing the rotatory strengths
*/
void OBElectronicTransitionData::SetRotatoryStrengthsLength(const std::vector<double> & vRotatoryStrengthsLength)
{
  this->_vRotatoryStrengthsLength = vRotatoryStrengthsLength;
}

//
//member functions for OBVibrationData class
//

/*!
**\brief Assign the data
**\param vLx Normal modes in 1/sqrt(a.u.)
**\param vFrequencies Harmonic frequencies in inverse centimeters
**\param vIntensities Infrared absorption intensities in KM/Mole
*/
void OBVibrationData::SetData(const std::vector< std::vector< vector3 > > & vLx,
                              const std::vector<double> & vFrequencies,
                              const std::vector<double> & vIntensities)
{
  this->_vLx = vLx;
  this->_vFrequencies = vFrequencies;
  this->_vIntensities = vIntensities;
}

/*!
**\brief Assign the data
**\param vLx Normal modes in 1/sqrt(a.u.)
**\param vFrequencies Harmonic frequencies in inverse centimeters
**\param vIntensities Infrared absorption intensities in KM/Mole
**\param vRamanActivities Raman activities
*/
void OBVibrationData::SetData(const std::vector< std::vector< vector3 > > & vLx,
                              const std::vector<double> & vFrequencies,
                              const std::vector<double> & vIntensities,
                              const std::vector<double> & vRamanActivities)
{
  this->_vLx = vLx;
  this->_vFrequencies = vFrequencies;
  this->_vIntensities = vIntensities;
  this->_vRamanActivities = vRamanActivities;
}


/*!
**\brief Get the number of frequencies
*/
unsigned int OBVibrationData::GetNumberOfFrequencies() const
{
  return !this->_vFrequencies.empty() ? this->_vFrequencies.size() : 0;
}

} //end namespace OpenBabel

//! \file generic.cpp
//! \brief Handle OBGenericData classes. Custom data for atoms, bonds, etc.
