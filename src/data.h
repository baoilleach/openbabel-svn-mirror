/**********************************************************************
data.h - Global data and resource file parsers.
 
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

#ifndef OB_DATA_H
#define OB_DATA_H

#include "babelconfig.h"

#include <stdio.h>

#include <fstream>
#include <vector>
#include <string>

namespace OpenBabel
{

  class OBAtom;
  class OBMol;
  class OBBitVec;

  //! \brief Base data table class, handles reading data files
  //!
  //! Base data table class--reads ASCII data files in various formats
  //! -# Checks for the environment variable _envvar (defaults to "BABEL_DATADIR")
  //!     - Tries the _subdir directory if defined (def. "data") and then the main directory
  //! -# Checks for the directory _dir (def. determined by the build environment)
  //!     - Tries the subdirectory corresponding to this version, then the main directory
  //! -# Reverts to the compiled-in default data
  class OBAPI OBGlobalDataBase
    {
    protected:
      bool         _init;		//!< Whether the data been read already
      const char  *_dataptr;//!< Default data table if file is unreadable
      std::string  _filename;//!< File to search for
      std::string  _dir;		//!< Data directory for file if _envvar fails
      std::string  _subdir;	//!< Subdirectory (if using environment variable)
      std::string  _envvar;	//!< Environment variable to check first

    public:
      //! Constructor
      OBGlobalDataBase(): _init(false), _dataptr(NULL) { }
      //! Destructor
      virtual ~OBGlobalDataBase()                  {}
      //! Read in the data file, falling back as needed
      void  Init();
      //! \return the size of the database (for error checking)
      virtual unsigned int GetSize()                 { return 0;}
      //! Set the directory before calling Init()
      void  SetReadDirectory(char *dir)            { _dir = dir;    }
      //! Set the environment variable to use before calling Init()
      void  SetEnvironmentVariable(char *var)      { _envvar = var; }
      //! Specified by particular table classes (parses an individual data line)
      virtual void ParseLine(const char*)          {}
    };

  //! \brief Individual element data type
  //!
  //! Stores a variety of data about an individual element
  class OBAPI OBElement
    {
      int _num;
      char _symbol[3];
      std::string _name;
      double _Rcov,_Rvdw,_mass,_elNeg,_ionize,_elAffinity;
      double _red, _green, _blue;
      int _maxbonds;
    public:
      OBElement()    {}
      OBElement(int num, const char *sym, double rcov, double rvdw,
                int maxbo, double mass, double elNeg, double ionize,
                double elAffin, double red, double green, double blue,
                std::string name) :
        _num(num), _name(name), _Rcov(rcov), _Rvdw(rvdw), _mass(mass), 
        _elNeg(elNeg), _ionize(ionize), _elAffinity(elAffin), 
        _red(red), _green(green), _blue(blue),
        _maxbonds(maxbo)
        {
          strncpy(_symbol, sym, 3);
        }

      //! \return the atomic number of this element
      int GetAtomicNum()         {       return(_num);    }
      //! \return the atomic symbol for this element
      char *GetSymbol()          {       return(_symbol); }
      //! \return the covalent radius of this element
      double GetCovalentRad()    {       return(_Rcov);   }
      //! \return the van der Waals radius of this element
      double GetVdwRad()         {       return(_Rvdw);   }
      //! \return the standard atomic mass for this element (in amu)
      double GetMass()           {       return(_mass);   }
      //! \return the maximum expected number of bonds to this element
      int GetMaxBonds()          {       return(_maxbonds);}
      //! \return the Pauling electronegativity for this element
      double GetElectroNeg()     {       return(_elNeg);  }
      //! \return the ionization potential (in eV) of this element
      double GetIonization()     {       return(_ionize);  }
      //! \return the electron affinity (in eV) of this element
      double GetElectronAffinity(){      return(_elAffinity);  }
      //! \return the name of this element (in English)
      std::string GetName()      {       return(_name);    }
      //! \return the red component of this element's default visualization color
      double GetRed()            {       return(_red);     }
      //! \return the green component of this element's default color
      double GetGreen()          {       return(_green);   }
      //! \return the blue component of this element's default color
      double GetBlue()           {       return(_blue);    }
    };

  // class introduction in data.cpp
  class OBAPI OBElementTable : public OBGlobalDataBase
    {
      std::vector<OBElement*> _element;

    public:

      OBElementTable(void);
      ~OBElementTable();

      void  ParseLine(const char*);

      //! \return the number of elements in the periodic table
      unsigned int		GetNumberOfElements();
      //! \return the number of elements in the periodic table
      unsigned int    GetSize() { return GetNumberOfElements(); }

      //! \deprecated Does not properly handle 'D' or 'T' hydrogen isotopes
      int   GetAtomicNum(const char *);
      //! \return the atomic number matching the element symbol passed
      //! or 0 if not defined. For 'D' or 'T' hydrogen isotopes, will return
      //! a value in the second argument
      int   GetAtomicNum(const char *, int &iso);
      //! \return the element symbol matching the atomic number passed
      char *GetSymbol(int);
      //! \return the van der Waals radius for this atomic number
      double GetVdwRad(int);
      //! \return the covalent radius for this atomic number
      double GetCovalentRad(int);
      //! \return the average atomic mass for this element.
      //! For exact isotope masses, use OpenBabel::OBIsotopeTable
      double GetMass(int);
      //! \return a "corrected" bonding radius based on the hybridization.
      //! Scales the covalent radius by 0.95 for sp2 and 0.90 for sp hybrids
      double CorrectedBondRad(int,int = 3); // atomic #, hybridization
      //! \return a "corrected" vdW radius based on the hybridization.
      //! Scales the van der Waals radius by 0.95 for sp2 and 0.90 for sp hybrids
      double CorrectedVdwRad(int,int = 3); // atomic #, hybridization
      //! \return the maximum expected number of bonds to this element
      int	GetMaxBonds(int);
      //! \return the Pauling electronegativity for this element
      double GetElectroNeg(int);
      //! \return the ionization potential (in eV) for this element
      double GetIonization(int);
      //! \return the electron affinity (in eV) for this element
      double GetElectronAffinity(int);
      //! \return a vector with red, green, blue color values for this element
      std::vector<double> GetRGB(int);
      //! \return the name of this element
      std::string GetName(int);
    };

  // class introduction in data.cpp
  class OBAPI OBIsotopeTable : public OBGlobalDataBase
    {
      std::vector<std::vector<std::pair <unsigned int, double> > > _isotopes;

    public:

      OBIsotopeTable(void);
      ~OBIsotopeTable()    {}

      //! \return the number of elements in the isotope table
      unsigned int GetSize() { return _isotopes.size(); }

      void	ParseLine(const char*);
      //! \return the exact masss of the isotope
      //!   (or by default (i.e. "isotope 0") the most abundant isotope)
      double	GetExactMass(const unsigned int atomicNum,
                           const unsigned int isotope = 0);
    };

  // class introduction in data.cpp
  class OBAPI OBTypeTable : public OBGlobalDataBase
    {
      int             _linecount;
      unsigned int    _ncols,_nrows;
      int             _from,_to;
      std::vector<std::string> _colnames;
      std::vector<std::vector<std::string> > _table;

    public:

      OBTypeTable(void);
      ~OBTypeTable() {}

      void ParseLine(const char*);

      //! \return the number of atom types in the translation table
      unsigned int GetSize() { return _table.size(); }

      //! Set the initial atom type to be translated
      bool SetFromType(const char*);
      //! Set the destination atom type for translation
      bool SetToType(const char*);
      //! Translate atom types
      bool Translate(char *to, const char *from); // to, from
      //! Translate atom types
      //! \return whether the translation was successful
      bool Translate(std::string &to, const std::string &from); // to, from
      //! Translate atom types
      //! \return the translated atom type, or an empty string if not possible
      std::string Translate(const std::string &from);

      //! \return the initial atom type to be translated
      std::string GetFromType();
      //! \return the destination atom type for translation
      std::string GetToType();
    };

  //! \brief Table of common biomolecule residues (for PDB or other files).
  //!   Can assign atom types and bond orders for arbitrary residues
  class OBAPI OBResidueData : public OBGlobalDataBase
    {
      int                                               _resnum;
      std::vector<std::string>                          _resname;
      std::vector<std::vector<std::string> >            _resatoms;
      std::vector<std::vector<std::pair<std::string,int> > > _resbonds;

      //variables used only temporarily for parsing resdata.txt
      std::vector<std::string>                          _vatmtmp;
      std::vector<std::pair<std::string,int> >          _vtmp;
    public:

      OBResidueData();
      void ParseLine(const char*);

      //! \return the number of residues in the table
      unsigned int GetSize() { return _resname.size(); }

      //! Sets the table to access the residue information for a specified
      //!  residue name
      //! \return whether this residue name is in the table
      bool SetResName(const std::string &);
      //! \return the bond order for the bond specified in the current residue
      //! \deprecated Easier to use the two-argument form
      int  LookupBO(const std::string &);
      //! \return the bond order for the bond specified between the two specified
      //! atom labels
      int  LookupBO(const std::string &, const std::string&);
      //! Look up the atom type and hybridization for the atom label specified
      //! in the first argument for the current residue
      //! \return whether the atom label specified is found in the current residue
      bool LookupType(const std::string &,std::string&,int&);
      //! Assign bond orders, atom types and residues for the supplied OBMol
      //! based on the residue information assigned to atoms
      //! \deprecated second OBBitVec argument is ignored
      bool AssignBonds(OBMol &,OBBitVec &);
    };

  //! Opens a datafile in a directory where OpenBabel expects to find it.
  OBAPI std::string OpenDatafile(std::ifstream& fs, const std::string& filename);

  // Used by other code for reading files
#ifdef WIN32
#define FILE_SEP_CHAR "\\"
#else
#define FILE_SEP_CHAR "/"
#endif

} // end namespace OpenBabel

#endif //DATA_H

//! \file data.h
//! \brief Global data and resource file parsers.
