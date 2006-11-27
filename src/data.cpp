/**********************************************************************
data.cpp - Global data and resource file parsers.
 
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

#ifdef WIN32
#pragma warning (disable : 4786)
#endif

#include "babelconfig.h"
#include "data.h"
#include "mol.h"

// data headers
#include "element.h"
#include "types.h"
#include "isotope.h"
#include "resdata.h"


#if !HAVE_STRNCASECMP
extern "C" int strncasecmp(const char *s1, const char *s2, size_t n);
#endif

using namespace std;

namespace OpenBabel
{

  OBElementTable   etab;
  OBTypeTable      ttab;
  OBIsotopeTable   isotab;
  OBResidueData    resdat;

  /** \class OBElementTable
      \brief Periodic Table of the Elements
 
      Translating element data is a common task given that many file
      formats give either element symbol or atomic number information, but
      not both. The OBElementTable class facilitates conversion between
      textual and numeric element information. An instance of the
      OBElementTable class (etab) is declared as external in data.cpp. Source
      files that include the header file mol.h automatically have an extern
      definition to etab. The following code sample demonstrates the use
      of the OBElementTable class:
      \code
      cout << "The symbol for element 6 is " << etab.GetSymbol(6) << endl;
      cout << "The atomic number for Sulfur is " << etab.GetAtomicNum(16) << endl;
      cout << "The van der Waal radius for Nitrogen is " << etab.GetVdwRad(7);
      \endcode
 
      Stored information in the OBElementTable includes elemental:
      - symbols
      - covalent radii
      - van der Waal radii
      - expected maximum bonding valence
      - molar mass (by IUPAC recommended atomic masses)
      - electronegativity
      - ionization potential
      - electron affinity
      - RGB colors for visualization programs
      - names (by IUPAC recommendation)
  */

  OBElementTable::OBElementTable()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "element.txt";
    _subdir = "data";
    _dataptr = ElementData;
  }

  OBElementTable::~OBElementTable()
  {
    vector<OBElement*>::iterator i;
    for (i = _element.begin();i != _element.end();++i)
      delete *i;
  }

  void OBElementTable::ParseLine(const char *buffer)
  {
    int num,maxbonds;
    char symbol[5];
    char name[256];
    double Rcov,Rvdw,mass, elNeg, ionize, elAffin;
    double red, green, blue;

    if (buffer[0] != '#') // skip comment line (at the top)
      {
        sscanf(buffer,"%d %5s %lf %*f %lf %d %lf %lf %lf %lf %lf %lf %lf %255s",
               &num,
               symbol,
               &Rcov,
               &Rvdw,
               &maxbonds,
               &mass,
               &elNeg,
               &ionize,
               &elAffin,
               &red,
               &green,
               &blue,
               name);

        OBElement *ele = new OBElement(num,symbol,Rcov,Rvdw,maxbonds,mass,elNeg,
                                       ionize, elAffin, red, green, blue, name);
        _element.push_back(ele);
      }
  }

  unsigned int OBElementTable::GetNumberOfElements()
  { 
    if (!_init)
      Init();
    
    return _element.size();
  }

  char *OBElementTable::GetSymbol(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return("\0");

    return(_element[atomicnum]->GetSymbol());
  }

  int OBElementTable::GetMaxBonds(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0);

    return(_element[atomicnum]->GetMaxBonds());
  }

  double OBElementTable::GetElectroNeg(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetElectroNeg());
  }

  double OBElementTable::GetIonization(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetIonization());
  }


  double OBElementTable::GetElectronAffinity(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetElectronAffinity());
  }
 
  vector<double> OBElementTable::GetRGB(int atomicnum)
  {
    if (!_init)
      Init();

    vector <double> colors;
    colors.reserve(3);
    
    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      {
        colors.push_back(0.0f);
        colors.push_back(0.0f);
        colors.push_back(0.0f);
        return(colors);
      }

    colors.push_back(_element[atomicnum]->GetRed());
    colors.push_back(_element[atomicnum]->GetGreen());
    colors.push_back(_element[atomicnum]->GetBlue());

    return (colors);
  }
 
  string OBElementTable::GetName(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return("Unknown");

    return(_element[atomicnum]->GetName());
  }

  double OBElementTable::GetVdwRad(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetVdwRad());
  }

  double OBElementTable::CorrectedBondRad(int atomicnum, int hyb)
  {
    double rad;
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(1.0);

    rad = _element[atomicnum]->GetCovalentRad();

    if (hyb == 2)
      rad *= 0.95;
    else if (hyb == 1)
      rad *= 0.90;

    return(rad);
  }

  double OBElementTable::CorrectedVdwRad(int atomicnum, int hyb)
  {
    double rad;
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(1.95);

    rad = _element[atomicnum]->GetVdwRad();

    if (hyb == 2)
      rad *= 0.95;
    else if (hyb == 1)
      rad *= 0.90;

    return(rad);
  }

  double OBElementTable::GetCovalentRad(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetCovalentRad());
  }

  double OBElementTable::GetMass(int atomicnum)
  {
    if (!_init)
      Init();

    if (atomicnum < 0 || atomicnum > static_cast<int>(_element.size()))
      return(0.0);

    return(_element[atomicnum]->GetMass());
  }

  int OBElementTable::GetAtomicNum(const char *sym)
  {
    int temp;
    return GetAtomicNum(sym, temp);
  }

  int OBElementTable::GetAtomicNum(const char *sym, int &iso)
  {
    if (!_init)
      Init();

    vector<OBElement*>::iterator i;
    for (i = _element.begin();i != _element.end();++i)
      if (!strncasecmp(sym,(*i)->GetSymbol(),2))
        return((*i)->GetAtomicNum());
    if (strcasecmp(sym, "D") == 0)
      {
        iso = 2;
        return(1);
      }
    else if (strcasecmp(sym, "T") == 0)
      {
        iso = 3;
        return(1);
      }
    else
      iso = 0;
    return(0);
  }

  /** \class OBIsotopeTable
      \brief Table of atomic isotope masses
 
  */

  OBIsotopeTable::OBIsotopeTable()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "isotope.txt";
    _subdir = "data";
    _dataptr = IsotopeData;
  }

  void OBIsotopeTable::ParseLine(const char *buffer)
  {
    unsigned int atomicNum;
    unsigned int i;
    vector<string> vs;

    pair <unsigned int, double> entry;
    vector <pair <unsigned int, double> > row;

    if (buffer[0] != '#') // skip comment line (at the top)
      {
        tokenize(vs,buffer);
        if (vs.size() > 3) // atomic number, 0, most abundant mass (...)
          {
            atomicNum = atoi(vs[0].c_str());
            for (i = 1; i < vs.size() - 1; i += 2) // make sure i+1 still exists
              {
                entry.first = atoi(vs[i].c_str()); // isotope
                entry.second = atof(vs[i + 1].c_str()); // exact mass
                row.push_back(entry);
              }
            _isotopes.push_back(row);
          }
        else
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse line in isotope table isotope.txt", obInfo);
      }
  }

  double        OBIsotopeTable::GetExactMass(const unsigned int ele,
                                             const unsigned int isotope)
  {
    if (!_init)
      Init();

    if (ele > _isotopes.size())
      return 0.0;

    unsigned int iso;
    for (iso = 0; iso < _isotopes[ele].size(); ++iso)
      if (isotope == _isotopes[ele][iso].first)
        return _isotopes[ele][iso].second;

    return 0.0;
  }

  /** \class OBTypeTable
      \brief Atom Type Translation Table
 
      Molecular file formats frequently store information about atoms in an
      atom type field. Some formats store only the element for each atom,
      while others include hybridization and local environments, such as the
      Sybyl mol2 atom type field. The OBTypeTable class acts as a translation
      table to convert atom types between a number of different molecular
      file formats. The constructor for OBTypeTable automatically reads the
      text file types.txt. Just as OBElementTable, an instance of
      OBTypeTable (ttab) is declared external in data.cpp and is referenced as
      extern OBTypeTable ttab in mol.h.  The following code demonstrates how
      to use the OBTypeTable class to translate the internal representation
      of atom types in an OBMol Internal to Sybyl Mol2 atom types.
 
      \code
      ttab.SetFromType("INT");
      ttab.SetToType("SYB");
      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      string src,dst;
      for (atom = mol.BeginAtom(i);atom;atom = mol.EndAtom(i))
      {
      src = atom->GetType();
      ttab.Translate(dst,src);
      cout << "atom number " << atom->GetIdx() << "has mol2 type " << dst << endl;
      }
      \endcode
 
      Current atom types include (defined in the top line of the data file types.txt):
      - INT (Open Babel internal codes)
      - ATN (atomic numbers)
      - HYB (hybridization)
      - MMD
      - MM2 (MM2 force field)
      - XYZ (element symbols from XYZ file format)
      - ALC (Alchemy file)
      - HAD
      - MCML
      - C3D (Chem3D)
      - SYB (Sybyl mol2)
      - MOL
      - MAP
      - DRE
      - XED (XED format)
      - DOK (Dock)
      - M3D
  */

  OBTypeTable::OBTypeTable()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "types.txt";
    _subdir = "data";
    _dataptr = TypesData;
    _linecount = 0;
    _from = _to = -1;
  }

  void OBTypeTable::ParseLine(const char *buffer)
  {
    if (buffer[0] == '#')
      return; // just a comment line

    if (_linecount == 0)
      sscanf(buffer,"%d%d",&_nrows,&_ncols);
    else if (_linecount == 1)
      tokenize(_colnames,buffer);
    else
      {
        vector<string> vc;
        tokenize(vc,buffer);
        if (vc.size() == (unsigned)_ncols)
          _table.push_back(vc);
        else
          {
            stringstream errorMsg;
            errorMsg << " Could not parse line in type translation table types.txt -- incorect number of columns";
            errorMsg << " found " << vc.size() << " expected " << _ncols << ".";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
          }
      }
    _linecount++;
  }

  bool OBTypeTable::SetFromType(const char* from)
  {
    if (!_init)
      Init();

    string tmp = from;

    unsigned int i;
    for (i = 0;i < _colnames.size();++i)
      if (tmp == _colnames[i])
        {
          _from = i;
          return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
  }

  bool OBTypeTable::SetToType(const char* to)
  {
    if (!_init)
      Init();

    string tmp = to;

    unsigned int i;
    for (i = 0;i < _colnames.size();++i)
      if (tmp == _colnames[i])
        {
          _to = i;
          return(true);
        }

    obErrorLog.ThrowError(__FUNCTION__, "Requested type column not found", obInfo);

    return(false);
  }

  //! Translates atom types (to, from), checking for size of destination
  //!  string and null-terminating as needed
  //! \deprecated Because there is no guarantee on the length of an atom type
  //!  you should consider using std::string instead
  bool OBTypeTable::Translate(char *to, const char *from)
  {
    if (!_init)
      Init();

    bool rval;
    string sto,sfrom;
    sfrom = from;
    rval = Translate(sto,sfrom);
    strncpy(to,(char*)sto.c_str(), sizeof(to) - 1);
    to[sizeof(to) - 1] = '\0';

    return(rval);
  }

  bool OBTypeTable::Translate(string &to, const string &from)
  {
    if (!_init)
      Init();

    if (from == "")
      return(false);

    if (_from >= 0 && _to >= 0 &&
        _from < _table.size() && _to < _table.size())
      {
        vector<vector<string> >::iterator i;
        for (i = _table.begin();i != _table.end();++i)
          if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
            {
              to = (*i)[_to];
              return(true);
            }
      }

    // Throw an error, copy the string and return false
    obErrorLog.ThrowError(__FUNCTION__, "Cannot perform atom type translation: table cannot find requested types.", obWarning);
    to = from;
    return(false);
  }

  std::string OBTypeTable::Translate(const string &from)
  {
    if (!_init)
      Init();

    if (!from.size())
      return("");

    if (_from >= 0 && _to >= 0 &&
        _from < _table.size() && _to < _table.size())
      {
        vector<vector<string> >::iterator i;
        for (i = _table.begin();i != _table.end();++i)
          if ((signed)(*i).size() > _from &&  (*i)[_from] == from)
            {
              return (*i)[_to];
            }
      }

    // Throw an error, copy the string and return false
    obErrorLog.ThrowError(__FUNCTION__, "Cannot perform atom type translation: table cannot find requested types.", obWarning);
    return("");
  }

  std::string OBTypeTable::GetFromType()
  {
    if (!_init)
      Init();

    if (_from > 0 && _from < _table.size())
      return( _colnames[_from] );
    else
      return( _colnames[0] );
  }

  std::string OBTypeTable::GetToType()
  {
    if (!_init)
      Init();

    if (_to > 0 && _to < _table.size())
      return( _colnames[_to] );
    else
      return( _colnames[0] );
  }

  void Toupper(string &s)
  {
    unsigned int i;
    for (i = 0;i < s.size();++i)
      s[i] = toupper(s[i]);
  }

  void Tolower(string &s)
  {
    unsigned int i;
    for (i = 0;i < s.size();++i)
      s[i] = tolower(s[i]);
  }

  ///////////////////////////////////////////////////////////////////////
  OBResidueData::OBResidueData()
  {
    _init = false;
    _dir = BABEL_DATADIR;
    _envvar = "BABEL_DATADIR";
    _filename = "resdata.txt";
    _subdir = "data";
    _dataptr = ResidueData;
  }

  bool OBResidueData::AssignBonds(OBMol &mol,OBBitVec &bv)
  {
    if (!_init)
      Init();

    OBAtom *a1,*a2;
    OBResidue *r1,*r2;
    vector<OBAtom*>::iterator i,j;
    vector3 v;

    int bo;
    unsigned int skipres=0;
    string rname = "";
    //assign residue bonds
    for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
      {
        r1 = a1->GetResidue();
        if (skipres && r1->GetNum() == skipres)
          continue;

        if (r1->GetName() != rname)
          {
            skipres = SetResName(r1->GetName()) ? 0 : r1->GetNum();
            rname = r1->GetName();
          }
        //assign bonds for each atom
        for (j=i,a2 = mol.NextAtom(j);a2;a2 = mol.NextAtom(j))
          {
            r2 = a2->GetResidue();
            if (r1->GetNum() != r2->GetNum())
              break;
            if (r1->GetName() != r2->GetName())
              break;

            if ((bo = LookupBO(r1->GetAtomID(a1),r2->GetAtomID(a2))))
              {
                v = a1->GetVector() - a2->GetVector();
                if (v.length_2() < 3.5) //check by distance
                  mol.AddBond(a1->GetIdx(),a2->GetIdx(),bo);
              }
          }
      }

    int hyb;
    string type;

    //types and hybridization
    rname = ""; // name of current residue
    skipres = 0; // don't skip any residues right now
    for (a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i))
      {
        if (a1->IsOxygen() && !a1->GetValence())
          {
            a1->SetType("O3");
            continue;
          }
        if (a1->IsHydrogen())
          {
            a1->SetType("H");
            continue;
          }

        //***valence rule for O-
        if (a1->IsOxygen() && a1->GetValence() == 1)
          {
            OBBond *bond;
            bond = (OBBond*)*(a1->BeginBonds());
            if (bond->GetBO() == 2)
              {
                a1->SetType("O2");
                a1->SetHyb(2);
              }
            else if (bond->GetBO() == 1)
              {
                a1->SetType("O-");
                a1->SetHyb(3);
                a1->SetFormalCharge(-1);
              }
            continue;
          }
	
        r1 = a1->GetResidue();
        if (skipres && r1->GetNum() == skipres)
          continue;

        if (r1->GetName() != rname)
          {
            // if SetResName fails, skip this residue
            skipres = SetResName(r1->GetName()) ? 0 : r1->GetNum();
            rname = r1->GetName();
          }

        if (LookupType(r1->GetAtomID(a1),type,hyb))
          {
            a1->SetType(type);
            a1->SetHyb(hyb);
          }
        else // try to figure it out by bond order ???
          {}
      }

    return(true);
  }

  void OBResidueData::ParseLine(const char *buffer)
  {
    int bo;
    string s;
    vector<string> vs;

    if (buffer[0] == '#')
      return;

    tokenize(vs,buffer);
    if (!vs.empty())
      {
        if (vs[0] == "BOND")
          {
            s = (vs[1] < vs[2]) ? vs[1] + " " + vs[2] :
              vs[2] + " " + vs[1];
            bo = atoi(vs[3].c_str());
            _vtmp.push_back(pair<string,int> (s,bo));
          }

        if (vs[0] == "ATOM" && vs.size() == 4)
          {
            _vatmtmp.push_back(vs[1]);
            _vatmtmp.push_back(vs[2]);
            _vatmtmp.push_back(vs[3]);
          }

        if (vs[0] == "RES")
          _resname.push_back(vs[1]);

        if (vs[0]== "END")
          {
            _resatoms.push_back(_vatmtmp);
            _resbonds.push_back(_vtmp);
            _vtmp.clear();
            _vatmtmp.clear();
          }
      }
  }

  bool OBResidueData::SetResName(const string &s)
  {
    if (!_init)
      Init();

    unsigned int i;

    for (i = 0;i < _resname.size();++i)
      if (_resname[i] == s)
        {
          _resnum = i;
          return(true);
        }

    _resnum = -1;
    return(false);
  }

  int OBResidueData::LookupBO(const string &s)
  {
    if (_resnum == -1)
      return(0);

    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();++i)
      if (_resbonds[_resnum][i].first == s)
        return(_resbonds[_resnum][i].second);

    return(0);
  }

  int OBResidueData::LookupBO(const string &s1, const string &s2)
  {
    if (_resnum == -1)
      return(0);
    string s;

    s = (s1 < s2) ? s1 + " " + s2 : s2 + " " + s1;

    unsigned int i;
    for (i = 0;i < _resbonds[_resnum].size();++i)
      if (_resbonds[_resnum][i].first == s)
        return(_resbonds[_resnum][i].second);

    return(0);
  }

  bool OBResidueData::LookupType(const string &atmid,string &type,int &hyb)
  {
    if (_resnum == -1)
      return(false);

    string s;
    vector<string>::iterator i;

    for (i = _resatoms[_resnum].begin();i != _resatoms[_resnum].end();i+=3)
      if (atmid == *i)
        {
          i++;
          type = *i;
          i++;
          hyb = atoi((*i).c_str());
          return(true);
        }

    return(false);
  }

  void OBGlobalDataBase::Init()
  {
    if (_init)
      return;
    _init = true;

    string buffer, subbuffer;
    ifstream ifs1, ifs2, ifs3, ifs4, *ifsP;
    // First, look for an environment variable
    if (getenv(_envvar.c_str()) != NULL)
      {
        buffer = getenv(_envvar.c_str());
        buffer += FILE_SEP_CHAR;

        if (!_subdir.empty())
          {
            subbuffer = buffer;
            subbuffer += _subdir;
            subbuffer += FILE_SEP_CHAR;
          }

        buffer += _filename;
        subbuffer += _filename;

        ifs1.open(subbuffer.c_str());
        ifsP= &ifs1;
        if (!(*ifsP))
          {
            ifs2.open(buffer.c_str());
            ifsP = &ifs2;
          }
      }
    // Then, check the configured data directory
    else // if (!(*ifsP))
      {
        buffer = _dir;
        buffer += FILE_SEP_CHAR;

        subbuffer = buffer;
        subbuffer += BABEL_VERSION;
        subbuffer += FILE_SEP_CHAR;

        subbuffer += _filename;
        buffer += _filename;

        ifs3.open(subbuffer.c_str());
        ifsP= &ifs3;
        if (!(*ifsP))
          {
            ifs4.open(buffer.c_str());
            ifsP = &ifs4;
          }
      }

    char charBuffer[BUFF_SIZE];
    if ((*ifsP))
      {
        while(ifsP->getline(charBuffer,BUFF_SIZE))
          ParseLine(charBuffer);
      }

    else
      // If all else fails, use the compiled in values
      if (_dataptr)
        {
          const char *p1,*p2;
          for (p1 = p2 = _dataptr;*p2 != '\0';++p2)
            if (*p2 == '\n')
              {
                strncpy(charBuffer, p1, (p2 - p1));
                charBuffer[(p2 - p1)] = '\0';
                ParseLine(charBuffer);
                p1 = ++p2;
              }
        }
      else
        {
          string s = "Unable to open data file '";
          s += _filename;
          s += "'";
          obErrorLog.ThrowError(__FUNCTION__, s, obWarning);
        }

    if (ifs1)
      ifs1.close();
    if (ifs2)
      ifs2.close();
    if (ifs3)
      ifs3.close();
    if (ifs4)
      ifs4.close();

    if (GetSize() == 0)
      {
        string s = "Cannot initialize database '";
        s += _filename;
        s += "' which may cause further errors.";
        obErrorLog.ThrowError(__FUNCTION__, "Cannot initialize database", obWarning);
      }
      
  }

  string OpenDatafile(ifstream& ifs, const string& filename)
  {
    /**Opens the filestream with the first file called filename found by looking 
       successively in the following directories:
       - the current directory
       - the directory in the environment variable BABEL_DATADIR
       or in the macro BABEL_DATADIR if the environment variable is not set
       - in a subdirectory of the BABEL_DATADIR directory with the version of 
       OpenBabel as its name
       Returns the full name of the file that was opened (i.e with path) except 
       if it is in current directory
    */
    ios_base::openmode imode = ios_base::in;
#ifdef ALL_READS_BINARY //Makes unix files compatible with VC++6
		imode = ios_base::in|ios_base::binary;
#endif
    ifs.close();
    ifs.clear();
    ifs.open(filename.c_str(),imode);
    if(ifs)
      return filename;

    string file;
    char* datadir = getenv("BABEL_DATADIR");
    if(!datadir)
      datadir = BABEL_DATADIR;
    if(datadir)
      {
        file = datadir;
        file += '/';
        file += filename;
      }

    ifs.clear();
    ifs.open(file.c_str(),imode);
    if(ifs)
      return file;;

    ifs.clear();
    ifs.open(file.c_str(),imode);
    if(ifs)
      return file;
    file = datadir;
    file += "/";
    file += BABEL_VERSION;
    file += "/" + filename;

    ifs.clear();
    ifs.open(file.c_str(),imode);
    return file;
  }

} // end namespace OpenBabel

//! \file data.cpp
//! \brief Global data and resource file parsers.
