/**********************************************************************
Copyright (C) 2005-2007 by Craig A. James, eMolecules Inc.
Some portions Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2008 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

// This code uses the old OpenEye SMILES parser
// but replaces the SMILES export with Craig James canonical smiles
// (For regular SMILES, the canonical order is not computed and ignored)

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/chiral.h>
#include <openbabel/atomclass.h>

#include <openbabel/canon.h>

#include <limits>

using namespace std;

namespace OpenBabel {

  //Base class for SMIFormat and CANSIFormat with most of the functionality
  class SMIBaseFormat : public OBMoleculeFormat
  {
  public:
    virtual const char* GetMIMEType() 
    { return "chemical/x-daylight-smiles"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ///////////////////////////////////////////////////////

    virtual const char* TargetClassDescription(){return OBMol::ClassDescription();};

    virtual const char* SpecificationURL()
    {return "http://www.daylight.com/smiles/";};

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      if(n==0) return 1; //already points after current line
      istream& ifs = *pConv->GetInStream();
      if (ifs.eof())
        return -1;
      
      int i=0;
      while(i<n && ifs.good())
        {
          if(!isNotSmiles(ifs.peek()))
            i++;
          ifs.ignore(numeric_limits<streamsize>::max(),'\n');
        }
      return ifs ? 1 : -1; 
    }

  private:
    bool isNotSmiles(char ch);

  };

  //**************************************************
  class SMIFormat : public SMIBaseFormat
  {
  public:
    //Register this format type ID
    SMIFormat()
    {
      OBConversion::RegisterFormat("smi",this, "chemical/x-daylight-smiles");
      OBConversion::RegisterFormat("smiles",this, "chemical/x-daylight-smiles");
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("t", this);
      OBConversion::RegisterOptionParam("r", this);
      OBConversion::RegisterOptionParam("a", this);
      OBConversion::RegisterOptionParam("h", this);
      OBConversion::RegisterOptionParam("x", this);
    }
    virtual const char* Description()
    {
      return
        "SMILES format\n"
        "A linear text format which can describe the connectivity\n"
        "and chirality of a molecule\n"
        "Write Options e.g. -xt\n"
        "  a  Output atomclass like [C:2], if available\n"
        "  c  Output in canonical form\n"
        "  h  explict hydrogen form, e.g. [CH3][CH3]\n"
        "  i  Do not include isotopic or chiral markings\n"
        "  n  No molecule name\n"
        "  r  Radicals lower case eg ethyl is Cc\n"
        "  t  Molecule name only\n"
        "  x  append X/Y coordinates in canonical-SMILES order\n"
        "\n";
    }


  };

  //Make an instance of the format class
  SMIFormat theSMIFormat;
  
  //**************************************************
  class CANSMIFormat : public SMIBaseFormat
  {
  public:
    //Register this format type ID
    CANSMIFormat()
    {
      OBConversion::RegisterFormat("can", this, "chemical/x-daylight-cansmiles");
    }

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
    {      
      //The "c" option sets us to use canonical ordering
      pConv->AddOption("c",OBConversion::OUTOPTIONS);
      return SMIBaseFormat::WriteMolecule(pOb, pConv);
    }

    ///////////////////////////////////////////////////////

    virtual const char* Description() {
      return
        "Canonical SMILES format.\n"
        "A linear text format which can describe the connectivity\n"
        "and chirality of a molecule, and has a single 'canonical'\n"
        "form for any particular molecule.\n"
        "Write Options e.g. -xt\n"
        "  a  Output atomclass like [C:2], if available\n"
        "  h  explict hydrogen form, e.g. [CH3][CH3]\n"
        "  i  Do not include isotopic or chiral markings\n"
        "  n  No molecule name\n"
        "  r  Radicals lower case eg ethyl is Cc\n"
        "  t  Molecule name only\n"
        "\n";
    };

  };

  // Make an instance of the format class
  CANSMIFormat theCANSMIFormat;

  //************************************************************

  class OBSmilesParser
  {
    int _bondflags;
    int _order;
    int _prev;
    char *_ptr;
    vector<int> _vprev;
    vector<vector<int> > _rclose;
    vector<vector<int> > _extbond;
    vector<int>          _path;
    vector<bool>         _avisit;
    vector<bool>         _bvisit;
    char _buffer[BUFF_SIZE];
    vector<int> PosDouble; //for extension: lc atoms as conjugated double bonds
    bool chiralWatch; // set when a chiral atom is read
    map<OBAtom*,OBChiralData*> _mapcd; // map of ChiralAtoms and their data
    OBAtomClassData _classdata; // to hold atom class data like [C:2]
  public:

    OBSmilesParser() { }
    ~OBSmilesParser() { }

    bool SmiToMol(OBMol&,string&);
    bool ParseSmiles(OBMol&);
    bool ParseSimple(OBMol&);
    bool ParseComplex(OBMol&);
    bool ParseRingBond(OBMol&);
    bool ParseExternalBond(OBMol&);
    bool CapExternalBonds(OBMol &mol);
    void FindAromaticBonds(OBMol &mol,OBAtom*,int);
    void FindAromaticBonds(OBMol&);
    void FindOrphanAromaticAtoms(OBMol &mol); //CM 18 Sept 2003
    void FixCisTransBonds(OBMol &);
    void CorrectUpDownMarks(OBBond *, OBAtom *);
  };

  /////////////////////////////////////////////////////////////////
  /*
    There is a set of characters which do not occur in SMILES strings.
    If the first char is one of this set, then line is read and discarded and
    parsing starts again with the next line.
    SMILES strings are terminated by characters from this set. If the
    terminating character is whitespace, the rest of the line is used as
    the title of the molecule and the input stream left at the start of the
    next line. If the terminating character is not whitespace, the input
    stream is left so that it will be the next character be read.
  */

  ///Returns true if character is not one used in a SMILES string.
  bool SMIBaseFormat::isNotSmiles(char ch)
  {
    static std::string notsmileschars(",<>\"\'!^&_|{}");
    return ch<=0x20 || notsmileschars.find(ch)!=string::npos;
  }

  //////////////////////////////////////////////////////////////////
  bool SMIBaseFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();

    istream &ifs = *pConv->GetInStream();
    //    const char* title = pConv->GetTitle();

    string ln;
    //Ignore lines that start with non-SMILES characters, including whitespace
    while(ifs && isNotSmiles(ifs.peek()))
      if(!getline(ifs, ln))
        return false;

    //Copy the input up to the first non-SMILES character
    string smiles;
    char ch=0;
    for(;ifs.good();ch=0) //exits with ch=0 at eof or failure
      {
        ifs.get(ch);
        if(isNotSmiles(ch))
          break;
        smiles.push_back(ch);
      }

    //when terminating char is...
    if(ch!=0 && ch!='\n') //...end of line: no title; stream ready for next line 
      {
        if(ch>0 && isspace(ch))
          {
            //...other whitespace: use the rest of the line as title
            getline(ifs, ln);
            Trim(ln);
            pmol->SetTitle(ln);
          }
        else
          //leave istream at any other terminating character
          ifs.unget();
      }

    pmol->SetDimension(0);
    OBSmilesParser sp;
    return sp.SmiToMol(*pmol, smiles);
  }


  //////////////////////////////////////////////

  bool OBSmilesParser::SmiToMol(OBMol &mol,string &s)
  {
    strncpy(_buffer,s.c_str(), BUFF_SIZE);
    _buffer[BUFF_SIZE - 1] = '\0';

    _vprev.clear();
    _rclose.clear();
    _prev=0;
    chiralWatch=false;

    if (!ParseSmiles(mol) || mol.NumAtoms() == 0)
      {
        mol.Clear();
        return(false);
      }

    mol.SetAutomaticFormalCharge(false);

    mol.SetChiralityPerceived(); //Avoid possibly buggy FindChiralCenters()

    return(true);
  }

  bool OBSmilesParser::ParseSmiles(OBMol &mol)
  {
    mol.BeginModify();
    
    for (_ptr=_buffer;*_ptr;_ptr++)
      {
        //        cerr << " parsing " << _ptr << endl;

        if (*_ptr<0 || isspace(*_ptr))
          continue;
        else if (isdigit(*_ptr) || *_ptr == '%') //ring open/close
          {
            if(!ParseRingBond(mol))
              return false;
            continue;
          }
        else if(*_ptr == '&') //external bond
          {
            ParseExternalBond(mol);
            continue;
          }
        else
          switch(*_ptr)
            {
            case '.':
              _prev=0;
              break;
            case '(':
              _vprev.push_back(_prev);
              break;
            case ')':
              if(_vprev.empty()) //CM
                return false;
              _prev = _vprev.back();
              _vprev.pop_back();
              break;
            case '[':
              if (!ParseComplex(mol))
                {
                  mol.EndModify();
                  mol.Clear();
                  return(false);
                }
              break;
            case '-':
              _order = 1;
              break;
            case '=':
              _order = 2;
              break;
            case '#':
              _order = 3;
              break;
            case ':':
              _order = 5;
              break;
            case '/':
              _bondflags |= OB_TORDOWN_BOND;   // initial mark, see FixCisTransBonds() below 
              break;
            case '\\':
              _bondflags |= OB_TORUP_BOND;     // initial mark, see FixCisTransBonds() below 
              break;
            default:
              if (!ParseSimple(mol))
                {
                  mol.EndModify();
                  mol.Clear();
                  return(false);
                }
            } // end switch
      } // end for _ptr

    // place dummy atoms for each unfilled external bond
    if(!_extbond.empty())
      CapExternalBonds(mol);

    //Save atom class values in OBGenericData object if there are any
    if(_classdata.size()>0)
      mol.SetData(new OBAtomClassData(_classdata));

    // Check to see if we've balanced out all ring closures
    // They are removed from _rclose when matched
    if ( _rclose.size() != 0) {
      mol.EndModify();
      mol.Clear();
      
      stringstream errorMsg;
      errorMsg << "Invalid SMILES string: " << _rclose.size() << " unmatched "
               << "ring bonds." << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false; // invalid SMILES since rings aren't properly closed
    }

    //set aromatic bond orders
    mol.SetAromaticPerceived();
    FindAromaticBonds(mol);
    FindOrphanAromaticAtoms(mol);// CM 18 Sept 2003
    mol.AssignSpinMultiplicity();
    mol.UnsetAromaticPerceived();

    FixCisTransBonds(mol);

    mol.EndModify();

    //Extension which interprets cccc with conjugated double bonds if niether
    //of its atoms is aromatic.
    vector<int>::iterator itr;
    for(itr=PosDouble.begin();itr!=PosDouble.end();++itr)
      {
        OBBond* bond = mol.GetBond(*itr);
        if(!bond->GetBeginAtom()->IsAromatic() && !bond->GetEndAtom()->IsAromatic())
          {
            bond->SetBO(2);
            mol.UnsetImplicitValencePerceived();
          }
      }
    
    //NE add the OBChiralData stored inside the _mapcd to the atoms now after end
    // modify so they don't get lost.
    if(_mapcd.size()>0)
      {
        OBAtom* atom;
        OBChiralData* cd;
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        for(ChiralSearch=_mapcd.begin();ChiralSearch!=_mapcd.end();ChiralSearch++)
          {
            atom=ChiralSearch->first;
            cd=ChiralSearch->second;
            atom->SetData(cd);
          }    
      }

    return(true);
  }

  void OBSmilesParser::FixCisTransBonds(OBMol &mol)
  {
    // OpenBabel's internal model for cis/trans uses an imaginary drawing,
    // in which the double bond is horizontal on the page, and an "up" bond
    // means, "above the double-bond on the page", and "down" means, "below
    // the double bond on the page".  Thus, a cis configuration can be
    // represented as either "up/up" or "down/down", and a trans
    // configuration can be either "up/down" or "down/up".
    //
    // When parsing a SMILES, '/' and '\' bonds are initially marked as
    // "down" and "up", respectively, but they don't mean "down" and "up"
    // yet -- they only mean "forward slash" and "backslash".  These bond
    // types have to be converted to an actual cis/trans specification.
    //
    // For example, consider the following trans ethene:
    //
    //   SMILES     Configuration  Initial parsing    Final designation
    //   ------     -------------  ---------------    -----------------
    //   C/C=C/C    	trans	     down/down           down/up
    //   C(\C)=C/C  	trans	     up/down             down/up
    //   C/C=C(/C)  	trans	     down/down           down/up
    //   C(=C/C)\C  	trans	     down/up             up/down
    //
    // The SMILES parsing rule is:
    //   Before double-bonded atom: '/' means down, '\' means up
    //   After  double-bonded atom: '/' means up,   '\' means down
    //
    // Note: This must be called *after* aromaticity detection.

    FOR_BONDS_OF_MOL(dbi, mol) {

      OBBond *dbl_bond = &(*dbi);

      // Not a double bond?
      if (!dbl_bond->IsDouble() || dbl_bond->IsAromatic())
        continue;

      // Find the single bonds around the atoms connected by the double bond.
      // While we're at it, note whether the pair of atoms on either end are
      // identical, in which case it's not cis/trans.

      OBAtom *a1 = dbl_bond->GetBeginAtom();
      OBAtom *a2 = dbl_bond->GetEndAtom();

      // Check that both atoms on the double bond have at least one
      // other neighbor, but not more than two other neighbors;
      int v1 = a1->GetValence();
      int v2 = a2->GetValence();
      if (v1 < 2 || v1 > 3 || v2 < 2 || v2 > 3) {
        continue;
      }

      // Get the bonds of neighbors of atom1 and atom2
      OBBond *a1_b1 = NULL, *a1_b2 = NULL, *a2_b1 = NULL, *a2_b2 = NULL;

      FOR_BONDS_OF_ATOM(bi, a1) {
        OBBond *b = &(*bi);
        if ((b) == (dbl_bond)) continue;  // skip the double bond we're working on
        if (NULL == a1_b1)
          a1_b1 = b;    // remember 1st bond of Atom1
        else
          a1_b2 = b;    // remember 2nd bond of Atom1
      }

      FOR_BONDS_OF_ATOM(bi, a2) {
        OBBond *b = &(*bi);
        if (b == dbl_bond) continue;
        if (NULL == a2_b1)
          a2_b1 = b;    // remember 1st bond of Atom2
        else
          a2_b2 = b;    // remember 2nd bond of Atom2
      }

      // Now check that at least two are marked up/down.
      int count = 0;
      if (a1_b1 && (a1_b1->IsUp() || a1_b1->IsDown())) count++;
      if (a1_b2 && (a1_b2->IsUp() || a1_b2->IsDown())) count++;
      if (a2_b1 && (a2_b1->IsUp() || a2_b1->IsDown())) count++;
      if (a2_b2 && (a2_b2->IsUp() || a2_b2->IsDown())) count++;
      if (count < 2) {
        continue;
      }

      // OK, now do what we're here for.  We have two, three or four
      // bonds, marked "up" or "down", but at this point it just
      // means '/' or '\', respectively.  In order to decide whether
      // a bond is "up" or "down", we examine the order of atoms in
      // the molecule.  See the comments at the top of this function.

      CorrectUpDownMarks(a1_b1, a1);
      CorrectUpDownMarks(a1_b2, a1);
      CorrectUpDownMarks(a2_b1, a2);
      CorrectUpDownMarks(a2_b2, a2);
    }
  }

  void OBSmilesParser::CorrectUpDownMarks(OBBond *b, OBAtom *a)
  {
    // This is an adjunct to FixCisTransBonds(), above.  See the comments
    // there.  In this function, atom a is one of the double-bonded atoms,
    // and bond b is a bond from atom a to one of the substituent atoms.

    if (!b || !a || !(b->IsUp() || b->IsDown())) return;
    
    OBAtom *ax = b->GetNbrAtom(a);

    // If the substituent comes before the double-bonded atom, then
    // the initial marks from the SMILES are correct.
    if (ax->GetIdx() < a->GetIdx())
      return;

    // The substituent comes after the double-bonded atom, so
    // the bond '/' means "up", and '\' means "down".
    if (b->IsUp()) {
      b->SetDown();
    } else {
      b->SetUp();
    }
  }


  void OBSmilesParser::FindOrphanAromaticAtoms(OBMol &mol)
  {
    //Facilitates the use lower case shorthand for radical entry
    //Atoms which are marked as aromatic but have no aromatic bonds
    //are taken to be radical centres
    OBAtom *atom;
    vector<OBAtom*>::iterator j;

    for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
      if(atom->IsAromatic())
        {
          if(atom->CountBondsOfOrder(5)<2) //bonds order 5 set in FindAromaticBonds()
            //not proper aromatic atoms - could be conjugated chain or radical centre
            atom->UnsetAromatic();
          else
            {
              //recognized as aromatic, so are not radicals
              atom->SetSpinMultiplicity(0);
            }
        }
  }

  void OBSmilesParser::FindAromaticBonds(OBMol &mol)
  {
    _path.clear();
    _avisit.clear();
    _bvisit.clear();
    _avisit.resize(mol.NumAtoms()+1);
    _bvisit.resize(mol.NumBonds());
    _path.resize(mol.NumAtoms()+1);

    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
      if (!bond->GetBeginAtom()->IsAromatic() ||
          !bond->GetEndAtom()->IsAromatic())
        _bvisit[bond->GetIdx()] = true;

    OBAtom *atom;
    vector<OBAtom*>::iterator j;

    for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
      if(!_avisit[atom->GetIdx()] && atom->IsAromatic())
        FindAromaticBonds(mol,atom,0);
  }

  void OBSmilesParser::FindAromaticBonds(OBMol &mol,OBAtom *atom,int depth )
  {
    OBBond *bond;
    vector<OBBond*>::iterator k;

    if (_avisit[atom->GetIdx()])
      {
        int j = depth-1;
        bond=mol.GetBond(_path[j--]);
        bond->SetBO(5);
        while( j >= 0 )
          {
            bond=mol.GetBond(_path[j--]);
            bond->SetBO(5);
            if(bond->GetBeginAtom() == atom || bond->GetEndAtom() == atom)
              break;
          }
      }
    else
      {
        _avisit[atom->GetIdx()] = true;
        for(bond = atom->BeginBond(k);bond;bond = atom->NextBond(k))
          if( !_bvisit[bond->GetIdx()])
            {
              _path[depth] = bond->GetIdx();
              _bvisit[bond->GetIdx()] = true;
              FindAromaticBonds(mol,bond->GetNbrAtom(atom),depth+1);
            }
      }
  }


  bool OBSmilesParser::ParseSimple(OBMol &mol)
  {
    char symbol[3];
    int element;
    bool arom=false;
    memset(symbol,'\0',sizeof(char)*3);

    if (isupper(*_ptr))
      switch(*_ptr)
        {
        case 'C':
          _ptr++;
          if (*_ptr == 'l')
            {
              strcpy(symbol,"Cl");
              element = 17;
            }
          else
            {
              symbol[0] = 'C';
              element = 6;
              _ptr--;
            }
          break;

        case 'N':
          element = 7;
          symbol[0] = 'N';
          break;
        case 'O':
          element = 8;
          symbol[0] = 'O';
          break;
        case 'S':
          element = 16;
          symbol[0] = 'S';
          break;
        case 'P':
          element = 15;
          symbol[0] = 'P';
          break;
        case 'F':
          element = 9;
          symbol[0] = 'F';
          break;
        case 'I':
          element = 53;
          symbol[0] = 'I';
          break;

        case 'B':
          _ptr++;
          if (*_ptr == 'r')
            {
              element = 35;
              strcpy(symbol,"Br");
            }
          else
            {
              element = 5;
              symbol[0] = 'B';
              _ptr--;
            }
          break;
        default:
          return(false);
        }
    else
      {
        arom = true;
        switch(*_ptr)
          {
          case 'c':
            element = 6;
            symbol[0] = 'C';
            break;
          case 'n':
            element = 7;
            symbol[0] = 'N';
            break;
          case 'o':
            element = 8;
            symbol[0] = 'O';
            break;
          case 'p':
            element = 15;
            symbol[0] = 'P';
            break;
          case 's':
            element = 16;
            symbol[0] = 'S';
            break;
          case '*':
            element = 0;
            strcpy(symbol,"Du");
            arom = false;
            break;
          case 'b':
            obErrorLog.ThrowError(__FUNCTION__, "Illegal aromatic element b", obWarning);
            element = 5;
            strcpy(symbol,"B");
            break;
          default:
            return(false);
          }
      }

    OBAtom *atom = mol.NewAtom();
    atom->SetAtomicNum(element);
    atom->SetType(symbol);

    if (arom)
      {
        atom->SetAromatic();
        atom->SetSpinMultiplicity(2); // CM 18 Sept 2003
      }
    else
      atom->ForceImplH();//ensures atom is never hydrogen deficient

    
    // Untrue, but necessary to avoid perception being called in OBAtom::IsAromatic()
    // on incomplete molecule. Undone at end of function. 
    mol.SetAromaticPerceived();
    
    if (_prev) //need to add bond
      {
        /* CM 18 Sept 2003
           An extension to the SMILES format has been added so that lower case c,n,o can 
           represent a radical centre: CcC is isopropyl radical;
           and cccc... a carbon chain bonded by conjugated double bonds.
           Fails sometimes when using c as both aromatic and as the extened form.
           For benzyl radical C6H5CH2. c1ccccc1c is ok; c1cc(c)ccc1 fails.
           Radical centres should not be involved in ring closure:
           for cyclohexyl radical C1cCCCC1 is ok, c1CCCCC1 is not.  

           Implementation
           Atoms c,n,o, etc initially added as a radical centre
           unless _prev is a radical centre when both are made a normal atoms
           connected by a double bond. But making this bond double is deferred until
           the molecule has been constructed, because it is not appropriate if
           either of the atoms is really part of an aromatic ring.

           Since they are still marked as aromatic, FindAromaticBonds() will
           replace the bonds by aromatic bonds if they are in a ring.
           FindOrphanAromand removes the aromatic tag from the atoms not found in this way
           and removes stray radical centres.
        */
        OBAtom* prevatom = mol.GetAtom(_prev);
        if(arom && prevatom->IsAromatic())
          {
            _order=5; //Potential aromatic bond
        
            if (prevatom->GetSpinMultiplicity())
              {
                //Previous atom had been marked, so bond is potentially a double bond
                //if it is not part of an aromatic ring. This will be decided when all
                //molecule has been constructed.
                PosDouble.push_back(mol.NumBonds()); //saves index of bond about to be added
                prevatom->SetSpinMultiplicity(0);
                atom->SetSpinMultiplicity(0);
              }
          }

        mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);
        
        //NE iterate through and see if atom is bonded to chiral atom
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        ChiralSearch=_mapcd.find(mol.GetAtom(_prev));
        if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
          {
            (ChiralSearch->second)->AddAtomRef(mol.NumAtoms(), input);
            // cerr << "Line 800: Adding "<<mol.NumAtoms()<<" to "<<ChiralSearch->second<<endl;
          }
      }

    //set values
    _prev = mol.NumAtoms();
    _order = 1;
    _bondflags = 0;

    mol.UnsetAromaticPerceived(); //undo 
    return(true);
  }

  bool OBSmilesParser::ParseComplex(OBMol &mol)
  {
    char symbol[7];
    int element=0;
    int isotope=0;
    int isoPtr=0;
    bool arom=false;
    memset(symbol,'\0',sizeof(char)*7);

    _ptr++;

    //grab isotope information
    for (;*_ptr && isdigit(*_ptr) && (isoPtr <= 6);_ptr++)
      {
        symbol[isoPtr] = *_ptr;
        isoPtr++;
      }
    if (isoPtr >= 6)
      return false;
    isotope = atoi(symbol);

    //parse element data
    if (isupper(*_ptr))
      switch(*_ptr)
        {
        case 'C':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element = 20;
              strcpy(symbol,"Ca");
              break;
            case 'd':
              element = 48;
              strcpy(symbol,"Cd");
              break;
            case 'e':
              element = 58;
              strcpy(symbol,"Ce");
              break;
            case 'f':
              element = 98;
              strcpy(symbol,"Cf");
              break;
            case 'l':
              element = 17;
              strcpy(symbol,"Cl");
              break;
            case 'm':
              element = 96;
              strcpy(symbol,"Cm");
              break;
            case 'o':
              element = 27;
              strcpy(symbol,"Co");
              break;
            case 'r':
              element = 24;
              strcpy(symbol,"Cr");
              break;
            case 's':
              element = 55;
              strcpy(symbol,"Cs");
              break;
            case 'u':
              element = 29;
              strcpy(symbol,"Cu");
              break;
            default:
              element =  6;
              symbol[0] = 'C';
              _ptr--;
            }
          break;

        case 'N':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element =  11;
              strcpy(symbol,"Na");
              break;
            case 'b':
              element =  41;
              strcpy(symbol,"Nb");
              break;
            case 'd':
              element =  60;
              strcpy(symbol,"Nd");
              break;
            case 'e':
              element =  10;
              strcpy(symbol,"Ne");
              break;
            case 'i':
              element =  28;
              strcpy(symbol,"Ni");
              break;
            case 'o':
              element = 102;
              strcpy(symbol,"No");
              break;
            case 'p':
              element =  93;
              strcpy(symbol,"Np");
              break;
            default:
              element =   7;
              symbol[0] = 'N';
              _ptr--;
            }
          break;

        case('O'):
          _ptr++;
          if(*_ptr == 's')
            {
              element = 76;
              strcpy(symbol,"Os");
            }
          else
            {
              element = 8;
              symbol[0] = 'O';
              _ptr--;
            }
          break;

        case 'P':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element = 91;
              strcpy(symbol,"Pa");
              break;
            case 'b':
              element = 82;
              strcpy(symbol,"Pb");
              break;
            case 'd':
              element = 46;
              strcpy(symbol,"Pd");
              break;
            case 'm':
              element = 61;
              strcpy(symbol,"Pm");
              break;
            case 'o':
              element = 84;
              strcpy(symbol,"Po");
              break;
            case 'r':
              element = 59;
              strcpy(symbol,"Pr");
              break;
            case 't':
              element = 78;
              strcpy(symbol,"Pt");
              break;
            case 'u':
              element = 94;
              strcpy(symbol,"Pu");
              break;
            default:
              element = 15;
              symbol[0] = 'P';
              _ptr--;
            }
          break;

        case('S'):
          _ptr++;
          switch(*_ptr)
            {
            case 'b':
              element = 51;
              strcpy(symbol,"Sb");
              break;
            case 'c':
              element = 21;
              strcpy(symbol,"Sc");
              break;
            case 'e':
              element = 34;
              strcpy(symbol,"Se");
              break;
            case 'i':
              element = 14;
              strcpy(symbol,"Si");
              break;
            case 'm':
              element = 62;
              strcpy(symbol,"Sm");
              break;
            case 'n':
              element = 50;
              strcpy(symbol,"Sn");
              break;
            case 'r':
              element = 38;
              strcpy(symbol,"Sr");
              break;
            default:
              element = 16;
              symbol[0] = 'S';
              _ptr--;
            }
          break;

        case 'B':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element = 56;
              strcpy(symbol,"Ba");
              break;
            case 'e':
              element =  4;
              strcpy(symbol,"Be");
              break;
            case 'i':
              element = 83;
              strcpy(symbol,"Bi");
              break;
            case 'k':
              element = 97;
              strcpy(symbol,"Bk");
              break;
            case 'r':
              element = 35;
              strcpy(symbol,"Br");
              break;
            default:
              element = 5;
              symbol[0] = 'B';
              _ptr--;
            }
          break;

        case 'F':
          _ptr++;
          switch(*_ptr)
            {
            case 'e':
              element = 26;
              strcpy(symbol,"Fe");
              break;
            case 'm':
              element = 100;
              strcpy(symbol,"Fm");
              break;
            case 'r':
              element = 87;
              strcpy(symbol,"Fr");
              break;
            default:
              element = 9;
              symbol[0] = 'F';
              _ptr--;
            }
          break;

        case 'I':
          _ptr++;
          switch(*_ptr)
            {
            case 'n':
              element = 49;
              strcpy(symbol,"In");
              break;
            case 'r':
              element = 77;
              strcpy(symbol,"Ir");
              break;
            default:
              element = 53;
              symbol[0] = 'I';
              _ptr--;
            }
          break;

        case 'A':
          _ptr++;
          switch(*_ptr)
            {
            case 'c':
              element = 89;
              strcpy(symbol,"Ac");
              break;
            case 'g':
              element = 47;
              strcpy(symbol,"Ag");
              break;
            case 'l':
              element = 13;
              strcpy(symbol,"Al");
              break;
            case 'm':
              element = 95;
              strcpy(symbol,"Am");
              break;
            case 'r':
              element = 18;
              strcpy(symbol,"Ar");
              break;
            case 's':
              element = 33;
              strcpy(symbol,"As");
              break;
            case 't':
              element = 85;
              strcpy(symbol,"At");
              break;
            case 'u':
              element = 79;
              strcpy(symbol,"Au");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'D':
          _ptr++;
          if (*_ptr == 'y')
            {
              element = 66;
              strcpy(symbol,"Dy");
            }
          else
            {
              _ptr--;
              return(false);
            }
          break;

        case 'E':
          _ptr++;
          switch(*_ptr)
            {
            case 'r':
              element = 68;
              strcpy(symbol,"Er");
              break;
            case 's':
              element = 99;
              strcpy(symbol,"Es");
              break;
            case 'u':
              element = 63;
              strcpy(symbol,"Eu");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'G':
          _ptr++;
          switch (*_ptr)
            {
            case 'a':
              element = 31;
              strcpy(symbol,"Ga");
              break;
            case 'd':
              element = 64;
              strcpy(symbol,"Gd");
              break;
            case 'e':
              element = 32;
              strcpy(symbol,"Ge");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'H':
          _ptr++;
          switch (*_ptr)
            {
            case 'e':
              element =  2;
              strcpy(symbol,"He");
              break;
            case 'f':
              element = 72;
              strcpy(symbol,"Hf");
              break;
            case 'g':
              element = 80;
              strcpy(symbol,"Hg");
              break;
            case 'o':
              element = 67;
              strcpy(symbol,"Ho");
              break;
            default:
              element = 1;
              symbol[0] = 'H';
              _ptr--;
            }
          break;

        case 'K':
          _ptr++;
          if(*_ptr == 'r')
            {
              element = 36;
              strcpy(symbol,"Kr");
            }
          else
            {
              element = 19;
              symbol[0] = 'K';
              _ptr--;
            }
          break;

        case 'L':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element =  57;
              strcpy(symbol,"La");
              break;
            case 'i':
              element =   3;
              strcpy(symbol,"Li");
              break;
            case 'r':
              element = 103;
              strcpy(symbol,"Lr");
              break;
            case 'u':
              element =  71;
              strcpy(symbol,"Lu");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'M':
          _ptr++;
          switch(*_ptr)
            {
            case 'd':
              element = 101;
              strcpy(symbol,"Md");
              break;
            case 'g':
              element =  12;
              strcpy(symbol,"Mg");
              break;
            case 'n':
              element =  25;
              strcpy(symbol,"Mn");
              break;
            case 'o':
              element =  42;
              strcpy(symbol,"Mo");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'R':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element = 88;
              strcpy(symbol,"Ra");
              break;
            case 'b':
              element = 37;
              strcpy(symbol,"Rb");
              break;
            case 'e':
              element = 75;
              strcpy(symbol,"Re");
              break;
            case 'h':
              element = 45;
              strcpy(symbol,"Rh");
              break;
            case 'n':
              element = 86;
              strcpy(symbol,"Rn");
              break;
            case 'u':
              element = 44;
              strcpy(symbol,"Ru");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case 'T':
          _ptr++;
          switch(*_ptr)
            {
            case 'a':
              element = 73;
              strcpy(symbol,"Ta");
              break;
            case 'b':
              element = 65;
              strcpy(symbol,"Tb");
              break;
            case 'c':
              element = 43;
              strcpy(symbol,"Tc");
              break;
            case 'e':
              element = 52;
              strcpy(symbol,"Te");
              break;
            case 'h':
              element = 90;
              strcpy(symbol,"Th");
              break;
            case 'i':
              element = 22;
              strcpy(symbol,"Ti");
              break;
            case 'l':
              element = 81;
              strcpy(symbol,"Tl");
              break;
            case 'm':
              element = 69;
              strcpy(symbol,"Tm");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;

        case('U'):  element = 92;
          symbol[0] = 'U';
          break;
        case('V'):  element = 23;
          symbol[0] = 'V';
          break;
        case('W'):  element = 74;
          symbol[0] = 'W';
          break;

        case('X'):
          _ptr++;
          if (*_ptr == 'e')
            {
              element = 54;
              strcpy(symbol,"Xe");
            }
          else
            {
              _ptr--;
              return(false);
            }
          break;

        case('Y'):
          _ptr++;
          if (*_ptr == 'b')
            {
              element = 70;
              strcpy(symbol,"Yb");
            }
          else
            {
              element = 39;
              symbol[0] = 'Y';
              _ptr--;
            }
          break;

        case('Z'):
          _ptr++;
          switch(*_ptr)
            {
            case 'n':
              element = 30;
              strcpy(symbol,"Zn");
              break;
            case 'r':
              element = 40;
              strcpy(symbol,"Zr");
              break;
            default:
              _ptr--;
              return(false);
            }
          break;
        }
    else
      {
        arom = true;
        switch(*_ptr)
          {
          case 'c':
            element = 6;
            symbol[0] = 'C';
            break;
          case 'n':
            element = 7;
            symbol[0] = 'N';
            break;
          case 'o':
            element = 8;
            symbol[0] = 'O';
            break;
          case 'p':
            element = 15;
            symbol[0] = 'P';
            break;
          case 'a':
            _ptr++;
            if (*_ptr == 's')
              {
                element = 33;
                strcpy(symbol,"As");
              }
            else
              return(false);
            break;
          case '*':
            element = 0;
            strcpy(symbol,"Du");
            arom = false;
            break;
          case 's': //note fall through
            _ptr++;
            if (*_ptr == 'e')
              {
                element = 34;
                strcpy(symbol,"Se");
                break;
              }
            else if (*_ptr == 'i' || *_ptr == 'n' || *_ptr == 'b')
              {
                _ptr--;
              }
            else
              {
                element = 16;
                symbol[0] = 'S';
                _ptr--;
                break;
              }
            //fall through
          default:
            strncpy(symbol, _ptr, 2);
            string symb(symbol);
            symbol[0] = toupper(symbol[0]);
            obErrorLog.ThrowError(__FUNCTION__, "Illegal aromatic element " + symb, obWarning);
            //But convert anyway
            ++_ptr;
            if(symb=="si")
              {
                element = 14;
                break;
              }
            else if(symb=="ge")
              {
                element = 32;
                break;
              }
            else if(symb=="sb")
              {
                element = 51;
                break;
              }
            else if(symb=="bi")
              {
                element = 83;
                break;
              }
            else if(symb=="te")
              {
                element = 52;
                break;
              }
            else if(symb=="sn")
              {
                element = 50;
                break;
              }
            else
              return(false);
          }
      }

    //handle hydrogen count, stereochemistry, and charge

    OBAtom *atom = mol.NewAtom();
    int hcount = 0;
    int charge=0;
    int rad=0;
    int clval=0;
    char tmpc[2];
    tmpc[1] = '\0';
    for (_ptr++;*_ptr && *_ptr != ']';_ptr++)
      {
        switch(*_ptr)
          {
          case '@':
            _ptr++;
            chiralWatch=true;
            _mapcd[atom]= new OBChiralData;
            if (*_ptr == '@')
              atom->SetClockwiseStereo();
            else
              {
                atom->SetAntiClockwiseStereo();
                _ptr--;
              }
            break;
          case '-':
            _ptr++;
            if (isdigit(*_ptr))
              {
                tmpc[0] = *_ptr;
                charge = -atoi(tmpc);
              }
            else
              {
                charge--;
                _ptr--;
              }
            break;
          case '+':
            _ptr++;
            if (isdigit(*_ptr))
              {
                tmpc[0] = *_ptr;
                charge = atoi(tmpc);
              }
            else
              {
                charge++;
                _ptr--;
              }
            break;

          case 'H':
            _ptr++;
            if (isdigit(*_ptr))
              {
                tmpc[0] = *_ptr;
                hcount = atoi(tmpc);
              }
            else
              {
                hcount = 1;
                _ptr--;
              }
            break;
          case '.': //CM Feb05
            rad=2;
            if(*(++_ptr)=='.')
              rad=3;
            else
              _ptr--;
            break;

          case ':':
            if(!isdigit(*(++_ptr)))
              {
                obErrorLog.ThrowError(__FUNCTION__,"The atom class following : must be a number", obError);
                return false;
              }
            while( isdigit(*_ptr) )
              clval = clval*10 + ((*_ptr++)-'0');
            --_ptr;
            _classdata.Add(atom->GetIdx(), clval);
            break;

          default:
            return(false);
          }
      }

    if (!*_ptr || *_ptr != ']')
      return(false); // we should have a trailing ']' now

    if (charge)
      atom->SetFormalCharge(charge);
    if (rad)
      atom->SetSpinMultiplicity(rad);
    atom->SetAtomicNum(element);
    atom->SetIsotope(isotope);
    atom->SetType(symbol);
    if (arom)
      atom->SetAromatic();

    if (_prev) //need to add bond
      {
        OBAtom* prevatom = mol.GetAtom(_prev);
        mol.SetAromaticPerceived();             // prevent aromaticity analysis
        if(arom && prevatom->IsAromatic())
          {
            _order=5; //Potential aromatic bond

            if (prevatom->GetSpinMultiplicity())
              {
                // Previous atom had been marked, so bond is potentially a double bond
                // if it is not part of an aromatic ring. This will be decided when all
                // molecule has been constructed.
                PosDouble.push_back(mol.NumBonds()); //saves index of bond about to be added
                prevatom->SetSpinMultiplicity(0);
                atom->SetSpinMultiplicity(0);
              }
          }
        mol.UnsetAromaticPerceived();
        mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);
        if(chiralWatch) // if chiral atom need to add its previous into atom4ref
          {
            if (_mapcd[atom] == NULL)
              _mapcd[atom]= new OBChiralData;
    
            (_mapcd[atom])->AddAtomRef((unsigned int)_prev,input);
            // cerr <<"line 1622: Added atom ref "<<_prev<<" to "<<_mapcd[atom]<<endl;
          }
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        ChiralSearch = _mapcd.find(mol.GetAtom(_prev));
        if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
          {
            (ChiralSearch->second)->AddAtomRef(mol.NumAtoms(), input);
            // cerr <<"line 1629: Added atom ref "<<mol.NumAtoms()<<" to "<<ChiralSearch->second<<endl;
          }
      }          

    //set values
    _prev = mol.NumAtoms();
    _order = 1;
    _bondflags = 0;

    //now add hydrogens
    if(hcount==0)
      atom->ForceNoH();//ensure AssignMultiplicity regards [C] as C atom

    for (int i = 0;i < hcount;i++)
      {
        atom = mol.NewAtom();
        atom->SetAtomicNum(1);
        atom->SetType("H");
        mol.AddBond(_prev,mol.NumAtoms(),1);
        if(chiralWatch)
          {
            if (_mapcd[mol.GetAtom(_prev)] != NULL)
              (_mapcd[mol.GetAtom(_prev)])->AddAtomRef(mol.NumAtoms(),input);
            // cerr << "line 1652: Added atom ref "<<mol.NumAtoms()<<" to "<<_mapcd[mol.GetAtom(_prev)]<<endl;
                       
          }
      }
    chiralWatch=false;
    return(true);
  }

  bool OBSmilesParser::CapExternalBonds(OBMol &mol)
  {

    if(_extbond.empty())
      return(true);

    OBAtom *atom;
    vector<vector<int> >::iterator bond;

    for(bond = _extbond.begin();bond != _extbond.end();bond++)
      {
        // create new dummy atom
        atom = mol.NewAtom();
        atom->SetAtomicNum(0);
        atom->SetType("*");

        // bond dummy atom to mol via external bond
        mol.AddBond((*bond)[1],atom->GetIdx(),(*bond)[2],(*bond)[3]);
        OBBond *refbond = atom->GetBond(mol.GetAtom((*bond)[1]));

        //record external bond information
        OBExternalBondData *xbd;
        if(mol.HasData(OBGenericDataType::ExternalBondData))
          xbd = (OBExternalBondData*)mol.GetData(OBGenericDataType::ExternalBondData);
        else
          {
            xbd = new OBExternalBondData;
            xbd->SetOrigin(fileformatInput);
            mol.SetData(xbd);
          }
        xbd->SetData(atom,refbond,(*bond)[0]);
        //this data gets cleaned up in mol.Clear.
      }

    return(true);
  }

  bool OBSmilesParser::ParseExternalBond(OBMol &mol)
  {
    int digit;
    char str[10];

    //*_ptr should == '&'
    _ptr++;

    switch (*_ptr) // check for bond order indicators CC&=1.C&1
      {
      case '-':
        _order = 1;
        _ptr++;
        break;
      case '=':
        _order = 2;
        _ptr++;
        break;
      case '#':
        _order = 3;
        _ptr++;
        break;
      case ';':
        _order = 5;
        _ptr++;
        break;
      case '/': //chiral, but _order still == 1
        _bondflags |= OB_TORDOWN_BOND;
        _ptr++;
        break;
        _ptr++;
      case '\\': // chiral, but _order still == 1
        _bondflags |= OB_TORUP_BOND;
        _ptr++;
        break;
      default: // no bond indicator just leave order = 1
        break;
      }

    if (*_ptr == '%') // external bond indicator > 10
      {
        _ptr++;
        str[0] = *_ptr;
        _ptr++;
        str[1] = *_ptr;
        str[2] = '\0';
      }
    else // simple single digit external bond indicator
      {
        str[0] = *_ptr;
        str[1] = '\0';
      }
    digit = atoi(str);  // convert indicator to digit

    //check for dot disconnect closures
    vector<vector<int> >::iterator j;
    int bondFlags,bondOrder;
    for(j = _extbond.begin();j != _extbond.end();j++)
      {
        if((*j)[0] == digit)
          {
            bondFlags = (_bondflags > (*j)[3]) ? _bondflags : (*j)[3];
            bondOrder = (_order > (*j)[2]) ? _order : (*j)[2];
            mol.AddBond((*j)[1],_prev,bondOrder,bondFlags);
            
            // after adding a bond to atom "_prev"
            // search to see if atom is bonded to a chiral atom
            map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
            ChiralSearch = _mapcd.find(mol.GetAtom(_prev));
            if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
              {
                (ChiralSearch->second)->AddAtomRef((*j)[1], input);
                // cerr << "Added external "<<(*j)[1]<<" to "<<ChiralSearch->second<<endl;
              }
            
            _extbond.erase(j);
            _bondflags = 0;
            _order = 0;
            return(true);
          }
      }

    //since no closures save another ext bond
    vector<int> vtmp(4);
    vtmp[0] = digit;
    vtmp[1] = _prev;
    vtmp[2] = _order;
    vtmp[3] = _bondflags;

    _extbond.push_back(vtmp);
    _order = 1;
    _bondflags = 0;

    return(true);

  }

  bool OBSmilesParser::ParseRingBond(OBMol &mol)
  {
    int digit;
    char str[10];

    if (*_ptr == '%')
      {
        _ptr++;
        str[0] = *_ptr;
        _ptr++;
        str[1] = *_ptr;
        str[2] = '\0';
      }
    else
      {
        str[0] = *_ptr;
        str[1] = '\0';
      }
    digit = atoi(str);

    int bf,ord;
    vector<vector<int> >::iterator j;
    for (j = _rclose.begin();j != _rclose.end();j++)
      if ((*j)[0] == digit)
        {
          bf = (_bondflags > (*j)[3]) ? _bondflags : (*j)[3];
          ord = (_order > (*j)[2]) ? _order : (*j)[2];
          // Check if this ring closure bond may be aromatic and set order accordingly
          if (ord == 1) {
            OBAtom *a1 = mol.GetAtom((*j)[1]);
            OBAtom *a2 = mol.GetAtom(_prev);
            mol.SetAromaticPerceived();                 // prevent aromaticity analysis
            if (a1->IsAromatic() && a2->IsAromatic())
              ord = 5;
            mol.UnsetAromaticPerceived();
          }
          
          mol.AddBond((*j)[1],_prev,ord,bf,(*j)[4]);
            
          // after adding a bond to atom "_prev"
          // search to see if atom is bonded to a chiral atom
          // need to check both _prev and (*j)[1] as closure is direction independent
          map<OBAtom*,OBChiralData*>::iterator ChiralSearch,cs2;
          ChiralSearch = _mapcd.find(mol.GetAtom(_prev));
          cs2=_mapcd.find(mol.GetAtom((*j)[1]));
          if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
            {
              (ChiralSearch->second)->AddAtomRef((*j)[1], input);
              // cerr << "Added ring closure "<<(*j)[1]<<" to "<<ChiralSearch->second << endl;
            }
          if (cs2!=_mapcd.end() && cs2->second != NULL)
            {
              //Ensure that the closure atom index is inserted at the position
              //decided when the ring closure digit was encountered.
              //The order needs to be SMILES atom order, not OB atom index order.
              vector<unsigned int> refs = (cs2->second)->GetAtom4Refs(input);
              // make sure the vector is large enough for the insert call
              refs.resize((*j)[4] + 1);
              refs.insert(refs.begin()+(*j)[4], _prev);
              (cs2->second)->SetAtom4Refs(refs, input);
              // cerr <<"Added ring opening "<<_prev<<" to "<<cs2->second<<endl;
            }
            
          //CM ensure neither atoms in ring closure is a radical centre
          OBAtom* patom = mol.GetAtom(_prev);
          patom->SetSpinMultiplicity(0);
          patom = mol.GetAtom((*j)[1]);
          patom->SetSpinMultiplicity(0);
          //CM end
          _rclose.erase(j);
          _bondflags = 0;
          _order = 1;
          return(true);
        }

    vector<int> vtmp(5);
    vtmp[0] = digit;
    vtmp[1] = _prev;
    vtmp[2] = _order;
    vtmp[3] = _bondflags;
    OBAtom* atom = mol.GetAtom(_prev);
    if(!atom)
      {
        obErrorLog.ThrowError(__FUNCTION__,"Number not parsed correctly as a ring bond", obError);
        return false;
      }

    vtmp[4] = atom->GetValence(); //store position to insert closure bond
    for (j = _rclose.begin();j != _rclose.end();j++) //correct for multiple closure bonds to a single atom
      if ((*j)[1] == _prev)
        vtmp[4]++;

    _rclose.push_back(vtmp);
    _order = 1;
    _bondflags = 0;

    return(true);
  }

  static bool IsCisOrTransH(OBAtom *atom)
  {
    if (!atom->IsHydrogen())
      return false;
    else
      FOR_BONDS_OF_ATOM(bond, atom)
        {
          if (bond->IsUp() || bond->IsDown())
            return true;
        }
    return false;
  }



  /*----------------------------------------------------------------------
   * CLASS: OBBondClosureInfo: For recording bond-closure digits as
   * work progresses on canonical SMILES.
   ----------------------------------------------------------------------*/

  class OBBondClosureInfo
  {
  public:
    OBAtom *toatom;       // second atom in SMILES order
    OBAtom *fromatom;     // first atom in SMILES order
    OBBond *bond;
    int    ringdigit;
    int    is_open;       // TRUE if SMILES processing hasn't reached 'toatom' yet

    OBBondClosureInfo(OBAtom *, OBAtom*, OBBond*, int, bool);
    ~OBBondClosureInfo();
  };

  OBBondClosureInfo::OBBondClosureInfo(OBAtom *a1, OBAtom *a2, OBBond *b, int rd, bool open)
  {
    toatom    = a1;
    fromatom  = a2;
    bond      = b;
    ringdigit = rd;
    is_open   = open;
  }

  OBBondClosureInfo::~OBBondClosureInfo()
  {
  }


  /*----------------------------------------------------------------------
   * CLASS: OBCanSmiNode: A Tree structure, each node of which is an atom in
   * the tree being built to write out the SMILES.
   ----------------------------------------------------------------------*/

  class OBCanSmiNode
  {
    OBAtom *_atom,*_parent;
    std::vector<OBCanSmiNode*> _child_nodes;
    std::vector<OBBond*> _child_bonds;

  public:
    OBCanSmiNode(OBAtom *atom);
    ~OBCanSmiNode();

    int Size()
    {
      return(_child_nodes.empty() ? 0 : _child_nodes.size());
    }

    void SetParent(OBAtom *a)
    {
      _parent = a;
    }

    void AddChildNode(OBCanSmiNode*,OBBond*);

    OBAtom *GetAtom()
    {
      return(_atom);
    }

    OBAtom *GetParent()
    {
      return(_parent);
    }

    OBAtom *GetChildAtom(int i)
    {
      return(_child_nodes[i]->GetAtom());
    }

    OBBond *GetChildBond(int i)
    {
      return(_child_bonds[i]);
    }

    OBCanSmiNode *GetChildNode(int i)
    {
      return(_child_nodes[i]);
    }
  };


  OBCanSmiNode::OBCanSmiNode(OBAtom *atom)
  {
    _atom = atom;
    _parent = NULL;
    _child_nodes.clear();
    _child_bonds.clear();
  }

  void OBCanSmiNode::AddChildNode(OBCanSmiNode *node,OBBond *bond)
  {
    _child_nodes.push_back(node);
    _child_bonds.push_back(bond);
  }

  OBCanSmiNode::~OBCanSmiNode()
  {
    vector<OBCanSmiNode*>::iterator i;
    for (i = _child_nodes.begin();i != _child_nodes.end();i++)
      delete (*i);
  }

  /*----------------------------------------------------------------------
   * CLASS OBMol2Cansmi - Declarations
   ----------------------------------------------------------------------*/

  class OBMol2Cansmi
  {
    std::vector<int> _atmorder;
    std::vector<bool> _aromNH;
    OBBitVec _uatoms,_ubonds;
    std::vector<OBBondClosureInfo> _vopen;
    std::string       _canorder;
    
    bool          _canonicalOutput; // regular or canonical SMILES

    OBConversion* _pconv;
    OBAtomClassData* _pac;

  public:
    OBMol2Cansmi()
    {
    }
    ~OBMol2Cansmi() {}

    void         Init(bool canonicalOutput = true, OBConversion* pconv=NULL);

    void         AssignCisTrans(OBMol*);
    char         GetCisTransBondSymbol(OBBond *, OBCanSmiNode *);
    void         AddHydrogenToChiralCenters(OBMol &mol, OBBitVec &frag_atoms);
    bool         AtomIsChiral(OBAtom *atom);
    bool         BuildCanonTree(OBMol &mol, OBBitVec &frag_atoms,
                                vector<unsigned int> &canonical_order,
                                OBCanSmiNode *node);
    void         CorrectAromaticAmineCharge(OBMol&);
    void         CreateFragCansmiString(OBMol&, OBBitVec&, bool, char *);
    bool         GetChiralStereo(OBCanSmiNode*,
                                 vector<OBAtom*>&chiral_neighbors,
                                 vector<unsigned int> &symmetry_classes,
                                 char*);
    bool         GetSmilesElement(OBCanSmiNode*,
                                  vector<OBAtom*>&chiral_neighbors,
                                  vector<unsigned int> &symmetry_classes,
                                  char*,
                                  bool isomeric);
    int          GetSmilesValence(OBAtom *atom);
    int          GetUnusedIndex();
    vector<OBBondClosureInfo>
    GetCanonClosureDigits(OBAtom *atom,
                          OBBitVec &frag_atoms,
                          vector<unsigned int> &canonical_order);
    bool         IsSuppressedHydrogen(OBAtom *atom);
    bool         SameChirality(vector<OBAtom*> &v1, vector<OBAtom*> &v2);
    void         ToCansmilesString(OBCanSmiNode *node,
                                   char *buffer,
                                   OBBitVec &frag_atoms,
                                   vector<unsigned int> &symmetry_classes,
                                   vector<unsigned int> &canonical_order,
                                   bool isomeric);

    std::string &GetOutputOrder()
    {
      return _canorder;
    }
  };


  /*----------------------------------------------------------------------
   * CLASS OBMol2Cansmi - implementation
   ----------------------------------------------------------------------*/

  /***************************************************************************
   * FUNCTION: Init
   *
   * DESCRIPTION:
   *       Initializes the OBMol2Cansmi writer object.
   ***************************************************************************/

  void OBMol2Cansmi::Init(bool canonical, OBConversion* pconv)
  {
    _atmorder.clear();
    _aromNH.clear();
    _uatoms.Clear();
    _ubonds.Clear();
    _vopen.clear();
    _canorder.clear();
    _pac = NULL;

    _pconv = pconv;
    _canonicalOutput = canonical;
  }


  /***************************************************************************
   * FUNCTION: GetUnusedIndex
   *
   * DESCRIPTION:
   *       Returns the next available bond-closure index for a SMILES.
   *
   *       You could just do this sequentially, not reusing bond-closure
   *       digits, thus:
   *
   *               c1cc2ccccc2cc1          napthalene
   *               c1ccccc1c2ccccc2        biphenyl
   *
   *       But molecules with more than ten rings, this requires the use of
   *       two-digit ring closures (like c1ccccc1C...c%11ccccc%11).  To help
   *       avoid digit reuse, this finds the lowest digit that's not currently
   *       "open", thus
   *
   *               c1cc2ccccc2cc1          napthalene (same)
   *               c1ccccc1c1ccccc1        biphenyl (reuses "1")
   *
   ***************************************************************************/


  int OBMol2Cansmi::GetUnusedIndex()
  {
    int idx=1;

    vector<OBBondClosureInfo>::iterator j;
    for (j = _vopen.begin();j != _vopen.end();)
      if (j->ringdigit == idx)
        {
          idx++; //increment idx and start over if digit is already used
          j = _vopen.begin();
        }
      else j++;

    return(idx);
  }

  /***************************************************************************
   * FUNCTION: CorrectAromaticAmineCharge
   *
   * DESCRIPTION:
   *       Finds all aromatic nitrogens, and updates the _aromNH vector to
   *       note which aromatic nitrogens require an H to balance the charge.
   ***************************************************************************/

  void OBMol2Cansmi::CorrectAromaticAmineCharge(OBMol &mol)
  {
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;

    _aromNH.clear();
    _aromNH.resize(mol.NumAtoms()+1);

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      if (atom->IsNitrogen() && atom->IsAromatic())
        if (atom->GetHvyValence() == 2)
          {
            if (atom->GetValence() == 3 || atom->GetImplicitValence() == 3)
              _aromNH[atom->GetIdx()] = true;
          }
  }

  /***************************************************************************
   * FUNCTION: OBMol2Cansmi
   *
   * DESCRIPTION:
   *       Traverse the tree searching for acyclic olefins. If an olefin is
   *       found with at least one heavy atom attachment on each end, assign
   *       stereochemistry.
   ***************************************************************************/

  void OBMol2Cansmi::AssignCisTrans(OBMol *pmol)
  {
    OBBond *bond;
    vector<OBEdgeBase*>::iterator j, k;

    FOR_BONDS_OF_MOL(dbi, pmol) {

      bond = &(*dbi);

      // Not double, or in a ring?  Skip it.
      if (bond->GetBO() != 2 || bond->IsInRing())
        continue;

      OBAtom *b = bond->GetBeginAtom();
      OBAtom *c = bond->GetEndAtom();

      // skip allenes
      if (b->GetHyb() == 1 || c->GetHyb() == 1)
        continue;

      // Skip if only hydrogen on either end (Note that GetHvyValence()
      // is counting the atom across the double bond, too, so the atom
      // must have at least two heavy atoms, i.e. at most one hydrogen.)
      if (b->GetHvyValence() < 2 || c->GetHvyValence() < 2)
        continue;

      // Ok, looks like a cis/trans double bond.

      OBAtom *a,*d;

      // Look for bond with assigned stereo as in poly-ene
      for (a = b->BeginNbrAtom(j); a; a = b->NextNbrAtom(j))
        if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown())
          break;

      if (!a) {
        for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
          if (a != c && !a->IsHydrogen())
            break;
      }
      for (d = c->BeginNbrAtom(k);d;d = c->NextNbrAtom(k)) {
        if (d != b && !d->IsHydrogen())
          break;
      }

      // Calculate the torsion angle between the "substituent" atoms.  This is an
      // odd use of the CalcTorsionAngle() function.  It measures how much you'd
      // have to twist around the double bond to bring both substituents to the
      // same side.  Cis bonds are already on the same side, so they'll have a
      // torsion angle of zero.  Trans bonds are opposite, so you'd have to twist
      // around the double bond by 180 degrees.  So small (near zero) means cis,
      // and large (near 180) means trans.  This is cool because it also works in
      // any 3D orientation.
      double angle = fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(),
                                           c->GetVector(),d->GetVector()));

      if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown()) //stereo already assigned
        {
          if (angle > 10.0) {
            // 180 degrees == trans configuration
            if (((OBBond*)*j)->IsUp())
              ((OBBond*)*k)->SetDown();	// set bonds "opposite sides" (up/down)
            else
              ((OBBond*)*k)->SetUp();		// set bonds "same side" (both up)
          }
          else {
            // small angle == cis configuration
            if (((OBBond*)*j)->IsUp())
              ((OBBond*)*k)->SetUp();
            else
              ((OBBond*)*k)->SetDown();
          }
        }
      else //assign stereo to both ends
        {
          ((OBBond*)*j)->SetUp();
          // See comments above re: angle between substituents
          if (angle > 10.0) {
            ((OBBond*)*k)->SetDown();	// trans configuration, set bonds "opposite sides" (up/down)
          } else {
            ((OBBond*)*k)->SetUp();	// cis configuration, set bonds "same side" (both up)
          }
        }
    }
  }

  char OBMol2Cansmi::GetCisTransBondSymbol(OBBond *bond, OBCanSmiNode *node)
  {
    // Given a cis/trans bond and the node in the SMILES tree, figures out
    // whether to write a '/' or '\' symbol.
    // See the comments smilesformat.cpp: FixCisTransBonds().
    // 
    // The OBCanSmiNode is the most-recently-written atom in the SMILES string
    // we're creating.  If it is the double-bonded atom, then the substituent
    // follows, so that "up" means '/' and "down" means '\'.  If the OBCanSmiNode
    // atom is the single-bonded atom then the double-bonded atom comes next,
    // in which case "up" means '\' and "down" means '/'.
    //
    // Note that there's an ambiguity: What if both ends of the bond are
    // double-bonded atoms?  In this case, the first one takes precedence.

    if (!bond || (!bond->IsUp() && !bond->IsDown()))
      return '\0';

    OBAtom *atom = node->GetAtom();

    if (atom->HasDoubleBond()) {		// double-bonded atom is first in the SMILES?
      if (bond->IsUp())
        return '/';
      else
        return '\\';
    }
    else {				// double-bonded atom is second in the SMILES
      if (bond->IsUp())
        return '\\';
      else
        return '/';
    }
  }


  // Helper function
  // Is this atom an oxygen in a water molecule
  // We know the oxygen is connected to one ion, but check for non-hydrogens
  // Returns: true if the atom is an oxygen and connected to two hydrogens + one coordinated atom
  bool isWaterOxygen(OBAtom *atom)
  {
    if (!atom->IsOxygen())
      return false;
     
    int nonHydrogenCount = 0;
    int hydrogenCount = 0;
    FOR_NBORS_OF_ATOM(neighbor, *atom) {
      if (!neighbor->IsHydrogen())
        nonHydrogenCount++;
      else
        hydrogenCount++;
    }
    
    return (hydrogenCount == 2 && nonHydrogenCount == 1);
  }


  /***************************************************************************
   * FUNCTION: GetSmilesElement
   *
   * DESCRIPTION:
   *       Writes the symbol for an atom, e.g. "C" or "[NH2]" or "[C@@H]".
   *
   * RETURNS: true (always)
   ***************************************************************************/


  bool OBMol2Cansmi::GetSmilesElement(OBCanSmiNode *node,
                                      vector<OBAtom*>&chiral_neighbors,
                                      vector<unsigned int> &symmetry_classes,
                                      char *buffer,
                                      bool isomeric)
  {
    char symbol[10];
    bool bracketElement = false;
    bool normalValence = true;
    bool writeExplicitHydrogen = false;

    OBAtom *atom = node->GetAtom();

    int bosum = atom->KBOSum();
    int maxBonds = etab.GetMaxBonds(atom->GetAtomicNum());
    // default -- bracket if we have more bonds than possible
    // we have some special cases below
    bracketElement = !(normalValence = (bosum <= maxBonds));

    switch (atom->GetAtomicNum()) {
    case 0: break;
    case 5: 
      bracketElement = !(normalValence = (bosum > 3));
      break;
    case 6: break;
    case 7:
      if (atom->IsAromatic() 
          && atom->GetHvyValence() == 2 
          && atom->GetImplicitValence() == 3) {
        bracketElement = !(normalValence = false);
        break;
      }
      else
        bracketElement = !(normalValence = (bosum == 3 || bosum == 5));
      break;
    case 8: break;
    case 9: break;
    case 15: break;
    case 16:
      bracketElement = !(normalValence = (bosum == 2 || bosum == 4 || bosum == 6));
      break;
    case 17: break;
    case 35: break;
    case 53: break;
    default: bracketElement = true;
    }

    if (atom->GetFormalCharge() != 0) //bracket charged elements
      bracketElement = true;

    if(isomeric && atom->GetIsotope())
      bracketElement = true;

    //If the molecule has Atom Class data and -xa option set and atom has data
    if(_pac && _pac->HasClass(atom->GetIdx()))
      bracketElement = true;

    char stereo[5] = "";
    if (GetSmilesValence(atom) > 2 && atom->IsChiral()) {
      if (GetChiralStereo(node, chiral_neighbors, symmetry_classes, stereo))
        strcat(buffer,stereo);
    }
    if (stereo[0] != '\0') 
      bracketElement = true;


    if (atom->GetSpinMultiplicity()) {
      //For radicals output bracket form anyway unless r option specified
      if(!(_pconv && _pconv->IsOption ("r")))
        bracketElement = true;
    }
    
    // Add brackets and explicit hydrogens for coordinated water molecules
    // PR#2505562
    if (isWaterOxygen(atom)) {
      bracketElement = true;
      writeExplicitHydrogen = true;
    }

    //Output as [CH3][CH3] rather than CC if -xh option has been specified
    if (!bracketElement && _pconv && _pconv->IsOption("h") && atom->ExplicitHydrogenCount() > 0) {
      bracketElement = true;
      writeExplicitHydrogen = true;
    }

    if (!bracketElement) {

      // Ordinary non-bracket element
      if (atom->GetAtomicNum()) {
        strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
        if (atom->IsAromatic())
          symbol[0] = tolower(symbol[0]);

        //Radical centres lc if r option set
        if(atom->GetSpinMultiplicity() && _pconv && _pconv->IsOption ("r"))
          symbol[0] = tolower(symbol[0]);
      }

      // Atomic number zero - either '*' or an external atom
      else {
        bool external = false;
        vector<pair<int,pair<OBAtom *,OBBond *> > > *externalBonds =
          (vector<pair<int,pair<OBAtom *,OBBond *> > > *)((OBMol*)atom->GetParent())->GetData("extBonds");
        vector<pair<int,pair<OBAtom *,OBBond *> > >::iterator externalBond;

        if (externalBonds)
          for(externalBond = externalBonds->begin();externalBond != externalBonds->end();externalBond++) {
            if (externalBond->second.first == atom) {
              external = true;
              strcpy(symbol,"&");
              OBBond *bond = externalBond->second.second;
              if (bond->IsUp()) {
                if ( (bond->GetBeginAtom())->HasDoubleBond() ||
                     (bond->GetEndAtom())->HasDoubleBond() )
                  strcat(symbol,"\\");
              }
              if (bond->IsDown()) {
                if ( (bond->GetBeginAtom())->HasDoubleBond() ||
                     (bond->GetEndAtom())->HasDoubleBond() )
                  strcat(symbol,"/");
              }
              if (bond->GetBO() == 2 && !bond->IsAromatic())
                strcat(symbol,"=");
              if (bond->GetBO() == 2 && bond->IsAromatic())
                strcat(symbol,":");
              if (bond->GetBO() == 3)
                strcat(symbol,"#");
              sprintf(symbol+strlen(symbol),"%d",externalBond->first);
              break;
            }
          }

        if(!external)
          strcpy(symbol,"*");
      }

      strcpy(buffer, symbol);
      return(true);
    }

    // Bracketed atoms, e.g. [Pb], [OH-], [C@]

    strcpy(buffer, "[");
    if (isomeric && atom->GetIsotope()) {
      char iso[4];
      sprintf(iso,"%d",atom->GetIsotope());
      strcat(buffer,iso);
    }
    if (!atom->GetAtomicNum())
      strcpy(symbol,"*");
    else {
      strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
      if (atom->IsAromatic())
        symbol[0] = tolower(symbol[0]);
    }
    strcat(buffer,symbol);

    // If chiral, append '@' or '@@'
    if (stereo[0] != '\0')
      strcat(buffer, stereo);

    // Add extra hydrogens.  If this is a bracket-atom *only* because the
    // "-xh" option was specified, then we're writing a SMARTS, so we
    // write the explicit hydrogen atom count only.  Otherwise we write
    // the proper total number of hydrogens.
    if (!atom->IsHydrogen()) {      
      int hcount;
      if (writeExplicitHydrogen) 
        hcount = atom->ExplicitHydrogenCount();
      else
        // if "isomeric", doesn't count isotopic H
        hcount = atom->ImplicitHydrogenCount() + atom->ExplicitHydrogenCount(isomeric);
      if (hcount != 0) {
        strcat(buffer,"H");
        if (hcount > 1) {
          char tcount[10];
          sprintf(tcount,"%d", hcount);
          strcat(buffer,tcount);
        }
      }
    }

    // Append charge to the end
    if (atom->GetFormalCharge() != 0) {
      if (atom->GetFormalCharge() > 0)
        strcat(buffer,"+");
      else
        strcat(buffer,"-");

      if (abs(atom->GetFormalCharge()) > 1)
        sprintf(buffer+strlen(buffer), "%d", abs(atom->GetFormalCharge()));
    }
    
    //atom class e.g. [C:2]
    if (_pac)
      strcat(buffer, _pac->GetClassString(atom->GetIdx()).c_str());

    strcat(buffer,"]");

    return(true);
  }

  /***************************************************************************
   * FUNCTION: OBMol2Cansmi::SameChirality
   *
   * DESCRIPTION:
   *       Given two atom vectors representing the chiral configuration around
   *       an atom, returns true/false indicating that they represent the same
   *       or opposite forms.
   *
   *       This is used when canonicalizing a SMILES.  During canonicalization,
   *       the order in which the atoms are printed is often changed, and we
   *       need to compare "before and after" to see if we've altered the
   *       chirality.
   *
   *               (NOTE: This should be integrated with OBChiralData, but
   *               isn't yet because that is a bigger project that requires
   *               rewriting this same section of smilesformat.cpp, as well as
   *               other code that uses the OBChiralData object.)
   *
   *       Throughout this code, we represent chirality as an ordered set of
   *       neighbor atoms, as follows.  Call the neighbors A, B, C, and D, and the
   *       central (chiral) atom X.  If the SMILES is A[X@](B)(C)D, then the
   *       vector could contain (A, B, C, D), in that order.  When "writing"
   *       down these vectors, we ALWAYS write them in anti-clockwise order,
   *       and we leave out the center, chiral atoms X.
   *
   *       However, there are many possible ways to write each chiral center 
   *       (hence the complexity of this function).  For example, the following
   *       all represent the exact same chirality:
   *
   *               A[X@](B)(C)D            "looking" from A to X
   *               B[X@](A)(D)C            "looking" from B to X
   *               C[X@](A)(B)D            "looking" from C to X
   *               D[X@](A)(C)B            "looking" from D to X
   *
   *       Furthermore, the choice of the second atom in the vector is arbitrary;
   *       you can "rotate" the last three atoms in the SMILES without altering
   *       the implied chirality, e.g. the following three represent the same
   *       chirality:
   *
   *               A[X@](B)(C)D
   *               A[X@](C)(D)B
   *               A[X@](D)(B)C
   *       
   *       These two sets of equalities (choice of first atom, choice of second
   *       atom) mean there are transformations of the vector that don't alter
   *       its meaning.  Using the first atom, we see that the following transformations
   *       are available:
   *       
   *               0 1 2 3                 original order
   *               1 0 3 2                 B A D C
   *               2 0 1 3                 C A B D
   *               3 0 2 1                 D A C B
   *
   *       Since the choice of the second atom is also arbitrary, by "rotatating" the
   *       last three atoms, the following transformations are available:
   *
   *               0 1 2 3                 A B C D         Original order
   *               0 2 3 1                 A C D B
   *               0 3 1 2                 A D B C
   *
   *       This function uses these transformations to determine whether two
   *       vectors represent the same or opposite chirality for a particular atom.
   *       Given two vectors, v1 and v2:
   *
   *               Transform v2 so that v2[0] == v1[0]
   *               Transform v2 so that v2[1] == v1[1]
   *
   *       After these two transformations, the third and fourth atoms of v1
   *       and v2 will either be the same or opposite, indication that the
   *       chirality represented by the vectors is the same or opposite.
   ***************************************************************************/

  bool OBMol2Cansmi::SameChirality(vector<OBAtom*> &v1, vector<OBAtom*> &v2)
  {
    vector<OBAtom*> vtmp;

    // First transform v2 so that the first atom matches v1
    if (v2[1] == v1[0]) {
      vtmp[0] = v2[1];
      vtmp[1] = v2[0];
      vtmp[2] = v2[3];
      vtmp[3] = v2[2];
      v2 = vtmp;
    }
    else if (v2[2] == v1[0]) {
      vtmp[0] = v2[2];
      vtmp[1] = v2[0];
      vtmp[2] = v2[1];
      vtmp[3] = v2[3];
      v2 = vtmp;
    }
    else if (v2[3] == v1[0]) {
      vtmp[0] = v2[3];
      vtmp[1] = v2[0];
      vtmp[2] = v2[2];
      vtmp[3] = v2[1];
      v2 = vtmp;
    }
    // else -- the first atoms already match.

    // Now rotate the last three atoms of v2 so that the
    // second atom matches v1's second atom

    if (v1[1] == v2[2]) {
      v2[2] = v2[3];
      v2[3] = v2[1];
      v2[1] = v1[1];      // use v1[1] rather than tmp var since it's got what we need
    }
    else if (v1[1] == v2[3]) {
      v2[3] = v2[2];
      v2[2] = v2[1];
      v2[1] = v1[1];      // ditto re tmp usage
    }

    // Now, the third and fourth atoms of v1/v2 are the same or opposite, indicating
    // the same or opposite chirality.
    return (v1[3] == v2[3]);
  }

  /***************************************************************************
   * FUNCTION: AtomIsChiral
   *
   * DESCRIPTION:
   *       Returns TRUE if the atom is genuinely chiral, that is, it meets
   *       the criteria from OBAtom::IsChiral, and additionally it actually
   *       has a connected hash or wedge bond.
   * 
   *       We arbitrarily reject chiral nitrogen because for our purposes there's
   *       no need to consider it.
   *
   *       NOTE: This is a simplistic test.  When the full SMILES canonicalization
   *       includes chiral markings, this should check the symmetry classes
   *       of the neighbors, not the hash/wedge bonds.
   ***************************************************************************/

  bool OBMol2Cansmi::AtomIsChiral(OBAtom *atom)
  {
    if (!atom->IsChiral())
      return false;
    if (atom->IsNitrogen())
      return false;
    // Added by ghutchis 2007-06-04 -- make sure to check for 3D molecules
    // Fixes PR#1699418
    if (atom->GetParent()->GetDimension() == 3)
      return true; // no hash/wedge for 3D molecules

    vector<int> symclass;
    FOR_BONDS_OF_ATOM(bond, atom) {
      if (bond->IsHash() || bond->IsWedge())
        return true;
    }
    return false;
  }

  /***************************************************************************
   * FUNCTION: GetChiralStereo
   *
   * DESCRIPTION:
   *       If the atom is chiral, fills in the string with either '@' or '@@',
   *       and returns true, otherwise returns false.
   ***************************************************************************/

  bool OBMol2Cansmi::GetChiralStereo(OBCanSmiNode *node,
                                     vector<OBAtom*> &chiral_neighbors,
                                     vector<unsigned int> &symmetry_classes,
                                     char *stereo)
  {
    double torsion;
    OBAtom *atom = node->GetAtom();
    OBMol *mol = (OBMol*) atom->GetParent();

    // If the molecule has no coordinates but DOES have chirality specified, it
    // must have come from a SMILES.  In this case, the atoms' GetIdx() values 
    // will be in the same order they appeared in the original SMILES, so we
    // can deduce the meaning of @ or @@ via "IsClockwise()" or "IsAnticlockwise()".
    // For example, if X is the center atom and A,B,C,D are the connected atoms,
    // appearing sequentially in the input SMILES, then A[X@](B)(C)D is
    //              
    //             B
    //            / 
    //      A -- X
    //           |\  
    //           C D (up wedge bond on D)
    //
    // and "@@" would be the opposite (with C and D switched).
    //

    if (!mol->HasNonZeroCoords()) {               // no coordinates?

      // NOTE: THIS SECTION IS WRONG.  IT'S JUST A COPY OF THE ORIGINAL OPENBABEL
      // CODE, AND DOESN'T ACCOUNT FOR THE FACT THAT THE CANONICAL SMILES IS REORDERED.
      // NEEDS TO BE REWRITTEN, BUT IN COORDINATION WITH A REWRITE OF CHIRALITY IN
      // THE smilesformat.cpp FILE.  -- CJ

      if (!atom->HasChiralitySpecified())         //   and no chirality on this atom?
        return(false);                            //   not a chiral atom -- all done.

      // Ok, it's a chiral atom, so we need to get the A, B, C, D atoms, in the order in
      // which they appeared in the original SMILES.  (NYI!!)
      if (atom->IsClockwise())
        strcpy(stereo,"@@");
      else if (atom->IsAntiClockwise())
        strcpy(stereo,"@");
      else
        return(false);
      return(true);
    }

    // If no chiral neighbors were passed in, we're done
    if (chiral_neighbors.size() < 4)
      return false;

    // If any of the neighbors have the same symmetry class, we're done.
    for (int i = 0; i < chiral_neighbors.size(); i++) {
      int idxI = chiral_neighbors[i]->GetIdx();
      int symclass = symmetry_classes[idxI-1];
      for (int j = i+1; j < chiral_neighbors.size(); j++) {
        int idxJ = chiral_neighbors[j]->GetIdx();
        if (symclass == symmetry_classes[idxJ-1])
          return false;
      }
    }

    // We have 3D coordinates for the four neighbor atoms of the chiral
    // center.  Use the "torsion angle" to deduce chirality.  If you're not
    // familiar with this, it helps to draw it on paper.  Imagine you have
    // A[X](B)(C)D.  Draw three vectors: A--X, X--B, B--C.  These three
    // vectors form a "torsion angle": If you imagine looking at X--B "end
    // on", the vectors A--X and B--C would form an angle.  If you flip the
    // chirality of X, that angle stays the same magnitude, but its sign
    // changes; thus, we can tell the chirality by whether the torsion angle
    // is positive or negative.  (Note: GetVector() is a bad name; it should
    // be called GetXYZ()).

    torsion = CalcTorsionAngle(chiral_neighbors[0]->GetVector(),
                               chiral_neighbors[1]->GetVector(),
                               chiral_neighbors[2]->GetVector(),
                               chiral_neighbors[3]->GetVector());

    strcpy(stereo,(torsion < 0.0) ? "@" : "@@");

    return(true);
  }


  /***************************************************************************
   * FUNCTION: BuildCanonTree
   *
   * DESCRIPTION:
   *       Builds the SMILES tree, in canonical order, for the specified
   *       molecular fragment.
   ***************************************************************************/

  bool OBMol2Cansmi::BuildCanonTree(OBMol &mol,
                                    OBBitVec &frag_atoms,
                                    vector<unsigned int> &canonical_order,
                                    OBCanSmiNode *node)
  {
    vector<OBEdgeBase*>::iterator i;
    OBAtom *nbr, *atom;
    vector<OBAtom *> sort_nbrs;
    vector<OBAtom *>::iterator ai;
    OBBond *bond;
    OBCanSmiNode *next;
    int idx;

    atom = node->GetAtom();

#if DEBUG
    cout << "BuildCanonTree: " << etab.GetSymbol(atom->GetAtomicNum()) << ", " << atom->GetIdx() << ", canorder " << canonical_order[atom->GetIdx()-1] << "\n";
#endif

    // Create a vector of neighbors sorted by canonical order, but favor
    // double and triple bonds over single and aromatic.  This causes
    // ring-closure digits to avoid double and triple bonds.
    //
    // Since there are typically just one to three neighbors, we just do a
    // ordered insertion rather than sorting.

    for (nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i)) {

      idx = nbr->GetIdx();
      if (nbr->IsHydrogen() && IsSuppressedHydrogen(nbr)) {
        _uatoms.SetBitOn(nbr->GetIdx());        // mark suppressed hydrogen, so it won't be considered
        continue;                               // later when looking for more fragments.
      }
      if (_uatoms[idx] || !frag_atoms.BitIsOn(idx))
        continue;

      OBBond *nbr_bond = atom->GetBond(nbr);
      int new_needs_bsymbol = nbr_bond->IsDouble() || nbr_bond->IsTriple();

      for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ai++) {
        bond = atom->GetBond(*ai);
        int sorted_needs_bsymbol = bond->IsDouble() || bond->IsTriple();
        if (new_needs_bsymbol && !sorted_needs_bsymbol) {
          sort_nbrs.insert(ai, nbr);
          ai = sort_nbrs.begin();//insert invalidated ai; set it to fail next test 
          break;
        }
        if (   new_needs_bsymbol == sorted_needs_bsymbol
               && canonical_order[idx-1] < canonical_order[(*ai)->GetIdx()-1]) {
          sort_nbrs.insert(ai, nbr);
          ai = sort_nbrs.begin();//insert invalidated ai; set it to fail next test 
          break;
        }
      }
      if (ai == sort_nbrs.end())
        sort_nbrs.push_back(nbr);
    }

    _uatoms.SetBitOn(atom->GetIdx());     //mark the atom as visited

    // Build the next layer of nodes, in canonical order
    for (ai = sort_nbrs.begin(); ai != sort_nbrs.end(); ai++) {
      nbr = *ai;
      idx = nbr->GetIdx();
      if (_uatoms[idx])   // depth-first search may have used this atom since
        continue;         // we sorted the bonds above
      bond = atom->GetBond(nbr);
      _ubonds.SetBitOn(bond->GetIdx());
      next = new OBCanSmiNode(nbr);
      next->SetParent(atom);
      node->AddChildNode(next, bond);
      BuildCanonTree(mol, frag_atoms, canonical_order, next);
    }

    return(true);
  }



  /***************************************************************************
   * FUNCTION: GetCanonClosureDigits
   *
   * DESCRIPTION:
   *       Given an atom, returns the ring-closure digits for that atom, in
   *       the form of a vector of digit/OBBond* pair.  Some of the digits may
   *       be for newly-opened rings (the matching digit occurs later in the
   *       SMILES string), and some may be for closing rings (the matching
   *       digit occured earlier in the string).
   *
   *       Canonicalization requires that atoms with more than one digit
   *       have the digits assigned in a canonical fashion.  For example,
   *       the SMILES  "CC12(NCCC2)CCC1" and "CC12(NCCC1)CCC2" are the
   *       same molecule; we need to assign the digits of the first "C12"
   *       such that it always comes out one way or the other.
   *
   *       This needs to find closing bonds (ring bonds already assigned a 
   *       digit) and opening bonds (ring bonds not encountered yet).
   *
   *    Closing Bonds:
   *       This is easy: open bonds are already stored in the _vopen vector, 
   *       in canonical order.  Just find open bonds to this atom and copy
   *       them from _vopen to our return vector.
   *
   *    Opening Bonds:
   *       This function looks through the bonds for this atoms and finds
   *       any that aren't on the _ubonds "used" list, (and also are non-H
   *       and are in this fragment).  Any such bonds must be ring-closure
   *       bonds.  If there is more than one, they are ordered by the
   *       canonical order of the bonds' neighbor atoms; that is, the bond 
   *       to the lowest canonical-ordered neighbor is assigned the first
   *       available number, and upwards in neighbor-atom canonical order.
   ***************************************************************************/

  vector<OBBondClosureInfo>
  OBMol2Cansmi::GetCanonClosureDigits(OBAtom *atom,
                                      OBBitVec &frag_atoms,
                                      vector<unsigned int> &canonical_order)
  {
    vector<OBBondClosureInfo> vp_closures;
    vector<OBBond*> vbonds;
    vector<OBBond*>::iterator bi;
    vector<OBEdgeBase*>::iterator i;
    OBBond *bond1, *bond2;
    OBAtom *nbr1, *nbr2;
    int nbr1_canorder, nbr2_canorder;

    vp_closures.clear();
    vbonds.clear();

    // Find new ring-closure bonds for this atom
    for (bond1 = atom->BeginBond(i); bond1; bond1 = atom->NextBond(i)) {

      // Is this a ring-closure neighbor?
      if (_ubonds.BitIsOn(bond1->GetIdx()))
        continue;
      nbr1 = bond1->GetNbrAtom(atom);
      // Skip hydrogens before checking canonical_order
      // PR#1999348
      if (   (nbr1->IsHydrogen() && IsSuppressedHydrogen(nbr1))
             || !frag_atoms.BitIsOn(nbr1->GetIdx()))
        continue;

      nbr1_canorder = canonical_order[nbr1->GetIdx()-1];

      // Insert into the bond-vector in canonical order (by neighbor atom order)
      for (bi = vbonds.begin(); bi != vbonds.end(); bi++) {
        bond2 = *bi;
        nbr2 = bond2->GetNbrAtom(atom);
        nbr2_canorder = canonical_order[nbr2->GetIdx()-1];
        if (nbr1_canorder < nbr2_canorder) {
          vbonds.insert(bi, bond1);
          bi = vbonds.begin();//insert invalidated bi; set it to fail next test
          break;
        }
      }
      if (bi == vbonds.end())     // highest one (or first one) - append to end
        vbonds.push_back(bond1);
    }

    // If we found new open bonds, assign a bond-closure digits to each one,
    // add it to _vopen, and add it to the return vector.
    for (bi = vbonds.begin(); bi != vbonds.end(); bi++) {
      bond1 = *bi;
      _ubonds.SetBitOn(bond1->GetIdx());
      int digit = GetUnusedIndex();
      int bo = (bond1->IsAromatic())? 1 : bond1->GetBO();  // CJ: why was this line added?  bo is never used?
      _vopen.push_back(OBBondClosureInfo(bond1->GetNbrAtom(atom), atom, bond1, digit, true));
      vp_closures.push_back(OBBondClosureInfo(bond1->GetNbrAtom(atom), atom, bond1, digit, true));
    }


    // Now look through the list of open closure-bonds and find any to this
    // atom (but watch out for the ones we just added).  For each one found,
    // add it to the return vector, and erase it from _vopen.

    if (!_vopen.empty()) {
      vector<OBBondClosureInfo>::iterator j;
      for (j = _vopen.begin(); j != _vopen.end(); ) {
        if (j->toatom == atom) {
          OBBondClosureInfo bci = *j;
          _vopen.erase(j);                // take bond off "open" list
          bci.is_open = false;            // mark it "closed"
          vp_closures.push_back(bci);     // and add it to this atom's list
          j = _vopen.begin();             // reset iterator
        }
        else
          j++;
      }
    }

    return(vp_closures);
  }


  /***************************************************************************
   * FUNCTION: IsSuppressedHydrogen
   *
   * DESCRIPTION:
   *       For a hydrogen atom, returns TRUE if the atom is not [2H] or [3H], only
   *       has one bond, and is not bonded to another hydrogen. 
   *
   *       NOTE: Return value is nonsensical if you pass it a non-hydrogen
   *       atom.  Presumably, you're calling this because you've found a 
   *       hydrogen and want to know if it goes in the SMILES.
   ***************************************************************************/

  bool OBMol2Cansmi::IsSuppressedHydrogen(OBAtom *atom)
  {
    if (atom->GetIsotope() != 0)          // Deuterium or Tritium
      return false;
    if (atom->GetValence() != 1)          // not exactly one bond
      return false;

    FOR_NBORS_OF_ATOM(nbr, atom) {
      if (nbr->GetAtomicNum() == 1)       // neighbor is hydrogen
        return false;
    }

    return true;
  }

  /***************************************************************************
   * FUNCTION: GetSmilesValence
   *
   * DESCRIPTION:
   *       This is like GetHvyValence(), but it returns the "valence" of an
   *       atom as it appears in the SMILES string.  In particular, hydrogens
   *       count if they will appear explicitly -- see IsSuppressedHydrogen()
   *       above.
   ***************************************************************************/

  int OBMol2Cansmi::GetSmilesValence(OBAtom *atom)
  {
    int count = 0;

    if (atom->IsHydrogen())
      return atom->GetValence();

    if (_pconv && _pconv->IsOption("h"))
      return atom->GetValence();

    FOR_NBORS_OF_ATOM(nbr, atom) {
      if (  !nbr->IsHydrogen()
            || nbr->GetIsotope() != 0
            || nbr->GetValence() != 1)
        count++;
    }
    return(count);
  }


  /***************************************************************************
   * FUNCTION: ToCansmilesString
   *
   * DESCRIPTION:
   *       Recursively writes the canonical SMILES string to a buffer.  Writes
   *       this node, then selects each of the child nodes (in canonical
   *       order) and writes them.
   *
   *       Chirality is the tricky bit here.  Before we can write out a chiral
   *       atom, we have to "look ahead" to determine the order in which the
   *       neighbor atoms are/will be written.
   *
   *       The SMILES language defines the order-of-appearance of a ring-closure
   *       bond as the position of the digit, in the SMILES, not the actual atom.
   *       For example, the fragments N[C@H](C)Br, and N[C@H]1(Br)CCCC1 have
   *       the same chiral center, because the "1" in the second one is a "stand
   *       in" for the "C" in the first, even though the actual carbon atom appears
   *       after the Bromine atom in the second string.
   ***************************************************************************/

  void OBMol2Cansmi::ToCansmilesString(OBCanSmiNode *node,
                                       char *buffer,
                                       OBBitVec &frag_atoms,
                                       vector<unsigned int> &symmetry_classes,
                                       vector<unsigned int> &canonical_order,
                                       bool isomeric)
  {
    OBAtom *atom = node->GetAtom();
    vector<OBAtom *> chiral_neighbors;

    // Get the ring-closure digits in canonical order.  We'll use these in
    // two places: First, for figuring out chirality, then later for writing
    // the actual ring-closure digits to the string.
    vector<OBBondClosureInfo> vclose_bonds = GetCanonClosureDigits(atom, frag_atoms, canonical_order);

    // First thing: Figure out chirality.  We start by creating a vector of the neighbors
    // in the order in which they'll appear in the canonical SMILES string.  This is more
    // complex than you'd guess because of implicit/explicit H and ring-closure digits.

    bool is_chiral = AtomIsChiral(atom);
    if (is_chiral) {

      // If there's a parent node, it's the first atom in the ordered neighbor-vector
      // used for chirality.
      if (node->GetParent()) {
        chiral_neighbors.push_back(node->GetParent());
      }

      // Next for chirality order will be hydrogen -- since it occurs
      // inside the atom's [] brackets, it's always before other neighbors.
      //
      // Note that we check the regular neighbor list, NOT the canonical
      // SMILES tree, because hydrogens normally aren't part of the canonical
      // SMILES, but we still need them to figure out chirality.
      //
      // There are two cases: it's explicit in the OBMol object but should be
      // written inside the brackets, i.e. "[C@H]", or it is explicit and
      // must be outside the brackets, such as for deuterium.  (A hydrogen
      // that will appear explicitly in the SMILES as a separate atom is
      // treated like any other atom when calculating the chirality.)

      FOR_NBORS_OF_ATOM(i_nbr, atom) {
        OBAtom *nbr = &(*i_nbr);
        if (nbr->IsHydrogen() && IsSuppressedHydrogen(nbr) ) {
          chiral_neighbors.push_back(nbr);
          break;        // quit loop: only be one H if atom is chiral
        }
      }

      // Ok, done with H.  Next in the SMILES will be the ring-closure characters.
      // So we need to find the corresponding atoms and add them to the list.
      // (We got the canonical ring-closure list earlier.)
      if (!vclose_bonds.empty()) {
        vector<OBBondClosureInfo>::iterator i;
        for (i = vclose_bonds.begin();i != vclose_bonds.end();i++) {
          OBBond *bond = i->bond;
          OBAtom *nbr = bond->GetNbrAtom(atom);
          chiral_neighbors.push_back(nbr);
        }
      }

      // Finally, add the "regular" neighbors, the "child" nodes in the
      // canonical-SMILES tree, to the chiral-neighbors list.
      for (int i = 0; i < node->Size(); i++) {
        OBAtom *nbr = node->GetChildAtom(i);
        chiral_neighbors.push_back(nbr);
      }
    }

    // Write the current atom to the string
    GetSmilesElement(node, chiral_neighbors, symmetry_classes, buffer+strlen(buffer), isomeric);

    _atmorder.push_back(atom->GetIdx());  //store the atom ordering

    // Write ring-closure digits
    if (!vclose_bonds.empty()) {
      vector<OBBondClosureInfo>::iterator bci;
      for (bci = vclose_bonds.begin();bci != vclose_bonds.end();bci++) {
        if (!bci->is_open) {
          char bs[2];
          bs[0] = GetCisTransBondSymbol(bci->bond, node);
          bs[1] = '\0';
          if (bs[0]) {
            strcat(buffer, bs);	// append "/" or "\"
          }
          else {
            if (bci->bond->GetBO() == 2 && !bci->bond->IsAromatic())  strcat(buffer,"=");
            if (bci->bond->GetBO() == 3)                              strcat(buffer,"#");
          }
        }
        if (bci->ringdigit > 9) strcat(buffer,"%");
        sprintf(buffer+strlen(buffer), "%d", bci->ringdigit); 
      }
    }

    // Write child bonds, then recursively follow paths to child nodes
    // to print the SMILES for each child branch.

    OBBond *bond;
    for (int i = 0;i < node->Size();i++) {
      bond = node->GetChildBond(i);
      if (i+1 < node->Size()) {
        strcat(buffer,"(");
      }
      if (bond->IsUp() || bond->IsDown()) {
        char cc[2];
        cc[0] = GetCisTransBondSymbol(bond, node);
        cc[1] = '\0';
        strcat(buffer, cc);
      }
      else if (bond->GetBO() == 2 && !bond->IsAromatic())
        strcat(buffer,"=");
      else if (bond->GetBO() == 3)
        strcat(buffer,"#");

      ToCansmilesString(node->GetChildNode(i),buffer, frag_atoms, symmetry_classes, canonical_order, isomeric);
      if (i+1 < node->Size()) strcat(buffer,")");
    }
  }

  /****************************************************************************
   * FUNCTION: StandardLabels
   * 
   * DESCRIPTION:
   *        Creates a set of non-canonical labels for the fragment atoms
   * ***************************************************************************/
  void StandardLabels(OBMol *pMol, OBBitVec &frag_atoms, 
                      vector<unsigned int> &symmetry_classes,
                      vector<unsigned int> &labels)
  {		
    FOR_ATOMS_OF_MOL(atom, *pMol) {
      if (frag_atoms.BitIsOn(atom->GetIdx())) {
        labels.push_back(atom->GetIdx() - 1);
        symmetry_classes.push_back(atom->GetIdx() - 1);
      }
      else{
        labels.push_back(2147483647); //to match situation when canonical ordering. Just a big number?
        symmetry_classes.push_back(2147483647);
      }
    }
  }

  /***************************************************************************
   * FUNCTION: CreateFragCansmiString
   *
   * DESCRIPTION:
   *       Selects the "root" atom, which will be first in the SMILES, then
   *       builds a tree in canonical order, and finally generates the SMILES.
   *       If there are then atoms that haven't been visited (i.e. a molecule
   *       with disconnected parts), selects a new root from the remaining
   *       atoms and repeats the process.
   ***************************************************************************/

  void OBMol2Cansmi::CreateFragCansmiString(OBMol &mol, OBBitVec &frag_atoms, bool isomeric, char *buffer)
  {
    OBAtom *atom;
    OBCanSmiNode *root;
    buffer[0] = '\0';
    vector<OBNodeBase*>::iterator ai;
    vector<unsigned int> symmetry_classes, canonical_order;

    //Pointer to Atom Class data set if -xa option and the molecule has any; NULL otherwise.
    if(_pconv->IsOption("a"))
      _pac = static_cast<OBAtomClassData*>(mol.GetData("Atom Class"));

    // First, create a canonical ordering vector for the atoms.  Canonical
    // labels are zero indexed, corresponding to "atom->GetIdx()-1".
    if (_canonicalOutput)
      CanonicalLabels(&mol, frag_atoms, symmetry_classes, canonical_order);
    else {
      StandardLabels(&mol, frag_atoms, symmetry_classes, canonical_order);
    }

    // OUTER LOOP: Handles dot-disconnected structures.  Finds the 
    // lowest unmarked canorder atom, and starts there to generate a SMILES.
    // Repeats until no atoms remain unmarked.

    while (1) {

      // It happens that the lowest canonically-numbered atom is usually 
      // a good place to start the canonical SMILES.
      OBAtom *root_atom;
      int lowest_canorder = 999999;
      root_atom = NULL;
      for (atom = mol.BeginAtom(ai); atom; atom = mol.NextAtom(ai)) {
        int idx = atom->GetIdx();
        if (!atom->IsHydrogen()       // don't start with a hydrogen
            && !_uatoms[idx]          // skip atoms already used (for fragments)
            && frag_atoms.BitIsOn(idx)// skip atoms not in this fragment
            //&& !atom->IsChiral()    // don't use chiral atoms as root node
            && canonical_order[idx-1] < lowest_canorder) {
          root_atom = atom;
          lowest_canorder = canonical_order[idx-1];
        }
      }

      // If we didn't pick an atom, it is because the fragment is made
      // entirely of hydrogen atoms (e.g. [H][H]).  Repeat the loop but
      // allow hydrogens this time.
      if (root_atom == NULL) {
        for (atom = mol.BeginAtom(ai); atom; atom = mol.NextAtom(ai)) {
          int idx = atom->GetIdx();
          if (!_uatoms[idx]           // skip atoms already used (for fragments)
              && frag_atoms.BitIsOn(idx)// skip atoms not in this fragment
              && canonical_order[idx-1] < lowest_canorder) {
            root_atom = atom;
            lowest_canorder = canonical_order[idx-1];
          }
        }
      }

      // No atom found?  We've done all fragments.
      if (root_atom == NULL) 
        break;

      // Clear out closures in case structure is dot disconnected
      //      _atmorder.clear();
      _vopen.clear();

      // Dot disconnected structure?
      if (strlen(buffer) > 0) strcat(buffer,"."); 
      root = new OBCanSmiNode (root_atom);

      BuildCanonTree(mol, frag_atoms, canonical_order, root);
      ToCansmilesString(root, buffer, frag_atoms, symmetry_classes, canonical_order, isomeric);
      delete root;
    }

    // save the canonical order as a space-separated string
    // which will be returned by GetOutputOrder() for incorporation
    // into an OBPairData keyed "canonical order"
    if (_atmorder.size()) {
      stringstream temp;
      vector<int>::iterator can_iter = _atmorder.begin();
      if (can_iter != _atmorder.end()) {
        temp << (*can_iter++);
      }

      for (; can_iter != _atmorder.end(); ++can_iter) {
        if (*can_iter <= mol.NumAtoms())
          temp << " " << (*can_iter);
      }
      _canorder = temp.str(); // returned by GetOutputOrder()
    }
  }

  /***************************************************************************
   * FUNCTION: OBMol2Cansmi::AddHydrogenToChiralCenters
   *
   * DESCRIPTION:
   *       Adds an explicit hydrogen to any chiral center that only has three
   *       atoms.  This makes analysis much easier since the algorithms can
   *       assume that all tetrahedral carbons have four neighbors.
   ***************************************************************************/

  void OBMol2Cansmi::AddHydrogenToChiralCenters(OBMol &mol, OBBitVec &frag_atoms)
  {
    bool is_modified = false;
    vector <OBAtom *> atomList;

    // Find all appropriate atoms to add hydrogens
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        if (!frag_atoms[atom->GetIdx()] || !AtomIsChiral(&*atom))
          continue;

        if (GetSmilesValence(&*atom) == 3 && atom->GetValence() == 3) {       // implicit H?
          atomList.push_back(&*atom);
        }
      }

    // Now add hyrdogens to the list
    if (atomList.size() > 0) {
      mol.BeginModify();

      vector<OBAtom*>::iterator i;
      //      OBAtom *atom;
      for (i = atomList.begin(); i != atomList.end(); ++i) {
        // Get the (x,y,z) coordinates where best to put the H
        vector3 v;
        (*i)->GetNewBondVector(v, 1.0);   // Returns (x,y,z) of the "empty" area, for a new bond

#if DEBUG
        cout << "AddHydrogenToChiralCenters: Adding H to atom " << atom->GetIdx() << "\n";
#endif

        // Add the H atom
        OBAtom *h = mol.NewAtom();
        h->SetAtomicNum(1);
        h->SetType("H");
        mol.AddBond((*i)->GetIdx(), h->GetIdx(), 1, 0, -1);

        // Set its (x,y,z) coordinates
        h->SetVector(v);

        frag_atoms.SetBitOn(h->GetIdx());
      }

      mol.EndModify();
    }
  }

  /*----------------------------------------------------------------------
   * END OF CLASS: OBMol2Cansmi
   ----------------------------------------------------------------------*/



  /***************************************************************************
   * FUNCTION: CreateCansmiString
   *
   * DESCRIPTION:
   *       Writes the canonical SMILES for a molecule or molecular fragment
   *       to the given buffer.
   *
   *       frag_atoms represents atoms in a fragment of the molecule; the
   *       SMILES will contain those atoms only.
   *
   *       (Note: This is an ordinary public C++ function, not a member
   *       of any class.)
   *
   ***************************************************************************/

  void CreateCansmiString(OBMol &mol, char *buffer, OBBitVec &frag_atoms, bool iso, OBConversion* pConv)
  {
    bool canonical = pConv->IsOption("c")!=NULL;

    // This is a hack to prevent recursion problems.
    //  we still need to fix the underlying problem -GRH
    if (mol.NumAtoms() > 1000) {
#ifdef HAVE_SSTREAM
      stringstream errorMsg;
#else
      strstream errorMsg;
#endif
      errorMsg <<
        "SMILES Conversion failed: Molecule is too large to convert."
        "Open Babel is currently limited to 1000 atoms." << endl;
      errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return;
    }

    // If we're doing isomeric (stereo), make a copy
    OBMol *pmol;
    if (iso) 
      pmol = new OBMol(mol);
    else
      pmol = &mol;

    OBMol2Cansmi m2s;
    m2s.Init(canonical, pConv);
    // GRH Added 208-06-05
    // This came from the WriteMolecule call below
    // It doesn't seem to have any effect
    m2s.CorrectAromaticAmineCharge(mol);

    // Figure out Cis/Trans 
    if (mol.Has2D()) // i.e. 2D or 3D
      m2s.AssignCisTrans(pmol);

    // If the molecule has 2D coordinate AND has hash/wedge bonds,
    // create pseudo-Z coordinates by "pushing" the up/down bonds
    // to +/-1 in the Z direction.  This will be used in the next
    // section when we deduce chirality from the coordinates.

    if (iso) {
      if (!pmol->Has3D()) {

        FOR_ATOMS_OF_MOL(iatom, *pmol) {
          OBAtom *atom = &(*iatom);

          if (!atom->IsChiral()) continue;
          if (m2s.GetSmilesValence(atom) < 3) continue;

          vector3 v;
          OBAtom *nbr;

          FOR_BONDS_OF_ATOM(bond, atom) {

            // The bond's "start atom" is the pointy end of the hash or wedge
            // bond, so we need to know whether the pointy end of the bond is
            // toward the center atom (normal case) or toward the neighbor atom
            // (poor drawing style, but it happens).  The amount to push up/down
            // is "z", and is normally 1.0, but is set to 0.5 for non-terminal
            // atoms.  This keeps adjacent chiral centers from screwing each other up.

            nbr = bond->GetNbrAtom(atom);
            double z = (nbr->GetHvyValence() > 1) ? 0.5 : 1.0;
            v = nbr->GetVector();
            if (bond->GetBeginAtom() == atom) {       // The pointy end is at the central atom
              if (bond->IsWedge())
                v.SetZ(z);
              else if (bond->IsHash())
                v.SetZ(-z);
            }
            else {                                    // The pointy end is at the neighbor atom
              if (bond->IsWedge())
                v.SetZ(-z);
              else if (bond->IsHash())
                v.SetZ(z);
            }
            nbr->SetVector(v);
          }
        }
      }

      m2s.AddHydrogenToChiralCenters(*pmol, frag_atoms);
    }

    else {
      // Not isomeric - be sure there are no Z coordinates, clear
      // all stereo-center and cis/trans information.
      OBBond *bond;
      OBAtom *atom;
      vector<OBEdgeBase*>::iterator bi;
      vector<OBNodeBase*>::iterator ai;
      for (bond = pmol->BeginBond(bi); bond; bond = pmol->NextBond(bi)) {
        bond->UnsetUp();
        bond->UnsetDown();
        bond->UnsetHash();
        bond->UnsetWedge();
      }
      for (atom = pmol->BeginAtom(ai); atom; atom = pmol->NextAtom(ai)) {
        atom->UnsetStereo();
        vector3 v = atom->GetVector();
        if (v[2] != 0.0) {
          v.SetZ(0.0);
          atom->SetVector(v);
        }
      }
    }

    // If the fragment includes ordinary hydrogens, get rid of them.
    // They won't appear in the SMILES anyway (unless they're attached to
    // a chiral center, or it's something like [H][H]).
    FOR_ATOMS_OF_MOL(iatom, *pmol) {
      OBAtom *atom = &(*iatom);
      if (frag_atoms.BitIsOn(atom->GetIdx()) && atom->IsHydrogen() 
          && (!iso || m2s.IsSuppressedHydrogen(atom))) {
        frag_atoms.SetBitOff(atom->GetIdx());
      }
    }

    m2s.CreateFragCansmiString(*pmol, frag_atoms, iso, buffer);
    if (iso) {
      pmol->Clear();
      delete pmol;
    }

    // Could also save canonical bond order if anyone desires
    if (canonical && !mol.HasData("Canonical Atom Order")) {
      OBPairData *canData = new OBPairData;
      canData->SetAttribute("Canonical Atom Order");
      canData->SetValue(m2s.GetOutputOrder());
      mol.SetData(canData);
    }
  }

  //////////////////////////////////////////////////
  bool SMIBaseFormat::WriteMolecule(OBBase* pOb,OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    // Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    // Title only option?
    if(pConv->IsOption("t")) {
      ofs << mol.GetTitle() <<endl;
      return true;
    }

    char buffer[BUFF_SIZE];
    *buffer = '\0'; // clear the buffer

    // This is a hack to prevent recursion problems.
    //  we still need to fix the underlying problem (mainly chiral centers) -GRH
    if (mol.NumAtoms() > 1000) {
#ifdef HAVE_SSTREAM
      stringstream errorMsg;
#else
      strstream errorMsg;
#endif
      errorMsg <<
        "SMILES Conversion failed: Molecule is too large to convert."
        "Open Babel is currently limited to 1000 atoms." << endl;
      errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return(false);
    }

    // If there is data attached called "SMILES_Fragment", then it's
    // an ascii OBBitVec, representing the atoms of a fragment.  The
    // SMILES generated will only include these fragment atoms.

    OBBitVec fragatoms(mol.NumAtoms());

    OBPairData *dp = (OBPairData *) mol.GetData("SMILES_Fragment");
    if (dp) {
      fragatoms.FromString(dp->GetValue(), mol.NumAtoms());
    }

    // If no "SMILES_Fragment" data, fill the entire OBBitVec
    // with 1's so that the SMILES will be for the whole molecule.
    else {
      FOR_ATOMS_OF_MOL(a, mol)
        {
          fragatoms.SetBitOn(a->GetIdx());
        }
    }

    if (mol.NumAtoms() > 0) {
      CreateCansmiString(mol, buffer, fragatoms, !pConv->IsOption("i"), pConv);
    }

    ofs << buffer;
    if(!pConv->IsOption("smilesonly")) {

      if(!pConv->IsOption("n"))
        ofs << '\t' <<  mol.GetTitle();

      if (pConv->IsOption("x") && mol.HasData("Canonical Atom Order")) {
        vector<string> vs;
        string canorder = mol.GetData("Canonical Atom Order")->GetValue();
        tokenize(vs, canorder);
        ofs << '\t';
        for (int i = 0; i < vs.size(); i++) {
          int idx = atoi(vs[i].c_str());
          OBAtom *atom = mol.GetAtom(idx);
          if (i > 0)
            ofs << ",";
          ofs << atom->GetX() << "," << atom->GetY();
        }
      }

      if(!pConv->IsOption("nonewline"))
        ofs << endl;
    }

    return true;
  }
  
  //********************************************************
  class FIXFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FIXFormat()
    {
      OBConversion::RegisterFormat("fix",this);
    }

    virtual const char* Description() //required
    {
      return
        "SMILES FIX format\n"
        "  No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  FIXFormat theFIXFormat;

  /////////////////////////////////////////////////////////////////

  bool FIXFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    OBMol2Cansmi m2s;

    // This is a hack to prevent recursion problems.
    //  we still need to fix the underlying problem -GRH
    if (mol.NumAtoms() > 1000)
      {
        stringstream errorMsg;
        errorMsg << "SMILES Conversion failed: Molecule is too large to convert. Open Babel is currently limited to 1000 atoms." << endl;
        errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        return(false);
      }

    // Write the SMILES in a FIX with canonical order
    m2s.Init(true, pConv);
    // From 2.1 code.
    m2s.CorrectAromaticAmineCharge(mol);

    // We're outputting a full molecule
    // so we pass a bitvec for all atoms
    OBBitVec allbits(mol.NumAtoms());
    FOR_ATOMS_OF_MOL(a, mol) {
      allbits.SetBitOn(a->GetIdx());
    }

    if (mol.NumAtoms() > 0) {
      CreateCansmiString(mol, buffer, allbits, !pConv->IsOption("i"), pConv);
    }
    ofs << buffer << endl;

    OBAtom *atom;
    vector<int>::iterator i;
    // Retrieve the canonical order of the molecule
    string orderString = m2s.GetOutputOrder();
    vector<string> canonical_order;
    tokenize(canonical_order, orderString);

    int j;
    int atomIdx;
    for (j = 0;j < mol.NumConformers();j++)
      {
        mol.SetConformer(j);
        for (unsigned int index = 0; index < canonical_order.size(); 
             ++index) {
          atomIdx = atoi(canonical_order[index].c_str());
          atom = mol.GetAtom(atomIdx);
          sprintf(buffer,"%9.3f %9.3f %9.3f",atom->GetX(),atom->GetY(),atom->GetZ());
          ofs << buffer<< endl;
        }
      }
    return(true);
  }

} // end namespace OpenBabel
