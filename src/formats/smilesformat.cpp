/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
//Contains SMIFormat and FIXFormat classes

#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"

using namespace std;

namespace OpenBabel
{
class SMIFormat : public OBMoleculeFormat
{
public:
	//Register this format type ID
	SMIFormat()
	{
		OBConversion::RegisterFormat("smi",this, "chemical/x-daylight-smiles");
		OBConversion::RegisterOptionParam("n", this);
		OBConversion::RegisterOptionParam("t", this);
	}

  virtual const char* GetMIMEType() 
  { return "chemical/x-daylight-smiles"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ///////////////////////////////////////////////////////

    virtual const char* Description()
    {
        return
            "SMILES format\n \
            A linear text format which can describe the connectivity\n \
            and chirality of a molecule\n \
            Write Options e.g. -xt\n \
            -n no molecule name\n \
            -t molecule name only\n \
	    -r radicals lower case eg ethyl is Cc\n\n";
    };

  virtual const char* SpecificationURL()
  {return "http://www.daylight.com/smiles/f_smiles.html";};

  virtual int SkipObjects(int n, OBConversion* pConv)
  {
    if(n==0) return 1; //already points after current line
    string temp;
    istream& ifs = *pConv->GetInStream();
    int i;
    for(i=0;i<n && ifs.good();i++)
      getline(ifs, temp);
    return ifs.good() ? 1 : -1;	
  };	
};

//Make an instance of the format class
SMIFormat theSMIFormat;

//////////////////////////////////////////////////////////////////
class OBSmiNode
{
    OBAtom *_atom,*_parent;
    std::vector<OBSmiNode*> _nextnode;
    std::vector<OBBond*> _nextbond;
public:
    OBSmiNode(OBAtom *atom);
    ~OBSmiNode();
    int        Size()
    {
        return((_nextnode.empty())?0:_nextnode.size());
    }
    void       SetParent(OBAtom *a)
    {
        _parent = a;
    }
    void       SetNextNode(OBSmiNode*,OBBond*);
    OBAtom    *GetAtom()
    {
        return(_atom);
    }
    OBAtom    *GetParent()
    {
        return(_parent);
    }
    OBAtom    *GetNextAtom(int i)
    {
        return(_nextnode[i]->GetAtom());
    }
    OBBond    *GetNextBond(int i)
    {
        return(_nextbond[i]);
    }
    OBSmiNode *GetNextNode(int i)
    {
        return(_nextnode[i]);
    }
};

class OBMol2Smi
{
    std::vector<int> _atmorder;
    std::vector<int> _storder;
    std::vector<bool> _aromNH;
    OBBitVec _uatoms,_ubonds;
    std::vector<OBEdgeBase*> _vclose;
    std::vector<std::pair<OBAtom*,std::pair<int,int> > > _vopen;
		OBConversion* _pconv;
public:
    OBMol2Smi()
    {
        _vclose.clear();
    }
    ~OBMol2Smi()
    {}
    int          GetUnusedIndex();
    void         Init(OBConversion* pconv=NULL);
    void         CreateSmiString(OBMol&,char*);
    void         GetClosureAtoms(OBAtom*,std::vector<OBNodeBase*>&);
    void         FindClosureBonds(OBMol&);
    void         ToSmilesString(OBSmiNode *node,char *buffer);
    void         RemoveUsedClosures();
    void         AssignCisTrans(OBSmiNode*);
    bool         BuildTree(OBSmiNode*);
    bool         GetSmilesElement(OBSmiNode*,char*);
    bool         GetChiralStereo(OBSmiNode*,char*);
    void         CorrectAromaticAmineCharge(OBMol&);
    std::vector<std::pair<int,OBBond*> >  GetClosureDigits(OBAtom*);
    std::vector<int> &GetOutputOrder()
    {
        return(_atmorder);
    }
};

bool WriteTheSmiles(OBMol & mol,char *out);

/////////////////////////////////////////////////////////////////
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
    bool chiralWatch; // set when a chiral atom is read
    map<OBAtom*,OBChiralData*> _mapcd; // map of ChiralAtoms and their data
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
};

/////////////////////////////////////////////////////////////////
bool SMIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    //Taken unchanged from ReadSmiles
    char buffer[BUFF_SIZE];

    if (!ifs.getline(buffer,BUFF_SIZE))
        return(false);
    vector<string> vs;
    tokenize(vs,buffer);

    // RWT 10/3/2000
    //
    // added the following to allow spaces in compound names (titles).
    // Essentially everything after the first space on a SMILES file line
    // is treated as the name.
    // Also had to change the condition a few lines below from:
    // if (vs.size() == 2) ...  to
    // if (vs.size() >= 2)

    if (vs.size() > 2)
    {
        for (unsigned int i=2;i<vs.size(); i++)
        {
            vs[1]=vs[1]+" "+vs[i];
        }
    }

    if (vs.empty())
      return false;
    mol.SetDimension(0);

    if (vs.size() >= 2)
      mol.SetTitle(vs[1].c_str());

    OBSmilesParser sp;
    return sp.SmiToMol(mol,vs[0]);
}

//////////////////////////////////////////////////
bool SMIFormat::WriteMolecule(OBBase* pOb,OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    if(pConv->IsOption("t")) //Title only option
    {
        ofs << mol.GetTitle() <<endl;
        return true;
    }
    char buffer[BUFF_SIZE];
		*buffer='\0'; //empty buffer

    // This is a hack to prevent recursion problems.
    //  we still need to fix the underlying problem (mainly chiral centers) -GRH
    if (mol.NumAtoms() > 1000)
      {
#ifdef HAVE_SSTREAM
	stringstream errorMsg;
#else
	strstream errorMsg;
#endif
	errorMsg << "SMILES Conversion failed: Molecule is too large to convert. Open Babel is currently limited to 1000 atoms." << endl;
	errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
	obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
	return(false);
      }

    if(mol.NumAtoms()!=0)
    {
      OBMol2Smi m2s;
      m2s.Init(pConv);
      m2s.CorrectAromaticAmineCharge(mol);
      m2s.CreateSmiString(mol,buffer);
    }

    ofs << buffer ;
    if(!pConv->IsOption("n"))
        ofs << '\t' <<  mol.GetTitle();
    ofs << endl;

    return true;
}

//////////////////////////////////////////////

//CM not needed extern OBAromaticTyper  aromtyper;
//OBAtomTyper atomtyperx; //CM

bool OBSmilesParser::SmiToMol(OBMol &mol,string &s)
{
    strcpy(_buffer,s.c_str());

    _vprev.clear();
    _rclose.clear();
    _prev=0;
		chiralWatch=false;

    if (!ParseSmiles(mol))
    {
        mol.Clear();
        return(false);
    }

//    mol.AddHydrogens(); // need to add implicit hydrogens

    return(true);
}

bool OBSmilesParser::ParseSmiles(OBMol &mol)
{
    mol.BeginModify();
    
    for (_ptr=_buffer;*_ptr;_ptr++)
    {
        if (isspace(*_ptr))
            continue;
        else if (isdigit(*_ptr) || *_ptr == '%') //ring open/close
        {
            ParseRingBond(mol);
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
                _bondflags |= OB_TORDOWN_BOND;
                break;
            case '\\':
                _bondflags |= OB_TORUP_BOND;
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

    //set aromatic bond orders
    mol.SetAromaticPerceived();
    FindAromaticBonds(mol);
    FindOrphanAromaticAtoms(mol);// CM 18 Sept 2003
    mol.AssignSpinMultiplicity();
    mol.UnsetAromaticPerceived();

    mol.EndModify();
    
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

// CM 18 Sept 2003
void OBSmilesParser::FindOrphanAromaticAtoms(OBMol &mol)
{
    //Facilitates the use lower case shorthand for radical entry
    //Atoms which are marked as aromatic but have no aromatic bonds
    //are taken to be radical centres
    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;

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
    vector<OBEdgeBase*>::iterator i;
    for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
        if (!bond->GetBeginAtom()->IsAromatic() ||
                !bond->GetEndAtom()->IsAromatic())
            _bvisit[bond->GetIdx()] = true;

    OBAtom *atom;
    vector<OBNodeBase*>::iterator j;

    for (atom = mol.BeginAtom(j);atom;atom = mol.NextAtom(j))
        if(!_avisit[atom->GetIdx()] && atom->IsAromatic())
            FindAromaticBonds(mol,atom,0);
}

void OBSmilesParser::FindAromaticBonds(OBMol &mol,OBAtom *atom,int depth )
{
    OBBond *bond;
    vector<OBEdgeBase*>::iterator k;

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
        connected by a double bond. 
        Since they are still marked as aromatic, FindAromaticBonds() will
        replace the bonds by aromatic bonds if they are in a ring.
        FindOrphanAromand removes the aromatic tag from the atoms not found in this way
        and removes stray radical centres in .

        To avoid difficulties in complex aromatics with 5 membered rings containing N and O,
        the above scheme modified to prevent recognition of aromatic structures is not confused.
        - the double bond is made single if it would exceed valence of atom (not all aromatics have conjugated bonds)
        - the radical centre is removed on both atoms when forming a ring (in ParseRingBond())
           and on the new atom if the valence of the prev atom is being exceeded.
        Note that the bond orders made here are reassigned in aromatic structures in FindAromaticBonds()
        */
        if(arom)
        {
            OBAtom* prevatom = mol.GetAtom(_prev);

            //Calculate available valency on prevatom
            //This is far more difficult than it should be!
            //Data not always updated during molecule constuction.
            int val=0;
            if(prevatom->IsCarbon())
                val=4;
            else if(prevatom->IsNitrogen())
                val=3;
            else if(prevatom->IsPhosphorus())
                val=3;
            else if(prevatom->IsOxygen())
                val=2;
            else if(prevatom->IsSulfur())
                val=2;

            /*			int sumBO=0;
            			vector<OBEdgeBase*>::iterator itr;
            			OBBond* bond = prevatom->BeginBond(itr);
            			while(bond)
            			{
            				sumBO +=bond->GetBO();
            				bond=prevatom->NextBond(itr);
            			}
            			if(prevatom->BOSum() != sumBO)
            				cerr << "BOSum != sumBO" << endl;
            */
            int AvailableValence = val + prevatom->GetFormalCharge() - prevatom->BOSum();//sumBO;

            if (prevatom->GetSpinMultiplicity())
            {
                prevatom->SetSpinMultiplicity(0);
                atom->SetSpinMultiplicity(0);

                //Make the new bond double unless not allowed by prevatom's valence
                _order = AvailableValence>=2  ?  2 : 1 ;
            }
            else
                if(AvailableValence<1) //Must be complex aromatic with O, N
                    atom->SetSpinMultiplicity(0); //radical centres not appropriate in complex aromatics
        }
        // CM end
        mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);
        
        //NE iterate through and see if atom is bonded to chiral atom
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        ChiralSearch=_mapcd.find(mol.GetAtom(_prev));
        if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
        {
	    (ChiralSearch->second)->AddAtomRef(mol.NumAtoms(), input);
          // cout << "Line 650: Adding "<<mol.NumAtoms()<<" to "<<ChiralSearch->second<<endl;
        }
    }
    //set values
    _prev = mol.NumAtoms();
    _order = 1;
    _bondflags = 0;

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
    for (;*_ptr && isdigit(*_ptr);_ptr++)
    {
        symbol[isoPtr] = *_ptr;
        isoPtr++;
    }
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
        case 's':
            _ptr++;
            if (*_ptr == 'e')
            {
                element = 34;
                strcpy(symbol,"Se");
            }
            else
            {
                element = 16;
                symbol[0] = 'S';
                _ptr--;
            }
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
        default:
            return(false);
        }
    }

    //handle hydrogen count, stereochemistry, and charge

    OBAtom *atom = mol.NewAtom();
    int hcount = 0;
    int charge=0;
		int rad=0;
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

        default:
            return(false);
        }
    }

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
       mol.AddBond(_prev,mol.NumAtoms(),_order,_bondflags);
        if(chiralWatch) // if chiral atom need to add its previous into atom4ref
        {
	  if (_mapcd[atom] == NULL)
            _mapcd[atom]= new OBChiralData;
	  
	  (_mapcd[atom])->AddAtomRef((unsigned int)_prev,input);
           // cout <<"line 1405: Added atom ref "<<_prev<<" to "<<_mapcd[atom]<<endl;
        }
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        ChiralSearch = _mapcd.find(mol.GetAtom(_prev));
        if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
        {
           (ChiralSearch->second)->AddAtomRef(mol.NumAtoms(), input);
          // cout <<"line 1431: Added atom ref "<<mol.NumAtoms()<<" to "<<ChiralSearch->second<<endl;
        }
       }          

    //set values
    _prev = mol.NumAtoms();
    _order = 1;
    _bondflags = 0;

    //now add hydrogens

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
          // cout << "line 1434: Added atom ref "<<mol.NumAtoms()<<" to "<<_mapcd[mol.GetAtom(_prev)]<<endl;
                       
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
            mol.SetData(xbd);
        }
        xbd->SetData(atom,refbond,(*bond)[0]);

        /* old code written by AGS -- mts
        {
        	externalbonds = (vector<pair<int,pair<OBAtom *,OBBond *> > > *)mol.GetData("extBonds");
        }
        else
        {
        	externalbonds = new vector<pair<int,pair<OBAtom *,OBBond *> > >;
        }

        //save data <external bond count, bond index>
        externalbonds->push_back(pair<int,pair<OBAtom *,OBBond *> > ((*bond)[0], pair<OBAtom *,OBBond *> (atom,mol.GetBond((*bond)[1],atom->GetIdx()))));
        mol.SetData("extBonds",externalbonds);
        */

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
    case '/':	//chiral, but _order still == 1
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
    digit = atoi(str);	// convert indicator to digit

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
                // cout << "Added external "<<(*j)[1]<<" to "<<ChiralSearch->second<<endl;
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
            mol.AddBond((*j)[1],_prev,ord,bf,(*j)[4]);
            
            // after adding a bond to atom "_prev"
            // search to see if atom is bonded to a chiral atom
            // need to check both _prev and (*j)[1] as closure is direction independant
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch,cs2;
        ChiralSearch = _mapcd.find(mol.GetAtom(_prev));
        cs2=_mapcd.find(mol.GetAtom((*j)[1]));
        if (ChiralSearch!=_mapcd.end() && ChiralSearch->second != NULL)
        {
	    (ChiralSearch->second)->AddAtomRef((*j)[1], input);
           //cout << "Added ring closure "<<(*j)[1]<<" to "<<ChiralSearch->second<<endl;
        }
        if (cs2!=_mapcd.end() && cs2->second != NULL)
        {
           (cs2->second)->AddAtomRef(_prev,input);
          // cout <<"Added ring opening "<<_prev<<" to "<<cs2->second<<endl;
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
    vtmp[4] = mol.GetAtom(_prev)->GetValence(); //store position to insert closure bond
    for (j = _rclose.begin();j != _rclose.end();j++) //correct for multiple closure bonds to a single atom
        if ((*j)[1] == _prev)
            vtmp[4]++;

    _rclose.push_back(vtmp);
    _order = 1;
    _bondflags = 0;

    return(true);
}


void OBMol2Smi::CreateSmiString(OBMol &mol,char *buffer)
{
    OBAtom *atom;
    OBSmiNode *root =NULL;
    buffer[0] = '\0';
    vector<OBNodeBase*>::iterator i;

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        //    if ((!atom->IsHydrogen() || atom->GetValence() == 0) && !_uatoms[atom->GetIdx()])
        if (!atom->IsHydrogen() && !_uatoms[atom->GetIdx()])
	  if (!atom->IsChiral()) //don't use chiral atoms as root node
	    {
	      //clear out closures in case structure is dot disconnected
	      _vclose.clear();
	      _atmorder.clear();
	      _storder.clear();
	      _vopen.clear();
	      //dot disconnected structure
	      if (strlen(buffer) > 0)
		strcat(buffer,".");
	      root = new OBSmiNode (atom);
	      BuildTree(root);
	      FindClosureBonds(mol);
	      if (mol.Has2D())
		AssignCisTrans(root);
	      ToSmilesString(root,buffer);
	      delete root;
	    }
    
    //If no starting node found e.g. [H][H] CM 21Mar05
    if(root==NULL)
      {
	root = new OBSmiNode(mol.GetFirstAtom());
	BuildTree(root);
	ToSmilesString(root,buffer);
	delete root;
      }
}

bool OBMol2Smi::BuildTree(OBSmiNode *node)
{
    vector<OBEdgeBase*>::iterator i;
    OBAtom *nbr,*atom = node->GetAtom();

    _uatoms.SetBitOn(atom->GetIdx()); //mark the atom as visited
    _atmorder.push_back(atom->GetIdx()); //store the atom ordering
    _storder.push_back(atom->GetIdx()); //store the atom ordering for stereo

    for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
        if ((!nbr->IsHydrogen()||nbr->GetIsotope()||atom->IsHydrogen()) //so D,T is explicit and H2 works!CM 21Mar05
					&& !_uatoms[nbr->GetIdx()]) 
        {
            _ubonds.SetBitOn((*i)->GetIdx());
            OBSmiNode *next = new OBSmiNode (nbr);
            next->SetParent(atom);
            node->SetNextNode(next,(OBBond*)*i);
            BuildTree(next);
        }

    return(true);
}

void OBMol2Smi::ToSmilesString(OBSmiNode *node,char *buffer)
{
    char tmpbuf[10];
    OBAtom *atom = node->GetAtom();

    //write the current atom to the string
    GetSmilesElement(node,tmpbuf);
    strcat(buffer,tmpbuf);

    //handle ring closures here
    vector<pair<int,OBBond*> > vc = GetClosureDigits(atom);
    if (!vc.empty())
    {
        vector<pair<int,OBBond*> >::iterator i;
        for (i = vc.begin();i != vc.end();i++)
        {
            if (i->second)
            {
                if (i->second->IsUp())
		  {
		    if ( (i->second->GetBeginAtom())->HasDoubleBond() ||
			 (i->second->GetEndAtom())->HasDoubleBond() )
		      strcat(buffer,"\\");
		  }
                if (i->second->IsDown())
		  {
		    if ( (i->second->GetBeginAtom())->HasDoubleBond() ||
			 (i->second->GetEndAtom())->HasDoubleBond() )
		      strcat(buffer,"/");
		  }
#ifndef KEKULE

                if (i->second->GetBO() == 2 && !i->second->IsAromatic())
                    strcat(buffer,"=");
#else

                if (i->second->GetBO() == 2)
                    strcat(buffer,"=");
#endif

                if (i->second->GetBO() == 3)
                    strcat(buffer,"#");
            }

            if (i->first > 9)
                strcat(buffer,"%");
            sprintf(tmpbuf,"%d",i->first);
            strcat(buffer,tmpbuf);
        }
    }

    //follow path to child atoms
    OBBond *bond;
    for (int i = 0;i < node->Size();i++)
    {
        bond = node->GetNextBond(i);
        if (i+1 < node->Size())
            strcat(buffer,"(");
        if (bond->IsUp())
	  {
	    if ( (bond->GetBeginAtom())->HasDoubleBond() ||
		 (bond->GetEndAtom())->HasDoubleBond() )
	      strcat(buffer,"\\");
	  }
        if (bond->IsDown())
	  {
	    if ( (bond->GetBeginAtom())->HasDoubleBond() ||
		 (bond->GetEndAtom())->HasDoubleBond() )
	      strcat(buffer,"/");
	  }
#ifndef KEKULE

        if (bond->GetBO() == 2 && !bond->IsAromatic())
            strcat(buffer,"=");
#else

        if (bond->GetBO() == 2)
            strcat(buffer,"=");
#endif

        if (bond->GetBO() == 3)
            strcat(buffer,"#");

        ToSmilesString(node->GetNextNode(i),buffer);
        if (i+1 < node->Size())
            strcat(buffer,")");
    }
}

void OBMol2Smi::GetClosureAtoms(OBAtom *atom,vector<OBNodeBase*> &va)
{

    //look through closure list for start atom
    vector<OBEdgeBase*>::iterator i;
    for (i = _vclose.begin();i != _vclose.end();i++)
        if (*i)
        {
            if (((OBBond*)*i)->GetBeginAtom() == atom)
                va.push_back(((OBBond*)*i)->GetEndAtom());
            if (((OBBond*)*i)->GetEndAtom() == atom)
                va.push_back(((OBBond*)*i)->GetBeginAtom());
        }

    OBAtom *nbr;
    vector<pair<OBAtom*,pair<int,int> > >::iterator j;
    for (j = _vopen.begin();j != _vopen.end();j++)
        for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
            if (nbr == j->first)
                va.push_back(nbr);
}

vector<pair<int,OBBond*> > OBMol2Smi::GetClosureDigits(OBAtom *atom)
{
    vector<pair<int,OBBond*> > vc;
    vc.clear();

    //look through closure list for start atom
    int idx,bo;
    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    for (i = _vclose.begin();i != _vclose.end();i++)
        if ((bond=(OBBond*)*i))
            if (bond->GetBeginAtom() == atom || bond->GetEndAtom() == atom)
            {
                idx = GetUnusedIndex();
                vc.push_back(pair<int,OBBond*> (idx,bond));
                bo = (bond->IsAromatic())? 1 : bond->GetBO();
                _vopen.push_back(pair<OBAtom*,pair<int,int> >
                                 (bond->GetNbrAtom(atom),pair<int,int>(idx,bo)));
                *i = NULL;//remove bond from closure list
            }

    //try to complete closures
    if (!_vopen.empty())
    {
        vector<pair<OBAtom*,pair<int,int> > >::iterator j;
        for (j = _vopen.begin();j != _vopen.end();)
            if (j->first == atom)
            {
                vc.push_back(pair<int,OBBond*> (j->second.first,(OBBond*)NULL));
                _vopen.erase(j);
                j = _vopen.begin();
            }
            else
                j++;
    }

    return(vc);
}

void OBMol2Smi::FindClosureBonds(OBMol &mol)
{
    //find closure bonds
    OBAtom *a1,*a2;
    OBBond *bond;
    vector<OBEdgeBase*>::iterator i;
    OBBitVec bv;
    bv.FromVecInt(_storder);

    for (bond = mol.BeginBond(i);bond;bond = mol.NextBond(i))
        if (!_ubonds[bond->GetIdx()] && bv[bond->GetBeginAtomIdx()])
        {
            a1 = bond->GetBeginAtom();
            a2 = bond->GetEndAtom();
            if (!a1->IsHydrogen() && !a2->IsHydrogen())
                _vclose.push_back(bond);
        }

    vector<OBEdgeBase*>::reverse_iterator j;
    vector<int>::iterator k;

    //modify _order to reflect ring closures
    for (j = _vclose.rbegin();j != _vclose.rend();j++)
    {
        bond = (OBBond*)*j;
        a1 = a2 = NULL;

        for (k = _storder.begin();k != _storder.end();k++)
            if (bond->GetBeginAtomIdx() == static_cast<unsigned int>(*k) ||
                    bond->
                    GetEndAtomIdx() == static_cast<unsigned int>(*k))
                if (!a1) a1 = mol.GetAtom(*k);
                else if (!a2)
                {
                    a2 = mol.GetAtom(*k)
                         ;
                    _storder.erase(k);
                    break;
                }

        for (k = _storder.begin()
                 ;
                k != _storder.end();
                k++)
            if (a1->GetIdx()
                    == static_cast<unsigned int>(*k))
            {
                k++;
                if (k != _storder.end())
                    _storder.insert(k,a2->GetIdx());
                else
                    _storder.push_back(a2->GetIdx());
                break;
            }
    }
}

int OBMol2Smi::GetUnusedIndex()
{
    int idx=1;

    vector<pair<OBAtom*,pair<int,int> > >::iterator j;
    for (j = _vopen.begin();j != _vopen.end();)
        if (j->second.first == idx)
        {
            idx++; //increment idx and start over if digit is already used
            j = _vopen.begin();
        }
        else
            j++;

    return(idx);
}

void OBMol2Smi::CorrectAromaticAmineCharge(OBMol &mol)
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

void OBMol2Smi::AssignCisTrans(OBSmiNode *node)
{
    //traverse the tree searching for acyclic olefins - if it
    //has at least one heavy atom attachment on each end assign stereochem

    OBBond *bond;
    for (int i = 0;i < node->Size();i++)
    {
        bond = node->GetNextBond(i);
        if (bond->GetBO() == 2 && !bond->IsInRing())
        {
            OBAtom *b = node->GetAtom();
            OBAtom *c = bond->GetNbrAtom(b);

            //skip allenes
            if (b->GetHyb() == 1 || c->GetHyb() == 1)
                continue;

            if (b->GetHvyValence() > 1 && c->GetHvyValence() > 1)
            {
                OBAtom *a,*d;
                vector<OBEdgeBase*>::iterator j,k;

                //look for bond with assigned stereo as in poly-ene
                for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
                    if (((OBBond*)*j)->IsUp() ||((OBBond*)*j)->IsDown())
                        break;

                if (!a)
                    for (a = b->BeginNbrAtom(j);a;a = b->NextNbrAtom(j))
                        if (a != c && !a->IsHydrogen())
                            break;
                for (d = c->BeginNbrAtom(k);d;d = c->NextNbrAtom(k))
                    if (d != b && !d->IsHydrogen())
                        break;
		//                obAssert(a);
		//                obAssert(d);

                if (((OBBond*)*j)->IsUp() || ((OBBond*)*j)->IsDown()) //stereo already assigned
                {
                    if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(),
                                              c->GetVector(),d->GetVector())) > 10.0)
                        if (((OBBond*)*j)->IsUp())
                            ((OBBond*)*k)->SetUp();
                        else
                            ((OBBond*)*k)->SetDown();
                    else
                        if (((OBBond*)*j)->IsUp())
                            ((OBBond*)*k)->SetDown();
                        else
                            ((OBBond*)*k)->SetUp();
                }
                else //assign stereo to both ends
                {
                    ((OBBond*)*j)->SetUp();
                    if (fabs(CalcTorsionAngle(a->GetVector(),b->GetVector(),
                                              c->GetVector(),d->GetVector())) > 10.0)
                        ((OBBond*)*k)->SetUp();
                    else
                        ((OBBond*)*k)->SetDown();
                }
            }
        }
        AssignCisTrans(node->GetNextNode(i));
    }
}

void OBMol2Smi::Init(OBConversion* pconv)
{
    _vclose.clear();
    _atmorder.clear();
    _storder.clear();
    _aromNH.clear();
    _uatoms.Clear();
    _ubonds.Clear();
    _vopen.clear();
		_pconv = pconv;
}

bool OBMol2Smi::GetSmilesElement(OBSmiNode *node,char *element)
{
    //***handle reference atom stuff here and return***
    char symbol[10];
    bool bracketElement = false;
    bool normalValence = true;

    OBAtom *atom = node->GetAtom();

    int bosum = atom->KBOSum();
    atom->BOSum(); //CM temp
    switch (atom->GetAtomicNum())
    {
    case 0:
        break;
    case 5: /*bracketElement = !(normalValence = (bosum == 3)); break;*/
        break;
    case 6:
        break;
    case 7:
        if (atom->IsAromatic() && atom->GetHvyValence() == 2 && atom->GetImplicitValence() == 3)
        {
            bracketElement = !(normalValence = false);
            break;
        }
        else
            bracketElement = !(normalValence = (bosum == 3 || bosum == 5));
        break;
    case 8:
        break;
    case 9:
        break;
    case 15:
        break;
    case 16:
        bracketElement = !(normalValence = (bosum == 2 || bosum == 4 || bosum == 6));
        break;
    case 17:
        break;
    case 35:
        break;
    case 53:
        break;

    default:
        bracketElement = true;
    }

    if (atom->GetHvyValence() > 2 && atom->IsChiral())
        if (((OBMol*)atom->GetParent())->HasNonZeroCoords() || atom->HasChiralitySpecified())
            bracketElement = true;

    if (atom->GetFormalCharge() != 0) //bracket charged elements
        bracketElement = true;

    if(atom->GetIsotope()) //CM 19Mar05
      bracketElement = true;

    //CM begin 18 Sept 2003

    //This outputs form [CH3][CH3] rather than CC if -h option has been specified
    if (((OBMol*)atom->GetParent())->HasHydrogensAdded())
        bracketElement = true;
		else
		{
			if (atom->GetSpinMultiplicity())
			{
				//For radicals output bracket form anyway unless r option specified
				if(!(_pconv && _pconv->IsOption ("r")))
					bracketElement = true;
			}
		}
    //CM end

    if (!bracketElement)
    {
        if (!atom->GetAtomicNum())
        {
            bool external = false;
            vector<pair<int,pair<OBAtom *,OBBond *> > > *externalBonds =
                (vector<pair<int,pair<OBAtom *,OBBond *> > > *)((OBMol*)atom->GetParent())->GetData("extBonds");
            vector<pair<int,pair<OBAtom *,OBBond *> > >::iterator externalBond;

            if (externalBonds)
                for(externalBond = externalBonds->begin();externalBond != externalBonds->end();externalBond++)
                {
                    if (externalBond->second.first == atom)
                    {
                        external = true;
                        strcpy(symbol,"&");
                        OBBond *bond = externalBond->second.second;
                        if (bond->IsUp())
			  {
			    if ( (bond->GetBeginAtom())->HasDoubleBond() ||
				 (bond->GetEndAtom())->HasDoubleBond() )
			      strcat(symbol,"\\");
			  }
                        if (bond->IsDown())
			  {
			    if ( (bond->GetBeginAtom())->HasDoubleBond() ||
				 (bond->GetEndAtom())->HasDoubleBond() )
			      strcat(symbol,"/");
			  }
#ifndef KEKULE

                        if (bond->GetBO() == 2 && !bond->IsAromatic())
                            strcat(symbol,"=");
                        if (bond->GetBO() == 2 && bond->IsAromatic())
                            strcat(symbol,";");
#else

                        if (bond->GetBO() == 2)
                            strcat(symbol,"=");
#endif

                        if (bond->GetBO() == 3)
                            strcat(symbol,"#");
                        sprintf(symbol,"%s%d",symbol,externalBond->first);
                        break;
                    }
                }

            if(!external)
                strcpy(symbol,"*");
        }
        else
        {
            strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
#ifndef KEKULE

            if (atom->IsAromatic())
                symbol[0] = tolower(symbol[0]);
#endif

					//Radical centres lc if r option set
				    if(atom->GetSpinMultiplicity() && _pconv && _pconv->IsOption ("r"))
				        symbol[0] = tolower(symbol[0]);
        }
        strcpy(element,symbol);

        return(true);
    }

    strcpy(element,"[");
		if(atom->GetIsotope()) //CM 19Mar05
		{
			char iso[4];
			sprintf(iso,"%d",atom->GetIsotope());
			strcat(element,iso);
		}
    if (!atom->GetAtomicNum())
        strcpy(symbol,"*");
    else
    {
        strcpy(symbol,etab.GetSymbol(atom->GetAtomicNum()));
#ifndef KEKULE

        if (atom->IsAromatic())
            symbol[0] = tolower(symbol[0]);
#endif

    }
    strcat(element,symbol);

    //if (atom->IsCarbon() && atom->GetHvyValence() > 2 && atom->IsChiral())
    if (atom->GetHvyValence() > 2 && atom->IsChiral())
    {
        char stereo[5];
        if (GetChiralStereo(node,stereo))
            strcat(element,stereo);
    }

    //add extra hydrogens
    //  if (!normalValence && atom->ImplicitHydrogenCount())
    if (atom->ImplicitHydrogenCount() && !atom->IsHydrogen()) //CM 21Mar05
    {
        strcat(element,"H");
        if (atom->ImplicitHydrogenCount() > 1)
        {
            char tcount[10];
            sprintf(tcount,"%d",atom->ImplicitHydrogenCount());
            strcat(element,tcount);
        }
    }

    //cat charge on the end
    if (atom->GetFormalCharge() != 0)
    {

        /*
        		if (atom->ImplicitHydrogenCount())
        		{
        			cerr << "imp = " << atom->GetAtomicNum() << ' ' << atom->GetImplicitValence() << endl;
              strcat(element,"H");
              if (atom->ImplicitHydrogenCount() > 1)
        			{
        				char tcount[10];
        				sprintf(tcount,"%d",atom->ImplicitHydrogenCount());
        				strcat(element,tcount);
        			}
        		}
        	*/
        if (atom->GetFormalCharge() > 0)
            strcat(element,"+");
        else
            strcat(element,"-");

        if (abs(atom->GetFormalCharge()) > 1)
        {
            char tcharge[10];
            sprintf(tcharge,"%d",abs(atom->GetFormalCharge()));
            strcat(element,tcharge);
        }
    }

    strcat(element,"]");

    return(true);
}

bool OBMol2Smi::GetChiralStereo(OBSmiNode *node,char *stereo)
{
    bool is2D=false;
    double torsion;
    OBAtom *a,*b,*c,*d,hydrogen;

    b = node->GetAtom();
    OBMol *mol = (OBMol*)b->GetParent();

    if (!mol->HasNonZeroCoords()) //must have come in from smiles string
    {
        if (!b->HasChiralitySpecified())
            return(false);
        if (b->IsClockwise())
            strcpy(stereo,"@@");
        else if (b->IsAntiClockwise())
            strcpy(stereo,"@");
        else
            return(false);
        //if (b->GetHvyValence() == 3) strcat(stereo,"H");
        return(true);
    }

    //give peudo Z coords if mol is 2D
    if (!mol->Has3D())
    {
        vector3 v,vz(0.0,0.0,1.0);
        is2D = true;
        OBAtom *nbr;
        OBBond *bond;
        vector<OBEdgeBase*>::iterator i;
        for (bond = b->BeginBond(i);bond;bond = b->NextBond(i))
        {
            nbr = bond->GetEndAtom();
            if (nbr != b)
            {
                v = nbr->GetVector();
                if (bond->IsWedge())
                    v += vz;
                else
                    if (bond->IsHash())
                        v -= vz;

                nbr->SetVector(v);
            }
            else
            {
                nbr = bond->GetBeginAtom();
                v = nbr->GetVector();
                if (bond->IsWedge())
                    v -= vz;
                else
                    if (bond->IsHash())
                        v += vz;

                nbr->SetVector(v);
            }
        }
    }

    c = d = NULL;
    a = node->GetParent();
    //    obAssert(a); //chiral atom can't be used as root node - must have parent

    if (b->GetHvyValence() == 3) //must have attached hydrogen
    {
        if (b->GetValence() == 4)//has explicit hydrogen
        {
            vector<OBEdgeBase*>::iterator i;
            for (c = b->BeginNbrAtom(i);c;c = b->NextNbrAtom(i))
                if (c->IsHydrogen())
                    break;
	    //            obAssert(c);
        }
        else  //implicit hydrogen
        {
            vector3 v;
            b->GetNewBondVector(v,1.0);
            hydrogen.SetVector(v);
            c = &hydrogen;
        }
    }

    //get connected atoms in order
    OBAtom *nbr;
    vector<int>::iterator j;

    //try to get neighbors that are closure atoms in the order they appear in the string
    vector<OBNodeBase*> va;
    GetClosureAtoms(b,va);
    if (!va.empty())
    {
        vector<OBNodeBase*>::iterator k;
        for (k = va.begin();k != va.end();k++)
            if (*k != a)
            {
                if (!c)
                    c = (OBAtom*)*k;
                else if (!d)
                    d = (OBAtom*)*k;
            }
    }

    for (j = _storder.begin();j != _storder.end();j++)
    {
        nbr = mol->GetAtom(*j);
        if (!b->IsConnected(nbr))
            continue;
        if (nbr == a || nbr == b || nbr == c)
            continue;
        if (!c)
            c = nbr;
        else if (!d)
            d = nbr;
    }

    torsion = CalcTorsionAngle(a->GetVector(),b->GetVector(),
                               c->GetVector(),d->GetVector());

    strcpy(stereo,(torsion<0.0)?"@":"@@");
    //if (b->GetHvyValence() == 3) strcat(stereo,"H");

    //re-zero psuedo-coords
    if (is2D)
    {
        vector3 v;
        OBAtom *atom;
        vector<OBNodeBase*>::iterator k;
        for (atom = mol->BeginAtom(k);atom;atom = mol->NextAtom(k))
        {
            v = atom->GetVector();
            v.SetZ(0.0);
            atom->SetVector(v);
        }
    }

    return(true);
}
//********************************************************
class FIXFormat : public OBFormat
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
            "SMILES FIX format\n \
            No comments yet\n \
            ";
    };

    virtual const char* SpecificationURL(){return
            "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
    /// The "Convert" interface functions
    virtual bool WriteChemObject(OBConversion* pConv)
    {
        //Retrieve the target OBMol
        OBBase* pOb = pConv->GetChemObject();
        OBMol* pmol = dynamic_cast<OBMol*> (pOb);
        bool ret=false;
        if(pmol)
            ret=WriteMolecule(pmol,pConv);
        delete pOb;
        return ret;
    };
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
    OBMol2Smi m2s;

    // This is a hack to prevent recursion problems.
    //  we still need to fix the underlying problem -GRH
    if (mol.NumAtoms() > 1000)
    {
#ifdef HAVE_SSTREAM
      stringstream errorMsg;
#else
      strstream errorMsg;
#endif
      errorMsg << "SMILES Conversion failed: Molecule is too large to convert. Open Babel is currently limited to 1000 atoms." << endl;
      errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms " << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
      return(false);
    }

    m2s.Init();
    //m2s.AssignCisTrans(mol);
    m2s.CorrectAromaticAmineCharge(mol);
    m2s.CreateSmiString(mol,buffer);

    OBAtom *atom;
    vector<int>::iterator i;
    vector<int> order = m2s.GetOutputOrder();
    ofs << buffer << endl;

    int j;
    for (j = 0;j < mol.NumConformers();j++)
    {
        mol.SetConformer(j);
        for (i = order.begin();i != order.end();i++)
        {
            atom = mol.GetAtom(*i);
            sprintf(buffer,"%9.3f %9.3f %9.3f",atom->GetX(),atom->GetY(),atom->GetZ());
            ofs << buffer<< endl;
        }
    }
    return(true);
}

OBSmiNode::OBSmiNode(OBAtom *atom)
{
    _atom = atom;
    _parent = NULL;
    _nextnode.clear();
    _nextbond.clear();
}

void OBSmiNode::SetNextNode(OBSmiNode *node,OBBond *bond)
{
    _nextnode.push_back(node);
    _nextbond.push_back(bond);
}

OBSmiNode::~OBSmiNode()
{
    vector<OBSmiNode*>::iterator i;
    for (i = _nextnode.begin();i != _nextnode.end();i++)
        delete (*i);
}


bool WriteTheSmiles(OBMol & mol,char *out)
{
    char buffer[2*BUFF_SIZE];

    OBMol2Smi m2s;

    m2s.Init();
    m2s.CorrectAromaticAmineCharge(mol);
    m2s.CreateSmiString(mol,buffer);

    strcpy(out,buffer);
    return(true);

}
}
