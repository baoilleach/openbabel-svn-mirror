/**********************************************************************
Copyright (C) 2008 Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <vector>
#include <map>

#include <sstream>

using namespace std;
namespace OpenBabel
{

  class PQRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PQRFormat()
    {
      OBConversion::RegisterFormat("pqr",this, "chemical/x-pqr");
    }

    virtual const char* Description() //required
    {
      return
        "PQR format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-pqr"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
  	virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PQRFormat thePQRFormat;

  ////////////////////////////////////////////////////////////////
  /// Utility functions
  static bool parseAtomRecord(char *buffer, OBMol &mol,int /*chainNum*/);
  static double parseAtomRadius(char *buffer, OBMol &mol);
  static double parseAtomCharge(char *buffer, OBMol &mol);

  /////////////////////////////////////////////////////////////////
 	int PQRFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
      ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          -- n;
      }
      
    return ifs.good() ? 1 : -1;       
  }
  /////////////////////////////////////////////////////////////////
  bool PQRFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int chainNum = 1;
    char buffer[BUFF_SIZE];
    OBBitVec bs;
    vector<double> charges, radii;
    string line, key, value;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDB-like file, we read all chain/res info.

    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          break;
        if (EQn(buffer,"END",3)) {
          // eat anything until the next ENDMDL
          while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
          break;
        }
        if (EQn(buffer,"TER",3)) {
          chainNum++;
          continue;
        }
        if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
          {
            if( ! parseAtomRecord(buffer,mol,chainNum))
              {
                stringstream errorMsg;
                errorMsg << "WARNING: Problems reading a PQR file\n"
                         << "  Problems reading a ATOM/HETATM record.\n";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
              }

            if (EQn(buffer,"ATOM",4))
              bs.SetBitOn(mol.NumAtoms());

            // Read in the partial charge and radius too
            charges.push_back( parseAtomCharge(buffer, mol) );
            radii.push_back( parseAtomRadius(buffer, mol) );
            continue;
          }
        }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    // Use residue definitions to assign bond orders
    resdat.AssignBonds(mol,bs);

    mol.EndModify();

    /*Now assign hetatm bonds based on distance*/
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    FOR_ATOMS_OF_MOL(a, mol) {
      // WARNING: Atom index issue here
      a->SetPartialCharge(charges[a->GetIdx() - 1]);

      cerr << " charge : " << charges[a->GetIdx() - 1] << endl;

      if (!a->HasData("Radius")) {
        std::ostringstream s;
        s << radii[ a->GetIdx()-1 ];
        OBPairData *p = new OBPairData;
        p->SetAttribute("Radius");
        p->SetValue( s.str() );
        a->SetData(p);
      }

      cerr << " radius : " << radii[a->GetIdx() - 1] << endl;

    }
    mol.SetPartialChargesPerceived();

    // clean out remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    return(true);
  }

  static double parseAtomCharge(char *buffer, OBMol &mol)
  // In PQR format, either:
  // Field name, atom number, atom name, residue name, residue number
  //    x y z charge radius
  // OR
  // Field, atom number, atom name, chain id, residue number, X, Y, Z, chg, rad
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 10)
      return atof(vs[8].c_str());
    else if (vs.size() == 11)
      return atof(vs[9].c_str());

    return 0.0;
  }

  static double parseAtomRadius(char *buffer, OBMol &mol)
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 10)
      return atof(vs[9].c_str());
    else if (vs.size() == 11)
      return atof(vs[10].c_str());

    return 0.0;
  }
  
  static bool parseAtomRecord(char *buffer, OBMol &mol,int /*chainNum*/)
  /* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */
  {
    string sbuf = &buffer[6];
    if (sbuf.size() < 48)
      return(false);

    bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;

    /* serial number */
    string serno = sbuf.substr(0,5);
    //SerialNum(the_atom) = atoi(tmp_str);

    /* atom name */
    string atmid = sbuf.substr(6,4);
    
    /* chain */
    char chain = sbuf.substr(15,1)[0];

    /* element */
    string element;
    if (sbuf.size() > 71)
      element = sbuf.substr(70,2);
    else
      element = "  ";

    //trim spaces on the right and left sides
    while (!atmid.empty() && atmid[0] == ' ')
      atmid = atmid.substr(1,atmid.size()-1);

    while (!atmid.empty() && atmid[atmid.size()-1] == ' ')
      atmid = atmid.substr(0,atmid.size()-1);

    /* residue name */

    string resname = sbuf.substr(11,3);
    if (resname == "   ")
      resname = "UNK";
    else
      {
        while (!resname.empty() && resname[0] == ' ')
          resname = resname.substr(1,resname.size()-1);

        while (!resname.empty() && resname[resname.size()-1] == ' ')
          resname = resname.substr(0,resname.size()-1);
      }

    string type;
    if (EQn(buffer,"ATOM",4))
      {
        type = atmid.substr(0,2);
        if (isdigit(type[0])) {
          // sometimes non-standard files have, e.g 11HH
          if (!isdigit(type[1])) type = atmid.substr(1,1);
          else type = atmid.substr(2,1); 
        } else if (sbuf[6] == ' ' &&
                   strncasecmp(type.c_str(), "Zn", 2) != 0 &&
                   strncasecmp(type.c_str(), "Fe", 2) != 0 ||
		               isdigit(type[1]))	//type[1] is digit in Platon
          type = atmid.substr(0,1);     // one-character element
        

        if (resname.substr(0,2) == "AS" || resname[0] == 'N')
          {
            if (atmid == "AD1")
              type = "O";
            if (atmid == "AD2")
              type = "N";
          }
        if (resname.substr(0,3) == "HIS" || resname[0] == 'H')
          {
            if (atmid == "AD1" || atmid == "AE2")
              type = "N";
            if (atmid == "AE1" || atmid == "AD2")
              type = "C";
          }
        if (resname.substr(0,2) == "GL" || resname[0] == 'Q')
          {
            if (atmid == "AE1")
              type = "O";
            if (atmid == "AE2")
              type = "N";
          }
        // fix: #2002557
        if (atmid[0] == 'H' && 
            (atmid[1] == 'D' || atmid[1] == 'E' || 
             atmid[1] == 'G' || atmid[1] == 'H')) // HD, HE, HG, HH, ..
          type = "H";
      }
    else //must be hetatm record
      {
        if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' ')))
          {
            if (isalpha(element[0]))
              type = element.substr(0,2);
            else
              type = element.substr(1,1);
            if (type.size() == 2)
              type[1] = tolower(type[1]);
          }
        else
          {
            if (isalpha(atmid[0])) 
              {
              if (atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' '))
                type = atmid.substr(0,2);
              else if (atmid[0] == 'A') // alpha prefix
                type = atmid.substr(1, atmid.size() - 1);
              else
                type = atmid.substr(0,1);
              }
            else if (atmid[0] == ' ')
              type = atmid.substr(1,1); // one char element
            else
              type = atmid.substr(1,2);

            // Some cleanup steps
            if (atmid == resname)
              {
                type = atmid;
                if (type.size() == 2)
                  type[1] = tolower(type[1]);
              }
            else
              if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
                  resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                  resname == "NDP" || resname == "ABA")
                {
                  if (type.size() > 1)
                    type = type.substr(0,1);
                  //type.erase(1,type.size()-1);
                }
              else
                if (isdigit(type[0]))
                  {
                    type = type.substr(1,1);
                  }
                else
                  if (type.size() > 1 && isdigit(type[1]))
                    type = type.substr(0,1);
                  else
                    if (type.size() > 1 && isalpha(type[1])) {
                      if (type[0] == 'O' && type[1] == 'H')
                        type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
                      else if(isupper(type[1]))
                        {
                          type[1] = tolower(type[1]);
                        }
                    }
          }
        
      }

    OBAtom atom;
    /* X, Y, Z */
    string xstr = sbuf.substr(24,8);
    string ystr = sbuf.substr(32,8);
    string zstr = sbuf.substr(40,8);
    vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
    atom.SetVector(v);
    atom.ForceImplH();

    // useful for debugging unknown atom types (e.g., PR#1577238)
    //    cout << mol.NumAtoms() + 1 << " " << atmid << " type: " << type << endl;
    atom.SetAtomicNum(etab.GetAtomicNum(type.c_str()));

    /* residue sequence number */
    string resnum = sbuf.substr(16,4);
    OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
    if (res == NULL || res->GetName() != resname 
        || res->GetNumString() != resnum)
      {
        vector<OBResidue*>::iterator ri;
        for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname 
              && res->GetNumString() == resnum
              && static_cast<int>(res->GetChain()) == chain)
            break;

        if (res == NULL)
          {
            res = mol.NewResidue();
            res->SetChain(chain);
            res->SetName(resname);
            res->SetNum(resnum);
          }
      }

    if (!mol.AddAtom(atom))
      return(false);
    else
      {
        OBAtom *atom = mol.GetAtom(mol.NumAtoms());

        res->AddAtom(atom);
        res->SetSerialNum(atom, atoi(serno.c_str()));
        res->SetAtomID(atom, sbuf.substr(6,4));
        res->SetHetAtom(atom, hetatm);

        return(true);
      }
  }

} //namespace OpenBabel
