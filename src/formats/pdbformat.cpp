/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003-2006 Geoffrey R. Hutchison
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

  class PDBFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PDBFormat()
    {
      OBConversion::RegisterFormat("pdb",this, "chemical/x-pdb");
      OBConversion::RegisterFormat("ent",this, "chemical/x-pdb");
    }

    virtual const char* Description() //required
    {
      return
        "Protein Data Bank format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.wwpdb.org/docs.html";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-pdb"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
  	virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PDBFormat thePDBFormat;

  /////////////////////////////////////////////////////////////////

  // Now in atomrecord.cpp (shared with PQR format)
  bool ParseAtomRecord(char *, OBMol &,int);
  bool ParseConectRecord(char *,OBMol &);

  //extern OBResidueData    resdat; now in mol.h

  /////////////////////////////////////////////////////////////////
 	int PDBFormat::SkipObjects(int n, OBConversion* pConv)
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
  bool PDBFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
    string line, key, value;
    OBPairData *dp;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDB file, we read all chain/res info.

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
            ParseAtomRecord(buffer,mol,chainNum);
            if (EQn(buffer,"ATOM",4))
              bs.SetBitOn(mol.NumAtoms());
            continue;
          }

        if (EQn(buffer,"CONECT",6)) {
          ParseConectRecord(buffer,mol);
          continue;
        }

        // another record type, add it as an OBPairData entry
        line = buffer;
        // if the file is valid, all lines should have more than 6 characters
        if (line.length() < 6)
		  {
            stringstream errorMsg;
            errorMsg << "ERROR: not a valid PDB file" << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
            return false;
          }
       key = line.substr(0,6); // the first 6 characters are the record name
        Trim(key);
        value = line.substr(6);
        
        // We haven't found this record yet
        if (!mol.HasData(key)) {
          dp = new OBPairData;
          dp->SetAttribute(key);
          dp->SetValue(value);
          dp->SetOrigin(fileformatInput);
          mol.SetData(dp);
        } 
        // Add on additional lines
        else {
          dp = static_cast<OBPairData*>(mol.GetData(key));
          line = dp->GetValue();
          line += '\n';
          line += value;
          dp->SetValue(line);
        }
      }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    resdat.AssignBonds(mol,bs);
    /*assign hetatm bonds based on distance*/

    mol.EndModify();
    // Clear all virtual bond data
    vector<OBGenericData*> vbonds = mol.GetAllData(OBGenericDataType::VirtualBondData);
    mol.DeleteData(vbonds);


    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // clean out remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    return(true);
  }

  /////////////////////////////////////////////////////////////////////////
  //! Utility function to read a 5-digit integer starting from a specified column
  /*! This function reads a 5-digit integer, starting from column
    columnAsSpecifiedInPDB from the buffer, converts it to a long
    integer, and returns either false or true, if the conversion was
    successful or not. If the conversion was not successful, the target
    is set to a random value.
 
    For instance, the PDB Format Description for a CONECT record specifies
 
    COLUMNS        DATA TYPE        FIELD           DEFINITION
    ---------------------------------------------------------------------------------
    1 -  6         Record name      "CONECT"
    7 - 11         Integer          serial          Atom serial number
    ...
 
    To read the Atom serial number, you would call
 
    long int target;
    if ( readIntegerFromRecord(buffer, 7, &target) == false ) {
    cerr << "Could not parse" << endl;
    }
  
    This function does not check the length of the buffer, or
    strlen(buffer). If the buffer is not long enough => SEGFAULT. 
  */
  static bool readIntegerFromRecord(char *buffer, unsigned int columnAsSpecifiedInPDB, long int *target)
  {
    char integerBuffer[6];
    integerBuffer[5] = '\0';

    strncpy(integerBuffer, buffer+columnAsSpecifiedInPDB-1, 5);

    char *errorCheckingEndPtr;
    *target = strtol(integerBuffer, &errorCheckingEndPtr, 10);
    if (integerBuffer == errorCheckingEndPtr)
      return(false);
    return(true);
  }

  //! Read a CONECT record
  /*! This function reads a CONECT record, as specified
    http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html,
    in short:
 
    COLUMNS         DATA TYPE        FIELD           DEFINITION
    ---------------------------------------------------------------------------------
    1 -  6         Record name      "CONECT"
    7 - 11         Integer          serial          Atom serial number
    12 - 16         Integer          serial          Serial number of bonded atom
    17 - 21         Integer          serial          Serial number of bonded atom
    22 - 26         Integer          serial          Serial number of bonded atom
    27 - 31         Integer          serial          Serial number of bonded atom
    32 - 36         Integer          serial          Serial number of hydrogen bonded atom
    37 - 41         Integer          serial          Serial number of hydrogen bonded atom
    42 - 46         Integer          serial          Serial number of salt bridged atom
    47 - 51         Integer          serial          Serial number of hydrogen bonded atom
    52 - 56         Integer          serial          Serial number of hydrogen bonded atom
    57 - 61         Integer          serial          Serial number of salt bridged atom
 
    Hydrogen bonds and salt bridges are ignored. --Stefan Kebekus.
  */

  bool ParseConectRecord(char *buffer,OBMol &mol)
  {
    stringstream errorMsg;
    string clearError;

    // Setup strings and string buffers
    vector<string> vs;
    buffer[70] = '\0';
    if (strlen(buffer) < 70)
      {
        errorMsg << "WARNING: Problems reading a PDB file\n"
                 << "  Problems reading a CONECT record.\n"
                 << "  According to the PDB specification,\n"
                 << "  the record should have 70 columns, but OpenBabel found "
                 << strlen(buffer) << " columns." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obInfo);
        errorMsg.str(clearError);
      }

    // Serial number of the first atom, read from column 7-11 of the
    // connect record, to which the other atoms connect to.
    long int startAtomSerialNumber;
    // A pointer to the first atom.
    OBAtom *firstAtom = NULL;
    // Serial numbers of the atoms which bind to firstAtom, read from
    // columns 12-16, 17-21, 22-27 and 27-31 of the connect record. Note
    // that we reserve space for 5 integers, but read only four of
    // them. This is to simplify the determination of the bond order;
    // see below.
    long int boundedAtomsSerialNumbers[5]  = {0,0,0,0,0};
    // Bools which tell us which of the serial numbers in
    // boundedAtomsSerialNumbers are read from the file, and which are
    // invalid
    bool boundedAtomsSerialNumbersValid[5] = {false, false, false, false, false};

    // Pragmatic approach -- too many non-standard PDB files out there
    // (including some old ones from us)
    // So if we have a small number of atoms, then try to break by spaces
    // Otherwise (i.e., NumAtoms() > 9,999 we need to go by position)
    // We'll switch back and forth a few times to save duplicating common code

    if (mol.NumAtoms() <= 9999)
      {
        // make sure we don't look at salt bridges or whatever, so cut the buffer short
        buffer[32] = '\0';
        tokenize(vs,buffer);
        if( vs.empty() || vs.size() < 2) 
          return false;
        vs.erase(vs.begin()); // remove "CONECT"

        startAtomSerialNumber = atoi(vs[0].c_str());
      }
    else
      {
        if (readIntegerFromRecord(buffer, 7, &startAtomSerialNumber) == false)
          {
            errorMsg << "WARNING: Problems reading a PDB file\n"
                     << "  Problems reading a CONECT record.\n"
                     << "  According to the PDB specification,\n"
                     << "  columns 7-11 should contain the serial number of an atom.\n"
                     << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
      }

    vector<OBAtom*>::iterator i;
    for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) {
      // atoms may not have residue information, but if they do,
      // check serial numbers
      if (a1->GetResidue() != NULL && 
          static_cast<long int>(a1->GetResidue()->
                                GetSerialNum(a1)) == startAtomSerialNumber)
        {
          firstAtom = a1;
          break;
        }
    }

    if (firstAtom == NULL)
      {
        errorMsg << "WARNING: Problems reading a PDB file:\n"
                 << "  Problems reading a CONECT record.\n"
                 << "  According to the PDB specification,\n"
                 << "  columns 7-11 should contain the serial number of an atom.\n"
                 << "  No atom was found with this serial number.\n"
                 << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
        return(false);
      }

    if (mol.NumAtoms() < 9999)
      {
        if (vs.size() > 1) boundedAtomsSerialNumbers[0] = atoi(vs[1].c_str());
        if (vs.size() > 2) boundedAtomsSerialNumbers[1] = atoi(vs[2].c_str());
        if (vs.size() > 3) boundedAtomsSerialNumbers[2] = atoi(vs[3].c_str());
        if (vs.size() > 4) boundedAtomsSerialNumbers[3] = atoi(vs[4].c_str());

        unsigned int limit = 4;
        if (vs.size() <= 4)
          limit = vs.size() - 1;

        for (unsigned int s = 0; s < limit; ++s)
          boundedAtomsSerialNumbersValid[s] = true;
      }
    else
      {
        // Now read the serial numbers. If the first serial number is not
        // present, this connect record probably contains only hydrogen
        // bonds and salt bridges, which we ignore. In that case, we just
        // exit gracefully.
        boundedAtomsSerialNumbersValid[0] = readIntegerFromRecord(buffer, 12, boundedAtomsSerialNumbers+0);
        if (boundedAtomsSerialNumbersValid[0] == false)
          return(true);
        boundedAtomsSerialNumbersValid[1] = readIntegerFromRecord(buffer, 17, boundedAtomsSerialNumbers+1);
        boundedAtomsSerialNumbersValid[2] = readIntegerFromRecord(buffer, 22, boundedAtomsSerialNumbers+2);
        boundedAtomsSerialNumbersValid[3] = readIntegerFromRecord(buffer, 27, boundedAtomsSerialNumbers+3);
      }

    // Now iterate over the VALID boundedAtomsSerialNumbers and connect
    // the atoms.
    for(unsigned int k=0; boundedAtomsSerialNumbersValid[k]; k++)
      {
        // Find atom that is connected to, write an error message
        OBAtom *connectedAtom = 0L;
        for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) {
          // again, atoms may not have residues, but if they do, check serials
          if (a1->GetResidue() != NULL &&
              static_cast<long int>(a1->GetResidue()->
                                    GetSerialNum(a1)) == boundedAtomsSerialNumbers[k])
            {
              connectedAtom = a1;
              break;
            }
        }
        if (connectedAtom == 0L)
          {
            errorMsg << "WARNING: Problems reading a PDB file:\n"
                     << "  Problems reading a CONECT record.\n"
                     << "  According to the PDB specification,\n"
                     << "  Atoms with serial #" << startAtomSerialNumber 
                     << " and #" << boundedAtomsSerialNumbers[k]
                     << " should be connected\n"
                     << "  However, an atom with serial #" << boundedAtomsSerialNumbers[k] << " was not found.\n"
                     << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }

        // Figure the bond order
        unsigned char order = 0;
        while(boundedAtomsSerialNumbersValid[k+order+1] && (boundedAtomsSerialNumbers[k+order]
                                                            == boundedAtomsSerialNumbers[k+order+1]))
          order++;
        k += order;
	
        // Generate the bond
        mol.AddBond(firstAtom->GetIdx(), connectedAtom->GetIdx(), order+1);
      }
    return(true);
  }

  //////////////////////////////////////////////////////////////////////////////
  bool PDBFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    char type_name[10], padded_name[10];
    char the_res[10];
    const char *element_name;
    int res_num;
    bool het=true;
    int model_num = 0;
    if (!pConv->IsLast() || pConv->GetOutputIndex() > 1)
      { // More than one molecule record
      model_num = pConv->GetOutputIndex(); // MODEL 1-based index
      snprintf(buffer, BUFF_SIZE, "MODEL %8d", model_num);
      ofs << buffer << endl;
      }

    if (strlen(mol.GetTitle()) > 0)
      snprintf(buffer, BUFF_SIZE, "COMPND    %s ",mol.GetTitle());
    else
      snprintf(buffer, BUFF_SIZE, "COMPND    UNNAMED");
    ofs << buffer << endl;

    snprintf(buffer, BUFF_SIZE, "AUTHOR    GENERATED BY OPEN BABEL %s",BABEL_VERSION);
    ofs << buffer << endl;

    // before we write any records, we should check to see if any coord < -1000
    // which will cause errors in the formatting

    double minX, minY, minZ;
    minX = minY = minZ = -999.0f;
    FOR_ATOMS_OF_MOL(a, mol)
      {
        if (a->GetX() < minX)
          minX = a->GetX();
        if (a->GetY() < minY)
          minY = a->GetY();
        if (a->GetZ() < minZ)
          minZ = a->GetZ();
      }
    vector3 transV = VZero;
    if (minX < -999.0)
      transV.SetX(-1.0*minX - 900.0);
    if (minY < -999.0)
      transV.SetY(-1.0*minY - 900.0);
    if (minZ < -999.0)
      transV.SetZ(-1.0*minZ - 900.0);

    // if minX, minY, or minZ was never changed, shift will be 0.0f
    // otherwise, move enough so that smallest coord is > -999.0f
    mol.Translate(transV);

    OBAtom *atom;
    OBResidue *res;
    for (i = 1; i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        strncpy(type_name, etab.GetSymbol(atom->GetAtomicNum()), sizeof(type_name));
        type_name[sizeof(type_name) - 1] = '\0';

        //two char. elements are on position 13 and 14 one char. start at 14
        if (strlen(type_name) > 1)
          type_name[1] = toupper(type_name[1]);
        else
          {
            char tmp[10];
            strncpy(tmp, type_name, 10);
            snprintf(type_name, sizeof(type_name), " %-3s", tmp);
          }

        if ( (res = atom->GetResidue()) != 0 )
          {
            het = res->IsHetAtom(atom);
            snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
            snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());

            //two char. elements are on position 13 and 14 one char. start at 14
            if (strlen(etab.GetSymbol(atom->GetAtomicNum())) == 1)
              {
                if (strlen(type_name) < 4)
                  {
                    char tmp[16];
                    strncpy(tmp, type_name, 16);
                    snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
                    strncpy(type_name,padded_name,4);
                    type_name[4] = '\0';
                  }
                else
                  {
                    type_name[4] = type_name[3];
                    type_name[3] = type_name[2];
                    type_name[2] = type_name[1];
                    type_name[1] = type_name[0];
                    type_name[0] = type_name[4];
                    type_name[4] = '\0';
                  }
              }
            res_num = res->GetNum();
          }
        else
          {
            strcpy(the_res,"UNK");
            snprintf(padded_name,sizeof(padded_name), "%s",type_name);
            strncpy(type_name,padded_name,4);
            type_name[4] = '\0';
            res_num = 1;
          }

        element_name = etab.GetSymbol(atom->GetAtomicNum());
        snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
                 het?"HETATM":"ATOM  ",
                 i,
                 type_name,
                 the_res,
                 res_num,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 element_name);
        ofs << buffer;
      }

    OBAtom *nbr;
    vector<OBBond*>::iterator k;
    for (i = 1; i <= mol.NumAtoms(); i ++)
      {
        atom = mol.GetAtom(i);
        if (atom->GetValence() == 0)
          continue; // no need to write a CONECT record -- no bonds

        snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
        ofs << buffer;
        // Write out up to 4 real bonds per line PR#1711154
        int currentValence = 0;
        for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
          {
            snprintf(buffer, BUFF_SIZE, "%5d", nbr->GetIdx());
            ofs << buffer;
            if (++currentValence % 4 == 0) {
              // Add the trailing space to finish this record
              ofs << "                                       \n";
              // write the start of a new CONECT record
              snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
              ofs << buffer;              
            }
          }

        // Add trailing spaces
        int remainingValence = atom->GetValence() % 4;
        for (int count = 0; count < (4 - remainingValence); count++) {
          snprintf(buffer, BUFF_SIZE, "     ");
          ofs << buffer;
        }
        ofs << "                                       \n";
      }

    snprintf(buffer, BUFF_SIZE, "MASTER        0    0    0    0    0    0    0    0 ");
    ofs << buffer;
    snprintf(buffer, BUFF_SIZE, "%4d    0 %4d    0\n",mol.NumAtoms(),mol.NumAtoms());
    ofs << buffer;

    ofs << "END\n";
    if (model_num) {
      ofs << "ENDMDL" << endl;
    }

    return(true);
  }


} //namespace OpenBabel
