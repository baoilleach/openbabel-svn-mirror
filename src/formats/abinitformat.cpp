/**********************************************************************
Copyright (C) 2011 by Geoffrey Hutchison

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

using namespace std;
namespace OpenBabel
{

#define BOHR_TO_ANGSTROM 0.5291772108
#define ANGSTROM_TO_BOHR 1.0 / BOHR_TO_ANGSTROM


  class ABINITFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ABINITFormat()
    {
      OBConversion::RegisterFormat("abinit",this);
    }

    virtual const char* Description() //required
    {
      return
        "ABINIT Output Format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://abinit.org/" ;}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  ABINITFormat theABINITFormat;

  /////////////////////////////////////////////////////////////////
  bool ABINITFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str;
    vector<string> vs;

    OBAtom *atom;
    int natom;
    vector<int> atomicNumbers, atomTypes;
    double x, y, z;
    vector<vector3> atomPositions;
    vector<double> energies;

    // Translation vectors (if present)
    // aka rprim
    vector3 translationVectors[3];
    double acell[3]; // scale of lattice vectors
    int numTranslationVectors = 0;
    int symmetryCode = 0;

    SpaceGroup *sg = NULL;

    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE))
      {
        // tokens are listed in alphabetical order for code clarity
        if (strstr(buffer, "acell")) {
          tokenize(vs, buffer);
          if (vs.size() < 4)
            continue; // invalid line

          // acell=  7.6967369631E+00  7.6967369631E+00  7.6967369631E+00
          for (int i = 0; i < 3; ++i) {
            acell[i] = atof(vs[i+1].c_str());
          }
        }
        // Sometimes Cartesian has lower-case letter
        else if (strstr(buffer, "artesian coordinates")) {
          double unit = BOHR_TO_ANGSTROM;
          if (strstr(buffer, "ngstrom"))
            unit = 1.0; // no conversion needed

          ifs.getline(buffer,BUFF_SIZE);
          tokenize(vs, buffer);
          while(vs.size() == 3) {
            x = atof(vs[1].c_str());
            y = atof(vs[2].c_str());
            z = atof(vs[3].c_str());

            atomPositions.push_back(vector3(x, y, z));
            // get next line
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs, buffer);
          }
        }
        else if (strstr(buffer, "natom")) {
          tokenize(vs, buffer);
          if (vs.size() != 2)
            continue;
          natom = atoi(vs[1].c_str());
        }
        else if (strstr(buffer, "rprim")) {
          numTranslationVectors = 0;
          int column;
          for (int i = 0; i < 3; ++i) {
            tokenize(vs, buffer);
            if (vs.size() < 3)
              break;

            // first line, rprim takes up a token
            if (i == 0)
              column = 1;
            else
              column = 0;

            x = atof((char*)vs[column].c_str()) * BOHR_TO_ANGSTROM;
            y = atof((char*)vs[column+1].c_str()) * BOHR_TO_ANGSTROM;
            z = atof((char*)vs[column+2].c_str()) * BOHR_TO_ANGSTROM;
            translationVectors[numTranslationVectors++].Set(x, y,z);
            ifs.getline(buffer,BUFF_SIZE);
          }
        }
        else if (strstr(buffer, "Symmetries")) {
          tokenize(vs, buffer, "()");
          // Should be something like (#160)
          symmetryCode = atoi(vs[1].substr(1).c_str());
        }
        else if (strstr(buffer, "typat")) {
          tokenize(vs, buffer);
          atomTypes.clear();
          for (int i = 1; i < vs.size(); ++i) {
            atomTypes.push_back(atoi(vs[i].c_str()));
          }
        }
        else if (strstr(buffer, "znucl")) {
          tokenize(vs, buffer);
          // make sure znucl is first token
          if (vs[0] != "znucl")
            continue;
          // push back the remaining tokens into atomicNumbers
          atomicNumbers.clear();
          atomicNumbers.push_back(0); // abinit starts typat with 1
          for (int i = 1; i < vs.size(); ++i)
            atomicNumbers.push_back(int(atof(vs[i].c_str())));
        }
        // xangst
        // forces
      }

    for (int i = 0; i < natom; ++i) {
      atom = mol.NewAtom();
      //set atomic number
      int idx = atom->GetIdx();
      int type = atomTypes[idx - 1];
      atom->SetAtomicNum(atomicNumbers[type]);
      // we set the coordinates by conformers in another loop
    }

    mol.EndModify();

    int numConformers = atomPositions.size() / natom;
    for (int i = 0; i < numConformers; ++i) {
      double *coordinates = new double[natom * 3];
      for (int j = 0; j < natom; ++j) {
        vector3 currentPosition = atomPositions[i*natom + j];
        cout << currentPosition << endl;
        coordinates[j*3] = currentPosition.x();
        coordinates[j*3 + 1] = currentPosition.y();
        coordinates[j*3 + 2] = currentPosition.z();
      }
      mol.AddConformer(coordinates);
    }
    mol.DeleteConformer(0); // first one is created by EndModify, bunch of 0s
    mol.SetConformer(mol.NumConformers() - 1);

    cout << atomPositions.size() << " " << mol.NumConformers() << endl;

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // Attach unit cell translation vectors if found
    if (numTranslationVectors > 0) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(acell[0] * translationVectors[0], acell[1] * translationVectors[1], acell[2] * translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      if (symmetryCode)
        uc->SetSpaceGroup(symmetryCode);
      mol.SetData(uc);
    }

    mol.SetTitle(title);
    return(true);
  }

} //namespace OpenBabel
