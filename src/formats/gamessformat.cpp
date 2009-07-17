/**********************************************************************
  Copyright (C) 2000 by OpenEye Scientific Software, Inc.
  Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
  Some portions Copyright (C) 2004 by Chris Morley
  Some portions Copyright (C) 2006 by Donald E. Curtis

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

#include <algorithm>

using namespace std;
namespace Gamess
{
} 

namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989
  class GAMESSOutputFormat : public OBMoleculeFormat
  {

  public:
    //Register this format type ID
    GAMESSOutputFormat()
    {
      OBConversion::RegisterFormat("gam", this, "chemical/x-gamess-output");
      OBConversion::RegisterFormat("gamout",this);
      OBConversion::RegisterFormat("gamess",this);
    }

    virtual const char* Description() //required
    {
      return
        "GAMESS Output\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n"
        "  c  Read multiple conformers\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

    virtual const char* GetMIMEType() 
    { return "chemical/x-gamess-output"; };

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

  private:
    //! \brief Parse GAMESS options section.
    //void ParseSection(char *tag, OBSetData *set, istream &ifs);

  };
  //***

  //Make an instance of the format class
  GAMESSOutputFormat theGAMESSOutputFormat;


  class GAMESSInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GAMESSInputFormat()
    {
      OBConversion::RegisterFormat("inp",this, "chemical/x-gamess-input");
      OBConversion::RegisterFormat("gamin",this);
      // Command-line keywords
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
    }


    virtual const char* Description() //required
    {
      return
        "GAMESS Input\n"
        "Write Options e.g. -xk\n"
        "  k  \"keywords\" Use the specified keywords for input\n"
        "  f    <file>     Read the file specified for input keywords\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

    virtual const char* GetMIMEType() 
    { return "chemical/x-gamess-input"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return WRITEONEONLY; // | NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  GAMESSInputFormat theGAMESSInputFormat;

  /////////////////////////////////////////////////////////////////
  /* this function is for parsing default options too.  it is decided that
   * we should only parse parameters that the user specified in the input
   * deck and not EVERY option which is defaulted to by GAMESS.
   void GAMESSOutputFormat::ParseSection(char *tag, OBSetData *set, istream &ifs)
   {
   char buffer[BUFF_SIZE];
   OBSetData *curset = (OBSetData *)set->GetData(tag);
   if(!curset)
   {
   curset = new OBSetData();
   curset->SetOrigin(fileformatInput);
   curset->SetAttribute(tag);
   set->AddData(curset);
   }

   string attr, value;
   char *ptr;

   for( ; ; )
   {
   ifs.getline(buffer,BUFF_SIZE);
   ptr = buffer;

   // trim initial line whitespace
   while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;
   // If this is it be done
   if(*ptr == '\0') break;

   // parse a line
   while(true)
   {
   attr.clear();
   value.clear();

   // Trim leading whitespace
   while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

   // Read the attribute name
   while(*ptr != ' ' && *ptr != '=' && *ptr != '\0') attr += toupper(*(ptr++));

   // If this is it, be done
   if(*ptr == '\0') break;

   // Read to next non-whitespace
   while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

   // Keywords are only one word.  So we must have extra data we don't want.
   // So in this case we just ignore it and go on like we're ready for the
   // next pair.
   if(*ptr != '=') continue;

   // Read to next non-whitespace
   while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

   while((*ptr == ' ' || *ptr == '\t' || *ptr == '=') && *ptr != '\0') ptr++;

   // Read the attribute value.
   while(*ptr != ' ' && *ptr != '\0') value += toupper(*(ptr++));


   if(attr == "IGAUSS") { attr = "NGAUSS"; }

   // cout << attr << "/" << value << endl;

   OBPairData *data = new OBPairData();
   data = new OBPairData();
   data->SetAttribute(attr);
   data->SetValue(value);
   data->SetOrigin(fileformatInput);

   curset->AddData(data);
   }
   }
   }
  */
  bool GAMESSOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    bool hasPartialCharges = false;
    int HOMO = 0;
    vector<double> orbitals;

    vector<double> frequencies, intensities;
    vector< vector<vector3> > displacements;
    int lowFreqModesBegin; // the number of the first low frequency mode
    int lowFreqModesEnd; // the number of the last low frequency mode
    int numFreq, numIntens, numDisp; // GAMESS prints rotations & transl., which we ignore
    numFreq = numIntens = numDisp = 0;

    // must build generic data while we parse then add at the end.
    OBSetData *gmsset = new OBSetData();
    gmsset->SetAttribute("gamess");
    gmsset->SetOrigin(fileformatInput);

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"ATOMIC                      COORDINATES (BOHR)") != NULL)
          {
            mol.Clear();
            mol.BeginModify();

            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
                x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
                y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
                z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
                atom->SetVector(x,y,z);
                vs[1].erase(vs[1].size() - 2, 2);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"MULTIPOLE COORDINATES, ELECTRONIC AND NUCLEAR CHARGES") != NULL)
          {
            /*This set of EFP coordinates belongs only to the
             * conformer directly above this (ATOMIC   COORDINATES (BOHR))
             */
            ifs.getline(buffer,BUFF_SIZE);      // column headings
            ifs.getline(buffer,BUFF_SIZE);
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while(vs.size() == 6)
              {
                int atomicNum;
                /* For the included EFP1 potentials,
                 * the atom name may start with "Z"
                 * or have a non-zero nuclear charge
                 */
                if (atof((char*)vs[5].c_str()) > 0.0) {
                  atom = mol.NewAtom();
                  atomicNum=etab.GetAtomicNum(vs[0].substr(0,1).c_str()); 
                  atom->SetAtomicNum(atomicNum);
                  x = atof((char*)vs[1].c_str())* BOHR_TO_ANGSTROM;
                  y = atof((char*)vs[2].c_str())* BOHR_TO_ANGSTROM;
                  z = atof((char*)vs[3].c_str())* BOHR_TO_ANGSTROM;
                  atom->SetVector(x,y,z);

                }
                else if (vs[0].substr(0,1) == "Z") {
                  atom = mol.NewAtom();
                  atomicNum=etab.GetAtomicNum(vs[0].substr(1,1).c_str());
                  atom->SetAtomicNum(atomicNum);
                  x = atof((char*)vs[1].c_str())* BOHR_TO_ANGSTROM;
                  y = atof((char*)vs[2].c_str())* BOHR_TO_ANGSTROM;
                  z = atof((char*)vs[3].c_str())* BOHR_TO_ANGSTROM;
                  atom->SetVector(x,y,z);

                }
                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }

        else if(strstr(buffer,"COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL)
          {
            mol.Clear();
            mol.BeginModify();

            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
                x = atof((char*)vs[2].c_str());
                y = atof((char*)vs[3].c_str());
                z = atof((char*)vs[4].c_str());
                atom->SetVector(x,y,z);
                vs[1].erase(vs[1].size() - 2, 2);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
            if(strstr(buffer,"COORDINATES OF FRAGMENT") != NULL)
              {
                ifs.getline(buffer,BUFF_SIZE);      // column headings
                ifs.getline(buffer,BUFF_SIZE);
                ifs.getline(buffer,BUFF_SIZE);
                //ifs.getline(buffer,BUFF_SIZE);    //FRAGNAME
                tokenize(vs,buffer);
                //while(vs.size() == 4)
                while(vs.size() > 0) {
                  if (vs.size() == 1) {
                    vector<string> vs2;
                    char delim[] = "=";
                    tokenize(vs2,buffer,delim);
                  }
                  else {
                    atom = mol.NewAtom();
                    /* For the included EFP1 potentials,
                     * the atom name may start with "Z"
                     */
                    int atomicNum;
                    if ( vs[0].substr(0,1) == "Z" ) 
                      atomicNum=etab.GetAtomicNum(vs[0].substr(1,1).c_str()); 
                    else 
                      atomicNum=etab.GetAtomicNum(vs[0].substr(0,1).c_str()); 
                    atom->SetAtomicNum(atomicNum);
                    x = atof((char*)vs[1].c_str());
                    y = atof((char*)vs[2].c_str());
                    z = atof((char*)vs[3].c_str());
                    atom->SetVector(x,y,z);
                  }
			    

                  if (!ifs.getline(buffer,BUFF_SIZE))
                    break;
                  tokenize(vs,buffer);
                }
              }
          }
        else if(strstr(buffer,"ELECTROSTATIC MOMENTS") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE); //-----
            ifs.getline(buffer,BUFF_SIZE);  // blank line
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE); // point charges TODO
            ifs.getline(buffer,BUFF_SIZE);	// column headings dipole moment
            ifs.getline(buffer,BUFF_SIZE);
            
            tokenize(vs, buffer);
            if (vs.size() == 4) {
              OBVectorData *dipoleMoment = new OBVectorData;
              dipoleMoment->SetAttribute("Dipole Moment");
              double x, y, z;
              x = atof(vs[0].c_str());
              y = atof(vs[1].c_str());
              z = atof(vs[2].c_str());
              dipoleMoment->SetData(x, y, z);
              dipoleMoment->SetOrigin(fileformatInput);
              mol.SetData(dipoleMoment);
            } 
          }
        else if(strstr(buffer,"MOPAC CHARGES") != NULL)
          {
            hasPartialCharges = true;
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 4)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"TOTAL MULLIKEN") != NULL)
          {
            hasPartialCharges = true;
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() >= 4)
              { // atom number, atomic symbol, mulliken pop, charge
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[3].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if (strstr(buffer,"NUMBER OF OCCUPIED ORBITALS") != NULL)
          {
            tokenize(vs, buffer);
            if (vs.size() == 7) // alpha
              HOMO = atoi(vs[6].c_str());
            else if (vs.size() == 8) //beta
              HOMO = atoi(vs[7].c_str());
          }
        else if (strstr(buffer, "TAKEN AS ROTATIONS AND TRANSLATIONS") != NULL)
          {
            tokenize(vs, buffer);
            if (vs.size() < 4)
              break;
            lowFreqModesBegin = atoi(vs[1].c_str());
            lowFreqModesEnd = atoi(vs[3].c_str());
          }
        else if (strstr(buffer,"FREQUENCY:") != NULL)
          {
            tokenize(vs, buffer);
            for (unsigned int i = 1; i < vs.size(); ++i) {
              if (vs[i] == "I") // artifact from previous imaginary frequency
                continue;
              ++numFreq;
              if (numFreq < lowFreqModesBegin) // imaginary frequency
                frequencies.push_back(-atof(vs[i].c_str()));
              if (numFreq > lowFreqModesEnd)
                frequencies.push_back(atof(vs[i].c_str()));
            }
            ifs.getline(buffer, BUFF_SIZE); // reduced mass
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            for (unsigned int i = 2; i < vs.size(); ++i) {
              ++numIntens;
              if (numIntens < lowFreqModesBegin || numIntens > lowFreqModesEnd)
                intensities.push_back(atof(vs[i].c_str()));
            }
            ifs.getline(buffer, BUFF_SIZE); // blank

            // now real work -- read displacements
            int prevModeCount = displacements.size();
            int newModes = frequencies.size() - displacements.size();
            vector<vector3> displacement;
            for (unsigned int i = 0; i < newModes; ++i) {
              displacements.push_back(displacement);
            }

            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            int modeCount = vs.size() - 3;
            double massNormalization;
            vector<double> x, y, z;
            while(modeCount >= 1) {
              // 1/sqrt(atomic mass)
              atom = mol.GetAtom(atoi(vs[0].c_str()));
              massNormalization = 1 / sqrt( atom->GetAtomicMass() );

              x.clear();
              // not a typo -- e.g., atom number, atom label, x, then data
              for (unsigned int i = 3; i < vs.size(); ++i) {
                x.push_back(massNormalization * atof(vs[i].c_str()));
              }
              y.clear();
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              for (unsigned int i = 1; i < vs.size(); ++i) {
                y.push_back(massNormalization * atof(vs[i].c_str()));
              }
              
              z.clear();
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              for (unsigned int i = 1; i < vs.size(); ++i) {
                z.push_back(massNormalization * atof(vs[i].c_str()));
              }

              // OK, now we have x, y, z for all new modes for one atom
              if (displacements.size()) {
                numDisp = prevModeCount;
                for (unsigned int i = 0; i < modeCount;  ++i) {
                  if (i >= modeCount - newModes){
                    displacements[numDisp++].push_back(vector3(x[i], y[i], z[i]));
		  }
                }
              }

              // Next set of atoms
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              modeCount = vs.size() - 3;
            }
          }
        else if (strstr(buffer,"EIGENVECTORS") != NULL ||
                 strstr(buffer,"MOLECULAR ORBITALS") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE); // ------ line
            ifs.getline(buffer,BUFF_SIZE); // blank
            orbitals.clear();

            while (strstr(buffer,"END OF RHF CALCULATION") == NULL &&
                   strstr(buffer,"-------") == NULL)
              {
                //loop
                ifs.getline(buffer,BUFF_SIZE); // orbitals!
                ifs.getline(buffer,BUFF_SIZE); // energies in hartree
                tokenize(vs, buffer);
                for (unsigned int i = 0; i < vs.size(); i++)
                  orbitals.push_back(27.21 * atof(vs[i].c_str()));

                ifs.getline(buffer,BUFF_SIZE); // symmetries
                // orbital coefficients
                while (ifs.getline(buffer,BUFF_SIZE) && strlen(buffer)
                       && strstr(buffer,"END") == NULL 
                       && strstr(buffer, "---") == NULL)
                  { }
                if (!ifs.good())
                  break;
              }
          }
        else if(strstr(buffer, "INPUT CARD> $"))
          {
            string attr, value;
            char *ptr;

            for( ; ; )
              {
                ptr = buffer + 14;
                tokenize(vs, ptr);

                if(vs.size() > 2)
                  {
                    OBSetData *curset = (OBSetData *)gmsset->GetData(vs[0]);
                    if(!curset)
                      {
                        curset = new OBSetData();
                        curset->SetAttribute(vs[0]);
                        curset->SetOrigin(fileformatInput);
                        gmsset->AddData(curset);
                      }

                    for(unsigned int i=1;i < vs.size() &&  vs[i].substr(0,4) != "$END"; i++) {
                      string::size_type loc = vs[i].find("=",0);
                      if(loc != string::npos)
                        {
                          // cout << vs[i].substr(0,loc) << " !!!!! " << vs[i].substr(loc+1) << endl;
                          OBPairData *data = new OBPairData();
                          data = new OBPairData();
                          data->SetAttribute(vs[i].substr(0,loc));
                          data->SetValue(vs[i].substr(loc+1));
                          data->SetOrigin(fileformatInput);
                          curset->AddData(data);
                        }
                    }
                  }

                break;

              }
          }
        /*
          else if(strstr(buffer, "$CONTRL OPTIONS"))
          {
          ParseSection("CONTRL", gmsset, ifs);
          }
          else if(strstr(buffer, "$SYSTEM OPTIONS"))
          {
          ParseSection("SYSTEM", gmsset, ifs);
          }
          else if(strstr(buffer, "BASIS OPTIONS"))
          {
          ParseSection("BASIS", gmsset, ifs);
          }
          else if(strstr(buffer, "GUESS OPTIONS"))
          {
          ParseSection("GUESS", gmsset, ifs);
          }
        */
      }
    //    cerr << " HOMO " << HOMO - 1 << " LUMO " << HOMO << " size: " << orbitals.size() << endl;
    //    cerr << title << " " << HOMO << " " << orbitals[HOMO - 1] << " " << orbitals[HOMO] << endl;

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    const char *keywordsEnable = pConv->IsOption("k",OBConversion::GENOPTIONS);

    if(keywordsEnable)
      {
        // add our gamess set
        pmol->SetData(gmsset);

        // if we have basis set data we should set our global pair data
        OBSetData *cset = (OBSetData *) gmsset->GetData("CONTRL");
        OBSetData *bset = (OBSetData *) gmsset->GetData("BASIS");

        string model = "b3lyp";
        string basis;
        string method;

        if(cset)
          {
            OBPairData *pd = NULL;

            pd = (OBPairData *) cset->GetData("SCFTYP");
            if(pd)
              {
                if(pd->GetValue() == "RHF")
                  {
                    model = "rhf";
                  }
              }

            pd = (OBPairData *) cset->GetData("DFTTYP");
            if(pd)
              {
                if(pd->GetValue() == "BLYP")
                  {
                    model = "b3lyp";
                  }
              }

            pd = (OBPairData *) cset->GetData("MPLEVL");
            if(pd)
              {
                if(pd->GetValue() == "2")
                  model = "mp2";
              }

            pd = (OBPairData *) cset->GetData("CCTYP");
            if(pd)
              {
                if(pd->GetValue() == "CCSD(T)")
                  model = "ccsd(t)";
              }

            pd = (OBPairData *) cset->GetData("RUNTYP");
            if(pd)
              {
                string value = pd->GetValue();
                if(value == "GRADIENT" || value == "HESSIAN" || value == "OPTIMIZE" || value == "SADPOINT")
                  {
                    method = pd->GetValue();
                    transform(method.begin(), method.end(), method.begin(), ::tolower);
                  }
              }

          }


        if(bset)
          {
            OBPairData *gbasis = (OBPairData *) bset->GetData("GBASIS");
            OBPairData *ngauss = (OBPairData *) bset->GetData("NGAUSS");

            if(gbasis)
              {
                string value = gbasis->GetValue();

                if( value == "am1" )
                  {
                    model = "am1";
                  }
                else if( value == "pm3" )
                  {
                    model = "pm3";
                  }
                else if(ngauss)
                  {
                    if(value == "STO")
                      {
                        basis.clear();
                        basis += "sto-";
                        basis += ngauss->GetValue();
                        basis += "g";
                      }
                    else if(ngauss->GetValue() == "3" || ngauss->GetValue() == "6")
                      {
                        basis.clear();
                        basis = ngauss->GetValue();
                        basis += "-";
                        basis += gbasis->GetValue().substr(1);
                        basis += "G(d)";
                      }
                  }
              }
          }
        OBPairData *nd = NULL;
        if(model != "")
          {
            nd = new OBPairData();
            nd->SetAttribute("model");
            nd->SetValue(model);
            nd->SetOrigin(fileformatInput);
            pmol->SetData(nd);
          }
        if(basis != "")
          {
            nd = new OBPairData();
            nd->SetAttribute("basis");
            nd->SetValue(basis);
            nd->SetOrigin(fileformatInput);
            pmol->SetData(nd);
          }
        if(method != "")
          {
            nd = new OBPairData();
            nd->SetAttribute("method");
            nd->SetValue(method);
            nd->SetOrigin(fileformatInput);
            pmol->SetData(nd);
          }

        /*
          cout << "model: " << model << endl;
          cout << "basis: " << basis << endl;
          cout << "method: " << method << endl;
        */
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges) {
      mol.SetPartialChargesPerceived();

      // Annotate that partial charges come from Mulliken
      OBPairData *dp = new OBPairData;
      dp->SetAttribute("PartialCharges");
      dp->SetValue("Mulliken");
      dp->SetOrigin(fileformatInput);
      mol.SetData(dp);
    }
    if (frequencies.size() != 0) { // we found some vibrations
      OBVibrationData *vd = new OBVibrationData;
      vd->SetData(displacements, frequencies, intensities);
      vd->SetOrigin(fileformatInput);
      mol.SetData(vd);
    }


    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////
  bool GAMESSInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    //const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    bool hasPartialCharges = false;
    string efragName; // used to save identifiers of EFRAG sections

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"$DATA") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE);	// title
            tokenize(vs,buffer);
            mol.SetTitle(buffer);
            ifs.getline(buffer,BUFF_SIZE);  // C1
            ifs.getline(buffer,BUFF_SIZE);
            while (strstr(buffer, "$END") == NULL)
              {
                tokenize(vs,buffer);
                if(vs.size() == 5) 
                  {
                    atom = mol.NewAtom();
                    atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
                    x = atof((char*)vs[2].c_str());
                    y = atof((char*)vs[3].c_str());
                    z = atof((char*)vs[4].c_str());
                    atom->SetVector(x,y,z);
                    vs[1].erase(vs[1].size() - 2, 2);
                  }

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
              }
          }
        if(strstr(buffer,"$EFRAG") != NULL)
          {
            while (strstr(buffer,"FRAGNAME") == NULL)
              {
                //read $EFRAG parameters
                tokenize(vs, buffer, "=");
                if (vs.size() > 1)
                  efragName = vs[1];
                if(!ifs.getline(buffer,BUFF_SIZE))
                  break;
              }
            while(strstr(buffer,"$END") == NULL)
              {
                tokenize(vs,buffer);
                if(vs.size() == 4)
                  {
                    atom = mol.NewAtom();
                    int atomicNum;
                    if( vs[0].substr(0,1) == "Z" || vs[0].substr(0,1) == "z" ) 
                      atomicNum=etab.GetAtomicNum(vs[0].substr(1,1).c_str());
                    else
                      atomicNum=etab.GetAtomicNum(vs[0].substr(0,1).c_str());
                    atom->SetAtomicNum(atomicNum);
                    x = atof((char*)vs[1].c_str());
                    y = atof((char*)vs[2].c_str());
                    z = atof((char*)vs[3].c_str());
                    atom->SetVector(x,y,z);

                    // Tag these atoms as part of a specific EFP fragment
                    OBPairData *dp = new OBPairData;
                    dp->SetAttribute("EFRAG");
                    dp->SetValue(efragName);
                    dp->SetOrigin(fileformatInput);
                    atom->SetData(dp);
                  }
                if(!ifs.getline(buffer,BUFF_SIZE))
                  break;
              }
          }
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    //mol.SetTitle(title);
    return(true);
  }


  bool GAMESSInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    //   unsigned int i;
    char buffer[BUFF_SIZE];

    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordsEnable = pConv->IsOption("k",OBConversion::GENOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);

    string defaultKeywords = " $CONTRL COORD=CART UNITS=ANGS $END";

    if(keywords)
      {
        defaultKeywords = keywords;
      }

    if (keywordsEnable)
      {

        OBSetData *gmsset = (OBSetData *)pmol->GetData("gamess");
        if(gmsset)
          {
            std::vector<OBGenericData*>::iterator i,j;

            for(i = gmsset->GetBegin(); i != gmsset->GetEnd(); i++)
              {
                OBSetData *cset = (OBSetData *)(*i);
                if(cset)
                  {
                    ofs << " $" << cset->GetAttribute();
                    for(j = cset->GetBegin(); j != cset->GetEnd(); j++)
                      {
                        OBPairData *pd = (OBPairData *) (*j);
                        if(pd)
                          {
                            ofs << " " << pd->GetAttribute() << "=" << pd->GetValue();
                          }
                      }
                    ofs << " $END"<<endl;
                  }
              }

          }
        else
          {
            ofs << "! Unable to translate keywords!" << endl;
            ofs << "! Defining default control keywords." << endl;
            ofs << defaultKeywords << endl;
          }

      }
    else if (keywordFile)
      {
        ifstream kfstream(keywordFile);
        string keyBuffer;
        if (kfstream)
          {
            while (getline(kfstream, keyBuffer))
              ofs << keyBuffer << endl;
          }
      }
    else
      {
        ofs << defaultKeywords << endl;
      }

    ofs << endl << " $DATA" << endl;
    ofs << mol.GetTitle() << endl;
    if (!mol.HasData(OBGenericDataType::SymmetryData))
      ofs << "C1" << endl;
    else
      {
        // \todo needs to be updated for point group symmetry recognition
        //   particularly for output of the symmetry elements
        //   and any necessary rotation for frame of reference for GAMESS
        ofs << "Put symmetry info here" << endl << endl;
      }

    //  OBAtom *atom;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%-3s %4d.0    %14.10f  %14.10f  %14.10f ",
                 etab.GetSymbol(atom->GetAtomicNum()),
                 atom->GetAtomicNum(),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << endl;
      }

    ofs << " $END" << endl << endl << endl;
    return(true);
  }

} //namespace OpenBabel
