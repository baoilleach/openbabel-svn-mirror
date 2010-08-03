/**********************************************************************
fastsearchformat.cpp: Preparation and searching of fingerprint-based index files
Copyright (C) 2005-2006 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel {

/// \brief Prepares and searches of fingerprint-based index files. See FastSearch class for details
class FastSearchFormat : public OBFormat
{
public:
//Register this format type ID
FastSearchFormat() : fsi(NULL) 
{
  OBConversion::RegisterFormat("fs",this);
  //Specify the number of option taken by options
  OBConversion::RegisterOptionParam("S", this, 1, OBConversion::GENOPTIONS);
  OBConversion::RegisterOptionParam("S", this, 1, OBConversion::INOPTIONS);
  OBConversion::RegisterOptionParam("f", this, 1);
  OBConversion::RegisterOptionParam("N", this, 1);
  OBConversion::RegisterOptionParam("u", this, 0);
  OBConversion::RegisterOptionParam("t", this, 1, OBConversion::INOPTIONS);
  OBConversion::RegisterOptionParam("l", this, 1, OBConversion::INOPTIONS);
  OBConversion::RegisterOptionParam("a", this, 0, OBConversion::INOPTIONS);
  OBConversion::RegisterOptionParam("e", this, 0, OBConversion::INOPTIONS);
}

virtual const char* Description() //required
{ return
  "Fastsearch format\n"
  "Fingerprint-aided substructure and similarity searching\n\n"

  "Writing to the fs format makes an index of a multi-molecule datafile:\n"
  "      babel dataset.sdf -ofs\n"
  "This prepares an index dataset.fs with default parameters, and is slow\n"
  "(~30 minutes for a 250,000 molecule file).\n\n"

  "Reading from the fs format searches and is much faster, a few seconds,\n"
  "and so can be done interactively.\n\n"
  "The search target is the parameter of the -s option and can be\n"
  "slightly extended SMILES (with [#n] atoms and ~ bonds) or\n"
  "the name of a file containing a molecule.\n\n"

  "Several types of searches are possible:\n"
  "- Identical molecule::\n"
  "      babel index.fs outfile.yyy -s SMILES exact\n"
  "- Substructure::\n"
  "      babel index.fs outfile.yyy  -s SMILES   or\n"
  "      babel index.fs outfile.yyy  -s filename.xxx\n"
  "   where xxx is a format id known to OpenBabel, e.g. sdf\n"
  "- Molecular similarity based on Tanimoto coefficient::\n"
  "      babel index.fs outfile.yyy -at15  -sSMILES (best 15 molecules)\n"
  "      babel index.fs outfile.yyy -at0.7 -sSMILES  (Tanimoto >0.7)\n"
  "      babel index.fs outfile.yyy -at0.7,0.9 -sSMILES\n"
  "   (Tanimoto >0.7 && Tanimoto < 0.9)\n"
  "The datafile plus the -ifs option can be used instead of the index file.\n\n"

  "Write Options (when making index) e.g. -xfFP3\n"
  " f# Fingerprint type\n"
  "     If not specified, the default fingerprint (currently FP2) is used\n"
  " N# Fold fingerprint to # bits\n"
  " u  Update an existing index\n\n"

  "Read Options (when searching) e.g. -at0.7\n"
  " t# Do similarity search: #mols or # as min Tanimoto\n"
  " a  Add Tanimoto coeff to title in similarity search\n"
  " l# Maximum number of candidates. Default<4000>\n"
  " e  Exact match\n"
  "     Alternative to using exact in -s parameter, see above\n"  
  " n  No further SMARTS filtering after fingerprint phase\n\n"
  ;
};

  virtual unsigned int Flags(){return READBINARY | READONEONLY | WRITEBINARY;};

  public:
    virtual bool ReadChemObject(OBConversion* pConv);
    virtual bool WriteChemObject(OBConversion* pConv);

  private:
    bool ObtainTarget(OBConversion* pConv, std::vector<OBMol>& patternMols, const std::string& indexname);
    void AddPattern(vector<OBMol>& patternMols, OBMol patternMol, int idx);

  private:
    ///big data structure which will remain in memory after it is loaded
    //until the program ends.
    FastSearch fs;
    FastSearchIndexer* fsi;
    streampos LastSeekpos; //used during update
    OBStopwatch sw; //used when preparing index
    int nmols; //number mols in data file
  };

  ///////////////////////////////////////////////////////////////
  //Make an instance of the format class
  FastSearchFormat theFastSearchFormat;

  ///////////////////////////////////////////////////////////////
  bool FastSearchFormat::ReadChemObject(OBConversion* pConv)
  {
    //Searches index file for structural matches
    //This function is called only once per search

    std::string auditMsg = "OpenBabel::Read fastsearch index ";
    std::string description(Description());
    auditMsg += description.substr(0,description.find('\n'));
    obErrorLog.ThrowError(__FUNCTION__,
                          auditMsg,
                          obAuditMsg);

    //Derive index name
    string indexname = pConv->GetInFilename();
    string::size_type pos=indexname.find_last_of('.');
    if(pos!=string::npos)
      {
        indexname.erase(pos);
        indexname += ".fs";
      }

    //Have to open input stream again because needs to be in binary mode
    ifstream ifs;
    stringstream errorMsg;
    if(!indexname.empty())
      ifs.open(indexname.c_str(),ios::binary);
    if(!ifs)
      {
        errorMsg << "Couldn't open " << indexname << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }

    string datafilename = fs.ReadIndex(&ifs);
    if(datafilename.empty())
      {
        errorMsg << "Difficulty reading from index " << indexname << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    
    vector<OBMol> patternMols;
    if(!ObtainTarget(pConv, patternMols, indexname))
      return false;

    bool exactmatch = pConv->IsOption("e",OBConversion::INOPTIONS)!=NULL;// -ae option
    
    //Open the datafile and put it in pConv
    //datafile name derived from index file probably won't have a file path
    //but indexname may. Derive a full datafile name
    string path;
    pos = indexname.find_last_of("/\\");
    if(pos==string::npos)
      path = datafilename;
    else
      path = indexname.substr(0,pos+1) + datafilename;
	
    ifstream datastream(path.c_str());
    if(!datastream)
      {
        errorMsg << "Difficulty opening " << path << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
    pConv->SetInStream(&datastream);
	
    //Input format is currently fs; set it appropriately
    if(!pConv->SetInAndOutFormats(pConv->FormatFromExt(datafilename.c_str()),pConv->GetOutFormat()))
			return false;

    //Now do searching
    const char* p = pConv->IsOption("t",OBConversion::INOPTIONS);
    if(p)
      {
        //Do a similarity search
        multimap<double, unsigned int> SeekposMap;
        string txt=p;
        if(txt.find('.')==string::npos)
          {
            //Finds n molecules with largest Tanimoto
            int n = atoi(p);
            fs.FindSimilar(&patternMols[0], SeekposMap, n);
          }
        else
          {
            //Finds molecules with Tanimoto > MinTani
            double MaxTani = 1.1;
            int pos = txt.find( ',' );
            if( pos != string::npos ) {
              MaxTani = atof( txt.substr( pos + 1 ).c_str() );
            }
            double MinTani = atof( txt.substr( 0, pos ).c_str() );
//            if(doSubset)
//              fs.FindSubset(SeekposMap, MinTani);
//            else
            fs.FindSimilar(&patternMols[0], SeekposMap, MinTani, MaxTani);
          }
		
        //Don't want to filter through SMARTS filter
        pConv->RemoveOption("s", OBConversion::GENOPTIONS);
		
        multimap<double, unsigned int>::reverse_iterator itr;
        for(itr=SeekposMap.rbegin();itr!=SeekposMap.rend();++itr)
          {
            datastream.seekg(itr->second);

            if(pConv->IsOption("a", OBConversion::INOPTIONS))
              {
                //Adds Tanimoto coeff to title
                //First remove any previous value
                pConv->RemoveOption("addtotitle", OBConversion::GENOPTIONS);
                stringstream ss;
                ss << " " << itr->first;
                pConv->AddOption("addtotitle",OBConversion::GENOPTIONS, ss.str().c_str());
			
              }
            pConv->SetOneObjectOnly();
            if(itr != --SeekposMap.rend())
              pConv->SetMoreFilesToCome();//so that not seen as last on output 
            pConv->Convert(NULL,NULL);
          }
      }

    else
    {
      //Structure search
      int MaxCandidates = 4000;
      p = pConv->IsOption("l",OBConversion::INOPTIONS);
      if(p && atoi(p))
        MaxCandidates = atoi(p);
	
      vector<unsigned int> SeekPositions;

      if(exactmatch)
      {
        //Find mols where all fingerprint bits are the same as the target
        fs.FindMatch(&patternMols[0], SeekPositions, MaxCandidates);
        // ensure that SMARTS filter in transform.cpp looks only for an exact match
        // by setting an option with the number of heavy atoms in the pattern mol included.
        stringstream ss;
        ss << patternMols[0].NumHvyAtoms();
        pConv->AddOption("exactmatch", OBConversion::GENOPTIONS, ss.str().c_str());
      }

      else
      {
        //Do a substructure search for each target
        vector<OBMol>::iterator iter;
        for(iter=patternMols.begin();iter!=patternMols.end();++iter)
          fs.Find(&*iter, SeekPositions, MaxCandidates);
        clog << SeekPositions.size() << " candidates from fingerprint search phase" << endl;
      }

      vector<unsigned int>::iterator seekitr, 
          begin = SeekPositions.begin(), end = SeekPositions.end();

      if(patternMols.size()>1)//only sort and elininate duplicates if necessary
      {
        sort(begin, end);
        end = unique(begin, end); //removed duplicates are after new end
      }

      //Output the candidate molecules, filtering through s filter, unless the
      //fingerprint type does not require it, or it was not requested
      if(fs.GetFingerprint()->Flags() & OBFingerprint::FPT_UNIQUEBITS
                    || pConv->IsOption("n",OBConversion::INOPTIONS) )
        pConv->RemoveOption("s",OBConversion::GENOPTIONS);

      pConv->SetLast(false);
      for(seekitr=begin; seekitr!=end; ++seekitr)
      {
        datastream.seekg(*seekitr);
        if(!pConv->GetInFormat()->ReadChemObject(pConv))
          return false;
        pConv->SetFirstInput(false); //needed for OpSort
      }
    }
    return false;	//To finish	
  }

  /////////////////////////////////////////////////////
  bool FastSearchFormat::WriteChemObject(OBConversion* pConv)
  {
    //Prepares or updates an index file. Called for each molecule indexed
    bool update = pConv->IsOption("u")!=NULL;

    static ostream* pOs;
    static bool NewOstreamUsed;
    if(fsi==NULL)
      {
        //First pass sets up FastSearchIndexer object
        pOs = pConv->GetOutStream();// with named index it is already open
        NewOstreamUsed=false;
        string mes("prepare an");
        if(update)
          mes = "update the";
        clog << "This will " << mes << " index of " << pConv->GetInFilename()
             <<  " and may take some time..." << flush;

        if(!pConv->IsLastFile())
        {
          obErrorLog.ThrowError(__FUNCTION__,
            "There should not be multiple input files. A .fs file is an index of a single datafile.",
            obError);
          return false;
        }

        std::string auditMsg = "OpenBabel::Write fastsearch index ";
        std::string description(Description());
        auditMsg += description.substr( 0, description.find('\n') );
        obErrorLog.ThrowError(__FUNCTION__,auditMsg,obAuditMsg);

        FptIndex* pidx; //used with update

        //if(pOs==&cout) did not work with GUI
        if(!dynamic_cast<ofstream*>(pOs))
          {
            //No index filename specified
            //Derive index name from datafile name
            string indexname=pConv->GetInFilename();
            string::size_type pos=indexname.find_last_of('.');
            if(pos!=string::npos)
              indexname.erase(pos);
            indexname += ".fs";

            bool idxok=true;
            if(update)
              {
                LastSeekpos = 0;

                //Read in existing index
                idxok=false;
                ifstream ifs(indexname.c_str(),ifstream::binary);
                if(ifs.good())
                  {
                    pidx = new FptIndex;
                    idxok = pidx->Read(&ifs);
                  }
              }//ifs closed here

            pOs = new ofstream(indexname.c_str(),ofstream::binary);

            if(!pOs->good() || !idxok)
              {
                stringstream errorMsg;
                errorMsg << "Trouble opening or reading " << indexname << endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
                static_cast<ofstream *>(pOs)->close(); // close the file before quitting
                delete pOs;
                return false;
              }
            NewOstreamUsed=true;
          }
        else // not cout
          {
            if(update)
              {
                obErrorLog.ThrowError(__FUNCTION__,
                  "Currently, updating is only done on index files that"
                  "have the same name as the datafile.\n"
                  "Do not specify an output file; use the form:\n"
	                "   babel datafile.xxx -ofs -xu", obError);
                return false;
              }
          }

        int nbits = 0;
        const char* p = pConv->IsOption("N");
        if(p)
          nbits = atoi(p);

        string fpid; //fingerprint type
        p=pConv->IsOption("f");
        if(p)
          fpid=p;

        //Prepare name without path
        string datafilename = pConv->GetInFilename();
        if(datafilename.empty())
          {
            obErrorLog.ThrowError(__FUNCTION__, "No datafile!", obError);
            return false;
          }
        string::size_type pos = datafilename.find_last_of("/\\");
        if(pos!=string::npos)
          datafilename=datafilename.substr(pos+1);

        nmols = pConv->NumInputObjects();
        if(nmols>0)
          clog << "\nIt contains " << nmols << " molecules" << flush;

        sw.Start();
		
        if(update)
          {
            fsi = new FastSearchIndexer(pidx, pOs, nmols);//using existing index

            //Seek to position in datafile of last of old objects
            LastSeekpos = *(pidx->seekdata.end()-1);
            pConv->GetInStream()->seekg(LastSeekpos);
          }
        else
          fsi = new FastSearchIndexer(datafilename, pOs, fpid, nbits, nmols);
		
        obErrorLog.StopLogging();
      }

    //All passes provide an object for indexing
    OBBase* pOb = pConv->GetChemObject();
    OBMol* pmol = dynamic_cast<OBMol*> (pOb);
    if(pmol)
      pmol->ConvertDativeBonds();//use standard form for dative bonds
	
    streampos seekpos = pConv->GetInPos();
    if(!update || seekpos>LastSeekpos) 
    {
      fsi->Add(pOb, seekpos );
      if(pConv->GetOutputIndex()==400 && nmols>1000)
      {
        clog << " Estimated completion time ";
        double secs = sw.Elapsed() * nmols / 400; //
        streamsize op = clog.precision(0);
        if(secs>150)
          clog << secs/60 << " minutes" << endl;
    else
          clog << secs << " seconds" << endl;
        clog.precision(op);
      }
    }
    else
      //Don't index old objects during update. Don't increment pConv->Index.
      pConv->SetOutputIndex(pConv->GetOutputIndex()-1);

    if(pConv->IsLast())
      {
        //Last pass 
        delete fsi; //saves index file
        if(NewOstreamUsed)
          delete pOs;

        //return to starting conditions
        fsi=NULL;

        obErrorLog.StartLogging();

        double secs = sw.Elapsed();
        if(secs>150)
          clog << "\n It took " << secs/60 << " minutes" << endl;
        else
          clog << "\n It took " << secs << " seconds" << endl;
      }
    delete pOb;
    return true;
  }

///////////////////////////////////////////////////////////////
  bool FastSearchFormat::ObtainTarget(OBConversion* pConv, vector<OBMol>& patternMols, const string& indexname)
  {
    //Obtains an OBMol from:
    // the filename in the -s option or
    // the SMARTS string in the -s option or
    // by converting the file in the -S or -aS options (deprecated).
    // If there is no -s -S or -aS option, information on the index file is displayed.

    OBMol patternMol;
    patternMol.SetIsPatternStructure();

    const char* p = pConv->IsOption("s",OBConversion::GENOPTIONS);
    
    bool OldSOption=false;
    //If no -s option, make OBMol from file in -S option or -aS option (both deprecated)
    if(!p)
    {
      p = pConv->IsOption("S",OBConversion::GENOPTIONS);
      if(!p)
        p = pConv->IsOption("S",OBConversion::INOPTIONS);//for GUI mainly
      OldSOption = true; 
    }
    if(p) 
    {
      vector<string> vec;
      tokenize(vec, p);
      
      //ignore leading ~ (not relevant to fastsearch)
      if(vec[0][0]=='~')
        vec[0].erase(0,1);

      if(vec.size()>1 && vec[1]=="exact")
        pConv->AddOption("e", OBConversion::INOPTIONS);

      OBConversion patternConv;
      OBFormat* pFormat;
      //Interpret as a filename if possible
      string filename(p);
      if( filename.empty() ||
          filename.find('.')==string::npos ||
          !(pFormat = patternConv.FormatFromExt(filename.c_str())) ||
          !patternConv.SetInFormat(pFormat) ||
          !patternConv.ReadFile(&patternMol, filename) ||
          patternMol.NumAtoms()==0)
        //if false, have a valid patternMol from a file
      {
        //is SMARTS/SMILES
        //Replace e.g. [#6] in SMARTS by C so that it can be converted as SMILES
        //for the fingerprint phase, but allow more generality in the SMARTS phase.
        string& txt =vec [0];
        for(;;)
        {
          string::size_type pos1, pos2;
          pos1 = txt.find("[#");
          if(pos1==string::npos)
            break;
          pos2 = txt.find(']');
          int atno;
          if(pos2!=string::npos &&  (atno = atoi(txt.substr(pos1+2, pos2-pos1-2).c_str())) && atno>0)
            txt.replace(pos1, pos2-pos1+1, etab.GetSymbol(atno));
          else
          {
            obErrorLog.ThrowError(__FUNCTION__,"Ill-formed [#n] atom in SMARTS", obError);
            return false;
          }
        }

        bool hasTildeBond;
        if(hasTildeBond = (txt.find('~')!=string::npos))
        {
          //Find ~ bonds and make versions of query molecule with a single and aromatic bonds
          //To avoid having to parse the SMILES here, replace ~ by $ (quadruple bond)
          //and then replace this in patternMol. Check first that there are no $ already
          //Sadly, isocynanides may have $ bonds.
          if(txt.find('$')!=string::npos)
          {
            obErrorLog.ThrowError(__FUNCTION__,
              "Cannot use ~ bonds in patterns with $ (quadruple) bonds.)", obError);
            return false;
          }
          replace(txt.begin(),txt.end(), '~' , '$');
        }

        //read as standard SMILES
        patternConv.SetInFormat("smi");
        patternConv.ReadString(&patternMol, vec[0]);

        if(hasTildeBond)
        {
          AddPattern(patternMols, patternMol, 0); //recursively add all combinations of tilde bond values
          return true;
        }
      }
    }

    if(OldSOption) //only when using deprecated -S and -aS options
    {
      //make -s option for later SMARTS test
      OBConversion conv;
      if(conv.SetOutFormat("smi"))
      {
        string optiontext = conv.WriteString(&patternMol, true);
        pConv->AddOption("s", OBConversion::GENOPTIONS, optiontext.c_str());
      }
    }

    if(!p)
    {
      //neither -s or -S options provided. Output info rather than doing search
      const FptIndexHeader& header = fs.GetIndexHeader();
      string id(header.fpid);
      if(id.empty())
        id = "default";
      clog << indexname << " is an index of\n " << header.datafilename 
           << ".\n It contains " << header.nEntries 
           << " molecules. The fingerprint type is " << id << " with "
           << OBFingerprint::Getbitsperint() * header.words << " bits.\n"
           << "Typical usage for a substructure search:\n"
           << "babel indexfile.fs -osmi -sSMILES\n"
           << "(-s option in GUI is 'Convert only molecules matching SMARTS')" << endl;
      return false;
    }

    patternMols.push_back(patternMol);
    return true;
  }

  void FastSearchFormat::AddPattern(vector<OBMol>& patternMols, OBMol patternMol, int idx)
  {
    //Recursive function to generate all combinations of aromatic/single bonds for each tilde bond
    //Copying an OBMol, which happens when adding it to a vector, kekulizes it,
    // changing aromatic (bo=5) bonds. So set order after adding. Should work here,
    // but is dangerous if the vector needs to be reallocated.

    if(idx>=patternMol.NumBonds())
      return;
    if(patternMol.GetBond(idx)->GetBO()==4)
    {
      patternMol.GetBond(idx)->SetBO(1);
      patternMols.push_back(patternMol);
      AddPattern(patternMols, patternMol,idx+1);

      patternMols.push_back(patternMol);
      patternMols.back().GetBond(idx)->SetBO(5);
    }
    AddPattern(patternMols, patternMol,idx+1);
  }

  /* Accept ~ bonds. Need to generate two PatternMols for each '~'
  i.e. '-' or nothing for single and ':' for aromatic
  or 2^n patterns for n bonds. SMILES format will accept : but cannot
  provide an isolated aromatic bond. So need to edit OBMol. Retain the atom
  indices of the bond's atoms and change the bond's order.
  ObtainTarget() will return a vector of OBMol and the Find() in L260 will be done
  for each. All fs matches will go into SeekPositions. At the end this
  will be sorted and duplicates removed with unique.  
  */

}//Openbabel

//! \file fastsearchformat.cpp
//! \brief Preparation and searching of fingerprint-based index files
