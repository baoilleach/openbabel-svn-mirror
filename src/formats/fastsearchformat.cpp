/**********************************************************************
fastsearchformat.cpp: Preparation and searching of fingerprint-based index files
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "babelconfig.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include "mol.h"
#include "obconversion.h"
#include "fingerprint.h"

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
	}
	
	virtual const char* Description() //required
	{ return
"FastSearching\n \
Uses molecular fingerprints in an index file.\n \
Writing to the fs format makes an index (a very slow process)\n \
  babel datafile.xxx index.fs\n \
Reading from the fs format does a fast search for:\n \
  Substructure\n \
    babel index.fs -sSMILES outfile.yyy   or\n \
    babel datafile.xxx -ifs -sSMILES outfile.yyy\n \
  Molecular similarity based on Tanimoto coefficient\n \
    babel index.fs -sSMILES outfile.yyy -t0.7  (Tanimoto >0.7)\n \
    babel index.fs -sSMILES outfile.yyy -t15   (best 15 molecules)\n \
  The structure spec can be a molecule from a file: -Spatternfile.zzz\n \
\n \
Write Options (when making index) e.g. -xfFP3 \n \
 f# Fingerprint type\n \
 N# Fold fingerprint to # bits\n \
 u  Update an existing index\n\n \
Read Options (when searching) e.g. -at0.7\n \
 t# Do similarity search: #mols or # as min Tanimoto\n \
 a  Add Tanimoto coeff to title in similarity search\n \
 l# Maximum number of candidates. Default<4000>\n \
 S\"filename\"  Structure spec in a file (rather than as SMARTS):\n\n \
";
	};

	virtual unsigned int Flags(){return READONEONLY | WRITEBINARY;};

public:
	virtual bool ReadChemObject(OBConversion* pConv);
	virtual bool WriteChemObject(OBConversion* pConv);

private:
	///big data structure which will remain in memory after it is loaded
	//until the program ends.
	FastSearch fs;
	FastSearchIndexer* fsi;
	streampos LastSeekpos; //used during update
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

	OBMol patternMol;
	stringstream smiles(stringstream::out);		
	ifstream patternstream;
	OBConversion PatternConv(&patternstream,&smiles);

	//Convert the SMARTS string to an OBMol
	const char* p = pConv->IsOption("s",OBConversion::GENOPTIONS);
	string txt;
	if(p) 
	{
		txt=p;
		stringstream smarts(txt, stringstream::in);		
		OBConversion Convsm(&smarts);
		if(!Convsm.SetInFormat("smi")) return false;
		Convsm.Read(&patternMol);
		
		//erase -s option in GeneralOptions since it will be rewritten
		pConv->RemoveOption("s",OBConversion::GENOPTIONS);
	}

	// or Make OBMol from file in -S option or -aS option	
	p = pConv->IsOption("S",OBConversion::GENOPTIONS);
	if(!p)
		p = pConv->IsOption("S",OBConversion::INOPTIONS);//for GUI mainly
	if(p && patternMol.Empty())
	{
		txt=p;
		string::size_type pos = txt.find_last_of('.');
		if(pos==string::npos)
		{
		  obErrorLog.ThrowError(__FUNCTION__, "Filename of pattern molecule in -S option must have an extension", obError);
			return false;
		}
		patternstream.open(txt.c_str());
		if(!patternstream)
		{
#ifdef HAVE_SSTREAM
		  stringstream errorMsg;
#else
		  strstream errorMsg;
#endif
		  
		  errorMsg << "Cannot open " << txt << endl;
		  obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
		  return false;
		}

		PatternConv.SetOneObjectOnly();
		if(PatternConv.SetInFormat(txt.substr(pos+1).c_str()))
			PatternConv.Read(&patternMol);
	}

	if(patternMol.Empty())
	{
	  obErrorLog.ThrowError(__FUNCTION__, "Cannot derive a molecule from the -s or -S options", obWarning);
		return false;
	}
	patternMol.ConvertDativeBonds();//use standard form for dative bonds

	//Convert to SMILES and generate a -s option for use in the final filtering
	if(!PatternConv.SetOutFormat("smi"))
		return false;
	PatternConv.Write(&patternMol);
	//remove name to leave smiles string
	string smilesstr(smiles.str());
	string::size_type pos = smilesstr.find(' ');
	if(pos!=string::npos)
		smilesstr = smilesstr.substr(0,pos);
	pConv->AddOption("s", OBConversion::GENOPTIONS, smilesstr.c_str());
	
	//Derive index name
	string indexname = pConv->GetInFilename();
	pos=indexname.find_last_of('.');
	if(pos!=string::npos)
	{
		indexname.erase(pos);
		indexname += ".fs";
	}

	//Have to open input stream again because needs to be in binary mode
	ifstream ifs;
#ifdef HAVE_SSTREAM
    stringstream errorMsg;
#else
    strstream errorMsg;
#endif
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
	pConv->AddOption("b",OBConversion::GENOPTIONS);


	//Now do searching
	p = pConv->IsOption("t",OBConversion::INOPTIONS);
	if(p)
	{
		//Do a similarity search
		multimap<double, unsigned int> SeekposMap;
		txt=p;
		if(txt.find('.')==string::npos)
		{
			//Finds n molecules with largest Tanimoto
			int n = atoi(p);
			fs.FindSimilar(&patternMol, SeekposMap, n);
		}
		else
		{
			//Finds molecules with Tanimoto > MinTani
			double MinTani = atof(txt.c_str());
			fs.FindSimilar(&patternMol, SeekposMap, MinTani);
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
		//Do a substructure search
		int MaxCandidates = 4000;
		p = pConv->IsOption("l",OBConversion::INOPTIONS);
		if(p && atoi(p))
			MaxCandidates = atoi(p);
		
		vector<unsigned int> SeekPositions;
		fs.Find(&patternMol, SeekPositions, MaxCandidates);
		clog << SeekPositions.size() << " candidates from fingerprint search phase" << endl;

		//Output the candidate molecules 
		//filtering through s filter, unless the fingerprint type does not require it
		if(fs.GetFingerprint()->Flags() & OBFingerprint::FPT_UNIQUEBITS)
			pConv->RemoveOption("s",OBConversion::GENOPTIONS);

		vector<unsigned int>::iterator itr;
		for(itr=SeekPositions.begin();itr!=SeekPositions.end();itr++)
		{
			datastream.seekg(*itr);
			//	datastream.seekg(*itr - datastream.tellg(), ios_base::cur); //Avoid retrieving start

			//debugging kludge to output all candidates directly
			if(pConv->IsOption("c",OBConversion::GENOPTIONS))
			{
				string ln;
				getline(datastream,ln);
				datastream.seekg(*itr);
				*pConv->GetOutStream() << "** " << ln << endl;
			}
			pConv->SetOneObjectOnly();
			if(itr+1 != SeekPositions.end())
				pConv->SetMoreFilesToCome();//so that not seen as last on output 
			pConv->Convert(NULL,NULL);
		}
	}
	return false;	//To finish	
}

/////////////////////////////////////////////////////
bool FastSearchFormat::WriteChemObject(OBConversion* pConv)
{
	//Prepares or updates an index file

	bool update = pConv->IsOption("u")!=NULL;
	string mes("prepare an");
	if(update)
		mes = "update the";
	if(fsi==NULL)
		clog << "This will " << mes << " index of " << pConv->GetInFilename()
		     <<  " and may take some time..." << flush;
			
	OBStopwatch sw;
//	sw.Start(); //seems stupid but makes gcc-4 happy
	
	std::string auditMsg = "OpenBabel::Write fastsearch index ";
	std::string description(Description());
        auditMsg += description.substr( 0, description.find('\n') );
        obErrorLog.ThrowError(__FUNCTION__,auditMsg,obAuditMsg);

	ostream* pOs = pConv->GetOutStream();// with named index it is already open
	bool NewOstreamUsed=false;
	if(fsi==NULL)
	{
		//First pass sets up FastSearchIndexer object
		sw.Start();
		
		FptIndex* pidx; //used with update

		if(pOs==&cout)
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
			#ifdef HAVE_SSTREAM
				stringstream errorMsg;
			#else
				strstream errorMsg;
			#endif
				errorMsg << "Trouble opening or reading " << indexname << endl;
				obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
				return false;
			}
			NewOstreamUsed=true;
		}
		else // not cout
		{
			if(update)
			{	obErrorLog.ThrowError(__FUNCTION__,
"Currently, updating	can only be done with index files that \
have the same name as the datafile. Use the form:\n \
	babel datafile.xxx -ofs -xu", obError);
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
		unsigned int pos = datafilename.find_last_of("/\\");
		if(pos!=string::npos)
			datafilename=datafilename.substr(pos+1);

		if(update)
		{
			fsi = new FastSearchIndexer(pidx,pOs);//using existing index

			//Seek to position in datafile of last of old objects
			LastSeekpos = *(pidx->seekdata.end()-1);
			pConv->GetInStream()->seekg(LastSeekpos);
		}
		else
			fsi = new FastSearchIndexer(datafilename, pOs, fpid, nbits);
		
		obErrorLog.StopLogging();
	}

	//All passes provide an object for indexing
	OBBase* pOb = pConv->GetChemObject();
	OBMol* pmol = dynamic_cast<OBMol*> (pOb);
	if(pmol)
		pmol->ConvertDativeBonds();//use standard form for dative bonds
	
	streampos seekpos = pConv->GetInPos();
	if(!update || seekpos>LastSeekpos) 
		fsi->Add(pOb, seekpos );
	else
		//Don't index old objects during update. Don't increment pConv->Index.
		pConv->SetOutputIndex(pConv->GetOutputIndex()-1);

	if(pConv->IsLast())
	{
		//Last pass 
		if(NewOstreamUsed)
			delete pOs;

		delete fsi; //saves index file
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

}//Openbabel

//! \file fastsearchformat.cpp
//! \brief Preparation and searching of fingerprint-based index files
