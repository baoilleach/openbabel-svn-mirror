/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"

namespace OpenBabel {
  bool SafeOpen(std::ifstream &fs, char *filename);
  bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

void GenerateSmartsReference();

int main(int argc,char *argv[])
{
  if (argc != 1) {
    if (strncmp(argv[1], "-g", 2))
      {
	cout << "Usage: smartstest" << endl;
	cout << "   Tests Open Babel SMILES/SMARTS pattern matching." << endl;
	return 0;
      } else {
	GenerateSmartsReference();
	return 0;
      }
  }

  cout << endl << "Testing SMARTS...  " << endl;

#ifdef TESTDATADIR
  string testdatadir = TESTDATADIR;
  string smarts_file = testdatadir + "smartstest.txt";
  string results_file = testdatadir + "smartsresults.txt";
  string smilestypes_file = testdatadir + "attype.00.smi";
#else
  string smarts_file = "smartstest.txt";
  string results_file = "smartsresults.txt";
  string smilestypes_file = "attype.00.smi";
#endif

  std::ifstream ifs;
  if (!SafeOpen(ifs, (char*)smarts_file.c_str())){
    return -1; // test failed
  }

  //read in the SMARTS test patterns
  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);) 
    {
      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer)) vsp.push_back(sp);
      else                  delete sp;
    }
  ifs.close();
  
  std::ifstream rifs;
  if (!SafeOpen(rifs, (char*)results_file.c_str())) {
    return -1; // test failed
  }
  unsigned int npats;
  rifs.getline(buffer,BUFF_SIZE);
  sscanf(buffer,"%d %*s",&npats);

  //make sure the number of SMARTS patterns is the same as in the 
  //reference data
  if (npats != vsp.size())
    {
      ThrowError("Correct number of patterns not read in");
      sprintf(buffer,"Read in %d, expected %d", vsp.size(), npats);
      ThrowError(buffer);
      return -1; // test failed
    }

  std::ifstream mifs;
  if (!SafeOpen(mifs, (char*)smilestypes_file.c_str())) {
    return -1; // test failed
  }

  unsigned int k;
  unsigned int res_line = 0;
  OBMol mol(SMI,SMI);
  vector<string> vs;
  vector<OBSmartsPattern*>::iterator i;
  vector<vector<int> > mlist;
  OBFileFormat ff;

  //read in molecules, match SMARTS, and compare results to reference data
  for (;mifs;)
    {
      mol.Clear();
      ff.ReadMolecule(mifs, mol);
      if (mol.Empty()) continue;
      
      for (i = vsp.begin();i != vsp.end();i++)
	{
	  if (!rifs.getline(buffer,BUFF_SIZE))
	    {
	      ThrowError("error reading reference data");
	      return -1; // test failed
	    }
	  res_line++;

	  tokenize(vs,buffer);
	  (*i)->Match(mol);
	  mlist = (*i)->GetMapList();
	  if (mlist.size() != vs.size())
	    {
	      ThrowError("number of matches different than reference");
	      cerr << "expected " << vs.size() << " matches, found " 
		   << mlist.size() << endl;
	      cerr << "error with molecule " << mol.GetTitle();
	      cerr << " on pattern " << (*i)->GetSMARTS() << endl;
	      if (mlist.size())
		cerr << "first match: atom #" << mlist[0][0] << endl;
	      return -1; // test failed
	    }
	  
	  if (mlist.size())
	    {
	      for (k = 0;k < vs.size();k++)
		if (atoi((char*)vs[k].c_str()) != mlist[k][0])
		{
		  ThrowError("matching atom numbers different than reference");
		  cerr << "expected " << vs[k] << " but found " 
		       << mlist[k][0] << endl;
		  ThrowError((char*)mol.GetTitle());
		  ThrowError((*i)->GetSMARTS());
		  return -1; // test failed
		}
	    }
	}
    }

  // Passed Test
  return 0;
}

void GenerateSmartsReference()
{
  std::ifstream ifs;
  if (!SafeOpen(ifs,"smartstest.txt")) return;

  char buffer[BUFF_SIZE];
  vector<OBSmartsPattern*> vsp;
  for (;ifs.getline(buffer,BUFF_SIZE);)
    {
      OBSmartsPattern *sp = new OBSmartsPattern;

      if (sp->Init(buffer)) vsp.push_back(sp);
      else                  delete sp;
    }

  std::ofstream ofs;
  if (!SafeOpen(ofs,"smartsresults.txt")) return;
    
  ofs << vsp.size() << " patterns" << endl;
  std::ifstream mifs;
  if (!SafeOpen(mifs,"attype.00.smi")) return;

  vector<int> vm;
  vector<vector<int> > mlist;
  vector<vector<int> >::iterator j;
  vector<OBSmartsPattern*>::iterator i;
  OBMol mol(SMI,SMI);
  OBFileFormat ff;

  for (;mifs;)
    {
      mol.Clear();
      ff.ReadMolecule(mifs, mol);

      if (mol.Empty()) continue;
      for (i = vsp.begin();i != vsp.end();i++)
	{
	  (*i)->Match(mol);
	  mlist = (*i)->GetMapList();
	  for (j = mlist.begin();j != mlist.end();j++)
	    {
	     sprintf(buffer,"%3d",*(j->begin()));
	     ofs << buffer;
	    }
	  ofs << endl; 
	}
    }
  

  ThrowError("SMARTS test results written successfully");
}

