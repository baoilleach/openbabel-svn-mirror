/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
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
#include "babelconfig.h"
#ifdef WIN32
#pragma warning (disable : 4786)
#pragma warning (disable : 4251) //
#endif
#include <string>
#include <iomanip>
#include "obmolecformat.h"
#include "kinetics.h"

using namespace std;
namespace OpenBabel
{
class ThermoFormat : public OBMoleculeFormat
{
public:
  ThermoFormat()
  {
      OBConversion::RegisterFormat("therm",this);
      OBConversion::RegisterFormat("tdd",this);
  }

  virtual const char* Description()
  {
      return
"Thermo format\n \
Reads and writes old-style NASA polynomials in original fixed format.\n";
  }

  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pReact, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pReact, OBConversion* pConv);
};

////////////////////////////////////////////
//Make an instance of the format class
ThermoFormat theThermoFormat;
////////////////////////////////////////////

bool ThermoFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(!pmol)
		return false;
	pmol->SetDimension(0);
	OBNasaThermoData* pND = new OBNasaThermoData; //to store rate constant data
	pmol->SetData(pND);

	istream &ifs = *pConv->GetInStream();

	double DefaultMidT = 1500;
	char ln[BUFF_SIZE];
	int i;

	//find line with 1 in col 80
	do
	{
		if(!ifs.getline(ln,BUFF_SIZE) || !strncasecmp(ln,"END",3))
			return false; 
	}while(ln[79]!='1');

	char phase, nam[25], dum[7], elname[3];
	elname[2]=0;
	int elnum;
	double Coeff[14];

	sscanf(ln,"%18s%6s",nam,dum);
	pmol->SetTitle(nam);
	char* p=ln+24;
	if(ln[80]=='&')
	{
		//Reaction Design extension
		p+=20;
		string line;
		if(!getline(ifs,line))return false; 
		vector<string> toks;
		tokenize(toks,line," \t\n\r");
		for(i=0;i<toks.size();i+=2)
		{
			OBAtom atom;
			atom.SetAtomicNum(etab.GetAtomicNum(toks[i].c_str()));
			elnum = atoi(toks[i+1].c_str());
			atom.ForceNoH();
			for(;elnum>0;--elnum)
				pmol->AddAtom(atom);
		}	
	}
	else
	{
		for(i=0;i<4;i++,p+=5)
		{
			char snum[4]={0,0,0,0};//Was problem with F   10   0 reading as ten
			sscanf(p,"%c%c%c%c%c",elname,elname+1,snum,snum+1,snum+2);
			elnum=atoi(snum);
			if(elname[0]!=' ' && elname[0]!='0')
			{
				if(elname[1]==' ')
					elname[1]=0;
				OBAtom atom;
				atom.SetAtomicNum(etab.GetAtomicNum(elname));
				atom.ForceNoH();
				for(;elnum>0;--elnum)
					pmol->AddAtom(atom);
			}
		}
	}
	double LoT, HiT, MidT=0;
	int nc = sscanf(p,"%c%10lf%10lf10%lf",&phase, &LoT, &HiT, &MidT);
	pND->SetPhase(phase);
	pND->SetLoT(LoT);
	pND->SetHiT(HiT);
	if(MidT>HiT || MidT<LoT)
		MidT=DefaultMidT;
	pND->SetMidT(MidT);
	if (!ifs.getline(ln, BUFF_SIZE)) return false;
	p=ln;
	for(i=0;i<5;i++,p+=15)
		sscanf(p,"%15lf",&Coeff[i]);
	if (!ifs.getline(ln, BUFF_SIZE)) return false;
	p=ln;
	for(i=5;i<10;i++,p+=15)
		sscanf(p,"%15lf",&Coeff[i]);
	if (!ifs.getline(ln, BUFF_SIZE)) return false;
	p=ln;
	for(i=10;i<14;i++,p+=15)
		sscanf(p,"%15lf",&Coeff[i]);
	
	for(i=0;i<14;++i)
		pND->SetCoeff(i, Coeff[i]);

	pmol->AssignSpinMultiplicity();
	return true;
}

////////////////////////////////////////
bool ThermoFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	string title(pmol->GetTitle());
	OBNasaThermoData* pND = static_cast<OBNasaThermoData*>(pmol->GetData(ThermoData));
	if(!pND)
	{
    obErrorLog.ThrowError(__FUNCTION__,"No thermo data in " + title, obWarning);
		return false;
	}
	ostream &ofs = *pConv->GetOutStream();
	int i;

	string formula = pmol->GetSpacedFormula();
	vector<string> toks;
	tokenize(toks,formula);

	ofs << left << setw(24) << title.substr(0,24);
	//Check that atom numbers are less than 999
	bool toobig=toks.size()>8;
	for(i=0;i<toks.size() && !toobig ;i+=2)
		if(atoi(toks[i+1].c_str())>999)
			toobig =true;
	if(toobig)
		//Reaction Design extension
		ofs << string(20,' ');
	else
	{
		toks.resize(8);
		for(i=0;i<8;i+=2)
		ofs << left << setw(2) << toks[i] << right << setw(3) << toks[i+1];
	}
	ofs << right << pND->GetPhase() << fixed << setprecision(3) << setw(10) << pND->GetLoT();
	ofs << setw(10) << pND->GetHiT() << setw(9) << pND->GetMidT() << "    01";
	
	if(toobig)
		ofs << "&\n" << formula << '\n'; 
	else
		ofs << '\n';

	ofs << scientific << setprecision(7);
	for(i=0;i<5;++i)
		ofs << setw(15) << pND->GetCoeff(i);
	ofs << "    2\n";
	for(i=5;i<10;++i)
		ofs << setw(15) << pND->GetCoeff(i);
	ofs << "    3\n";
	for(i=10;i<14;++i)
		ofs << setw(15) << pND->GetCoeff(i);
	ofs << "                   4\n";
	return true;
}

}//OpenBabel namespace

