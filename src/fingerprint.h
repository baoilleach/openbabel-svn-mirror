/**********************************************************************
fingerprint.h - Base class for fingerprints and fast searching 
 
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

#ifndef OB_FINGERPRINT_H
#define OB_FINGERPRINT_H

#include <list>
#include <map>
#include <set>
#include <vector>
#include <string>

namespace OpenBabel
{
	class OBBase; //Forward declaration; used only as pointer.

/// \brief The base class for fingerprints
class OBAPI OBFingerprint
{
//see end of cpp file for detailed documentation
public:
	/// Sets the nth bit
	void SetBit(std::vector<unsigned int>& vec, unsigned int n);	

	/// Repeatedly ORs the top half with the bottom half until no smaller than nbits 
	void Fold(std::vector<unsigned int>& vec, unsigned int nbits); 

	/// Returns fingerprint in vector, which may be resized, folded to nbits (if nbits!=0)
	virtual bool GetFingerprint(OBBase* pOb, std::vector<unsigned int>& fp, int nbits=0)=0;

	/// Required short description of the fingerprint type.
	virtual std::string Description()=0;

	/// Optional flags
	enum FptFlag{FPT_UNIQUEBITS=1};
	virtual unsigned int Flags() { return 0;}; 

	/// Obtain info on available fingerprints
	static bool GetNextFPrt(std::string& id, OBFingerprint*& pFPrt);

	/// Returns a pointer to a fingerprint (the default if ID is empty), or NULL if not available
	static OBFingerprint* FindFingerprint(std::string& ID);

	/// Returns the Tanimoto coefficient between two vectors (vector<unsigned int>& SeekPositions)
	static double Tanimoto(const std::vector<unsigned int>& vec1, const std::vector<unsigned int>& vec2);
	
	/// Inline version of Tanimoto() taking a pointer for the second vector
	static double OBFingerprint::Tanimoto(const std::vector<unsigned int>& vec1, const unsigned int* p2) 
	{
		///If used for two vectors, vec1 and vec2, call as Tanimoto(vec1, &vec2[0]);
		int andbits=0, orbits=0;
		unsigned int i;
		for (i=0;i<vec1.size();++i)
		{
			int andfp = vec1[i] & p2[i];
			int orfp = vec1[i] | p2[i];
			//Count bits
			for(;andfp;andfp=andfp<<1)
				if(andfp<0) ++andbits;
			for(;orfp;orfp=orfp<<1)
				if(orfp<0) ++orbits;
		}
			return((double)andbits/(double)orbits);
	};
	
	static unsigned int Getbitsperint(){ return bitsperint; }

private:
	///Function object to set bits
	struct bit_or
	{
		unsigned int operator()(const unsigned int a, const unsigned int b)
		{
			return a | b;	
		}
	};
	
	typedef std::map<std::string, OBFingerprint*> FPMapType;
	typedef FPMapType::iterator Fptpos;

protected:
	///This static function returns a reference to the FPtsMap
	///which, because it is a static local variable is constructed only once.
	///This fiddle is to avoid the "static initialization order fiasco"
	///See Marshall Cline's C++ FAQ Lite document, www.parashift.com/c++-faq-lite/". 
	static FPMapType& FPtsMap()
	{
		static FPMapType* fptm = new FPMapType;
		return *fptm;
	};

	OBFingerprint(std::string ID, bool IsDefault=false)
	{
		FPtsMap()[ID] = this; //registers the derived fingerprint class
		if(IsDefault || FPtsMap().empty())
			_pDefault=this;
	};
	
private:
	static OBFingerprint* _pDefault;
	static const unsigned int bitsperint;// = 8 * sizeof(unsigned int);
	static int rubbish;
};




//*************************************************************
//Fast search routines
///Header for fastsearch index file
struct OBAPI FptIndexHeader
{
	unsigned int headerlength;///<offset to data: sizeof(FptIndexHeader)
	unsigned int nEntries;    ///<number of fingerprints
	unsigned int words;				///<number 32bit words per fingerprint
	char fpid[16];            ///<ID of the fingerprint type
	char datafilename[256];   ///<the data that this is an index to
};
/// Structure of fastsearch index files
struct OBAPI FptIndex
{
	FptIndexHeader header;
	std::vector<unsigned int> fptdata;
	std::vector<unsigned int> seekdata;
	bool Read(std::istream* pIndexstream);
	///\brief Returns pointer to FP used or NULL and an error message
	OBFingerprint* CheckFP();
};

/// \brief Class to search fingerprint index files
class OBAPI FastSearch
{
//see end of cpp file for detailed documentation
public:
  std::string ReadIndex(std::istream* pIndexstream);
	virtual ~FastSearch(){};

	/// \brief Does substructure search and returns vector of the file positions of matches 
	bool    Find(OBBase* pOb, std::vector<unsigned int>& SeekPositions, unsigned int MaxCandidates);

	/// \brief Returns multimap containing objects whose Tanimoto coefficients with the target
	///	is greater than the value specified.
	bool    FindSimilar(OBBase* pOb, std::multimap<double, unsigned int>& SeekposMap,
		double MinTani);

	/// \brief Returns multimap containing the nCandidates objects with largest Tanimoto
	///  coefficients with the target.
	bool    FindSimilar(OBBase* pOb, std::multimap<double, unsigned int>& SeekposMap,
		int nCandidates=0);

	/// \brief Returns a pointer to the fingerprint type used to constuct the index
	OBFingerprint* GetFingerprint() const{ return _pFP;};

private:
	FptIndex   _index;
	OBFingerprint* _pFP;
};

//**********************************************
/// \brief Class to prepare fingerprint index files See FastSearch class for details
class OBAPI FastSearchIndexer
{
//see end of cpp file for detailed documentation
public:
	///\brief Constructor with a new index
	FastSearchIndexer(std::string& datafilename, std::ostream* os, std::string& fpid,
			int FptBits=0);

	///\brief Constructor using existing index
	FastSearchIndexer(FptIndex* pindex, std::ostream* os);
	
	~FastSearchIndexer();

	///\brief Called for each object
	bool Add(OBBase* pOb, std::streampos seekpos);

private:
	std::ostream* _indexstream;
	FptIndex*		_pindex;
	OBFingerprint* _pFP;
	int _nbits;
};

} //namespace OpenBabel
#endif

//! \file fingerprint.h
//! \brief Declaration of OBFingerprint base class and fastsearch classes
