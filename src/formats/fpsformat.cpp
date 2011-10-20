/**********************************************************************
Copyright (C) 2011 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <vector>
#include <string>
#include <iomanip>

#include <openbabel/obmolecformat.h>
#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel
{

  /// \brief Outputs collections of fingerprints
  class FPSFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FPSFormat() {OBConversion::RegisterFormat("fps",this);}

    virtual const char* Description() //required
    { return
    "FPS text fingerprint format (Dalke)\n"
      "Any molecule without a title is given one of the form #n.\n"
      "Write Options e.g. -xf FP3 -xN 128\n"
      " f<id> fingerprint type\n"
      " N# fold to specified number of bits, 32, 64, 128, etc.\n"
      " p use full input path as source, not just filename\n"
      " t <text> use <text> as source in header\n\n";
    }
  virtual const char* SpecificationURL()
  { return "http://code.google.com/p/chem-fingerprints/wiki/FPS"; }

    virtual unsigned int Flags(){return NOTREADABLE;};
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  private:
    string getTimeStr();
  private:
    int _nbits;
    OBFingerprint* _pFP;
  };

  ////////////////////////////////////////////////////
  //Make an instance of the format class
  FPSFormat theFPSFormat;

//*******************************************************************
static inline unsigned short bswap_16(unsigned short x) {
  return (x>>8) | (x<<8);
}
static inline unsigned int bswap_32(unsigned int x) {
  return (bswap_16(x&0xffff)<<16) | (bswap_16(x>>16));
}
/////////////////////////////////////////////////////////////////////

bool FPSFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  ostream &ofs = *pConv->GetOutStream();
  vector<unsigned int> fptvec;

  if(pConv->GetOutputIndex()==1)
  {
    string fpid;
    const char* p=pConv->IsOption("f");
    if(p)
    {
      fpid=p;
      fpid = fpid.substr(0,fpid.find('"'));
    }

    _pFP = OBFingerprint::FindFingerprint(fpid.c_str());
    if(!_pFP)
    {
      stringstream errorMsg;
      errorMsg << "Fingerprint type '" << fpid << "' not available" << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    p=pConv->IsOption("N");
    _nbits=0;
    if(p)
      _nbits = atoi(p);
    if(_nbits<0)
      obErrorLog.ThrowError(__FUNCTION__,
      "The number of bits to fold to, in the-xN option, should be >=0", obWarning);

    if(_nbits==0) //if not folded, use first number on second line of the description
    {
      // A quirk is that the fingerprint has to be used before the
      // desciption contains the number of bits and the version.
      _pFP->GetFingerprint(pOb, fptvec, _nbits);
      _nbits = atoi(strchr(_pFP->Description(),'\n')+1); 
    }

    //Write metadata to header
    const char* txt = pConv->IsOption("t");
    string source = txt ? txt : pConv->GetInFilename();
    if(!pConv->IsOption("p"))
    {
      string::size_type pos = source.find_last_of("/\\");
      if(pos!=string::npos)
        source.erase(0,pos+1);
    }
    ofs << "#FPS1\n"
        << "#num_bits=" << _nbits << '\n'
        << "#type=OpenBabel-" << _pFP->GetID() << "/1" << '\n'
        << "#software=OpenBabel/" << BABEL_VERSION << '\n'
        << "#source=" << source << '\n'
        << "#date=" << getTimeStr() << endl;
  }

  stringstream molID;
  if(strlen(pOb->GetTitle())==0)
    molID << '#' << pConv->GetOutputIndex() << ' ';
  else
    molID << pOb->GetTitle() << ' ';

  if(!_pFP->GetFingerprint(pOb, fptvec, _nbits))
    return false;

  stringstream ss;
  for(unsigned i=0;i<(_nbits+31)/32;++i)
  {
    ss << hex << setw(8) << setfill('0');
#ifdef WORDS_BIGENDIAN
    ss <<fptvec[i];
#else
    ss << bswap_32(fptvec[i]);
#endif
  }
  // truncate to hex from whole number of bytes (seems to be the way fps does it)
  ofs << dec << ss.str().erase(2*((_nbits+7)/8));
  ofs << '\t' << Trim(molID.str()) << endl;

  return true;
}

/////////////////////////////////////
string FPSFormat::getTimeStr()
{
  //e.g. 2011-09-25T09:56:19
  const int TIME_STR_SIZE = 64;
  time_t akttime;                              /* Systemtime                        */
  char timestr[TIME_STR_SIZE + 1] = "";        /* Timestring                        */
  size_t time_res;                             /* Result of strftime                */

  /* ---- Get the system-time ---- */
  akttime = time((time_t *) NULL);
  time_res = strftime(timestr,
                      TIME_STR_SIZE,
                      "%Y-%m-%dT%H:%M:%S",
                      gmtime((time_t *) &akttime)
                      );
  return string(timestr);
}
}//namespace
 