/**********************************************************************
xout.cpp - OpExtraOut to write a OBBase object additionally to a file

Copyright (C) 2009 by Chris Morley

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
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/descriptor.h>
#include <iostream>

namespace OpenBabel
{

//This format is not registered and is called from OpExtraOut::Do()
class ExtraFormat : public OBFormat
{
public:
  ExtraFormat(OBConversion* pOrigConv, OBConversion* pExtraConv)
    : _pOrigConv(pOrigConv), _pExtraConv(pExtraConv){ }

  virtual const char* Description() { return "Outputs to two formats"; }

  virtual bool WriteChemObject(OBConversion* pConv)
  {
    OBBase* pOb = pConv->GetChemObject();

    OBMol* pMolCopy;
    if(_pOrigConv)
    {
      //Need to copy pOb. But currently OBBase does not have a virtual Clone() function.
      //So do it just for OBMol.
      OBMol* pmol = dynamic_cast<OBMol*>(pOb);
      if(pmol)
      {
        pMolCopy = new OBMol(*pmol);

        //need to set output index manually
        //Output to the original output.
        _pOrigConv->SetOutputIndex(pConv->GetOutputIndex()-2);
        if(!_pOrigConv->AddChemObject(pMolCopy))
          pConv->SetLast(true); //error on main output stops conversion
        _pOrigConv->SetLast(pConv->IsLast());
      }
    }

    //Output to extra output
    if(_pExtraConv)
    {
      _pExtraConv->SetOutputIndex(pConv->GetOutputIndex()-2);
      if(!_pExtraConv->AddChemObject(pOb))
        _pExtraConv = NULL; //error on extra output stops only it
      _pExtraConv->SetLast(pConv->IsLast());
    }

    if(pConv->IsLast())
    {
      _pOrigConv->AddChemObject(pMolCopy); //dummy add to empty queue
      if(_pExtraConv)
      {
        _pExtraConv->AddChemObject(pOb);  //dummy add to empty queue
        delete _pExtraConv->GetOutStream(); //filestream
      }
      pConv->SetOutFormat(_pOrigConv->GetOutFormat()); //ReportNumberConverted() uses this at end
      delete _pOrigConv;
      delete _pExtraConv;
      _pOrigConv=NULL;
      _pExtraConv=NULL;
      delete this;//self destruction; was made in new in an OBOp
    }
    return true;
  }

private:
  OBConversion *_pOrigConv, *_pExtraConv;
};

//*****************************************************************
/**
  * This op is mainly intended to be used to generate svg files for display
  * from the GUI. But it also can be used to generate two output files with
  * different formats:
  *   babel infile.xxx outfile.yyy --xout extra.zzz
 **/
 /*
  This and other ops doing simple things are rather more complicated and
  difficult to code correctly than is desirable. The reasons are:
   - code segregation. The ops and formats are plugins and so can not use
     be called or controlled by any specific code in the main program. Also
     the commandline, the GUI and OBConversion know nothing of the chemistry
     (OBMol, etc.).
   - the conversion process is designed not to be store the OBBase objects,
     making it potentially applicable to any size of dataset. However, some
     operations, such as sorting and uniqueness checking obviosly require this,
     requiring special handling.
   - The conversion process uses a queue of two so that the output formats can
     do special things with the last object, such as adding </cml> to close a
     top level element in XML. This requires special handling if the normal
     conversion process is subverted.
 */
class OpExtraOut : public OBOp
{
public:
  OpExtraOut(const char* ID) : OBOp(ID, false){};
  const char* Description(){ return "<file.xxx> Additional file output\n"
       "Mainly intended to be used to generate svg files for display from the GUI,\n"
       "but can also be used to output to two different formats:\n"
       "      obabel infile.sdf  -osmi  -0xout secondout.svg"; }

  virtual bool WorksWith(OBBase* pOb)const{ return dynamic_cast<OBMol*>(pOb)!=NULL; }
  virtual bool Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv);
};

/////////////////////////////////////////////////////////////////
OpExtraOut theOpExtraOut("0xout"); //Global instance

/////////////////////////////////////////////////////////////////
bool OpExtraOut::Do(OBBase* pOb, const char* OptionText, OpMap* pmap, OBConversion* pConv)
{
  /*
    OptionText contains an output filename with a format extension.
    Make an OBConversion object with this as output destination.
    Make a copy the current OBConversion and replace the output format by
    an instance of ExtraFormat. This then does all the subsequent work.
  */
  if(!pConv || !OptionText || *OptionText=='\0')
    return true; //silent no-op. false would prevent the main output

  if(pConv->IsFirstInput())
  {
    OBConversion* pExtraConv = new OBConversion(*pConv); //copy ensures OBConversion::Index>-1
    std::ofstream* ofs;
    if( (ofs = new std::ofstream(OptionText)) ) // extra parens to indicate truth value
      pExtraConv->SetOutStream(ofs);
    if(!ofs || !pExtraConv->SetOutFormat(OBConversion::FormatFromExt(OptionText)))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Error setting up extra output file", obError);
      return true;
    }
    OBConversion* pOrigConv = new OBConversion(*pConv);

    //Make an instance of ExtraFormat and divert the output to it. It will delete itself.
    pConv->SetOutFormat(new ExtraFormat(pOrigConv, pExtraConv));
  }
  return true;
}

}//namespace
