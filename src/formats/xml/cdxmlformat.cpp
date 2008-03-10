/**********************************************************************
Copyright (C) 2006 by Geoff Hutchison
 
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
#include <openbabel/xml.h>

#ifdef WIN32
#pragma warning (disable : 4800)
#endif

using namespace std;
namespace OpenBabel
{

class ChemDrawXMLFormat : public XMLMoleculeFormat
{
public:
	ChemDrawXMLFormat(): Order (-1)
	{
		OBConversion::RegisterFormat("cdxml", this, "chemical/x-cdxml");
		XMLConversion::RegisterXMLFormat(this, false, "http://www.camsoft.com/xml/cdxml.dtd");
		XMLConversion::RegisterXMLFormat(this);
	}
	virtual const char* NamespaceURI()const{return "http://www.cambridgesoft.com/xml/cdxml.dtd";}
  virtual const char* Description()
  {
      return " \
ChemDraw CDXML format \n \
Minimal support of chemical structure information only.\n \
\n";
};

  virtual const char* GetMIMEType() 
  { return "chemical/x-cdxml"; };

  virtual const char* SpecificationURL()
  {return "http://www.cambridgesoft.com/services/documentation/sdk/chemdraw/cdx/";}


  virtual unsigned int Flags()
  {
      return READXML;
  };

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
	virtual bool DoElement(const string& name);
	virtual bool EndElement(const string& name);
 
	// EndTag is used so that the stream buffer is is filled with the XML from
	// complete objects, as far as possible. 
	virtual const char* EndTag(){ return "/fragment>"; };

    //atoms and bonds might have no content, so EndElement is not always called
    // that's why we need to ensure that atoms and bonds are really added.
    void EnsureEndElement(void);

private:
  OBAtom _tempAtom; //!< A temporary atom as the atom tag is read
  int Begin, End, Order, Flag; // Data for current bond
  map <int, int> atoms; // maps chemdraw atom id to openbabel idx.
  int _offset; // used to ensure that atoms have different ids.
  double _scale; // current scale
};

////////////////////////////////////////////////////////////////////

ChemDrawXMLFormat theChemDrawXMLFormat;

////////////////////////////////////////////////////////////////////

bool ChemDrawXMLFormat::DoElement(const string& name)
{
  string buf;
  if(name=="fragment") 
  {
    //This is the start of the molecule we are extracting and it will
    //be put into the OBMol* _pmol declared in the parent class.
    //initialise everything
    _tempAtom.Clear();
    atoms.clear();

    _pmol->SetDimension(2);
    _pmol->BeginModify();
  }
  else if(name=="n")
  {
    EnsureEndElement();
    buf = _pxmlConv->GetAttribute("Type");
    if (buf.length())
    {
      if (buf != "Unspecified" && buf != "Element")
      {
        cerr << "CDXML Format: Node type \"" << buf <<
                  "\" is not currently supported." << endl;
        return false; // FIXME: use as many types as possible
	  }
    }
    _tempAtom.SetAtomicNum(6); // default is carbon
    buf = _pxmlConv->GetAttribute("id");
    if (buf.length())
      _tempAtom.SetIdx(atoi(buf.c_str()));
    buf = _pxmlConv->GetAttribute("Element");
    if (buf.length())
      _tempAtom.SetAtomicNum(atoi(buf.c_str()));

    buf = _pxmlConv->GetAttribute("p"); // coords
    if (buf.length())
    {
      double x = 0., y = 0.;
      sscanf(buf.c_str(), "%lf %lf", &x, &y);
      _tempAtom.SetVector(x, y, 0.);
    }
    buf = _pxmlConv->GetAttribute("Charge");
    if (buf.length())
      _tempAtom.SetFormalCharge(atoi(buf.c_str()));
  }
  else if(name=="b")
  {
    EnsureEndElement();
	bool invert_ends = false;
	Begin = End = Flag = 0;
    buf = _pxmlConv->GetAttribute("Order");
    if (buf.length())
      Order = atoi(buf.c_str());
    else
      Order = 1; //default value
    buf = _pxmlConv->GetAttribute("Display");
    if (buf.length())
    {
	  if (buf == "WedgeEnd")
      {
        invert_ends = true;
        Flag = OB_WEDGE_BOND;
	  }
      else if (buf == "WedgeBegin")
        Flag = OB_WEDGE_BOND;
      else if (buf == "WedgedHashBegin")
      {
        invert_ends = true;
        Flag = OB_HASH_BOND;
	  }
      else if (buf == "Hash" || buf == "WedgedHashEnd")
        Flag = OB_HASH_BOND;
    }
    buf = _pxmlConv->GetAttribute("B");
    if (buf.length())
    {
      if (invert_ends)
        End = atoms[atoi(buf.c_str())];
      else
        Begin = atoms[atoi(buf.c_str())];
    }
    buf = _pxmlConv->GetAttribute("E");
    if (buf.length())
    {
      if (invert_ends)
        Begin = atoms[atoi(buf.c_str())];
      else
        End = atoms[atoi(buf.c_str())];
    }
  }

  return true;
}

bool ChemDrawXMLFormat::EndElement(const string& name)
{
  unsigned int i;
  if(name=="n")
  {
    _pmol->AddAtom(_tempAtom);
    atoms[_tempAtom.GetIdx()] = _pmol->NumAtoms();
    _tempAtom.Clear();
  }
  else if(name=="b")
  {
    _pmol->AddBond(Begin, End, Order, Flag);
	Order = -1;
  }
  else if(name=="fragment") //this is the end of the molecule we are extracting
  {
    EnsureEndElement();
    _pmol->EndModify();
    atoms.clear();
    return false;//means stop parsing
  }
  return true;
}	

void ChemDrawXMLFormat::EnsureEndElement(void)
{
  if (_tempAtom.GetAtomicNum() != 0)
  {
    _pmol->AddAtom(_tempAtom);
    atoms[_tempAtom.GetIdx()] = _pmol->NumAtoms();
    _tempAtom.Clear();
  }
  else if (Order >= 0)
  {
    _pmol->AddBond(Begin, End, Order, Flag);
	Order = -1;
  }
}

bool ChemDrawXMLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  static const xmlChar C_MOLECULE[]         = "fragment";
  static const xmlChar C_CDXML[]            = "CDXML";
  static const xmlChar C_BONDLENGTH[]       = "BondLength";
  static const xmlChar C_PAGE[]             = "page";
  static const xmlChar C_ATOM[]             = "n";
  static const xmlChar C_BOND[]             = "b";
  static const xmlChar C_ID[]               = "id";

  static const xmlChar C_CHARGE[]           = "Charge";
  static const xmlChar C_COORDS[]           = "p";
  static const xmlChar C_ELEMENT[]          = "Element";
  static const xmlChar C_ORDER[]            = "Order";
  static const xmlChar C_BEGIN[]            = "B";
  static const xmlChar C_END[]              = "E";
  static const xmlChar C_DISPLAY[]          = "Display";

  _pxmlConv = XMLConversion::GetDerived(pConv,false);
  if(!_pxmlConv)
    return false;

  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
	return false;
  OBMol &mol = *pmol;

  OBBond *pbond;
  vector<OBBond*>::iterator j;
  if(_pxmlConv->GetOutputIndex() == 1)
  {
    xmlTextWriterStartDocument(writer(), NULL, NULL, NULL);
    xmlTextWriterWriteDTD(writer(), BAD_CAST "CDXML", NULL, BAD_CAST "http://www.camsoft.com/xml/cdxml.dtd", NULL);
    xmlTextWriterStartElement(writer(), C_CDXML);
    xmlTextWriterWriteFormatAttribute(writer(), C_BONDLENGTH , "30");
    xmlTextWriterStartElement(writer(), C_PAGE); // put everything on one page
    // now guess the average bond size for the first molecule and scale to 30.
    _scale = 0.;
    if (mol.NumBonds())
    {
      for (pbond = mol.BeginBond(j);pbond;pbond = mol.NextBond(j))
        _scale += pbond->GetLength();
      _scale /= mol.NumBonds();
	}
    else
      _scale = 1.; // FIXME: what happens if the molecule has no bond?
    _scale = 30. / _scale;
    _offset = 0;
  }
	
  xmlTextWriterStartElement(writer(), C_MOLECULE);

  OBAtom *patom;
  vector<OBAtom*>::iterator i;
  int n;
  for (patom = mol.BeginAtom(i);patom;patom = mol.NextAtom(i))
  {
    xmlTextWriterStartElement(writer(), C_ATOM);

    xmlTextWriterWriteFormatAttribute(writer(), C_ID , "%d", patom->GetIdx() + _offset);
    xmlTextWriterWriteFormatAttribute(writer(), C_COORDS , "%f %f", patom->GetX() * _scale, patom->GetY() * _scale);
	n = patom->GetAtomicNum();
    if (n != 6)
    {
      xmlTextWriterWriteFormatAttribute(writer(), C_ELEMENT , "%d", n);
    }
    n = patom->GetFormalCharge();
    if (n != 0)
    {
      xmlTextWriterWriteFormatAttribute(writer(), C_CHARGE , "%d", n);
    }
    xmlTextWriterEndElement(writer());
  }

  for (pbond = mol.BeginBond(j);pbond;pbond = mol.NextBond(j))
  {
    xmlTextWriterStartElement(writer(), C_BOND);
	patom = pbond->GetBeginAtom();
    xmlTextWriterWriteFormatAttribute(writer(), C_BEGIN , "%d", patom->GetIdx() + _offset);
	patom = pbond->GetEndAtom();
    xmlTextWriterWriteFormatAttribute(writer(), C_END , "%d", patom->GetIdx() + _offset);
	n = pbond->GetBO();
    if (n != 1)
    {
      xmlTextWriterWriteFormatAttribute(writer(), C_ORDER , "%d", n);
    }
    if (pbond->IsWedge())
      xmlTextWriterWriteFormatAttribute(writer(), C_DISPLAY , "WedgeBegin");
    else if (pbond->IsHash())
      xmlTextWriterWriteFormatAttribute(writer(), C_DISPLAY , "WedgedHashEnd");
    xmlTextWriterEndElement(writer());
  }
  _offset += mol.NumAtoms ();

  xmlTextWriterEndElement(writer());//molecule

  if(_pxmlConv->IsLast())
  {
    xmlTextWriterEndDocument(writer()); // page
    xmlTextWriterEndDocument(writer()); //document
    OutputToStream();
  }
  return true;
}

}//namespace
