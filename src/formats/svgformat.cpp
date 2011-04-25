/*
svgformat.cpp  Format for rendering multiple molecules by SVG
Copyright (C) 2009 by Chris Morley

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
#include <openbabel/op.h>
#include <openbabel/text.h>
#include <openbabel/depict/svgpainter.h>
#include <openbabel/depict/depict.h>
#include <openbabel/alias.h>

using namespace std;
namespace OpenBabel
{

class SVGFormat : public OBFormat
{
public:
  SVGFormat() : _ncols(0), _nrows(0), _nmax(0)
  {
    OBConversion::RegisterFormat("svg",this);
    OBConversion::RegisterOptionParam("N", this, 1, OBConversion::OUTOPTIONS);
    OBConversion::RegisterOptionParam("rows", this, 1, OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("cols", this, 1, OBConversion::GENOPTIONS);
    OBConversion::RegisterOptionParam("px", this, 1, OBConversion::GENOPTIONS);
 }

  virtual const char* NamespaceURI()const{return "http://www.w3.org/2000/svg";}
  virtual const char* Description()
  {
    return
      "SVG 2D depiction\n"
      "Scalable Vector Graphics 2D rendering of molecular structure.\n\n"

      "Single molecules are displayed at a fixed scale, as in normal diagrams,\n"
      "but multiple molecules are displayed in a table which expands to fill\n"
      "the containing element, such as a browser window.\n\n"

      "Multiple molecules are displayed in a grid of dimensions specified by\n"
      "the ``-xr`` and ``-xc`` options (number of rows and columns respectively\n"
      "and ``--rows``, ``--cols`` with babel).\n"
      "When displayed in an appropriate program, e.g. Firefox, there is\n"
      "javascript support for zooming (with the mouse wheel)\n"
      "and panning (by dragging with the left mouse button).\n\n"

      "If both ``-xr`` and ``-xc`` are specified, they define the maximum number of\n"
      "molecules that are displayed.\n"
      "If only one of them is displayed, then the other is calculated so that\n"
      "ALL the molecules are displayed.\n"
      "If neither are specified, all the molecules are output in an\n"
      "approximately square table.\n\n"

      "By default, 2D atom coordinates are generated (using gen2D) unless they\n"
      "are already present. This can be slow with a large number of molecules.\n"
      "(3D coordinates are ignored.) Include ``--gen2D`` explicitly if you wish\n"
      "any existing 2D coordinates to be recalculated.\n\n"

      "Write Options e.g. -xu\n"
      " u no element-specific atom coloring\n"
      "    Use this option to produce a black and white diagram\n"
      " U do not use internally-specified color\n"
      "    e.g. atom color read from cml or generated by internal code\n"
      " b black background\n"
      "    The default is white. The atom colors work with both.\n"
      " C do not draw terminal C (and attached H) explicitly\n"
      "    The default is to draw all hetero atoms and terminal C explicitly,\n"
      "    together with their attched hydrogens.\n"
      " a draw all carbon atoms\n"
      "    So propane would display as H3C-CH2-CH3\n"
      " d do not display molecule name\n"
      " s use asymmetric double bonds\n"
      " t use thicker lines\n"
      " e embed molecule as CML\n"
      "    OpenBabel can read the resulting svg file as a cml file.\n"
      " p# scale to bondlength in pixels(single mol only)\n"
      " px# scale to bondlength in pixels(single mol only)(not displayed in GUI)\n"
      " c# number of columns in table\n"
      " cols# number of columns in table(not displayed in GUI)\n"
      " r# number of rows in table\n"
      " rows# number of rows in table(not displayed in GUI)\n"
      " N# max number objects to be output\n"
      " l draw grid lines\n"
      " i add index to each atom\n"
      "    These indices are those in sd or mol files and correspond to the\n"
      "    order of atoms in a SMILES string.\n"
      " j do not embed javascript\n"
      "    Javascript is not usually embedded if there is only one molecule,\n"
      "    but it is if the rows and columns have been specified as 1: ``-xr1 -xc1``\n"
      " w generate wedge/hash bonds(experimental)\n"
      " x omit XML declaration (not displayed in GUI)\n"
      "    Useful if the output is to be embedded in another xml file.\n"
      " A display aliases, if present\n"
      "    This applies to structures which have an alternative, usually\n"
      "    shorter, representation already present. This might have been input\n"
      "    from an A or S superatom entry in an sd or mol file, or can be\n"
      "    generated using the --genalias option. For example::\n \n"

      "      obabel -:\"c1cc(C=O)ccc1C(=O)O\" -O out.svg\n"
      "             --genalias -xA\n \n"

      "    would add a aliases COOH and CHO to represent the carboxyl and\n"
      "    aldehyde groups and would display them as such in the svg diagram.\n"
      "    The aliases which are recognized are in data/superatom.txt, which\n"
      "    can be edited.\n\n"


      "If the input molecule(s) contain explicit hydrogen, you could consider\n"
      "improving the appearance of the diagram by adding an option ``-d`` to make\n"
      "it implicit. Hydrogen on hetero atoms and on explicitly drawn C is\n"
      "always shown.\n"

      "For example, if input.smi had 10 molecules::\n\n"

      "      obabel input.smi -O out.svg -xb -xC -xe\n\n"

      "would produce a svg file with a black background, with no explict\n"
      "terminal carbon, and with an embedded cml representation of each\n"
      "molecule. The structures would be in two rows of four and one row\n"
      "of two.\n\n"
    ;
  }

  virtual unsigned int Flags()
  {
      return NOTREADABLE | ZEROATOMSOK;
  }

  bool WriteChemObject(OBConversion* pConv);
  bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
  bool EmbedCML(OBMol* pmol, OBConversion* pConv);
  bool EmbedScript(ostream& ofs);
private:
  int _ncols, _nrows, _nmax;
  vector<OBBase*> _objects;
  OBText* _ptext;
  string::size_type _textpos;
};
/////////////////////////////////////////////////////////////////
SVGFormat theSVGFormat;

/////////////////////////////////////////////////////////////////
bool SVGFormat::WriteChemObject(OBConversion* pConv)
{
  //Molecules are stored here as pointers to OBOb objects, which are not deleted as usual.
  //When there are no more they are sent to WriteMolecule.
  //This allows their number to be determined whatever their source
  //(they may also have been filtered), so that the table can be properly dimensioned.

  //NOT CURRENTLY IMPLEMENTED
  //If the first object is OBText, the part of it before each insertion point (if it exists)
  //is output before every molecule. This allows molecule structures to be displayed
  //in a template. The x option to omit the XML header is set

  OBBase* pOb = pConv->GetChemObject();

  if(pConv->GetOutputIndex()==1)
  {
    _objects.clear();
    _nmax=0;

    const char* pc = pConv->IsOption("c");
    //alternative for babel because -xc cannot take a parameter, because some other format uses it
    //similarly for -xr -xp
    if(!pc)
      pc = pConv->IsOption("cols", OBConversion::GENOPTIONS);
    const char* pr = pConv->IsOption("r");
    if(!pr)
      pr = pConv->IsOption("rows", OBConversion::GENOPTIONS);
    if(pr)
      _nrows = atoi(pr);
    if(pc)
      _ncols = atoi(pc);
    if(pr && pc) // both specified: fixes maximum number objects to be output
      _nmax = _nrows * _ncols;

    //explicit max number of objects
    const char* pmax =pConv->IsOption("N");
    if(pmax)
      _nmax = atoi(pmax);

/*
    _ptext = dynamic_cast<OBText*>(pOb);
    if(_ptext)
    {
      pConv->AddOption("x");//omit XML header
      _textpos = 0;
      return true;
    }
*/
  }

  OBMoleculeFormat::DoOutputOptions(pOb, pConv);

  //save molecule
  _objects.push_back(pOb);

  bool ret=true;
  //Finish if no more input or if the number of molecules has reached the allowed maximum(if specified)
  bool nomore = _nmax && (_objects.size()==_nmax);
  if((pConv->IsLast() || nomore))
  {
    int nmols = _objects.size();
    //Set table properties according to the options and the number of molecules to be output
    if(!(nmols==0 ||                      //ignore this block if there is no input or
         (_nrows && _ncols) ||            //if the user has specified both rows and columns or
         (!_nrows && !_ncols) && nmols==1)//if neither is specified and there is one output molecule
      )
    {
      if(!_nrows && !_ncols ) //neither specified
      {
        //assign cols/rows in square
        _ncols = (int)ceil(sqrt(((double)nmols)));
      }

      if(_nrows)
        _ncols = (nmols-1) / _nrows + 1; //rounds up
      else if(_ncols)
        _nrows = (nmols-1) / _ncols + 1;
    }

    //output all collected molecules
    int n=0;
/*
    if(_ptext)
    {
      _textpos =0;
      //Output the text up to the first insertion point, or all of it if there is no insertion point.
      *pConv->GetOutStream() << _ptext->GetText(_textpos);
    }
*/
    vector<OBBase*>::iterator iter;
    for(iter=_objects.begin(); ret && iter!=_objects.end(); ++iter)
    {
      //need to manually set these to mimic normal conversion
      pConv->SetOutputIndex(++n);
      pConv->SetLast(n==_objects.size());

      ret=WriteMolecule(*iter, pConv);
/*
      //If there is a subsequent insertion point, output text up to it and update _textpos.
      //If there is not, do nothing.
      if(_ptext)
        *pConv->GetOutStream() << _ptext->GetText(_textpos, true);
*/
    }
/*
    //Output remaining text
    if(_ptext)
      *pConv->GetOutStream() << _ptext->GetText(_textpos);
*/
    //delete all the molecules
    for(iter=_objects.begin();iter!=_objects.end(); ++iter)
      delete *iter;
    delete _ptext;//delete text, NULL or not

    _objects.clear();
    _ptext = NULL;
    _nmax = _ncols = _nrows = 0;
  }
  //OBConversion decrements OutputIndex when returns false because it thinks it is an error
  //So we compensate.
  if(nomore)
    pConv->SetOutputIndex(pConv->GetOutputIndex()+1);
  return ret && !nomore;
}
////////////////////////////////////////////////////////////////
bool SVGFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  ostream &ofs = *pConv->GetOutStream();

  //*** Coordinate generation ***
  //Generate coordinates only if no existing 2D coordinates
  if( (pConv->IsOption("y") || !pmol->Has2D(true)) && !pConv->IsOption("n") )
  {
    OBOp* pOp = OBOp::FindType("gen2D");
    if(!pOp)
    {
      obErrorLog.ThrowError("SVGFormat", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(pmol))
    {
      obErrorLog.ThrowError("SVGFormat", string(pmol->GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!pmol->Has2D() && pmol->NumAtoms()>1)//allows 3D coordinates (if passed by -xn above)
  {
    string mes("Molecule ");
    mes += pmol->GetTitle();
    mes += " needs 2D coordinates to display in SVGformat";
    obErrorLog.ThrowError("SVGFormat", mes, obError);
    return false;
  }
  
  bool hasTable = (_nrows || _ncols);

  string background = pConv->IsOption("b") ? "black" : "white";
  string bondcolor  = pConv->IsOption("b") ? "white" : "black";

  if(pConv->GetOutputIndex()==1)
  {
    //For the first molecule...
    if(hasTable)
    {
      //multiple molecules - use a table
      //Outer svg has viewbox for 0 0 100 100 or adjusted for table shape,
      //and no width or height - it uses the whole of its containing element.
      //Inner svg with width, height, x, y of table cell,
      //and viewbox to match molecule min and max x and y
      if(!pConv->IsOption("x"))
        ofs << "<?xml version=\"1.0\"?>\n";

      ofs << "<svg version=\"1.1\" id=\"topsvg\"\n"
             "xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
             "xmlns:cml=\"http://www.xml-cml.org/schema\" ";

      //*** Outer viewbox ***
      double vbwidth=100, vbheight=100;
      if (_nrows>_ncols)
        vbwidth = (100*_ncols)/_nrows;
      else if(_ncols>_nrows)
        vbheight = (100*_nrows)/_ncols;

      ofs << "x=\"0\" y=\"0\" width=\"100%\" height=\"100%\" ";
      ofs << "viewBox=\"0 0 " << vbwidth << ' ' << vbheight << "\">\n";

      ofs << "<title>OBDepict</title>\n";
      // Draw the background
      ofs << "<rect x=\"0\" y=\"0\" width=\"" << vbwidth << "\" height=\"" << vbheight
          << "\" fill=\"" << background << "\"/>\n";
    }
  }

  //All mols
  double cellsize;
  if(hasTable)
  {
    //*** Parameter for inner svg ***
    int nc = _ncols ? _ncols : 1;
    int nr = (_nrows ? _nrows : 1);
    cellsize = 100. / std::max(nc, nr);
    int indx = pConv->GetOutputIndex() - 1;
    double innerX  = (indx % nc) * cellsize;
    double innerY  = (indx / nc) * cellsize;

    //*** Write molecule name ***
    if(!pConv->IsOption("d"))
      ofs << "<text text-anchor=\"middle\" font-size=\"" << 0.06*cellsize << "\""
      << " fill =\"" << bondcolor << "\" font-family=\"sans-serif\"\n"
      << "x=\"" << innerX + cellsize * 0.5 << "\" y=\"" << innerY + cellsize - 2.0/nr << "\" >"
      << pmol->GetTitle() << "</text>\n";

    SVGPainter painter(*pConv->GetOutStream(), true, cellsize,cellsize,innerX,innerY);
    OBDepict depictor(&painter);

    if(pConv->IsOption("w"))
      depictor.SetOption(OBDepict::genWedgeHash);
    if(!pConv->IsOption("C"))
      depictor.SetOption(OBDepict::drawTermC);// on by default
    if(pConv->IsOption("a"))
      depictor.SetOption(OBDepict::drawAllC);

    if(pConv->IsOption("A"))
    {
      AliasData::RevertToAliasForm(*pmol);
      depictor.SetAliasMode();
    }
    painter.SetFontFamily("sans-serif");
    painter.SetPenColor(OBColor(bondcolor));
    depictor.SetBondColor(bondcolor);
    if(pConv->IsOption("t"))
      painter.SetPenWidth(4);
    else
      painter.SetPenWidth(2);

    //No element-specific atom coloring if requested
    if(pConv->IsOption("u"))
      depictor.SetOption(OBDepict::bwAtoms);
    if(!pConv->IsOption("U"))
      depictor.SetOption(OBDepict::internalColor);
    if(pConv->IsOption("s"))
      depictor.SetOption(OBDepict::asymmetricDoubleBond);

    depictor.DrawMolecule(pmol);

    //Draw atom indices if requested
    if(pConv->IsOption("i"))
      depictor.AddAtomLabels(OBDepict::AtomIndex);

    //Embed CML of molecule if requested
    if(pConv->IsOption("e"))
      EmbedCML(pmol,pConv);
  }
  else //single molecule
  {
    //Nothing written until DrawMolecule call
    //Final </svg> written at the end of this block (painter destructor)
    //This leads to some code duplication.
    double factor = 1.0;
    SVGPainter painter(*pConv->GetOutStream(), false);
    OBDepict depictor(&painter);

    //Scale image by specifying the average bond length in pixels.
    const char* ppx = pConv->IsOption("p");
    if(!ppx)
      ppx= pConv->IsOption("px", OBConversion::GENOPTIONS);
    if(ppx)
    {
      double oldblen = depictor.GetBondLength();
      double newblen = atof(ppx);
      depictor.SetBondLength(newblen);
      factor = newblen / oldblen;
      //Scale bondspacing and font size by same factor
      depictor.SetBondSpacing(depictor.GetBondSpacing() * factor);
      depictor.SetFontSize((int)(depictor.GetFontSize() * factor));
    }

    if(pConv->IsOption("w"))
      depictor.SetOption(OBDepict::genWedgeHash);
    if(!pConv->IsOption("C"))
      depictor.SetOption(OBDepict::drawTermC);// on by default
    if(pConv->IsOption("a"))
      depictor.SetOption(OBDepict::drawAllC);

    if(pConv->IsOption("A"))
    {
      AliasData::RevertToAliasForm(*pmol);
      depictor.SetAliasMode();
    }

    painter.SetFontFamily("sans-serif");
    painter.SetPenColor(OBColor(bondcolor));
    depictor.SetBondColor(bondcolor);
    painter.SetFillColor(OBColor(background));
    if(pConv->IsOption("t"))
      painter.SetPenWidth(4);
    else
      painter.SetPenWidth(1);

    //No element-specific atom coloring if requested
    if(pConv->IsOption("u"))
      depictor.SetOption(OBDepict::bwAtoms);
    if(!pConv->IsOption("U"))
      depictor.SetOption(OBDepict::internalColor);
    if(pConv->IsOption("s"))
      depictor.SetOption(OBDepict::asymmetricDoubleBond);

    depictor.DrawMolecule(pmol);

    //Draw atom indices if requested
    if(pConv->IsOption("i"))
      depictor.AddAtomLabels(OBDepict::AtomIndex);


    //*** Write molecule name ***
    if(!pConv->IsOption("d"))
      ofs << "<text font-size=\"" << 18 * factor  << "\""
      << " fill =\"" << bondcolor << "\" font-family=\"sans-serif\"\n"
      << "x=\"" << 140 * factor << "\" y=\"" << 20 * factor << "\" >"
      << pmol->GetTitle() << "</text>\n";

    //*** Write page title name ***
    ofs << "<title>" << pmol->GetTitle() << " - OBDepict</title>\n";

    //Embed CML of molecule if requested
    if(pConv->IsOption("e"))
      EmbedCML(pmol,pConv);
  }

  if(hasTable && pConv->IsLast())
  {
    //Draw grid lines
    if(_nrows && _ncols && pConv->IsOption("l"))
    {
      for(int i=1; i<_nrows; ++i)
        ofs << " <line  stroke=\"gray\" stroke-width=\"0.1\" x1=\"0\" x2=\"100\""
            << " y1=\""  << i*cellsize << "\" y2=\""  << i*cellsize << "\"/>\n";
      for(int i=1; i<_ncols; ++i)
        ofs << " <line  stroke=\"gray\" stroke-width=\"0.1\" y1=\"0\" y2=\"100\""
            << " x1=\""  << i*cellsize << "\" x2=\""  << i*cellsize << "\"/>\n";
    }

    //Insert javascript for zooming and panning
    if(!pConv->IsOption("j"))
      EmbedScript(ofs);

    ofs << "</svg>\n" << endl;//Outer svg
  }

  return true;
}


/////////////////////////////////////////////////////////////
//returns true if the file "svgformat.script" was inserted into the output
bool SVGFormat::EmbedScript(ostream& ofs)
{
  ifstream ifs;
  if(!ifs || OpenDatafile(ifs, "svgformat.script").empty())
    return false;
  ofs << ifs.rdbuf(); //copy whole file
  return true;
}

///////////////////////////////////////////////////////////////////////////////////
bool SVGFormat::EmbedCML(OBMol* pmol, OBConversion* pConv)
{
  OBConversion CMLConv(*pConv);
  if(!CMLConv.SetOutFormat("cml"))
  {
    obErrorLog.ThrowError(__FUNCTION__, "CML format was not found\n",obError);
    return false;
  }
  CMLConv.AddOption("MolsNotStandalone",OBConversion::OUTOPTIONS);
  CMLConv.AddOption("N",OBConversion::OUTOPTIONS,"cml");
  CMLConv.AddOption("p",OBConversion::OUTOPTIONS); //include properties
//  CMLConv.AddOption("x",OBConversion::OUTOPTIONS);
  return CMLConv.Write(pmol);
}

/*
The script below was originally (and still could be) in data/svgformat.script,
the whole of which is embedded into the output.
It works adequately in Firefox 3 to zoom with the mouse wheel and pan by dragging,
but may need modification for other SVG viewers. (It works in Opera.)

<script type="text/ecmascript">
  <![CDATA[
    addEventListener('DOMMouseScroll', wheel, false);
    onmousewheel = wheel;
    var svgEl = document.getElementById("topsvg");
    function wheel(evt){
      var vb = new Array(4);
      var vbtext = svgEl.getAttributeNS(null,"viewBox");
      vb = vbtext.split(" ");
      var zoom = (evt.detail>0)? 1.41 : 0.71;
      //var dwidth = parseFloat(Math.max(vb[2],vb[3])) * (1-zoom);
      vb[0] = parseFloat(vb[0]) + parseFloat(vb[2])*(1-zoom) * evt.clientX/innerWidth;
      vb[1] = parseFloat(vb[1]) + parseFloat(vb[3])*(1-zoom) * evt.clientY/innerHeight;
      vb[2] = parseFloat(vb[2]) * zoom;
      vb[3] = parseFloat(vb[3]) * zoom;
      svgEl.setAttributeNS(null, "viewBox", vb.join(" "));
    }
    var startx=0;
    var starty=0;
    onmousedown = function(evt) {
      startx = evt.clientX;
      starty = evt.clientY;
    }
    onmousemove=function(evt) {
      if(startx!=0 && starty!=0
        && ((evt.clientX - startx)*(evt.clientX - startx)+(evt.clientY - starty)*(evt.clientY - starty)>100))
      {
        var vbtext = svgEl.getAttributeNS(null,"viewBox");
        vb = vbtext.split(" ");
        var maxwh = Math.max(parseFloat(vb[2]),parseFloat(vb[3]));
        vb[0] = parseFloat(vb[0]) - (evt.clientX - startx)*maxwh/innerWidth;
        vb[1] = parseFloat(vb[1]) - (evt.clientY - starty)*maxwh/innerHeight;
        svgEl.setAttributeNS(null, "viewBox", vb.join(" "));
        startx = evt.clientX;
        starty = evt.clientY;
      }
    }
    onmouseup=function() {
      startx=0;
      starty=0;
   }
  ]]>
</script>

Alternatively, svgformat.script could contain:

<script type="text/ecmascript" xlink:href="morescript.js" />

with the real script in morescript.js.
*/

}//namespace

