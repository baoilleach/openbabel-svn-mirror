/**********************************************************************
op.h - plugin options or operations
 
Copyright (C) 2007 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
***********************************************************************/

#ifndef OB_OP_H
#define OB_OP_H

#include <openbabel/babelconfig.h>
#include <string>
#include <map>
#include <openbabel/plugin.h>
#include <openbabel/base.h>

namespace OpenBabel
{

class OBAPI OBOp : public OBPlugin
{
  MAKE_PLUGIN(OBOp);

public:
  typedef const std::map<std::string, std::string> OpMap ;

  ///Provides the name of this kind of plugin. Use -L "ops" to list from commandline.
  virtual const char* TypeID(){ return "ops"; }

  ///Required function that does the work. Normally return true, unless object is not to be output.
  virtual bool Do(OBBase* pOb, OpMap* pOptions=NULL, const char* OptionText=NULL)=0;

  ///Returns true if this op is designed to work with the class of pOb, e.g. OBMol
  virtual bool WorksWith(OBBase* pOb)const=0;

  /// Returns string describing options, for display with -H and to make checkboxes in GUI
  static std::string OpOptions(OBBase* pOb)
  {
    std::string s;
    OBPlugin::PluginIterator itr;
    for(itr=OBPlugin::Begin("ops");itr!=OBPlugin::End("ops");++itr)
    {
      if(*(itr->first)=='_')//ignore ops with IDs that begin with '_'
        continue;
      OBOp* pOp = dynamic_cast<OBOp*>(itr->second);
      if(pOp && pOp->WorksWith(pOb))
      {
        s += "--";
        s += itr->first; //ID
        s += ' ';
        s += OBPlugin::FirstLine(pOp->Description()) + '\n';
      }
    }
    s += '\n';
    return s;
  }

  ///Call Do() of all the OBOps whose ID is a key in the map.
  ///Called from Transform(). The map has general options like -x or --multicharoption
  ///The key is the option name and the value, if any, is text which follows the option name.
  /// In some cases, there may be several parameters, space separated)
  ///Returns false indicating object should not be output, if any Do() returns false
  static bool DoOps(OBBase* pOb, OpMap* pOptions)
  {
    OpMap::const_iterator itr;
    for(itr=pOptions->begin();itr!=pOptions->end();++itr)
    {
      OBOp* pOp = FindType(itr->first.c_str());
      if(pOp)
        if(!pOp->Do(pOb, pOptions, itr->second.c_str()))
          return false; //Op has decided molecule should not be output
    }
    return true;
  }
};
}//namespace

/*
Classes derived from OBOp implement options for the babel program (for both
its commandline and GUI interfaces). It is intended for options that carry out some
modification on the molecule(or reaction) after it has been input, but before
it is output. An example is the --center option implemented in the OpCenter class
in ops.cpp, which is a duplicate of the built in -c option for centering coordinates.

The advantage of plugin classes is that no existing code has to be modified 
when a new class is added. You can list those that are present by 
babel -L ops 
or from a menu item in the GUI.

Any OBOp derived class has to have a constructor, a function returning a short description,
and a Do() function which does the work. It also needs a WorksWith() function
which is always the same when operating on OBMol objects. (It is not made a default
to reducecode dependencies.) A single global instance of the class needs to be
instantiated to define the ID, by which the class is subsequently accessed.

OBOp works by two of its static functions being called from code in transform.cpp:
 - OpOptions(OBBase* pOb, OpMap* pOptions) returns a string describing each of the
derivated classes relevant to objects of the class of the OBBase parameter,
for use in the help text and to set checkboxes in the GUI;
 - DoOps(OBBase* pOb) applies each option whose ID is listed in the  Opmap parameter
to the object (ususally an OBMol) in the OBBase parameter.

Options which need parameters are passed these (space delimited) in the text parameter
of the Do() function. They can also access other general options specified on the
command line by examining the the OpMap parameter.

To use an OBOp class from the API it is necessary to use an extra step in case it isn't
present. So to apply the OBOp clas with ID gen3D to your mol

OBOp* pOp = OBOp::FindType("gen3D");
if(!pOp)
  ...report error
pOp->Do(mol);

  */
#endif
