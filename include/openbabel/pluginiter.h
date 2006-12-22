/**********************************************************************
pluginiter.h - facilitates construction of plugin classes
 
Copyright (C) 2006 by Chris Morley
 
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

#ifndef OB_PLUGINITER_H
#define OB_PLUGINITER_H

#include <openbabel/babelconfig.h>

namespace OpenBabel
{
template<typename T>
class OBAPI PluginIter
{
private:
  typedef std::map<const std::string, T*> Maptype;
  Maptype _map;
  typename Maptype::iterator _itr;
  T* _default;

public:
  void Register(T* pType, const std::string ID, bool IsDefault)
  {
    _map[ID] = pType;
    if(IsDefault || _map.empty())
      _default=pType;
  }

  T* FindType(const std::string& ID)
  {
    if(ID.empty())
      return _default;
    _itr = _map.find(ID);
    if(_itr==_map.end())
      return NULL;
    else
      return _itr->second;
  }

  T* FindDefaultType() const { return _default; }

  std::string ID() const { return _itr->first; }

  void ToStart() { _itr = _map.begin(); }
  
  PluginIter& operator++()
  {
    ++_itr;
    return *this;
  }
  operator bool() const { return _itr != _map.end(); }
  T* operator->() const   { return _itr->second;  }
  T& operator*() const    { return *(_itr->second); }
};

//Macro to iterate over all sub-types
#define FOR_EACH(plugintype, f)  for(PluginIter<plugintype>& f=plugintype::Iter(); f; ++f )

//Macro to be added to definition of the base class
#define MAKE_PLUGIN(BaseClass)\
public:\
BaseClass(std::string ID, bool IsDefault=false)\
{Iter().Register(this, ID, IsDefault);}\
static PluginIter<BaseClass>& Iter()\
{static PluginIter<BaseClass>* p = NULL;\
if(!p) p = new PluginIter<BaseClass>;\
p->ToStart();\
return *p;}\
static BaseClass* FindDefaultType(){ return Iter().FindDefaultType();}\
  static BaseClass* FindType(const std::string& ID){ return Iter().FindType(ID);}

} //namespace
#endif //OB_PLUGINITER_H

/**
The code in this file makes it easy to make 'plugin' classes. These classes are
derived from an abstract base class, like OBFingerprint. The derived classes 
('sub-types' like fingerprint2) usually have a single instance. Plugin classes 
are only discovered at runtime, so no existing code needs to be changed when 
adding a new derived class. In some builds the new code can be added or removed 
by just moving a DLL or so file.

-----------------------------------------------------------------------------
<strong>Step-by-step Instructions for use with a fictitious base class, YourBaseClass.</strong>

1) In the header file for YourBaseClass.
#include "pluginiter.h" and in the definition of YourBaseClass add the
MAKE_PLUGIN macro and a pure virtual function Description().
\code
class YourBaseClass
{
  MAKE_PLUGIN(YourBaseClass)
  virtual string Description()=0;

  ...rest of implementation, probably involving virtual functions redefined
  in the sub-type classes
};
\endcode
See below for what this macro contains.

2) Declare each sub-type in a class derived from the base class
and give it a constructor which calls the base class constructor as shown:
\code
class YourSubType1 : public YourBaseClass
{
public:
  YourSubtype1(string ID, bool IsDefault=false) 
    : YourBaseClass(ID, IsDefault){}

  virtual string Description()
  { return "A short useful description";};

  ...rest ofimplementation
};
\endcode

3) Declare a global instance of the sub-type class which specifies its ID.
and, optionally, whether it is to be regarded as the default type of YourBaseClass.
\code
YourSubType1 theType1("FP2",true);
\endcode

4) The following functions are available:
YourBaseClass* YourBaseClass::FindType(const string& ID);
YourBaseClass* YourBaseClass::FindDefaultType();

PluginIter<YourBaseClass>& YourBaseClass::Iter();
This returns an object which looks like a pointer to YourBaseClass when used
with * or -> . It initially points to the first sub-type and can be subsequently
made to point to all the rest using the prefix ++ operator. When tested as a bool
it returns false when there are no more subtypes. 

But the easiest way to access all the subtypes is to use the macro FOR_EACH.
For example to print out all the subtypes with their descriptions:
\code
FOR_EACH(YourBaseClass, iter)
{
   cout << iter.ID() << ' ' << iter.Description() << endl;
}
\endcode

---------------------------------------------------------------------
H<strong>How it works</strong>

MAKE_PLUGIN(YourBaseClass) inserts the following code into YourBaseClass:\code

//This static function returns a reference to the PluginIter object,
//which contains the map of sub-types and the default sub-type.
//Because it is a static local variable it is constructed only once.
//This avoids the "static initialization order fiasco",
//see www.parashift.com/c++-faq-lite/.
//Every time this function is used it sets the iterator to the start of the map.
static PluginIter<YourBaseClass>& Iter()
{
  static PluginIter<YourBaseClass>* p = NULL;
  if(!p)
    p = new PluginIter<YourBaseClass>*;
  p->ToStart();
  return *p;
}

Each sub-type is registered from its constructor as it is loaded, which 
could be at program start up, or later.
YourBaseClass(std::string ID, bool IsDefault=false)
{
  Iter().Register(this, ID, IsDefault);
}

//The following just pass on the work to the PluginIter object.
YourBaseClass* FindDefaultType(){ return Iter().FindDefault();}
YourBaseClass* FindType(const string& ID){ return Iter().FindType(const string& ID);}
\endcode
*/