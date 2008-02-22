/**********************************************************************
descriptor.cpp - Implementation of the base class for molecular descriptors
 
Copyright (C) 2007 by Chris Morley
 
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
#include <openbabel/oberror.h>
#include <openbabel/generic.h>
#include <openbabel/base.h>
#include <openbabel/descriptor.h>

using namespace std;
namespace OpenBabel
{
/** 
     Compare() is a virtual function and can be overridden to allow different
     comparison behaviour. 
     The default implementation here is suitable for OBDescriptor classes
     which return a double value. The stringstream is parsed to retrieve a
     comparison operator, one of     > < >= <= = == != , and a numerical value.
     The function compares this the value returned by Predict() and returns
     the result. The stringstream is left after the number, and its state 
     reflects whether any errors have occurred.
     If noEval is true, the parsing is as normal but Predict is not called
     and the function returns false.
 **/
bool OBDescriptor::Compare(OBBase* pOb, istream& optionText, bool noEval)
{
  //Get comparison operator
  char ch1=0, ch2=0;
  while (optionText && !ispunctU(ch1))
    optionText >> ch1;
  if(ispunctU(optionText.peek()))
    optionText >> ch2;
  
  //Get number
  double filterval, val;
  optionText >> filterval;
  if (optionText)
  {
    if(noEval)
      return false;
    
    val = Predict(pOb);

    return DoComparison(ch1, ch2, val, filterval);
  }
  optionText.setstate(std::ios::badbit); //shows error
  obErrorLog.ThrowError(__FUNCTION__, "Error in filter string" , obError);
  return false;
}

/// Interprets the --filter option string and returns the combined result of all the comparisons it contains  
/**
    The string has the form:
    PropertyID1 predicate1 [booleanOp] PropertyID2 predicate2 ...
    The propertyIDs are the ID of instances of a OBDescriptor class or
    the Attributes of OBPairData, and contain only letters, numbers and underscores.
    The predicates must start with a punctuation character and are interpreted by
    the Compare function of the OBDescriptor class. The default implementation expects
    a comparison operator and a number, e.g. >=1.3  Whitespace is optional and is ignored.
    Each predicate and this OBBase object (usually OBMol) is passed to 
    the Compare function of a OBDescriptor. The result of each comparison 
    is combined in a boolean expression (which can include parentheses)
    in the normal way. The AND operator can be & or &&, the OR operator can be 
    | or ||, and a unitary NOT is !  The expected operator precedence
    is achieved using recursive calls of the function. If there is no boolean Op, all
    the tests have to return true for the function to return true, i.e. the default is AND.
    If the first operand of an AND is 0, or of an OR is 1, the parsing of the second operand
    continues but no comparisons are done since the result does not matter.
 **/
bool OBDescriptor::FilterCompare(OBBase* pOb, std::istream& optionText, bool noEval)
{
  for(;;)
  {
    bool negate=false, retFromCompare, ret=false;
    char ch=0;
    optionText >> ch; //skips whitespace
    if(!optionText)
      return false;

    if(ch=='!')
    {
      negate=true;
      optionText >> ch;
    }

    if(ch=='(')
    {
      //bracketted expression
      retFromCompare = FilterCompare(pOb, optionText, noEval);//noEval persists in subsidiary calls
      optionText >> ch;
      if(ch!=')')
      {
        obErrorLog.ThrowError(__FUNCTION__, "Missing ')' in filter string", obError);
        return retFromCompare;
      }
    }
    else //unbracketted expression
    {
      if(!ispunctU(ch))
        optionText.unget(); //must be start of ID
      else
      {
        string mes("Filter string has erroneous character : ");
        obErrorLog.ThrowError(__FUNCTION__, mes + ch, obError);
        optionText.setstate(std::ios::badbit); //shows error
        return false;
      }

      string descID = GetIdentifier(optionText);
      if(descID.empty())
      {
        optionText.setstate(std::ios::badbit); //shows error
        return false;
      }
      
      //If there is existing OBPairData use that
      if(MatchPairData(pOb, descID))
      {
        string value = pOb->GetData(descID)->GetValue();
        retFromCompare = CompareStringWithFilter(optionText, value, noEval, true);
      }
      else
      {
        //if no existing data see if it is an OBDescriptor
        OBDescriptor* pDesc = OBDescriptor::FindType(descID.c_str());
        if(pDesc && !noEval)
          retFromCompare = pDesc->Compare(pOb, optionText, noEval);
        else
        {
          //just parse
          char ch1,ch2=0;
          string svalue;
          ParsePredicate(optionText, ch1, ch2, svalue);
          //no existing data, not a descriptor result is false meaning "does not exist"
          retFromCompare = false;
        }
      }
    }

    if(negate)
      retFromCompare=!retFromCompare;

    if(!noEval)
      ret = retFromCompare;

    //Look for boolean operator
    ch=0;
    if(!(optionText >> ch))
      return ret; //end of filterstring

    if(ch==')')
    {
      optionText.unget();
      return ret;
    }

    if(!ispunctU(ch))
      optionText.unget(); //start of next ID or )
    else
      if(optionText.peek()==ch) //treat && and || as & and |
        optionText.ignore();
        
    if(ch=='|')
    {
      retFromCompare = FilterCompare(pOb, optionText, ret || noEval);
      return !noEval && (ret || retFromCompare); //always return false if noEval=true;
    }
    else //includes & and , and ;
      noEval=!ret; //if ret is false keep parsing but don't bother to evaluate
  }//go for next conditional expression
  return false;//never come here
}

//////////////////////////////////////////////////////////////
string OBDescriptor::GetIdentifier(istream& optionText)
{
  string descID;
  descID.clear();
  char ch;
  optionText >> ch; //ignore leading white space
  optionText.unsetf(ios::skipws);
  for(;;)
  {
    if(!optionText || isspace(ch))
      break;
    if(ispunctU(ch))
    {
      optionText.unget(); //put back char after ID
      break;
    }
    descID.push_back(ch);
    optionText >> ch;
  }
  optionText.setf(ios::skipws);
  return descID;
}


///Reads comparison operator and the following string. Return its value if possible else NaN
//The comparison operator characters in ch1 and ch2 if found, 0 otherwise.
double OBDescriptor::ParsePredicate(istream& optionText, char& ch1, char& ch2, string& svalue)
{
  double val = std::numeric_limits<double>::quiet_NaN();
  ch2=0;
  ch1=0;
  //Get comparison operator
  optionText >> ch1;
  if(!ch1 || isalnum(ch1) || ch1=='&' || ch1=='|' || ch1==')')
  {
    //no comparison operator
    optionText.unget();
    optionText.clear(); //not an error to reach eof
    ch1=0;
    return val;
  }
  else
  {
    if(optionText.peek()=='=')
    optionText >> ch2;
  }

  //Try to read a double. Rewind and read as a string
  streampos spos = optionText.tellg();
  optionText >> val;
  optionText.clear();
  optionText.seekg(spos);
  ReadStringFromFilter(optionText, svalue);
  return val;
}


/// Reads a string from the filter string  optionally preceded by = or !=
/** On entry the stringstream position should be just after the ID. On exit it is after the string.
    If there is an error, the stringstream badbit is set.
    Returns false if != found, to indicate negation.
    Can be of any of the following forms:
    mystring  =mystring ==mystring [must be terminated by a space or tab]
    "mystring" 'mystring'  ="mystring" ='mystring' [mystring can contain spaces or tabs]
    !=mystring !="mystring" [Returns false indicating negate]
    There can be spaces or tabs after the operator = == !=  
 **/
bool OBDescriptor::ReadStringFromFilter(istream& optionText, string& result)
{
  bool error=false;
  bool ret=true;
  char ch;
  
  if(optionText >> ch)
  {
    if(ch=='=' || ch=='!')
    {
      if(optionText.get()!='=')
        optionText.unget();
      if(ch=='!')
        ret=false; //to indicate negation
    }
    else  //no operator
      optionText.unget();
  
    optionText >> ch;
    if(ch=='\"' || ch=='\'')
    {
      getline(optionText, result, ch); //get quoted text
    }
    else // not quoted; get string up to next space or ')'
    {
      optionText.unget();
      char ch;
      result.clear();
      optionText >> ch; //ignore leading white space
      optionText.unsetf(ios::skipws);
      for(;;)
      {
        if(!optionText || isspace(ch) || ch==')')
        {
          optionText.unget();
          optionText.clear();
          break;
        }
        result.push_back(ch);
        optionText >> ch;
      }
      optionText.setf(ios::skipws);
    }
  }
  
  if(optionText.fail())
    obErrorLog.ThrowError(__FUNCTION__, "Error reading string from filter", obError);
  
  return ret;
}

double OBDescriptor::PredictAndSave(OBBase* pOb)
{
  string attr = GetID();
  string svalue;
  double val = GetStringValue(pOb, svalue);

  OBPairData *dp = static_cast<OBPairData *> (pOb->GetData(attr));
  if (dp == NULL) // already set
    dp = new OBPairData;
  dp->SetAttribute(attr);
  dp->SetValue( svalue );
  dp->SetOrigin(perceived);
  pOb->SetData(dp);
  return val;
}

/// This default version provides a string representation of the numeric value
double OBDescriptor::GetStringValue(OBBase* pOb, string& svalue)
{
  double val = Predict(pOb);
  if(val!=std::numeric_limits<double>::quiet_NaN())
  {
    stringstream ss;
    ss << val;
    svalue = ss.str();
  }
  return val;
}

bool OBDescriptor::CompareStringWithFilter(istream& optionText, string& sval, bool noEval, bool NoCompOK)
{
  char ch1=0, ch2=0;
  string sfilterval;
  double filterval = ParsePredicate(optionText, ch1, ch2, sfilterval);
  if(ch1==0 && NoCompOK)
  {
    // there is no comparison operator
    return true; // means that the identifier exists
  }
  
  stringstream ss(sval);
  double val;
  if((ss >> val) && !IsNan(filterval))
    //Do a numerical comparison if both values are numbers
    return DoComparison(ch1, ch2, val, filterval);
  else
  {
    //Do a string comparison if either the filter or the OBPair value is not a number
    string::size_type pos = sfilterval.find('*');
    if(pos!=string::npos)
    {
      //filter string contains an asterisk
      //cannot match if string is too small
      if(sval.size() < sfilterval.size())
        return false;

      if(pos==sfilterval.size()-1)
      {
        // '*' is last char; delete it and subsequent chars
        sfilterval.erase(pos);
        sval.erase(pos);
      }
      else if(pos==0)
      {
        // '*' is first char; delete up to it
        sfilterval.erase(0,1);
        sval.erase(0,sval.size()-sfilterval.size());
      }
      else
      {
        // '*' is some other character. Not not currently supported.
        obErrorLog.ThrowError("--filter option", "Wild card * can only be the first or last character.", obError);
        return false;
      }
    }
    return DoComparison(ch1, ch2, sval, sfilterval);
  }
}

void OBDescriptor::AddProperties(OBBase* pOb, const string& DescrList)
{
  vector<string> vs;
  tokenize(vs, DescrList.c_str(), " \t\r\n,/-*&;:|%+");
  vector<string>::iterator itr;
  for(itr=vs.begin();itr!=vs.end();++itr)
  {
    OBDescriptor* pDescr = OBDescriptor::FindType(itr->c_str());
    if(pDescr)
      pDescr->PredictAndSave(pOb);
    else
      obErrorLog.ThrowError(__FUNCTION__, *itr + " not recognized as a descriptor", obError);
  }
}

void OBDescriptor::DeleteProperties(OBBase* pOb, const string& DescrList)
{
  vector<string> vs;
  tokenize(vs, DescrList.c_str(), " \t\r\n,/-*&;:|%+");
  vector<string>::iterator itr;
  for(itr=vs.begin();itr!=vs.end();++itr)
  {
    if(MatchPairData(pOb, *itr))
      pOb->DeleteData(*itr);
  }
}

  //Reads list of descriptor IDs and OBPairData names and returns a list of values
  //each precede by a space or the first character in the list if it is whitespace or punctuation.
  //Used in OBMol::Transform() to append to title , but that is not done here to avoid
  //having to #include mol.h in this file.
  string OBDescriptor::GetValues(OBBase* pOb, const std::string& DescrList)
  {
    vector<string> vs;
    char delim = DescrList[0];
    delim = isspace(delim)||ispunctU(delim) ? delim : ' ';

    tokenize(vs, DescrList.c_str(), " \t\r\n,/-*&;:|%+");
    vector<string>::iterator itr;
    string values;
    for(itr=vs.begin();itr!=vs.end();++itr)
    {
      string thisvalue;
      //If there is existing OBPairData use that
      if(MatchPairData(pOb,*itr))
        thisvalue = pOb->GetData(*itr)->GetValue();
      else
      {
        //if it is an OBDescriptor
        OBDescriptor* pDesc = OBDescriptor::FindType(itr->c_str());
        if(!pDesc) 
        {
          obErrorLog.ThrowError(__FUNCTION__, 
            *itr + " not recognized as either a property or a descriptor", obError);
          return values;
        }
        pDesc->GetStringValue(pOb, thisvalue);
      }
      values += delim + thisvalue;
    }
    return values;
  }

  bool OBDescriptor::MatchPairData(OBBase* pOb, string& s)
  {
    //If s matches a PairData attribute return true
    //else if s with all '_' replaced by spaces matches return true and s is now the form with spaces
    //else return false.
    if(pOb->HasData(s))
      return true;
    if(s.find('_')==string::npos)
      return false;
    string temp(s);
    string::size_type pos=-1;
    //replace all underscores by spaces
    while((pos=temp.find('_', ++pos))!=string::npos)
      temp[pos]=' ';
    if(pOb->HasData(temp))
    {
      s = temp;
      return true;
    }
    return false;
  }

bool OBDescriptor::Display(std::string&txt, const char* param, const char* ID)
{
  //Use the base class version except when the parameter is a descriptor ID.
  //For a paramater which is the matching descriptor set verbose.
  //No display for other descriptors.
  //Allows babel descriptors HBA1
  if(param  && FindType(param))
  {
    if(strcmp(ID, param))
      return false;
    param = "verbose";
  }
  return OBPlugin::Display(txt,param,ID);
}

}//namespace

/* 
OBDescriptor and Filtering
           
On the command line, using the option --filter filter-string converts only
those molecules which meet the criteria specified in the filter-string. This
is useful to select particular molecules from a set.
It is used like:
babel dataset.xxx outfile.smi --filter "MW>200 SMARTS!=c1ccccc1 PUBCHEM_CACTVS_ROTATABLE_BOND<5"

The identifier , "PUBCHEM_CACTVS_ROTATABLE_BOND" is the name of an attribute
of an OBPairData which has probably been imported from a property in a SDF
or CML file. The identifier names are (currently) case dependent. A comparison
is made with the value in the OBPairData. This is a numeric comparison if both
operands can be converted to numbers (as in the example). If the 5 had been
enclosed in single or double quotes the comparison would have been a string
comparison, which gives a different result in some cases. OBPairData is searched
first to match an identifier.

If there are no OBPair attributes that match, the identifier is taken to be the
ID of an OBDescriptor class object. The class OBDescriptor is the base class
for classes that wrap molecular properties, descriptors or features. In the example
"MW" and "SMARTS" are OBDescriptor IDs and are case independent. They are plugin
classes, like fingerprints, forcefields and formats, so that new molecular features
can be added or old ones removed (to prevent code bloat) without altering old code.
A list of available descriptors is available from the commandline:
babel -L descriptors
or from the functions OBPlugin::List, OBPlugin::ListAsString and OBPlugin::ListAsVector.

The filter-string is interpreted by a static function of OBDescriptor,
FilterCompare(). This identifies the descriptor IDs and then calls a virtual 
function, Compare(), of each OBDescriptor class to interpret the rest of relational
expression, for example, ">200", or "=c1ccccc1". The default version of Compare()
is suitable for descriptors, like MW or logP, which return a double from
their Predict() method. Classes like SMARTS which need different semantics
provide their own.

By default, as in the example, OBDescriptor::FilterCompare() would AND each
comparison so that all the comparisons must be true for the test to succeed.
However filter-string could also be a full boolean expression, with &, |, !,
and parenthases allowing any combination of features to be selected. 
FilterCompareAs calls itself recursively to give AND precidence over OR and
evaluation is not carried out if not needed.

The aim has been to make interpretation of the filter-string as liberal as
possible, so that AND can be &&, there can be spaces or commas in places
that are reasonable.

The base class, OBDescriptor, uses pointers to OBBase in its functions,
like OBFormat, to improve extendability - reactions could have
features too. It does mean that a dynamic_cast is needed at the start of the
Predict(OBBase* pOb) functions.

To use a particular descriptor, like logP, when programming with the API, use
code like the following:
\code
  OBDescriptor* pDescr = OBDecriptor::FindType("logP");
  if(pDescr)
    double val = pDescr->Predict(mol);
\endcode
To add the descriptor ID and the predicted data to OBPairData attached to
the object, use PredictAndSave().

This facility can be called from the command line.Use the option
--add "descriptor list", which will add the requested descriptors to the 
molecule.  They are then visible as properties in SDF and CML formats. 
The IDs in the list can be separated by spaces or punctuation characters.
All Descriptors will provide an output value as a string through a  virtual
function GetStringValue((OBBase* pOb, string& svalue)) which
assigns the value of a string descriptor(like inchi) to svalue or a string
representation of a numerical property like logP.

The classes MWFilter and TitleFilter illustrate the code that has to be 
provided for numerical and non-numerical descriptors.

*/ 
