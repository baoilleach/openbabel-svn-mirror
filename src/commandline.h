/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef __COMMANDLINE_H_
#define __COMMANDLINE_H_

#include <stdlib.h>

#include <string>
#include <vector>
#include <map>

class ArgumentInfo
{
 public:
  int      argtyp;
  int      numargs;
  void    *argvar;
  bool     found;
  bool     list;
  bool     hasArgument;
  string   helpString;
  vector<string> arguments;

  ArgumentInfo() {
    numargs = 1;
    argvar  = NULL;
    argtyp  = -2;
    found   = false;
    list    = false;
    hasArgument = false;
    helpString= "";
  }

  // constructor with Help string 
  ArgumentInfo( const string & help) : argtyp(-2), numargs(1), argvar(NULL), 
                                found(false), list(false), hasArgument(false),
                                helpString(help) {}
};

typedef map<string, ArgumentInfo, less<string> > maptype;

class CommandLine
{
 protected:
  
  int argc, switches;
  char **argv;
  string arglist;

  maptype cmdline;
  void (*usage)(void);
  
  void SetVariable(ArgumentInfo &, char[]);
  bool _exitOnError;
 public:
  
  CommandLine() { switches = 0; usage = NULL; arglist = ""; _exitOnError= false; }
  CommandLine(bool exitOnError) { switches = 0; usage = NULL; arglist = ""; _exitOnError= exitOnError; }
  virtual ~CommandLine() {}

  // Switch/Flag with no Arguments
  
  void AddFlag( const char *, bool &var, bool list = false );
  void AddFlag( const char *, char *hlp, bool list = false );

  // Switch/Flag followed by Multiple Arguments

  void AddSwitch( const char *, int numargs,          bool list = false );
  void AddSwitch( const char *, vector<int>    &,     bool list = false );
  void AddSwitch( const char *, vector<float>  &,     bool list = false );
  void AddSwitch( const char *, vector<string> &,     bool list = false );
  void AddSwitch( const char *, void (*)(void *),     bool list = false );

  // Switch/Flag followed by Single Argument (with default value)

  void AddSwitch( const char *, bool   &, bool   def, bool list = false );
  void AddSwitch( const char *, int    &, int    def, bool list = false );
  void AddSwitch( const char *, float  &, float  def, bool list = false );
  void AddSwitch( const char *, char   *, char  *def, bool list = false );
  void AddSwitch( const char *, string &, string def, bool list = false );

  int ProcessArguments( char *argv[], int argc );
  int ProcessArguments( int argc, char *argv[] ) { return ProcessArguments(argv, argc); }

  bool WasCalledWith(const char *);
  bool HasArgument(const char *arg);

  char *GetArgument(const char *,int i = 0);
  char *GetArgList(void) { return (char *)arglist.c_str(); }

  void SetUsageFunction(void (*fxn)(void));
  virtual void Usage();

  // Add switched with help string
  void AddSwitch(const char *arg, string &var, string def, const char *help, bool list= false);
  void AddSwitch(const char *arg, int &var, int def, const char *help, bool list= false);
  void AddSwitch(const char *arg, float &var, float def, const char *help, bool list= false);

  void printHelp();

  bool getValue(const char *s, int &f);     
  bool getValue(const char *s, float &f);   
  bool getValue(const char *s, string &ret);
};

#endif

