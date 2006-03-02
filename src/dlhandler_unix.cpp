/**********************************************************************
dlhandler_unix.cpp - Dynamic loader for UNIX (handles file format shared obj.)

Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2004-2005 by Geoffrey R. Hutchison
 
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

#include "dlhandler.h"
#include "babelconfig.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <dlfcn.h>

#include <iostream>

using namespace std;

namespace OpenBabel {
OBAPI bool tokenize(vector<string> &, const char *, const char *);
}

#ifndef BUFF_SIZE
#define BUFF_SIZE 32768
#endif

//Globals for scandir()

int matchFiles (SCANDIR_CONST struct dirent *entry_p)
{
  return (strstr(entry_p->d_name, DLHandler::getFormatFilePattern()) != 0);
}

bool DLHandler::getConvDirectory(string& convPath)
{
    //Need to provide the directory from which this shared library was loaded.
    //This is the default directory for format shared library files.
  
  string testPath;
  testPath += OB_MODULE_PATH;
  convPath = testPath;

  return true;
}

int DLHandler::findFiles (std::vector <std::string>& file_list,const std::string& pattern, const std::string& path)
{
  string currentPath;
  vector<string> paths, vs;
  char buffer[BUFF_SIZE];

  if (!path.empty())
    paths.push_back(path);

  if (getenv("BABEL_LIBDIR") != NULL)
    {
        strncpy(buffer,getenv("BABEL_LIBDIR"), BUFF_SIZE - 1);
	// add a trailing NULL just in case
	buffer[BUFF_SIZE] = '\0';

	OpenBabel::tokenize(vs, buffer, "\r\n\t :");

	if (vs.size() > 0)
	  {
	    for (unsigned int i = 0; i < vs.size(); i++)
	      paths.push_back(vs[i]);
	  }
    }

  if (paths.size() == 0)
    paths.push_back("./"); // defaults to current directory

  struct dirent **entries_pp;
  int count;

  for (unsigned int i = 0; i < paths.size(); i++)
    {
      currentPath = paths[i];
      count = scandir (currentPath.c_str(), &entries_pp, SCANDIR_T matchFiles, NULL);

      for(int i=0; i<count; i++)
	{
	  file_list.push_back(currentPath + getSeparator() + (entries_pp[i])->d_name);
	  free(entries_pp[i]);
	}
    }

  if (entries_pp)
    free(entries_pp);
  return count;
}

bool DLHandler::openLib(const string& lib_name)
{
    return dlopen(lib_name.c_str(), RTLD_LAZY) != 0;
}

const char* DLHandler::getFormatFilePattern()
{
    return MODULE_EXTENSION;
}

char DLHandler::getSeparator()
{
    return '/';
}

//! \file dlhandler_unix.cpp
//! \brief Dynamic loader for UNIX (handles file format shared obj.)
