/**********************************************************************
spacegroup.h - Handle Space Groups.
 
Copyright (C) 2007 by Jean Bréfort

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or 
modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_SPACE_GROUP_H
#define OB_SPACE_GROUP_H

#include <openbabel/math/transform3d.h>
#include <string>
#include <list>

namespace OpenBabel
{

  class OBAPI SpaceGroup
    {
      public:
        SpaceGroup();
        ~SpaceGroup();

        void SetHMName(const char *name)
          { m_HM = name; }
        void SetHMName(const std::string &name)
          { m_HM = name; }
        void SetHallName(const char *name)
          { m_Hall = name; }
        void SetHallName(const std::string &name)
          { m_Hall = name; }
        void SetId(unsigned n)
          { m_id = n; }
        void AddTransform(const std::string &s);

        const std::string & GetHMName() const
          { return m_HM;}
        const std::string & GetHallName()const
          { return m_Hall;} 
        unsigned GetId() const
          { return m_id; }
        std::list<vector3> Transform(const vector3 &v) const;

        transform3d const * BeginTransform(transform3dIterator &i) const;
        transform3d const * NextTransform(transform3dIterator &i) const;

        // static methods
        /* The name might be either a HM or Hall name */
        static const SpaceGroup * GetSpaceGroup (char const *name);
        static const SpaceGroup * GetSpaceGroup (const std::string &name);
        static const SpaceGroup * GetSpaceGroup (unsigned id);
        static const SpaceGroup * Find (SpaceGroup* group);
        /* Use it if the space group is unknown (might happen if no database has
         been loaded or if the HM name is not usual. */
        // Unfortunately, this seems to confuse the SWIG parser
        // Fortunately, it's not needed for the scripting bindings, 
        // since this is internal code
#ifndef SWIG
        void RegisterSpaceGroup (int nb = 0, ...);
#endif

        bool operator ==(const SpaceGroup &) const;
        int operator!=(const SpaceGroup &other) const
          {
            return !((*this) == other);
          }
        bool IsValid() const;

      private:
        std::string m_HM, m_Hall;
        unsigned int m_id;
        std::list<transform3d*> m_transforms;
    };

}

#endif // OB_SPACE_GROUP_H

