/**********************************************************************
parsmart.h - SMART parser.

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (c) 2001-2003 by Geoffrey R. Hutchison

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

#ifndef OB_PARSMART_H
#define OB_PARSMART_H

#include <string>
#include <vector>

#include "mol.h"

/*==========================*/
/*  SMARTS Data Structures  */
/*==========================*/

#define AE_LEAF      0x01
#define AE_RECUR     0x02
#define AE_NOT       0x03
#define AE_ANDHI     0x04
#define AE_OR        0x05
#define AE_ANDLO     0x06

#define AL_CONST     0x01
#define AL_MASS      0x02
#define AL_AROM      0x03
#define AL_ELEM      0x04
#define AL_HCOUNT    0x05
#define AL_NEGATIVE  0x06
#define AL_POSITIVE  0x07
#define AL_CONNECT   0x08
#define AL_DEGREE    0x09
#define AL_IMPLICIT  0x0a
#define AL_RINGS     0x0b
#define AL_SIZE      0x0c
#define AL_VALENCE   0x0d
#define AL_CHIRAL    0x0e
#define AL_HYB       0x0f
#define AL_CLOCKWISE     1
#define AL_ANTICLOCKWISE 2

namespace OpenBabel
{

typedef union _AtomExpr {
        int type;
        struct {
            int type;
            int prop;
            int value;
        } leaf;
        struct {
            int type;
            void *recur;
        } recur;
        struct {
            int type;
            union _AtomExpr *arg;
        } mon;
        struct {
            int type;
            union _AtomExpr *lft;
            union _AtomExpr *rgt;
        } bin;
    } AtomExpr;

#define BE_LEAF      0x01
#define BE_ANDHI     0x02
#define BE_ANDLO     0x03
#define BE_NOT       0x04
#define BE_OR        0x05

#define BL_CONST     0x01
#define BL_TYPE      0x02

#define BT_SINGLE     0x01
#define BT_DOUBLE     0x02
#define BT_TRIPLE     0x03
#define BT_AROM       0x04
#define BT_UP         0x05
#define BT_DOWN       0x06
#define BT_UPUNSPEC   0x07
#define BT_DOWNUNSPEC 0x08
#define BT_RING       0x09

typedef union _BondExpr {
        int type;
        struct {
            int type;
            int prop;
            int value;
        } leaf;
        struct {
            int type;
            union _BondExpr *arg;
        } mon;
        struct {
            int type;
            union _BondExpr *lft;
            union _BondExpr *rgt;
        } bin;
    } BondExpr;

typedef struct {
        BondExpr *expr;
        int src,dst;
        int visit;
        bool grow;
    } BondSpec;

typedef struct {
  AtomExpr *expr;
  int visit;
  int part;
  int chiral_flag;
  int vb;
    } AtomSpec;

typedef struct {
  int aalloc,acount;
  int balloc,bcount;
  bool ischiral;
  AtomSpec *atom;
  BondSpec *bond;
  int parts;
} Pattern;

// class introduction in parsmart.cpp
class OBSmartsPattern
{
protected:
    std::vector<bool>          		_growbond;
    std::vector<std::vector<int> >	_mlist;
    Pattern				*_pat;
    std::string				_str;

public:
    OBSmartsPattern()  { _pat=NULL; }
    virtual ~OBSmartsPattern();

    OBSmartsPattern(const OBSmartsPattern& cp) { _pat = NULL; *this = cp; }
    OBSmartsPattern& operator=(const OBSmartsPattern& cp) 
    {
        if (_pat) 
            delete [] _pat; 
        _pat = NULL; 
        std::string s = cp._str; 
        Init(s); 
        return (*this);
    }

    unsigned int NumMatches() const { return (unsigned int)_mlist.size(); }
    unsigned int NumAtoms()   const { return _pat ? _pat->acount : 0;     }
    unsigned int NumBonds()   const { return _pat ? _pat->bcount : 0;     }

    int          GetAtomicNum(int);
    void         GetBond(int&,int&,int&,int);
    int          GetCharge(int);
    const std::string &GetSMARTS() const         {return _str;}
    std::string  &GetSMARTS()               {return _str;}
    int          GetVectorBinding(int idx) const { return(_pat->atom[idx].vb); }
    bool         Empty()                   const { return(_pat == NULL);       }
    bool         IsValid()                 const { return(_pat != NULL);       }
    bool         Init(const char*);
    bool         Init(const std::string&);
    void         WriteMapList(std::ostream&);

    bool Match(OBMol &mol, bool single=false);
    bool RestrictedMatch(OBMol &mol, std::vector<std::pair<int,int> > &pairs, bool single=false);
    bool RestrictedMatch(OBMol &mol, OBBitVec &bv, bool single=false);

    std::vector<std::vector<int> > &GetMapList() {return(_mlist);}
    std::vector<std::vector<int> > &GetUMapList();
    std::vector<std::vector<int> >::iterator BeginMList() {return(_mlist.begin());}
    std::vector<std::vector<int> >::iterator EndMList() {return(_mlist.end());}
};


class OBSSMatch //used for fast exhaustive matching
{
protected:
    bool        *_uatoms;
    OBMol       *_mol;
    Pattern     *_pat;
    std::vector<int>  _map;

public:
    OBSSMatch(OBMol&,Pattern*);
    ~OBSSMatch();
    void Match(std::vector<std::vector<int> > &v, int bidx=-1);
};

void SmartsLexReplace(std::string &,
		      std::vector<std::pair<std::string,std::string> > &);

} // end namespace OpenBabel

#endif // OB_PARSMART_H

