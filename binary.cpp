/**********************************************************************
Copyright (C) 1998-2000 by OpenEye Scientific Software, Inc.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "mol.h"
#include "rotor.h"
#include "binary.h"
#include "oeutil.h"

#define OE_TITLE_SIZE 254
#define OE_BINARY_SETWORD 32
namespace OpenEye 
{
//test byte ordering
static int SINT = 0x00000001;
static unsigned char *STPTR = (unsigned char*)&SINT;
bool SwabInt = (STPTR[0]!=0);

void SetRotorToAngle(float *c,OEAtom **ref,float ang,vector<int> atoms);

int Swab(int i)
{
  unsigned char tmp[4],c;
  memcpy(tmp,(char*)&i,sizeof(int));
  c = tmp[0]; tmp[0] = tmp[3]; tmp[3] = c;
  c = tmp[1]; tmp[1] = tmp[2]; tmp[2] = c;
  memcpy((char*)&i,tmp,sizeof(int));

  return(i);
}

OERotamerList::~OERotamerList()
{
  vector<unsigned char*>::iterator i;
  for (i = _vrotamer.begin();i != _vrotamer.end();i++)
    delete [] *i;

  vector<pair<OEAtom**,vector<int> > >::iterator j;
  for (j = _vrotor.begin();j != _vrotor.end();j++)
    delete [] j->first;

  //Delete the interal base coordinate list
  unsigned int k;
  for (k=0 ; k<_c.size() ; k++) delete [] _c[k];
}

void OERotamerList::GetReferenceArray(unsigned char *ref)
{
  int j;
  vector<pair<OEAtom**,vector<int> > >::iterator i;
  for (j=0,i = _vrotor.begin();i != _vrotor.end();i++)
    {
      ref[j++] = (unsigned char)(i->first[0])->GetIdx();
      ref[j++] = (unsigned char)(i->first[1])->GetIdx();
      ref[j++] = (unsigned char)(i->first[2])->GetIdx();
      ref[j++] = (unsigned char)(i->first[3])->GetIdx();
    }
}

void OERotamerList::Setup(OEMol &mol,OERotorList &rl)
{
  //clear the old stuff out if necessary
  _vres.clear();
  vector<unsigned char*>::iterator j;
  for (j = _vrotamer.begin();j != _vrotamer.end();j++) delete [] *j;
  _vrotamer.clear();

  vector<pair<OEAtom**,vector<int> > >::iterator k;
  for (k = _vrotor.begin();k != _vrotor.end();k++)
    delete [] k->first;
  _vrotor.clear();

  //create the new list
  OERotor *rotor;
  vector<OERotor*>::iterator i;
  vector<int> children;

  int ref[4];
  OEAtom **atomlist;
  for (rotor = rl.BeginRotor(i);rotor;rotor = rl.NextRotor(i))
    {
      atomlist = new OEAtom* [4];
      rotor->GetDihedralAtoms(ref);
      atomlist[0] = mol.GetAtom(ref[0]); atomlist[1] = mol.GetAtom(ref[1]);
      atomlist[2] = mol.GetAtom(ref[2]); atomlist[3] = mol.GetAtom(ref[3]);
      mol.FindChildren(children,ref[1],ref[2]);
      _vrotor.push_back(pair<OEAtom**,vector<int> > (atomlist,children));
      _vres.push_back(rotor->GetResolution());
    }

  vector<float>::iterator n;
  vector<vector<float> >::iterator m;
  for (m = _vres.begin();m != _vres.end();m++)
    for (n = m->begin();n != m->end();n++)
      *n *= RAD_TO_DEG;
}

void OERotamerList::Setup(OEMol &mol,unsigned char *ref,int nrotors)
{
  //clear the old stuff out if necessary
  _vres.clear();
  vector<unsigned char*>::iterator j;
  for (j = _vrotamer.begin();j != _vrotamer.end();j++) delete [] *j;
  _vrotamer.clear();

  vector<pair<OEAtom**,vector<int> > >::iterator k;
  for (k = _vrotor.begin();k != _vrotor.end();k++)
    delete [] k->first;
  _vrotor.clear();

  //create the new list
  int i;
  vector<int> children;

  int refatoms[4];
  OEAtom **atomlist;
  for (i = 0;i < nrotors;i++)
    {
      atomlist = new OEAtom* [4];
      refatoms[0] = (int)ref[i*4  ];
      refatoms[1] = (int)ref[i*4+1];
      refatoms[2] = (int)ref[i*4+2];
      refatoms[3] = (int)ref[i*4+3];
      mol.FindChildren(children,refatoms[1],refatoms[2]);
      atomlist[0] = mol.GetAtom(refatoms[0]); 
      atomlist[1] = mol.GetAtom(refatoms[1]);
      atomlist[2] = mol.GetAtom(refatoms[2]);
      atomlist[3] = mol.GetAtom(refatoms[3]);
      _vrotor.push_back(pair<OEAtom**,vector<int> > (atomlist,children));
    }

}

void OERotamerList::AddRotamer(float *c)
{
  int idx,size;
  float angle,res=255.0f/360.0f;
  Vector v1,v2,v3,v4;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (char) 0;

  vector<pair<OEAtom**,vector<int> > >::iterator i;
  for (size=1,i = _vrotor.begin();i != _vrotor.end();i++,size++)
    {
      idx = (i->first[0])->GetCIdx(); v1.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[1])->GetCIdx(); v2.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[2])->GetCIdx(); v3.Set(c[idx],c[idx+1],c[idx+2]);
      idx = (i->first[3])->GetCIdx(); v4.Set(c[idx],c[idx+1],c[idx+2]);

      angle = CalcTorsionAngle(v1,v2,v3,v4);
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[size] = (unsigned char)rint(angle*res);
    }

  _vrotamer.push_back(rot);
}

void OERotamerList::AddRotamer(int *arr)
{
  unsigned int i;
  float angle,res=255.0f/360.0f;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (unsigned char)arr[0];

  for (i = 0;i < _vrotor.size();i++)
    {
      angle = _vres[i][arr[i+1]];
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[i+1] = (unsigned char)rint(angle*res);
    }
  _vrotamer.push_back(rot);
}

void OERotamerList::AddRotamer(unsigned char *arr)
{
  unsigned int i;
  float angle,res=255.0f/360.0f;

  unsigned char *rot = new unsigned char [_vrotor.size()+1];
  rot[0] = (unsigned char)arr[0];

  for (i = 0;i < _vrotor.size();i++)
    {
      angle = _vres[i][(int)arr[i+1]];
      while (angle < 0.0f)   angle += 360.0f;
      while (angle > 360.0f) angle -= 360.0f;
      rot[i+1] = (unsigned char)rint(angle*res);
    }
  _vrotamer.push_back(rot);
}

void OERotamerList::AddRotamers(unsigned char *arr,int nrotamers)
{
  int i,size=_vrotor.size()+1;

  for (i = 0;i < nrotamers;i++)
    {
      unsigned char *rot = new unsigned char [size];
      memcpy(rot,&arr[i*size],sizeof(char)*size);
      _vrotamer.push_back(rot);
    }
}

void OERotamerList::ExpandConformerList(OEMol &mol,vector<float*> &clist)
{
  unsigned int j;
  float angle,invres=360.0f/255.0f;
  unsigned char *conf;
  vector<float*> tmpclist;
  vector<unsigned char*>::iterator i;

  for (i = _vrotamer.begin();i != _vrotamer.end();i++)
    {
      conf = *i;
      float *c = new float [mol.NumAtoms()*3];
      memcpy(c,clist[(int)conf[0]],sizeof(float)*mol.NumAtoms()*3);

      for (j = 0;j < _vrotor.size();j++)
	{
	  angle = invres*((float)conf[j+1]);
	  if (angle > 180.0) angle -= 360.0;
	  SetRotorToAngle(c,_vrotor[j].first,angle,_vrotor[j].second);
	}
      tmpclist.push_back(c);
    }

  //transfer the conf list
  vector<float*>::iterator k;
  for (k = clist.begin();k != clist.end();k++)
    delete [] *k;
  clist = tmpclist;
}

//Create a conformer list using the internal base set of coordinates
vector<float*> OERotamerList::CreateConformerList(OEMol& mol)
  {
    unsigned int j;
    float angle,invres=360.0f/255.0f;
    unsigned char *conf;
    vector<float*> tmpclist;
    vector<unsigned char*>::iterator i;
   
    for (i = _vrotamer.begin();i != _vrotamer.end();i++)
      {
        conf = *i;
        float *c = new float [mol.NumAtoms()*3];
        memcpy(c,_c[(int)conf[0]],sizeof(float)*mol.NumAtoms()*3);
   
        for (j = 0;j < _vrotor.size();j++)
          {
            angle = invres*((float)conf[j+1]);
            if (angle > 180.0) angle -= 360.0;
            SetRotorToAngle(c,_vrotor[j].first,angle,_vrotor[j].second);
          }
        tmpclist.push_back(c);
      }
   
    return tmpclist;
  }

//Copies the coordinates in bc, NOT the pointers, into the object
void OERotamerList::SetBaseCoordinateSets(vector<float*> bc, unsigned int N)
  {
    unsigned int i,j;

    //Clear out old data
    for (i=0 ; i<_c.size() ; i++) delete [] _c[i];
    _c.clear();

    //Copy new data
    float *c = NULL;
    float *cc= NULL;
    for (i=0 ; i<bc.size() ; i++) {
        c = new float [3*N];
        cc = bc[i];
        for (j=0 ; j<3*N ; j++) c[j] = cc[j];
        _c.push_back(c);
      }
    _NBaseCoords = N;
  }



int PackCoordinate(float c[3],float max[3])
{
  int tmp;
  tmp  = ((int)(c[0]*max[0])) << 20;
  tmp |= ((int)(c[1]*max[1])) << 10;
  tmp |= ((int)(c[2]*max[2]));
  return(tmp);
}

void UnpackCoordinate(float c[3],float max[3],int tmp)
{
  c[0] = (float)(tmp>>20);            c[0] *= max[0];
  c[1] = (float)((tmp&0xffc00)>>10);  c[1] *= max[1];
  c[2] = (float)(tmp&0x3ff);          c[2] *= max[2];
}

bool WriteBinary(ostream &ofs,OEMol &mol)
{
	/*
  if (mol.NumAtoms() >= 255 || mol.NumBonds() >= 255)
    {
      string s = "Unable to write molecule '";
      s += mol.GetTitle();
      s += "' to binary file - too many atoms or bonds";
      ThrowError("Unable to write");
      return(false);
    }
	*/
  int tmp,size;
  unsigned char buf[1000000];
  mol.SetOutputType(OEBINARY);
  WriteBinary(buf,size,mol);
  tmp = size;
  if (SwabInt) tmp = Swab(tmp);
  ofs.write((char*)&tmp,sizeof(int));
  ofs.write((char*)buf,size);

  return(true);
}

bool WriteBinary(unsigned char *buf,int &size,OEMol &mol)
{
  int m,tmp,idx;
  unsigned int k;
  OEAtom *atom;
  vector<OEAtom*>::iterator i;
  vector<float*>::iterator j;
  idx=0;

  //read title first
  int len = strlen(mol.GetTitle());
  if (len > OE_TITLE_SIZE) len = OE_TITLE_SIZE;
  if (len > 0)
    {
      buf[idx] = (char)len;
      idx += sizeof(char);
      memcpy(&buf[idx],mol.GetTitle(),sizeof(char)*len);
      idx += len;
    }
  else
    {
      buf[idx]=(char)0;
      idx += sizeof(char);
    }
  
  unsigned char c;
  tmp = (mol.NumAtoms() << 16) | mol.NumBonds();
  if (SwabInt) tmp = Swab(tmp);
  memcpy(&buf[idx],(const char*)&tmp,sizeof(int));idx += sizeof(int);

  for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
      c = (unsigned char)atom->GetAtomicNum();
      memcpy(&buf[idx],(const unsigned char*)&c,sizeof(unsigned char));
      idx += sizeof(unsigned char);
    }

  OEBond *bond;
  vector<OEBond*>::iterator bi;
  unsigned char bc[3];
  for (bond = mol.BeginBond(bi);bond;bond = mol.NextBond(bi))
    {
      bc[0] = (unsigned char)bond->GetBeginAtomIdx();
      bc[1] = (unsigned char)bond->GetEndAtomIdx();
      bc[2] = (unsigned char)bond->GetBO();
      memcpy(&buf[idx],(const unsigned char*)bc,sizeof(unsigned char)*3); 
      idx += sizeof(unsigned char)*3;
    }

  //Write out conformers and coordinates
  OERotamerList *rml = (OERotamerList *)mol.GetData(oeRotamerList);
  //find min and max
  int imin[3],imax[3];
  float min[3] = {10E10f,10E10f,10E10f};
  float max[3] = {-10E10f,-10E10f,-10E10f};
  vector<float*> clist;

  //If we have a rotamer list with internally stored base coordinates
  //have clist point to them rather than the molecules conformer coordinates
  if (rml && ((rml) ? rml->NumRotamers() : 0) && rml->NumBaseCoordinateSets()) {
      if (rml->NumAtoms() == mol.NumAtoms()) {//Error check, these should match
          for (k=0 ; k<rml->NumBaseCoordinateSets() ; k++) 
              clist.push_back(rml->GetBaseCoordinateSet(k));
        }
      else clist = mol.GetConformers();
    }
  else clist = mol.GetConformers();
 
  for (j = clist.begin();j != clist.end();j++)
    for (k = 0;k < mol.NumAtoms();k++)
      {
        if ((*j)[k*3  ] < min[0]) min[0] = (*j)[k*3  ];
        if ((*j)[k*3+1] < min[1]) min[1] = (*j)[k*3+1];
        if ((*j)[k*3+2] < min[2]) min[2] = (*j)[k*3+2];
        if ((*j)[k*3  ] > max[0]) max[0] = (*j)[k*3  ];
        if ((*j)[k*3+1] > max[1]) max[1] = (*j)[k*3+1];
        if ((*j)[k*3+2] > max[2]) max[2] = (*j)[k*3+2];
      }
 
  //store integer versions of min and max
  for (k = 0;k < 3;k++)
    {
      max[k] -= min[k];
      imin[k] = (int) (1000000.0f*min[k]);
      imax[k] = (int) (1000000.0f*max[k]);
      if (SwabInt)
        {
          imin[k] = Swab(imin[k]);
          imax[k] = Swab(imax[k]);
        }
    }
 
  //write the min and max
  memcpy((unsigned char *)&buf[idx], (const char*)imin,sizeof(int)*3); idx += sizeof(int)*3;
  memcpy((unsigned char *)&buf[idx], (const char*)imax,sizeof(int)*3); idx += sizeof(int)*3;
 
  //quantize max for packing coordinates
  for (k = 0;k < 3;k++) max[k] = (fabs(max[k])> 0.01) ? 1023.0f/max[k]:0.0;

  //write the number of confs and rotamers
  tmp = clist.size(); if (SwabInt) tmp = Swab(tmp);
  memcpy((unsigned char *)&buf[idx],
         (char*)&tmp,sizeof(int)); idx += sizeof(int);
  tmp = (rml) ? rml->NumRotamers() : 0; if (SwabInt) tmp = Swab(tmp);
  memcpy((unsigned char *)&buf[idx],
         (char*)&tmp,sizeof(int)); idx += sizeof(int);

  //The third boolean expression in the next if statement is an error check to make
  //sure than the number of atoms in the OERotamerLists internal coordinates are the
  //same as the molecules IF we are using the OERotamerLists internal coordinates.
  //If we are using the OERotamerList's internal coordinates but the number of atoms
  //in those coordinates don't match the number of atoms in the molecule then  the
  //rotamer list is incorrect for the molecule and we just default to writing the
  //molecules conformers.  I put this in because it strikes me as a very easy error
  //to make if the molecule is modified in any way and the user is not aware that
  //he/she is responsible for correctly updating the OERotamerList MM 4/20/01
  if (rml && ((rml) ? rml->NumRotamers() : 0) && (rml->NumBaseCoordinateSets()==0 || rml->NumAtoms() == mol.NumAtoms()) ) {//Store conformers as torsion list
      //Write base coordinates
      float tc[3];
      for (j = clist.begin();j != clist.end();j++)
        for (k = 0;k < mol.NumAtoms();k++)
          {
            for (m = 0;m < 3;m++) tc[m] = (*j)[k*3+m]-min[m];
            tmp = PackCoordinate(tc,max);
            if (SwabInt) tmp = Swab(tmp);
            memcpy((unsigned char *)&buf[idx],
                   (const char*)&tmp,sizeof(int)); idx += sizeof(int);
          }

      //Write rotors
      tmp = rml->NumRotors(); 
      if (SwabInt) tmp = Swab(tmp);
      memcpy((unsigned char *)&buf[idx],
             (const char*)&tmp,sizeof(int)); idx += sizeof(int);
      unsigned char *ref = new unsigned char [rml->NumRotors()*4];
      rml->GetReferenceArray(ref);
      memcpy((unsigned char *)&buf[idx],(const unsigned char*)ref,
             sizeof(unsigned char)*rml->NumRotors()*4);
      idx += sizeof(unsigned char)*rml->NumRotors()*4;
      delete [] ref;

      //Write Rotamers
      vector<unsigned char*>::iterator k;
      for (k = rml->BeginRotamer();k != rml->EndRotamer();k++)
        {
	  memcpy((unsigned char *)&buf[idx],(const unsigned char*)*k,
                 sizeof(unsigned char)*(rml->NumRotors()+1));
          idx += sizeof(char)*(rml->NumRotors()+1);
        }
    }
  else if (mol.NumConformers() > 1) {//Store conformers as coordinates
      //Write coordinate conformers 
      float tc[3];
      for (j = clist.begin();j != clist.end();j++)
        for (k = 0;k < mol.NumAtoms();k++)
          {
            for (m = 0;m < 3;m++) tc[m] = (*j)[k*3+m]-min[m];
            tmp = PackCoordinate(tc,max);
            if (SwabInt) tmp = Swab(tmp);
            memcpy((unsigned char *)&buf[idx],
                   (const char*)&tmp,sizeof(int)); idx += sizeof(int);
          }
    }
  else //must be storing single-conformer structure
    {
      //write the coordinates
      float coord[3];
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
	{
	  (atom->GetVector()).Get(coord);
	  for (k = 0;k < 3;k++) coord[k] -= min[k];
	  tmp = PackCoordinate(coord,max);
	  if (SwabInt) tmp = Swab(tmp);
	  memcpy((unsigned char *)&buf[idx],
		 (const char*)&tmp,sizeof(int)); idx += sizeof(int);
	}  
    }

  int nwords,word,bit;
  unsigned int *arobits;

  if (mol.NumAtoms()) //set bits on for aromatic atoms
    {
      nwords = mol.NumAtoms()/OE_BINARY_SETWORD;
      if (mol.NumAtoms()%OE_BINARY_SETWORD) nwords++;
      arobits = new unsigned int [nwords];
      memset((char*)arobits,'\0',sizeof(int)*nwords);

      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
	if (atom->IsAromatic())
	{
	  word = (atom->GetIdx()-1)/OE_BINARY_SETWORD;
	  bit = (atom->GetIdx()-1)%OE_BINARY_SETWORD;
	  arobits[word] |= (1<<bit);
	}
      
      if (SwabInt)
	for (m = 0;m < nwords;m++)
	  arobits[m] =  Swab(arobits[m]);
      memcpy(&buf[idx],(const char*)arobits,sizeof(int)*nwords);idx += sizeof(int)*nwords;
      delete [] arobits;
    }

  if (mol.NumBonds()) //set bits on for aromatic bonds
    {
      nwords = mol.NumBonds()/OE_BINARY_SETWORD;
      if (mol.NumBonds()%OE_BINARY_SETWORD) nwords++;

      arobits = new unsigned int [nwords];
      memset((char*)arobits,'\0',sizeof(int)*nwords);

      for (bond = mol.BeginBond(bi);bond;bond = mol.NextBond(bi))
	if (bond->IsAromatic())
	  {
	    word = (bond->GetIdx())/OE_BINARY_SETWORD;
	    bit = (bond->GetIdx())%OE_BINARY_SETWORD;
	    arobits[word] |= (1<<bit);
	  }
      
      if (SwabInt)
	for (m = 0;m < nwords;m++)
	  arobits[m] =  Swab(arobits[m]);
      memcpy(&buf[idx],(const char*)arobits,sizeof(int)*nwords);idx += sizeof(int)*nwords;
      delete [] arobits;
    }

  //Write pose information
  
  //Number of poses
  unsigned int numposes = mol.NumPoses();
  idx += OE_io_write_binary((char*)&buf[idx],(char*)&numposes, sizeof(unsigned int), 1); 

  //Specify a version number for the poses
  unsigned short int pose_version=0;
  idx += OE_io_write_binary((char*)&buf[idx],(char*)&pose_version,sizeof(unsigned short int), 1);
  
  for (k=0 ; k<mol.NumPoses() ; k++) 
	  idx += mol.GetPose(k).WriteBinary((char*)&buf[idx]); //Each pose

  size = idx;
  return(true);
}

bool ReadBinary(istream &ifs,OEMol &mol)
{
  int size = 0;
  unsigned char buf[1000000];
  if (!ifs.read((char*)&size,sizeof(int))) return(false);
  if (SwabInt) size = Swab(size);
  if (!ifs.read((char*)buf,sizeof(char)*size)) return(false);
  ReadBinary(buf,mol,size);

  return(true);
}

bool ReadBinary(istream &ifs, unsigned char **bin)
{
  int size = 0;
  unsigned char buf[100000];

  oeAssert(bin != NULL);

#ifdef __sgi

  if (!ifs.read((char*)&size,sizeof(int))) return(false);
  if (SwabInt) size = Swab(size);
  if (!ifs.read((char*)buf,sizeof(char)*size)) return(false);

  *bin = new unsigned char[sizeof(int) + (sizeof(char) * size)];
 
  memcpy(*bin, &size, sizeof(int));
  memcpy(*bin + sizeof(int), &buf[0], (sizeof(char) * size));

#else

  int temp = 0;

  if (!ifs.read((char*)&temp,sizeof(int))) return(false);
  if (SwabInt) size = Swab(temp);
  if (!ifs.read((char*)buf,sizeof(char)*size)) return(false);

  *bin = new unsigned char[sizeof(int) + (sizeof(char) * size)];
 
  memcpy(*bin, &temp, sizeof(int));
  memcpy(*bin + sizeof(int), &buf[0], (sizeof(char) * size));

#endif

  return(true);
}

bool ReadBinary(unsigned char *buf,OEMol &mol,int size)
{
  int i,j,k,idx,natoms,nbonds,tmp;
  char title[OE_TITLE_SIZE+1]; 
  idx = 0;

  //read title
  i = (int)buf[0];
  idx += sizeof(char);
  if (i > 0)
    {
      memcpy(title,&buf[idx],sizeof(char)*i);
      title[i] = '\0';
      idx += i;
    }
  else
    strcpy(title,"****");

  //readnumber of atoms and bonds
  memcpy(&tmp,(unsigned char*)&buf[idx],sizeof(int)); idx += sizeof(int);
  if (SwabInt) tmp = Swab(tmp);

  natoms = (tmp >> 16);
  nbonds = tmp & 0xffff;

  unsigned char *anum = new unsigned char [natoms];
  memcpy((unsigned char*)anum,
	 &buf[idx],sizeof(unsigned char)*natoms); 
  idx += sizeof(unsigned char)*natoms;
  
  mol.BeginModify();
  //read atom data
  OEAtom atom;
  for (i = 0;i < natoms;i++)
    {
      atom.SetAtomicNum((int)anum[i]);
      atom.SetType(etab.GetSymbol((int)anum[i]));
      mol.AddAtom(atom);
      atom.Clear();
    }
  delete [] anum;

  //read bond data
  int start,end,order;
  unsigned char *bnd = new unsigned char [nbonds*3];
  memcpy((unsigned char*)bnd,
	 &buf[idx],sizeof(unsigned char)*3*nbonds); 
  idx += sizeof(unsigned char)*3*nbonds;
  for (i = 0;i < nbonds;i++)
    {
      start = bnd[i*3  ];
      end   = bnd[i*3+1];
      order = bnd[i*3+2];
      mol.AddBond(start,end,order);
    }
  delete [] bnd;

  mol.EndModify();

  //read the min and max
  int imin[3],imax[3];
  float min[3],max[3];
  memcpy((char*)imin,&buf[idx],sizeof(int)*3); idx += sizeof(int)*3;
  memcpy((char*)imax,&buf[idx],sizeof(int)*3); idx += sizeof(int)*3;

  //unpack min and max
  for (i = 0;i < 3;i++)
    {
      if (SwabInt)
	{
	  imin[i] = Swab(imin[i]);
	  imax[i] = Swab(imax[i]);
	}
      min[i] = (float) imin[i]/1000000.0f;
      max[i] = (float) imax[i]/1000000.0f;
      max[i] = (fabs(max[i])> 0.01) ? max[i]/1023.0f:0.0; 
    }

  //read conformer information if available
  int nconfs,rotmrs;
  memcpy((char*)&nconfs,&buf[idx],sizeof(int)); idx += sizeof(int);
  memcpy((char*)&rotmrs,&buf[idx],sizeof(int)); idx += sizeof(int);
  if (SwabInt)
    {
      nconfs = Swab(nconfs);
      rotmrs = Swab(rotmrs);
    }

  if (nconfs == 1 && !rotmrs) //only a single conformer
    {
      Vector v;
      int *tmpi = new int [natoms];
      memcpy((char*)tmpi,&buf[idx],sizeof(int)*natoms); idx += sizeof(int)*natoms;
      float coord[3];
      for (i = 0;i < natoms;i++)
	{
	  if (SwabInt) tmpi[i] = Swab(tmpi[i]);
	  UnpackCoordinate(coord,max,tmpi[i]);
	  for (j = 0;j < 3;j++) coord[j] += min[j];
	  v.Set(coord);
	  (mol.GetAtom(i+1))->SetVector(v);
	}
      delete [] tmpi;
    }
  else
    {
      int *tmpi = new int [natoms];
      vector<float*> cltmp;
      for (i = 0;i < nconfs;i++)
	{
	  memcpy((char*)tmpi,(unsigned char*)&buf[idx],sizeof(int)*natoms); 
	  idx += sizeof(int)*natoms;
	  float *coord = new float [mol.NumAtoms()*3];
	  for (j = 0;j < natoms;j++) 
	    {
	      if (SwabInt) tmpi[j] = Swab(tmpi[j]);
	      UnpackCoordinate(&coord[j*3],max,tmpi[j]);
	      for (k = 0;k < 3;k++) coord[j*3+k] += min[k];
	    }
	  cltmp.push_back(coord);
	}

      delete [] tmpi;

      if (!rotmrs)
	mol.SetConformers(cltmp);
      else
	{
       
	  int nrotors;
	  OERotamerList *rml = new OERotamerList;
	  memcpy((char*)&nrotors,&buf[idx],sizeof(int)); idx += sizeof(int);
	  if (SwabInt) nrotors = Swab(nrotors);

	  unsigned char *ref = new unsigned char [nrotors*4];
	  memcpy((unsigned char*)ref,&buf[idx],sizeof(unsigned char)*nrotors*4);
	  idx += sizeof(unsigned char)*nrotors*4;
	  rml->Setup(mol,ref,nrotors);
	  delete [] ref;
	  
	  unsigned char *rotamers = new unsigned char [(nrotors+1)*rotmrs];
	  memcpy((unsigned char*)rotamers,
		 &buf[idx],sizeof(unsigned char)*(nrotors+1)*rotmrs);
	  idx += sizeof(unsigned char)*(nrotors+1)*rotmrs;
	  rml->AddRotamers(rotamers,rotmrs);
	  delete [] rotamers;
	 
          //Copy the base coordinate list into the OERotamerList object
          rml->SetBaseCoordinateSets(cltmp,mol.NumAtoms());
 
	  //expand rotamer information to a conformer list
	  rml->ExpandConformerList(mol,cltmp);
	  mol.SetConformers(cltmp);

          //Add the OERotamerList to the molecule as user data
          mol.SetData(rml);
	} // end else !rotmrs
		}  // end else nconf==1 && !rotmrs

  mol.SetTitle(title);

  if (idx >= size) return(true);
  
  int nwords;
  unsigned int *arobits;

  if (mol.NumAtoms()) //set bits on for aromatic atoms
    {
      nwords = mol.NumAtoms()/OE_BINARY_SETWORD;
      if (mol.NumAtoms()%OE_BINARY_SETWORD) nwords++;
      arobits = new unsigned int [nwords];

      memcpy((unsigned char*)arobits,&buf[idx],sizeof(int)*nwords);
      idx += sizeof(int)*nwords;

      if (SwabInt)
	for (i = 0;i < nwords;i++) arobits[i] =  Swab(arobits[i]);

      for (i = 0;i < (signed)mol.NumAtoms();i++)
	if ((arobits[i/OE_BINARY_SETWORD]>>(i%OE_BINARY_SETWORD))&1)
	  mol.GetAtom(i+1)->SetAromatic();

      delete [] arobits;
    }

  if (mol.NumBonds()) //set bits on for aromatic atoms
    {
      nwords = (mol.NumBonds()/OE_BINARY_SETWORD);
      if (mol.NumBonds()%OE_BINARY_SETWORD) nwords++;
      arobits = new unsigned int [nwords];

      memcpy((unsigned char*)arobits,&buf[idx],sizeof(int)*nwords);
      idx += sizeof(int)*nwords;

      if (SwabInt)
	for (i = 0;i < nwords;i++) arobits[i] =  Swab(arobits[i]);

      for (i = 0;i < (signed)mol.NumBonds();i++)
	if ((arobits[i/OE_BINARY_SETWORD]>>(i%OE_BINARY_SETWORD))&1)
	  mol.GetBond(i)->SetAromatic();

      delete [] arobits;
    }

  mol.SetAromaticPerceived();

  //Read in poses if present
  if (idx>=size) return(true);  //Backwards compatibility 
  unsigned int kk;
  mol.DeletePoses();
  unsigned int Nposes=0;
  OEPose pose;
  idx += OE_io_read_binary((char*)&buf[idx],(char*)&Nposes,sizeof(unsigned int), 1); //Read number of poses
  unsigned short int pose_version;
  idx += OE_io_read_binary((char*)&buf[idx],(char*)&pose_version,sizeof(unsigned short int), 1); //Read the version number
  if (pose_version == 0) {
      for (kk=0 ; kk<Nposes ; kk++) { //Read in the poses
          idx += pose.ReadBinary((char*)&buf[idx]);
          mol.AddPose(pose);
        }
    }
  else {
      cerr << "ERROR! in OEMol binary reader, pose version not supported" << endl;
      return false;
    }

  return(true);
}

void SetRotorToAngle(float *c,OEAtom **ref,float ang,vector<int> atoms)
     //this function will rotate the coordinates of 'atoms'
     //such that tor == ang - atoms in 'tor' should be ordered such 
     //that the 3rd atom is the pivot around which atoms rotate
     //ang is in degrees
{
  float v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  float c1x,c1y,c1z,c2x,c2y,c2z,c3x,c3y,c3z;
  float c1mag,c2mag,radang,costheta,m[9];
  float x,y,z,mag,rotang,sn,cs,t,tx,ty,tz;

  int tor[4];
  tor[0] = ref[0]->GetCIdx();
  tor[1] = ref[1]->GetCIdx();
  tor[2] = ref[2]->GetCIdx();
  tor[3] = ref[3]->GetCIdx();

  //
  //calculate the torsion angle
  //
  v1x = c[tor[0]]   - c[tor[1]];   v2x = c[tor[1]]   - c[tor[2]];
  v1y = c[tor[0]+1] - c[tor[1]+1]; v2y = c[tor[1]+1] - c[tor[2]+1];
  v1z = c[tor[0]+2] - c[tor[1]+2]; v2z = c[tor[1]+2] - c[tor[2]+2];
  v3x = c[tor[2]]   - c[tor[3]];
  v3y = c[tor[2]+1] - c[tor[3]+1];
  v3z = c[tor[2]+2] - c[tor[3]+2];

  c1x = v1y*v2z - v1z*v2y;   c2x = v2y*v3z - v2z*v3y;
  c1y = -v1x*v2z + v1z*v2x;  c2y = -v2x*v3z + v2z*v3x;
  c1z = v1x*v2y - v1y*v2x;   c2z = v2x*v3y - v2y*v3x;
  c3x = c1y*c2z - c1z*c2y;
  c3y = -c1x*c2z + c1z*c2x;
  c3z = c1x*c2y - c1y*c2x; 
  
  c1mag = SQUARE(c1x)+SQUARE(c1y)+SQUARE(c1z);
  c2mag = SQUARE(c2x)+SQUARE(c2y)+SQUARE(c2z);
  if (c1mag*c2mag < 0.01) costheta = 1.0; //avoid div by zero error
  else costheta = (c1x*c2x + c1y*c2y + c1z*c2z)/(sqrt(c1mag*c2mag));

  if (costheta < -0.999999) costheta = -0.999999f;
  if (costheta >  0.999999) costheta =  0.999999f;
			      
  if ((v2x*c3x + v2y*c3y + v2z*c3z) > 0.0) radang = -acos(costheta);
  else                                     radang = acos(costheta);

  //
  // now we have the torsion angle (radang) - set up the rot matrix
  //

  //find the difference between current and requested
  rotang = (DEG_TO_RAD*ang) - radang; 

  sn = sin(rotang); cs = cos(rotang);t = 1 - cs;
  //normalize the rotation vector
  mag = sqrt(SQUARE(v2x)+SQUARE(v2y)+SQUARE(v2z));
  x = v2x/mag; y = v2y/mag; z = v2z/mag;
  
  //set up the rotation matrix
  m[0]= t*x*x + cs;     m[1] = t*x*y + sn*z;  m[2] = t*x*z - sn*y;
  m[3] = t*x*y - sn*z;  m[4] = t*y*y + cs;    m[5] = t*y*z + sn*x;
  m[6] = t*x*z + sn*y;  m[7] = t*y*z - sn*x;  m[8] = t*z*z + cs;

  //
  //now the matrix is set - time to rotate the atoms
  //
  tx = c[tor[1]];ty = c[tor[1]+1];tz = c[tor[1]+2];
  vector<int>::iterator i;int j;
  for (i = atoms.begin();i != atoms.end();i++)
    {
      j = ((*i)-1)*3;
      c[j] -= tx;c[j+1] -= ty;c[j+2]-= tz;
      x = c[j]*m[0] + c[j+1]*m[1] + c[j+2]*m[2];
      y = c[j]*m[3] + c[j+1]*m[4] + c[j+2]*m[5];
      z = c[j]*m[6] + c[j+1]*m[7] + c[j+2]*m[8];
      c[j] = x; c[j+1] = y; c[j+2] = z;
      c[j] += tx;c[j+1] += ty;c[j+2] += tz;
    }
}

//OEBinaryDBase class - facilitates random access to OEBinary files

OEBinaryDBase::OEBinaryDBase(char *fname)
{
  int size;
  streampos pos;
  unsigned char buf[100000];

  if (!SafeOpen(_ifs,fname)) exit(0);

  for (;;)
    {
      pos = _ifs.tellg();
      if (!_ifs.read((char*)&size,sizeof(int))) break;
      if (SwabInt) size = Swab(size);
      if (!_ifs.read((char*)buf,sizeof(char)*size)) break;
      _vpos.push_back(pos);
    }
  _ifs.close();

  if (!SafeOpen(_ifs,fname)) exit(0);
}

OEBinaryDBase::OEBinaryDBase(string &fname)
{
  int size;
  streampos pos;
  unsigned char buf[100000];

  if (!SafeOpen(_ifs,(char*)fname.c_str())) exit(0);

  for (;;)
    {
      pos = _ifs.tellg();
      if (!_ifs.read((char*)&size,sizeof(int))) break;
      if (SwabInt) size = Swab(size);
      if (!_ifs.read((char*)buf,sizeof(char)*size)) break;
      _vpos.push_back(pos);
    }
  _ifs.close();

  if (!SafeOpen(_ifs,(char*)fname.c_str())) exit(0);
}

int OEBinaryDBase::Size()
{
  return(_vpos.size());
}

void OEBinaryDBase::GetMolecule(OEMol &mol,int idx)
{
  OEFileFormat ff;
  mol.Clear();
  mol.SetInputType(OEBINARY);
  _ifs.seekg(_vpos[idx]);
  ff.ReadMolecule(_ifs,mol);

}


} //namespace OpenEye

