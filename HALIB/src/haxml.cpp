/*! \file haxml.cpp

   Classes for to work with XML files in HARLEM

   \author Nikolay Simakov, Igor Kurnikov
   \date 2006-
*/

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <zlib.h>
#include <tinyxml.h>

#include "haobject.h"
#include "haxml.h"


//#include "expatxmlparser.h"

//ExpatXMLParser* HaObject::GetExpatXMLParser()
//{
//	ExpatXMLParser* p_xml_parser = new ExpatXMLParser();
//	return p_xml_parser;
//}

HaAutoTreeObject::HaAutoTreeObject(HaAutoTreeObject* Parent,const char* tag, int Position)
{
  ParentObject=Parent;
  PreviousSibling=NULL;
  NextSibling=NULL;
  TagName=tag;
  if(ParentObject!=NULL)Parent->AddChild(this,Position);
}

HaAutoTreeObject* HaAutoTreeObject::Create(HaAutoTreeObject* Parent,const char* tag, int Position)
{
  return new HaAutoTreeObject(Parent,tag,Position);
}

HaAutoTreeObject* HaAutoTreeObject::CreateChild(const char* Tag, int Position)
{
  CreateHaAutoTreeObject CHO=ChildrenBuildMap[Tag];
  if(CHO!=NULL) return (*CHO)(this,Tag,Position);
  else return NULL;
}

HaAutoTreeObject::~HaAutoTreeObject()
{
  Clear();
}
void HaAutoTreeObject::Clear()
{
  //std::cout<<"MyName"<<TagName<<"\n";
  //delete children
  int i;
  while(Children.size()!=0)
  {
    delete Children[0];
  }
  /*for(i=0;i<Children.size();i++)
  {
    //printf("%d %s Children.size=%d %p\n",i,TagName.c_str(),Children.size(),Children[i]);
    //HaAutoTreeObject* obj=Children[i];
    //delete obj;
    delete Children[i];
    i--;
  }*/
  //connect Siblings
  if(NextSibling!=NULL)
  {
    if(PreviousSibling!=NULL)
    {
      NextSibling->SetPreviousSibling(PreviousSibling);
      PreviousSibling->SetNextSibling(NextSibling);
    }
    else
      NextSibling->SetPreviousSibling(NULL);
  }
  else
  {
    if(PreviousSibling!=NULL)
    {
      PreviousSibling->SetNextSibling(NULL);
    }
  }
  //Remove from parent
  if(ParentObject!=NULL)
  {
    for(i=0;i<ParentObject->Children.size();i++)
    {
      if(ParentObject->Children[i]==this)
      {
        //std::vector<HaAutoTreeObject*>::iterator removeIter;
        ParentObject->Children.erase(ParentObject->Children.begin()+i);
        break;
      }
    }
  }
}

void HaAutoTreeObject::Save(const char* filename)
{
  TiXmlDocument Doc(filename);
  
  TiXmlDeclaration Decl("1.0","UTF-8","yes");
  Doc.InsertEndChild(Decl);
  
  TiXmlElement* RootElt=GetObject();
  
  Doc.InsertEndChild(*RootElt);
  
  Doc.SaveFile();
  
}
void HaAutoTreeObject::Load(const char* filename)
{
  TiXmlDocument Doc(filename);
  Doc.LoadFile();
  SetObject(Doc.RootElement());
}
void HaAutoTreeObject::SetObject(TiXmlElement* Elt)
{
  if(Elt!=NULL)
  {
    TagName=Elt->Value();
    TiXmlElement* ChildElt;
    for(ChildElt=Elt->FirstChildElement();ChildElt!=NULL;ChildElt=ChildElt->NextSiblingElement())
    {
      HaAutoTreeObject* Obj=CreateChild(ChildElt->Value());
      if(Obj!=NULL)Obj->SetObject(ChildElt);
    }
  }
}
TiXmlElement* HaAutoTreeObject::GetObject()
{
  int i;
  TiXmlElement* Elt= new TiXmlElement(TagName.c_str());
  for(i=0;i<Children.size();i++)
  {
    TiXmlElement* ChildElt=Children[i]->GetObject();
    Elt->LinkEndChild(ChildElt);
    //delete ChildElt;
  }
  return Elt;
}

void HaAutoTreeObject::AddChild(HaAutoTreeObject* Obj, int Position)
{
  //AddChild
  if(Position==-1)//push_back I gess most popular
  {
    Children.push_back(Obj);
    if(Children.size()>1)
    {
      Children[Children.size()-2]->SetNextSibling(Children[Children.size()-1]);
      Children[Children.size()-1]->SetPreviousSibling(Children[Children.size()-2]);
    }
  }
  else 
  {
    std::vector<HaAutoTreeObject*>::iterator iter;
    if(Position>-1)
    {
      iter=Children.insert(Children.begin()+Position,Obj);
    }
    else
    {
      iter=Children.insert(Children.end()-Position+1,Obj);
    }
    //Set Siblings of this Child
    if(iter==Children.begin())
    {
      (*iter)->SetNextSibling(*(iter+1));
      (*(iter+1))->SetPreviousSibling(*iter);
    }
    else if(iter==Children.end())
    {
      (*iter)->SetPreviousSibling(*(iter-1));
      (*(iter-1))->SetNextSibling(*iter);
    }
    else
    {
      (*iter)->SetNextSibling(*(iter+1));
      (*iter)->SetPreviousSibling(*(iter-1));
      (*(iter-1))->SetNextSibling(*iter);
      (*(iter+1))->SetPreviousSibling(*iter);
    }
  }
}
void HaAutoTreeObject::AppendChild(TiXmlElement* Elt)
{
  HaAutoTreeObject* Obj=CreateChild(Elt->Value());
  if(Obj!=NULL)Obj->SetObject(Elt);
}
void HaAutoTreeObject::Print(FILE *out)
{
  TiXmlElement* Elt=GetObject();
  if(Elt!=NULL)Elt->Print(out,0);
  delete Elt;
}
HaAutoTreeObject* HaAutoTreeObject::GetChild(const char* cTag,bool recursevly)
{
  int i;
  HaAutoTreeObject* Child;
  std::string Tag(cTag);
  for(i=0;i<Children.size();i++)
  {
    if(Tag==Children[i]->GetTagName())return Children[i];
    if(recursevly)if((Child=Children[i]->GetChild(Tag.c_str(),recursevly))!=NULL)return Child;
  }
  return NULL;
}
HaAutoTreeObject* HaAutoTreeObject::GetFirstChild(const char* Tag)
{
  if(Children.size()>0)return Children[0];
  return NULL;
}
std::vector<HaAutoTreeObject*> HaAutoTreeObject::GetChildList(bool recursevly)
{
  std::vector<HaAutoTreeObject*> ChildList;
  int i;
  for(i=0;i<Children.size();i++)
  {
    ChildList.push_back(Children[i]);
    if(recursevly)
    {
      std::vector<HaAutoTreeObject*> ChildList2=Children[i]->GetChildList(recursevly);
      if(ChildList2.size()>0)ChildList.insert(ChildList.end(), ChildList2.begin(), ChildList2.end());
    }
  }
  return ChildList;
}
std::vector<HaAutoTreeObject*> HaAutoTreeObject::GetChildList(const char* cTag,bool recursevly)
{
  std::vector<HaAutoTreeObject*> ChildList;
  std::string Tag(cTag);
  int i;
  if(std::string("")==Tag)
  {
    for(i=0;i<Children.size();i++)
    {
      ChildList.push_back(Children[i]);
      if(recursevly)
      {
        std::vector<HaAutoTreeObject*> ChildList2=Children[i]->GetChildList(Tag.c_str(),recursevly);
        if(ChildList2.size()>0)ChildList.insert(ChildList.end(), ChildList2.begin(), ChildList2.end());
      }
    }
  }
  else
  {
	for(i=0;i<Children.size();i++)
    {
      if(Tag==Children[i]->GetTagName())ChildList.push_back(Children[i]);
      if(recursevly)
	  {
        std::vector<HaAutoTreeObject*> ChildList2=Children[i]->GetChildList(Tag.c_str(),recursevly);
        ChildList.insert(ChildList.end(), ChildList2.begin(), ChildList2.end());
	  }
	}
  }
  return ChildList;
}
void HaXML::ReplaceCommaWithSpace(std::string * str)
{
  int i;
  for(i=0;i<str->size();i++)
  {
    if((*str)[i]==',')(*str)[i]=' ';
  }
}
//Single Variables
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,int i)
{
  char buffer[33];
  sprintf(buffer,"%d\0",i);
  Elt->SetAttribute(AtrName,buffer);
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,unsigned int i)
{
  char buffer[33];
  sprintf(buffer,"%d\0",i);
  Elt->SetAttribute(AtrName,buffer);
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,float f,const char* Format)
{
  char buffer[33];
  sprintf(buffer,Format,f);
  Elt->SetAttribute(AtrName,buffer);
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,double d,const char* Format)
{
  char buffer[33];
  sprintf(buffer,Format,d);
  Elt->SetAttribute(AtrName,buffer);
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,bool b)
{
  if(b==true)Elt->SetAttribute(AtrName,"true");
  else Elt->SetAttribute(AtrName,"false");
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,std::string str)
{
  Elt->SetAttribute(AtrName,str.c_str());
  return EXIT_SUCCESS;
}
int HaXML::SetAtribute(TiXmlElement* Elt,const char* AtrName,const char* str)
{
  Elt->SetAttribute(AtrName,str);
  return EXIT_SUCCESS;
}
//Set of Var
int HaXML::SetAtributeV(TiXmlElement* Elt,const char* AtrName,int* ivar,int n)
{
  std::ostringstream outs;
  int j;
  if(n<0)n=sizeof (ivar)/sizeof(int);
  for(j=0;j<n;j++)
  {
    outs<<ivar[j];
    if(j!=n-1)outs<<",";
  }
  Elt->SetAttribute(AtrName,outs.str().c_str());
  return EXIT_SUCCESS;
}
int HaXML::SetAtributeV(TiXmlElement* Elt,const char* AtrName,unsigned int* ivar,int n)
{
  std::ostringstream outs;
  int j;
  if(n<0)n=sizeof (ivar)/sizeof(unsigned int);
  for(j=0;j<n;j++)
  {
    outs<<ivar[j];
    if(j!=n-1)outs<<",";
  }
  Elt->SetAttribute(AtrName,outs.str().c_str());
  return EXIT_SUCCESS;
}
int HaXML::SetAtributeV(TiXmlElement* Elt,const char* AtrName,float* f,int n,const char* Format)
{
  char buffer[33];
  std::ostringstream outs;
  int j;
  if(n<0)n=sizeof (f)/sizeof(float);
  for(j=0;j<n;j++)
  {
    sprintf(buffer,Format,f[j]);
    outs<<buffer;
    if(j!=n-1)outs<<",";
  }
  Elt->SetAttribute(AtrName,outs.str().c_str());
  return EXIT_SUCCESS;
}
int HaXML::SetAtributeV(TiXmlElement* Elt,const char* AtrName,double* d,int n,const char* Format)
{
  char buffer[33];
  std::ostringstream outs;
  int j;
  if(n<0)n=sizeof (d)/sizeof(double);
  for(j=0;j<n;j++)
  {
    sprintf(buffer,Format,d[j]);
    outs<<buffer;
    if(j!=n-1)outs<<",";
  }
  Elt->SetAttribute(AtrName,outs.str().c_str());
  return EXIT_SUCCESS;
}
int HaXML::SetAtributeV(TiXmlElement* Elt,const char* AtrName,bool* b,int n)
{
  std::ostringstream outs;
  int j;
  if(n<0)n=sizeof (b)/sizeof(bool);
  for(j=0;j<n;j++)
  {
    if(b[j]==true)outs<<"true";
    else outs<<"false";
    if(j!=n-1)outs<<",";
  }
  Elt->SetAttribute(AtrName,outs.str().c_str());
  //
  return EXIT_SUCCESS;
}

//set vectors in Element
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector<int>* Vec)
{
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)outs<<(*Vec)[i]<<" ";
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector<string>* Vec)
{
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)outs<<(*Vec)[i]<<" ";
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,int** Vec,int n)
{
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i;
  if(n<0)n=sizeof (Vec)/sizeof(int);
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)outs<<(*Vec)[i]<<" ";
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector<float>* Vec,const char* Format)
{
  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,(*Vec)[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,float** Vec,int n,const char* Format)
{
  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i;
  if(n<0)n=sizeof (Vec)/sizeof(float);
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,(*Vec)[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector<double>* Vec,const char* Format)
{
  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,(*Vec)[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,double* Vec,int n,const char* Format)
{
  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i;
  if(n<0)n=sizeof (Vec)/sizeof(double);
  if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,Vec[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector< vector<float> >* VVec,const char* Format)
{
  char buffer[33];
  int j;
  for(j=0;j<VVec->size();j++)
  {
    TiXmlElement TextElt(TagName);
    std::ostringstream outs;
    int i,n=(*VVec)[j].size();
    
    if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
    {
      for(i=0;i<n;i++)
      {
        sprintf(buffer,Format,(*VVec)[j][i]);
        outs<<buffer<<" ";
      }
      
      TiXmlText Text(outs.str().c_str());
      
      TextElt.InsertEndChild(Text);
      Elt->InsertEndChild(TextElt);
    }
    else return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,std::vector< vector<double> >* VVec,const char* Format)
{
  char buffer[33];
  int j;
  for(j=0;j<VVec->size();j++)
  {
    TiXmlElement TextElt(TagName);
    std::ostringstream outs;
    int i,n=(*VVec)[j].size();
    
    if(SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
    {
      for(i=0;i<n;i++)
      {
        sprintf(buffer,Format,(*VVec)[j][i]);
        outs<<buffer<<" ";
      }
      
      TiXmlText Text(outs.str().c_str());
      
      TextElt.InsertEndChild(Text);
      Elt->InsertEndChild(TextElt);
    }
    else return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
int HaXML::SetTableElement(TiXmlElement* Elt,const char* TagName,const char* cFormat,const char* note,int n,...)
{
  char buffer[1024];
  TiXmlElement TableElt(TagName);
  SetAtribute(&TableElt,"Note",note);
  SetAtribute(&TableElt,"Format",cFormat);
  SetAtribute(&TableElt,"N",n);
  
  
  
  
  
  //up to 32 rows
  int i=1,NRows=0;
  char typev[32];
  std::string Format(cFormat);
  while(Format[i]!='\0')
  {
    switch(Format[i])
    {
      case 'd':
        //if(Format[i-1]=='l')typev[NRows]='D';
        //else 
        typev[NRows]='d';
        NRows++;
        break;
      case 'f':
        if(Format[i-1]=='l')typev[NRows]='F';
        else typev[NRows]='f';
        NRows++;
        break;
      case 'e':
        if(Format[i-1]=='l')typev[NRows]='F';
        else typev[NRows]='f';
        NRows++;
        break;
      case 'g':
        if(Format[i-1]=='l')typev[NRows]='F';
        else typev[NRows]='f';
        NRows++;
        break;
      case 's':
        typev[NRows]='s';
        NRows++;
        break;
      case 'S':
        typev[NRows]='S';
        Format[i]='s';
        NRows++;
        break;
    }
    i++;
  }
  //built row's format vector
  std::vector<std::string> FormatV;
  std::istringstream istr(Format,std::istringstream::in);
  while(!istr.eof())
  {
    std::string str;
    istr>>str;
    FormatV.push_back(str);
  }
  
  std::vector< std::vector<int>* > IntP;
  std::vector< std::vector<float>* > FloatP;
  std::vector< std::vector<double>* > DoubleP;
  std::vector< std::vector<std::string>* > StringP;
  IntP.resize(NRows);
  FloatP.resize(NRows);
  DoubleP.resize(NRows);
  StringP.resize(NRows);
  
  for(i=0;i<NRows;i++)
  {
    IntP[i]=NULL;
    FloatP[i]=NULL;
    DoubleP[i]=NULL;
  }
  va_list arg_list;
  va_start(arg_list,n);
  for(i=0;i<NRows;i++)
  {
    switch(typev[i])
    {
      case 'd':
        IntP[i]=va_arg(arg_list,std::vector<int>*);
        break;
      case 'f':
        FloatP[i]=va_arg(arg_list,std::vector<float>*);
        break;
      case 'F':
        DoubleP[i]=va_arg(arg_list,std::vector<double>*);
        break;
      case 'S':
        StringP[i]=va_arg(arg_list,std::vector<std::string>*);
        break;
    }
  }
  va_end(arg_list);
  std::ostringstream outs;
  int j;
  outs<<"\n";
  for(j=0;j<n;j++)
  {
    for(i=0;i<NRows;i++)
    {
      switch(typev[i])
      {
        case 'd':
          sprintf(buffer,FormatV[i].c_str(),(*IntP[i])[j]);
          //arg_list[i]=;
          break;
        case 'f':
          //arg_list[i]=(*FloatP[i])[j];
          sprintf(buffer,FormatV[i].c_str(),(*FloatP[i])[j]);
          break;
        case 'F':
          //arg_list[i]=(*DoubleP[i])[j];
          sprintf(buffer,FormatV[i].c_str(),(*DoubleP[i])[j]);
          break;
        case 'S':
          //arg_list[i]=(*DoubleP[i])[j];
          sprintf(buffer,FormatV[i].c_str(),(*StringP[i])[j].c_str());
          break;
      }
      outs<<buffer<<" ";
    }
    outs<<"\n";
  }
  TiXmlText Text(outs.str().c_str());
  TableElt.InsertEndChild(Text);
  Elt->InsertEndChild(TableElt);
  return EXIT_SUCCESS;
}
int HaXML::SetText(TiXmlElement* Elt,std::vector<int>* Vec)
{
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(HaXML::SetAtribute(Elt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)outs<<(*Vec)[i]<<" ";
    
    TiXmlText Text(outs.str().c_str());
    
    Elt->InsertEndChild(Text);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
//Get Single Var
std::string HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)return std::string(c_str);
  else return std::string("");
}
//int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,int* i)
//{
//  return HaXML::GetAtributeI(Elt,AtrName,i);
//}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,int* i)
{
  const char* Str=Elt->Attribute(AtrName);
  if(Str==NULL)
  {
    return EXIT_FAILURE;
  }
  else
  {
    *i=atoi(Str);
    return EXIT_SUCCESS;
  }
}
//int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,float* f)
//{
 // return HaXML::GetAtributeF(Elt,AtrName,f);
//}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,float* f)
{
  const char* Str=Elt->Attribute(AtrName);
  if(Str==NULL)
  {
    return EXIT_FAILURE;
  }
  else
  {
    *f=(float)atof(Str);
    return EXIT_SUCCESS;
  }
}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,double* d)
{
  const char* Str=Elt->Attribute(AtrName);
  if(Str==NULL)
  {
    return EXIT_FAILURE;
  }
  else
  {
    *d=atof(Str);
    return EXIT_SUCCESS;
  }
}
int HaXML::GetBool(const char* cstr,bool* b)
{
  std::string str(cstr);
  if(str=="true")
  {
    *b=true;
    return EXIT_SUCCESS;
  }
  else if (str=="false")
  {
    *b=false;
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,bool* b)
{
  const char* Str=Elt->Attribute(AtrName);
  if(Str!=NULL)
  {
    if(HaXML::GetBool(Str,b))return EXIT_SUCCESS;
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE;
}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,std::string* str)
{
  
  const char* Str=Elt->Attribute(AtrName);
  if(Str==NULL)
  {
    return EXIT_FAILURE;
  }
  else
  {
    *str=Str;
    return EXIT_SUCCESS;
  }
}
int HaXML::GetAtribute(TiXmlElement* Elt,const char* AtrName,char* str)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)
  {
    strcpy(str,c_str);
  }
  return EXIT_FAILURE;
}
//HaXML::Get of Vector of Var
int HaXML::GetAtributeV(TiXmlElement* Elt,const char* AtrName,int* ivar,int n)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)
  {
    std::string Str(c_str);
    ReplaceCommaWithSpace(&Str);
    
    std::istringstream ins(Str);
    
    int i;
    if(n<0)n=sizeof(ivar)/sizeof(int);
    for(i=0;i<n;i++)
    {
      ins>>ivar[i];
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetAtributeV(TiXmlElement* Elt,const char* AtrName,float* ivar,int n)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)
  {
    std::string Str(c_str);
    ReplaceCommaWithSpace(&Str);
    
    std::istringstream ins(Str);
    
    int i;
    if(n<0)n=sizeof(ivar)/sizeof(int);
    for(i=0;i<n;i++)
    {
      ins>>ivar[i];
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetAtributeV(TiXmlElement* Elt,const char* AtrName,double* ivar,int n)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)
  {
    std::string Str(c_str);
    ReplaceCommaWithSpace(&Str);
    
    std::istringstream ins(Str);
    
    int i;
    if(n<0)n=sizeof(ivar)/sizeof(int);
    for(i=0;i<n;i++)
    {
      ins>>ivar[i];
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetAtributeV(TiXmlElement* Elt,const char* AtrName,bool* bvar,int n)
{
  const char* c_str=Elt->Attribute(AtrName);
  if(c_str!=NULL)
  {
    std::string Str(c_str);
    char buff[32];
    ReplaceCommaWithSpace(&Str);
    
    std::istringstream ins(Str);
    int i;
    if(n<0)n=sizeof(bvar)/sizeof(int);
    for(i=0;i<n;i++)
    {
      ins>>buff;
      if(HaXML::GetBool(buff,&bvar[i])==EXIT_FAILURE)return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,int** Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    int n;
    if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
    {
	  if(n>0)
	  {
		TiXmlText* Text;
		TiXmlNode* XmlNode=TextElt->FirstChild();
		while(XmlNode!=NULL)
		{
			if((Text=XmlNode->ToText())!=NULL)break;
			TextElt->IterateChildren(XmlNode);
		}
		if(Text!=NULL)
		{
			std::istringstream ins(Text->Value());
			//if((*Vec).size()<n)(*Vec).resize(n);
			int *v=new int[n];
			int i=0;
			for(i=0;i<n;i++)
				ins>>v[i];
			
			*Vec=v;
        }
		else
		{
			return EXIT_FAILURE;
		}
	  }
	}
  }
  else
  {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,std::vector<double>* Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*Vec).size()<n)(*Vec).resize(n);
        int i=0;
        while(!ins.eof())
        {
          ins>>(*Vec)[i];
          i++;
        }
      }
      else
      {
        int i=0;
        while(!ins.eof())
        {
          int j;
          ins>>j;
          if((*Vec).size()<i+1)(*Vec).resize((*Vec).size()+1);
          (*Vec)[i]=j;
          i++;
        }
      }
    return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,std::vector< std::vector<float> >* VVec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  int k=0;
  while(TextElt!=NULL)
  {
    VVec->resize(VVec->size()+1);
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*VVec)[k].size()<n)(*VVec)[k].resize(n);
        int i=0;
        while(!ins.eof())
        {
          ins>>(*VVec)[k][i];
          i++;
        }
      }
      else
      {
        int i=0;
        while(!ins.eof())
        {
          int j;
          ins>>j;
          if((*VVec)[k].size()<i+1)(*VVec)[k].resize((*VVec)[k].size()+1);
          (*VVec)[k][i]=j;
          i++;
        }
      }
    }
    TextElt=TextElt->NextSiblingElement(TagName);
    k++;
  }
  return EXIT_SUCCESS;
}

int HaXML::GetVVIntElement(const TiXmlElement* Elt,const char* TagName,std::vector< std::vector<int> >* VVec)
{
	TiXmlElement* TextElt=(TiXmlElement*)Elt->FirstChildElement(TagName);
	int k=0;
	while(TextElt!=NULL)
	{
		VVec->resize(VVec->size()+1);
		TiXmlText* Text;
		TiXmlNode* XmlNode=TextElt->FirstChild();
		while(XmlNode!=NULL)
		{
			if((Text=XmlNode->ToText())!=NULL)break;
			TextElt->IterateChildren(XmlNode);
		}
		if(Text!=NULL)
		{
			std::istringstream ins(Text->Value());
			int n;
			if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
			{
				if((*VVec)[k].size()<n)(*VVec)[k].resize(n);
				int i=0;
				while(!ins.eof())
				{
					ins>>(*VVec)[k][i];
					i++;
				}
			}
			else
			{
				int i=0;
				while(!ins.eof())
				{
					int j;
					ins>>j;
					if((*VVec)[k].size()<i+1)(*VVec)[k].resize((*VVec)[k].size()+1);
					(*VVec)[k][i]=j;
					i++;
				}
			}
		}
		TextElt=TextElt->NextSiblingElement(TagName);
		k++;
	}
	return EXIT_SUCCESS;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,std::vector< std::vector<double> >* VVec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  int k=0;
  while(TextElt!=NULL)
  {
    VVec->resize(VVec->size()+1);
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*VVec)[k].size()<n)(*VVec)[k].resize(n);
        int i=0;
        while(!ins.eof())
        {
          ins>>(*VVec)[k][i];
          i++;
        }
      }
      else
      {
        int i=0;
        while(!ins.eof())
        {
          int j;
          ins>>j;
          if((*VVec)[k].size()<i+1)(*VVec)[k].resize((*VVec)[k].size()+1);
          (*VVec)[k][i]=j;
          i++;
        }
      }
    }
    TextElt=TextElt->NextSiblingElement(TagName);
    k++;
  }
  return EXIT_SUCCESS;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,std::vector<float>* Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*Vec).size()<n)(*Vec).resize(n);
        int i=0;
        while(!ins.eof())
        {
          ins>>(*Vec)[i];
          i++;
        }
      }
      else
      {
        int i=0;
        while(!ins.eof())
        {
          int j;
          ins>>j;
          if((*Vec).size()<i+1)(*Vec).resize((*Vec).size()+1);
          (*Vec)[i]=j;
          i++;
        }
      }
    return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,vector<string>* Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*Vec).size()<n)(*Vec).resize(n);
        int i=0;
        for(i=0;i<n;i++)
        {
          ins>>(*Vec)[i];
        }
      }
      else
      {
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,float** Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        //if((*Vec).size()<n)(*Vec).resize(n);
        float *v=new float[n];
        int i=0;
        for(i=0;i<n;i++)
          ins>>v[i];
         
        *Vec=v;
      }
      else
        return EXIT_FAILURE;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
  return EXIT_FAILURE;
}

int HaXML::GetTableElement(TiXmlElement* Elt,const char* TagName,...)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::string note, Format;
      int n;
      HaXML::GetAtribute(TextElt,"Note",&note);
      HaXML::GetAtribute(TextElt,"Format",&Format);
      HaXML::GetAtribute(TextElt,"N",&n);
  
  //built row's format vector
      std::vector<std::string> FormatV;
      std::istringstream istr(Format,std::istringstream::in);
      while(!istr.eof())
      {
        std::string str;
        istr>>str;
        FormatV.push_back(str);
      }
  //up to 32 rows
      int i=1,NRows=0;
      char typev[32];
  
      while(Format[i]!='\0')
      {
        switch(Format[i])
        {
          case 'd':
        //if(Format[i-1]=='l')typev[NRows]='D';
        //else 
            typev[NRows]='d';
            NRows++;
            break;
          case 'f':
            if(Format[i-1]=='l')typev[NRows]='F';
            else typev[NRows]='f';
            NRows++;
            break;
          case 'e':
            if(Format[i-1]=='l')typev[NRows]='F';
            else typev[NRows]='f';
            NRows++;
            break;
          case 'g':
            if(Format[i-1]=='l')typev[NRows]='F';
            else typev[NRows]='f';
            NRows++;
            break;
          case 's':
            typev[NRows]='s';
            NRows++;
            break;
          case 'S':
            typev[NRows]='S';
            Format[i]='s';
            NRows++;
            break;
        }
        i++;
      }
      std::vector< std::vector<int>* > IntP;
      std::vector< std::vector<float>* > FloatP;
      std::vector< std::vector<double>* > DoubleP;
      std::vector< std::vector<std::string>* > StringP;
      IntP.resize(NRows);
      FloatP.resize(NRows);
      DoubleP.resize(NRows);
      StringP.resize(NRows);
      for(i=0;i<NRows;i++)
      {
        IntP[i]=NULL;
        FloatP[i]=NULL;
        DoubleP[i]=NULL;
        StringP[i]=NULL;
      }
      va_list arg_list;
      va_start(arg_list,TagName);
      for(i=0;i<NRows;i++)
      {
        switch(typev[i])
        {
          case 'd':
            IntP[i]=va_arg(arg_list,std::vector<int>*);
            (*IntP[i]).resize(n);
            break;
          case 'f':
            FloatP[i]=va_arg(arg_list,std::vector<float>*);
            (*FloatP[i]).resize(n);
            break;
          case 'F':
            DoubleP[i]=va_arg(arg_list,std::vector<double>*);
            (*DoubleP[i]).resize(n);
            break;
          case 'S':
            StringP[i]=va_arg(arg_list,std::vector<std::string>*);
            (*StringP[i]).resize(n);
            break;
        }
      }
      va_end(arg_list);
	  if(n>0)
	  {
       std::istringstream ins(Text->Value());
  
       int j;
			 double dtmp;
       for(j=0;j<n;j++)
       {
        for(i=0;i<NRows;i++)
        {
          switch(typev[i])
          {
            case 'd':
              ins>>(*IntP[i])[j];
              break;
            case 'f':
							ins>>dtmp;
							(*FloatP[i])[j]=dtmp;
              break;
            case 'F':
              ins>>(*DoubleP[i])[j];
              break;
            case 'S':
              ins>>(*StringP[i])[j];
              break;
          }
        }
       }
	  }
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::GetText(TiXmlElement* Elt,std::vector<int>* Vec)
{
  std::istringstream ins(Elt->GetText());
  int n;
  if(HaXML::GetAtribute(Elt,"N",&n)==EXIT_SUCCESS)
  {
    if((*Vec).size()<n)(*Vec).resize(n);
    int i=0;
    while(!ins.eof())
    {
      ins>>(*Vec)[i];
      i++;
    }
  }
  else
  {
    int i=0;
    while(!ins.eof())
    {
      int j;
      ins>>j;
      if((*Vec).size()<i+1)(*Vec).resize((*Vec).size()+1);
      (*Vec)[i]=j;
      i++;
    }
  }
  return EXIT_SUCCESS;
}

int HaXML::GetStrIndex(std::vector<std::string> StrList,std::string StrToComp)
{
  int i;
  for(i=0;i<StrList.size();i++)
    if(StrList[i]==StrToComp)return i;
  return -1;
}

std::string HaXML::GetStrByIndex(std::vector<std::string> StrList,int Index)
{
  return StrList[Index];
}

#ifndef HALIB_BASE
int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,HaVec_float* Vec,const char* Format)
{

  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(HaXML::SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,(*Vec)[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,HaVec_float* Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        if((*Vec).size()<n)(*Vec).newsize(n);
        int i=0;
        for(i=0;i<n;i++)
        {
          ins>>(*Vec)[i];
        }
      }
      else
      {
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,HaVec_double* Vec,const char* Format)
{

  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,n=Vec->size();
  
  if(HaXML::SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    for(i=0;i<n;i++)
    {
      sprintf(buffer,Format,(*Vec)[i]);
      outs<<buffer<<" ";
    }
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,HaVec_double* Vec)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS)
      {
        //if((*Vec).size()<n)
        Vec->newsize(n);
        int i=0;
        for(i=0;i<n;i++)
        {
          ins>>(*Vec)[i];
        }
      }
      else
      {
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

int HaXML::SetElement(TiXmlElement* Elt,const char* TagName,HaMat_double* Mat,const char* Format)
{
  char buffer[33];
  TiXmlElement TextElt(TagName);
  std::ostringstream outs;
  int i,j,n=Mat->num_cols(),m=Mat->num_rows();
  
  if(HaXML::SetAtribute(&TextElt,"M",m)==EXIT_SUCCESS&&HaXML::SetAtribute(&TextElt,"N",n)==EXIT_SUCCESS)
  {
    outs<<"\n";
    for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
      {
        sprintf(buffer,Format,Mat->r0(i,j));
        outs<<buffer<<" ";
      }
      outs<<"\n";
    }
    outs<<"\n";
    
    TiXmlText Text(outs.str().c_str());
    
    TextElt.InsertEndChild(Text);
    Elt->InsertEndChild(TextElt);
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int HaXML::GetElement(TiXmlElement* Elt,const char* TagName,HaMat_double* Mat)
{
  TiXmlElement* TextElt=Elt->FirstChildElement(TagName);
  
  if(TextElt!=NULL)
  {
    TiXmlText* Text;
    TiXmlNode* XmlNode=TextElt->FirstChild();
    while(XmlNode!=NULL)
    {
      if((Text=XmlNode->ToText())!=NULL)break;
      TextElt->IterateChildren(XmlNode);
    }
    if(Text!=NULL)
    {
      std::istringstream ins(Text->Value());
      int n,m;
      if(HaXML::GetAtribute(TextElt,"N",&n)==EXIT_SUCCESS&&HaXML::GetAtribute(TextElt,"M",&m)==EXIT_SUCCESS)
      {
        Mat->newsize(m,n);
        int i,j;
        for(i=0;i<m;i++)
        {
          for(j=0;j<n;j++)
            ins>>Mat->r0(i,j);
        }
        /*
        while(!ins.eof())
        {
          ins>>(*Vec)[i];
          i++;
      }*/
      }
      else
      {
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }
  }
  else
  {
    return EXIT_FAILURE;
  }
}

TiXmlElement* BldXmlEltFromCstr(const char* Cstr)
{
  TiXmlElement* Elt=new TiXmlElement("TiXmlElement");
  istringstream ins(Cstr);
  ins>>*Elt;
  return Elt;
}

#endif


#ifndef MAP_IO_STRING_LENGTH
#define MAP_IO_STRING_LENGTH 512
#endif
#ifndef MAP_IO_FORCE
#define MAP_IO_FORCE 1
#endif
#ifndef MAP_IO_ADD
#define MAP_IO_ADD 2
#endif

//! Write float array in two column format, i.e. value number_of_repetition
int HaWriteMapGZinTwoColumns(gzFile file, float * fmap, int GSXYZ, float coef)
{
	int count;
	float current;
	int gridPoint;

	//Write
	count = 0;
	current = fmap[0];
	for (gridPoint = 0; gridPoint < GSXYZ; gridPoint++) {
		if (current == fmap[gridPoint])
			count++;
		else {
			if (!gzprintf(file, "%.7e %d\n", current*coef, count)) {
				fprintf(stderr, "ERROR 106: Problems writing to file\n");
				return EXIT_FAILURE;
			}
			count = 1;
			current = fmap[gridPoint];
		}
	}
	if (!gzprintf(file, "%.7e %d\n", current*coef, count)) {
		fprintf(stderr, "ERROR 106: Problems writing to file\n");
		exit(106);
	}
	return EXIT_SUCCESS;
}

int HaReadMapGZOneColumns(gzFile file, float * nmap, unsigned int N, float coef)
{
	float current;
	unsigned int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	unsigned int i;

	gridPoint = 0;
	while (gzgets(file, str, MAP_IO_STRING_LENGTH) != Z_NULL && gridPoint < N)
	{
		nmap[gridPoint] = coef * atof(str);
		gridPoint++;
	}

	return EXIT_SUCCESS;
}

int HaWriteMapGZ(const char *filename, TiXmlElement *header, float * fmap, int* gridsize, float coef, char* Comments)
{
	//Build XML Header
	header->SetArrOfIntAttribute("GridSize", gridsize, 3);
	header->SetAttribute("Comments", Comments);
	header->SetIntAttribute("Columns", 2);

	gzFile file;
	file = gzopen(filename, "wb");
	if (file == NULL) {
		fprintf(stderr, "ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	string StrHeader;
	StrHeader << *header;
	gzprintf(file, "%s\n", StrHeader.c_str());
	//header.Print(file);
	int GSXYZ = gridsize[0] * gridsize[1] * gridsize[2];
	HaWriteMapGZinTwoColumns(file, fmap, GSXYZ, coef);
	gzclose(file);
	return EXIT_SUCCESS;
}


int HaReadMapGZ(gzFile file, float *fmap, int GSXYZ, float coef, int Columns)
{
	int count;
	float current;
	int gridPoint;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int coll = MAP_IO_FORCE;

	if (Columns == 1)
	{
		HaReadMapGZOneColumns(file, fmap, GSXYZ, coef);
	}
	else
	{
		gridPoint = 0;
		while (gridPoint < GSXYZ)
		{
			if (Z_NULL == gzgets(file, str, MAP_IO_STRING_LENGTH))
			{
				fprintf(stderr, "readMap: problem reading from file\n");
				return EXIT_FAILURE;
			}
			if (2 != sscanf(str, "%f %i\n", &current, &count))
			{
				fprintf(stderr, "readMap: unkown format\n");
				return EXIT_FAILURE;
			}
			//gzscanf(file,"%f %i\n",&current,&count);
			if (gridPoint + count > GSXYZ) {
				fprintf(stderr, "readMap: too many points file\n");
				return EXIT_FAILURE;
			}
			for (i = gridPoint; i < gridPoint + count; i++) {
				if (coll == MAP_IO_ADD)
					fmap[i] += current * coef;
				else if (coll == MAP_IO_FORCE)
					fmap[i] = current * coef;
				else {
					if (fmap[i] == -1)
						fmap[i] = current * coef;
				}
			}
			gridPoint += count;
		}
	}
	return EXIT_SUCCESS;
}

int HaReadMapGZ(const char *filename, TiXmlElement** pheader, float **pfmap, int* gridsize, float coef)
{
	//! Read Map from file and return values *fmap
	gzFile file;
	char str[MAP_IO_STRING_LENGTH];
	int i;
	int GS[3], GSXYZ;
	int Columns;
	float *vfield = *pfmap;
	file = gzopen(filename, "rb");
	if (file == NULL) {
		fprintf(stderr, "ERROR 102: Can not open file %s\n", filename);
		return EXIT_FAILURE;
	}
	//Header
	if (Z_NULL == gzgets(file, str, MAP_IO_STRING_LENGTH)) {
		fprintf(stderr, "readMap: problem reading from file: %s\n", filename);
		return EXIT_FAILURE;
	}
	//string StrHeader(string);
	istringstream ins(str);
	TiXmlElement *header = new TiXmlElement("Field");
	*pheader = header;
	ins >> *header;
	header->GetArrOfIntAttribute("GridSize", GS, 3);
	if (header->GetIntAttribute("Columns", &Columns) == EXIT_FAILURE)Columns = 2;
	if (gridsize[0] >= 0)
	{
		if (gridsize[0] != GS[0] || gridsize[1] != GS[1] || gridsize[2] != GS[2])
		{
			fprintf(stderr, "readMap: Grid size of map in file do not coinside with request\n");
			return EXIT_FAILURE;
		}
	}
	else
	{
		gridsize[0] = GS[0];
		gridsize[1] = GS[1];
		gridsize[2] = GS[2];
	}
	GSXYZ = GS[0] * GS[1] * GS[2];

	if (vfield == NULL)vfield = new float[GSXYZ];
	HaReadMapGZ(file, vfield, GSXYZ, coef, Columns);
	*pfmap = vfield;
	gzclose(file);
	return EXIT_SUCCESS;
}