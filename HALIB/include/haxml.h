/*!  \file haxml.h

   Classes for to work with XML files in HARLEM
   
   \author Nikolay Simakov, Igor Kurnikov  
   \date 2006-

*/
#ifndef HAXML_H
#define HAXML_H

#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <cstdlib>

#include <stdarg.h>
#ifndef HALIB_BASE
#ifndef SWIG
#include "halinalg.h"
#endif
//class HaMat_double;
//class HaVec_float;
//class HaVec_double;
#endif

class TiXmlElement;
using namespace std;

class HaAutoTreeObject;
//!pointer on HaObjects creators
typedef HaAutoTreeObject* (*CreateHaAutoTreeObject)(HaAutoTreeObject*,const char*,int);


#include "haobject.h"


TiXmlElement* BldXmlEltFromCstr(const char* Cstr);
//!@brief like HaTreeObject plus more automatic stuff
//!ability to load whole tree by constractin objects using std::map<std::string, CreateHaObject > ChildrenBuildMap;
class HaAutoTreeObject
{
  public:
    //!@brief Default constractor
    /*!
    @param Parent Parent of the object, will call AddChild in the Parent object.
    In case if Parent is NULL then supposed that it is a root object
    @param tag TagName of the object
    @param Position look AddChild
    To be able to create children automaticly developer should built ChildrenBuildMap, which is using in automatic children creation, this function executed in HaObject constractor (via HaObject* Create), and the developer should impliment this function in its own code. Just thouse objects will be created which tag names will be present in ChildrenBuildMap.
    */
    HaAutoTreeObject(HaAutoTreeObject* Parent=NULL,const char* tag="HaObject",int Position=-1);
    //!Object Deconstractor, will remove the object from tree correctly
    //!developer should write his own destructor to be sure that object is deleted
    virtual ~HaAutoTreeObject();
    //!Delete all the children of this object
    void Clear();
    //!This function create object of this type look HaObject(), developer should impliment it
    //!@todo impliment me
    static HaAutoTreeObject* Create(HaAutoTreeObject* Parent=NULL,const char* tag="HaObject",int Position=-1);
    //!Save XML represantation of the object to file
    void Save(const char* filename);
    //!Load object from XML file
    void Load(const char* filename);
    //!@brief Set object from XML representation
    //!(developer should write its own implemintation of this function, look MyObject for an example)
    //!@todo impliment me
    virtual void SetObject(TiXmlElement* Elt);
    //!@brief Set up XML representation of the object
    //! (developer should write its own implemintation of this function, look MyObject for an example)
    //!@todo impliment me
    virtual TiXmlElement* GetObject();
    //!@brief Add Child to Object
    /*!
    @param Obj Child to add
    @param Position Position to where to locate (0,1,2,...). If Positon is negative then locate starting to count from the back: -1 the last, -2 the last but one, etc.
    */
    void AddChild(HaAutoTreeObject* Obj, int Position=-1);
    void AppendChild(TiXmlElement* Elt);
    //!Print the object content in xml format
    void Print(FILE *out=stdout);
    
    void SetPreviousSibling(HaAutoTreeObject* Obj){PreviousSibling=Obj;}
    void SetNextSibling(HaAutoTreeObject* Obj){NextSibling=Obj;}
    HaAutoTreeObject* GetPreviousSibling(){return PreviousSibling;}
    HaAutoTreeObject* GetNextSibling(){return NextSibling;}
    HaAutoTreeObject* GetParentObject(){return ParentObject;}
    virtual HaAutoTreeObject* CreateChild(const char* Tag, int Position=-1);
    
    HaAutoTreeObject* GetChild(const char* Tag,bool recursevly=false);
    HaAutoTreeObject* GetFirstChild(const char* Tag="");
    std::vector<HaAutoTreeObject*> GetChildList(bool recursevly=false);
    std::vector<HaAutoTreeObject*> GetChildList(const char* Tag,bool recursevly=false);
    const char* GetTagName(){return TagName.c_str();}
  protected:
    HaAutoTreeObject* ParentObject;
    HaAutoTreeObject* PreviousSibling;
    HaAutoTreeObject* NextSibling;
    //!Tag Name of the object, shouldn't be the same for whole object
    std::string TagName;
    
  public:
    //! list of pointers to the children
    std::vector<HaAutoTreeObject*> Children;
  protected:
    //!Map of HaObject::Create(...) functions (look CreateChildrenBuildMap()) left is a tag name and right is a pointer to HaObject::Create(...) which is refer to this pointer
    std::map<std::string, CreateHaAutoTreeObject > ChildrenBuildMap;
};

//!@brief functions to work with XML
class HaXML
{
  public:
    HaXML(){}
    virtual ~HaXML(){};
    static void ReplaceCommaWithSpace(std::string * str);
    static int GetBool(const char* cstr,bool* b);
    //Single Variables
    //!Set integer variable as atribute for element Elt with name AtrName
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,int i);
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,unsigned int i);
    //int SetAtributeI(TiXmlElement* Elt,const char* AtrName,int i);
    //!Set float variable as atribute for element Elt with name AtrName and formate Format
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,float f,const char* Format="%.8g\0");
    //int SetAtributeF(TiXmlElement* Elt,const char* AtrName,float f,const char* Format="%.8g\0");
    //!Set double variable as atribute for element Elt with name AtrName and formate Format
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,double d,const char* Format="%.10lg\0");
    //!Set bool variable as atribute for element Elt with name AtrName
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,bool b);
    //!Set string variable as atribute for element Elt  with name AtrName
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,std::string str);
    //!Set char* variable as atribute for element Elt  with name AtrName
    static int SetAtribute(TiXmlElement* Elt,const char* AtrName,const char* str);
    //Set of Variables
    //!Set a vector of integer variables as atribute for element Elt with name AtrName
    //!Size of vector is n if n=-1 then try to use sizeof
    static int SetAtributeV(TiXmlElement* Elt,const char* AtrName,int* i,int n=-1);
    static int SetAtributeV(TiXmlElement* Elt,const char* AtrName,unsigned int* i,int n=-1);
    //!Set a vector float variables as atribute for element Elt with name AtrName and formate Format
    static int SetAtributeV(TiXmlElement* Elt,const char* AtrName,float* f,int n=-1,const char* Format="%.8g\0");
    //!Set a vector double variable as atribute for element Elt with name AtrName and formate Format
    static int SetAtributeV(TiXmlElement* Elt,const char* AtrName,double* d,int n=-1,const char* Format="%.10lg\0");
    //!Set a vector bool atribute as variables for element Elt with name AtrName
    static int SetAtributeV(TiXmlElement* Elt,const char* AtrName,bool* b,int n=-1);
    //Save as element
    //!@brief Set vector on integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector<int>* Vec);
    //!@brief Set array of integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N,Size of vector is n if n=-1 then try to use sizeof
    
                  
    static int SetElement(TiXmlElement* Elt,const char* TagName,int** Vec,int n=-1);
    //!@brief Set vector on integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector<float>* Vec,const char* Format="%.8g\0");
    //!@brief Set array of integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N,Size of vector is n if n=-1 then try to use sizeof
    static int SetElement(TiXmlElement* Elt,const char* TagName,float** Vec,int n=-1,const char* Format="%.8g\0");
    //!@brief Set vector on integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector<double>* Vec,const char* Format="%.10lg\0");
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector<string>* Vec);
    //!@brief Set array of integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N,Size of vector is n if n=-1 then try to use sizeof
    static int SetElement(TiXmlElement* Elt,const char* TagName,double* Vec,int n=-1,const char* Format="%.10lg\0");
    //!@brief Set vector on integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector< vector<float> >* Vec,const char* Format="%.8g\0");
    //!@brief Set vector on integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N
    static int SetElement(TiXmlElement* Elt,const char* TagName,std::vector< vector<double> >* Vec,const char* Format="%.10lg\0");
    //Diffrent stuff
    //!Save Data in table format ,...=vector<T1>,vector<T2>,vector<T3>, 'S' means sdt::string
    static int SetTableElement(TiXmlElement* Elt,const char* TagName,const char* Format,const char* note,int n,...);
    //!Write variables in Text part of Element Elt
    static int SetText(TiXmlElement* Elt,std::vector<int>* Vec);
    //Single Variables
    static std::string GetAtribute(TiXmlElement* Elt,const char* AtrName);
    //!Get integer variable from atribute of element Elt, with name AtrName
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,int* i);
    //int GetAtributeI(TiXmlElement* Elt,const char* AtrName,int* i);
    //!Get float variable from atribute of element Elt, with name AtrName
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,float* f);
    //int GetAtributeF(TiXmlElement* Elt,const char* AtrName,float* f);
    //!Get double variable from atribute of element Elt, with name AtrName
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,double* d);
    //!Get bool variable from atribute of element Elt, with name AtrName
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,bool* b);
    //!Get string variable from atribute of element Elt, with name AtrName
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,std::string* str);
    //!Get char* variable from atribute of element Elt, with name AtrName, don't care about char size, yet...
    static int GetAtribute(TiXmlElement* Elt,const char* AtrName,char* str);
    //Set of Variables
    //!Get a vector of integer variables from atribute of element Elt with name AtrName
    //!Size of vector is n if n=-1 then try to use sizeof
    static int GetAtributeV(TiXmlElement* Elt,const char* AtrName,int* i,int n=-1);
    //!Get a vector of float variables from atribute of element Elt with name AtrName
    static int GetAtributeV(TiXmlElement* Elt,const char* AtrName,float* f,int n=-1);
    //!Get a vector of double variables from atribute of element Elt with name AtrName
    static int GetAtributeV(TiXmlElement* Elt,const char* AtrName,double* d,int n=-1);
    //!Get a vector of bool variables from atribute of element Elt with name AtrName
    static int GetAtributeV(TiXmlElement* Elt,const char* AtrName,bool* b,int n=-1);
    //!@brief Get vector of integer/float/double variables from element with name TagName, which is a child of element Elt
    //!
    //!In case if atribute "N" is present will resize vector,
    //!if not present will read all values if size of readable vector is bigger when Vec will encrease size of Vec,
    //!if smaller will do nothing
    static int GetElement(TiXmlElement* Elt,const char* AtrName,std::vector<double>* Vec);
    //!@brief Get array of integer variables as element with name TagName, which is a child of element Elt
    //!
    static int GetElement(TiXmlElement* Elt,const char* AtrName,std::vector<float>* Vec);
    static int GetElement(TiXmlElement* Elt,const char* AtrName,std::vector<string>* Vec); 
    
    //!Will read what have the smallest array win (reading one or in object)
    static int GetElement(TiXmlElement* Elt,const char* AtrName,int** Vec);
    static int GetElement(TiXmlElement* Elt,const char* TagName,std::vector< std::vector<float> >* VVec);
    static int GetElement(TiXmlElement* Elt,const char* TagName,std::vector< std::vector<double> >* VVec);
    //!@brief Set array of integer variables as element with name TagName, as child of element Elt
    //!
    //!Will write size of array as atribute N,Size of vector is n if n=-1 then try to use sizeof
    static int GetElement(TiXmlElement* Elt,const char* TagName,float** Vec);
		static int GetVVIntElement(const TiXmlElement* Elt,const char* TagName,std::vector< std::vector<int> >* VVec);
				
    #ifndef HALIB_BASE
    static int SetElement(TiXmlElement* Elt,const char* AtrName,HaVec_float* Vec,const char* Format="%.8g\0");
    static int GetElement(TiXmlElement* Elt,const char* AtrName,HaVec_float* Vec);
    static int SetElement(TiXmlElement* Elt,const char* AtrName,HaVec_double* Vec,const char* Format="%.10lg\0");
    static int GetElement(TiXmlElement* Elt,const char* AtrName,HaVec_double* Vec);
    static int SetElement(TiXmlElement* Elt,const char* AtrName,HaMat_double* Mat,const char* Format="%.10lg\0");
    static  int GetElement(TiXmlElement* Elt,const char* AtrName,HaMat_double* Mat);
    #endif
    //Diffrent stuff
    //!Load Data in table format ,...=vector<T1>,vector<T2>,vector<T3>
    static int GetTableElement(TiXmlElement* Elt,const char* TagName,...);
    static int GetText(TiXmlElement* Elt,std::vector<int>* Vec);
    
    static int GetStrIndex(std::vector<std::string> StrList,std::string StrToComp);
    
    static std::string GetStrByIndex(std::vector<std::string> StrList,int Index);
};


#endif // End of define HAXML_H